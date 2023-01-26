###############gene select####################
library(Seurat)
library(tidyverse)
library(car) #package to calculate Variance Inflation Factor
library(corrplot) #correlation plots
library(leaps) #best subsets regression
library(glmnet) #allows ridge regression, LASSO and elastic net
library(reshape2)
library(pROC)
library(rpart) #tissueification and regression trees
library(xgboost) #gradient boosting 
library(caret)
velo<-read.csv("veloGenes.csv",header = T)
velo<-velo[velo$fit_likelihood>0,]
slingshot<-read.table("slingshot.txt",header = T)
monocle3<-read.table("monocle3.gene_noCCL5_all.txt",header = T)
monocle3<-na.omit(monocle3)
Monocle2<-read.table("monocle2.gene_noCCL5.txt",header = T)
Monocle2<-na.omit(Monocle2)
destiny<-read.table("destiny_genes_noCCL5.txt",header = T)
destiny<-na.omit(destiny)

#macro<-read.table("marcophage genes.txt",header = T)

index<-intersect(monocle3$gene_short_name,slingshot$x)
Genes<-union(monocle3$gene_short_name,slingshot$x)
Genes<-union(Genes,Monocle2$x)%>%unique()
Genes<-union(Genes,velo$Gene)%>%unique()
Genes<-union(Genes,destiny$x)%>%unique()

load("D:/Bioinfrolf/data/tmpRdata/05.Mono_Macr-tmp.Rdata")
matrix=as.matrix(Mono_Macr@assays$RNA@data)
matrix[1:5,1:5]
g<-rownames(matrix)[match(Genes,rownames(matrix))]%>%na.omit()
matrix<-matrix[g,]
matrix<-t(matrix)
data<-as.data.frame(matrix)
results<-Mono_Macr@meta.data$tissue
data$tissue<-results
data$tissue<-as.factor(data$tissue)
rm(matrix,Mono_Macr)

library(Boruta)
set.seed(100)
feature.selection <- Boruta(tissue ~ ., data = data, doTrace = 1)
feature.selection$timeTaken
save(feature.selection,file = "Boruta-sub.Rdata")
table(feature.selection$finalDecision)
write.table(feature.selection$finalDecision,
            file = 'FinalDecision-sub.txt',sep="\t",quote=F,col.names=F)
fNames <- getSelectedAttributes(feature.selection) #withTentative = TRUE
fNames<-gsub("`","",fNames)
fNamesAd<-names(feature.selection$finalDecision[feature.selection$finalDecision!="Rejected"])
fNamesAd<-gsub("`","",fNamesAd)
table(feature.selection$finalDecision)
#orignal<-data
data<-data[,Genes]
data<-data[,fNamesAd]
set.seed(201)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1, ] #the training data set
test <- data[ind==2, ] #the test data set
table(train$tissue)
table(test$tissue)
fc<-ncol(data)

x <- as.matrix(train[, 1:(ncol(data)-1)])
y <- train[, ncol(data)]
rm(matrix,Mono_Macr,monocle3,feature.selection,Genes,slingshot)
#############xgbt#######################
cntrl = trainControl(
  method = "cv",
  number = 10,
  verboseIter = TRUE,
  returnData = FALSE,
  returnResamp = "final"                                                        
)
grid = expand.grid(
  nrounds = c(100),
  colsample_bytree = 1,
  min_child_weight = 1,
  eta = c(0.1,0.2), #0.3 is default,
  gamma = c(0.5),
  subsample = 0.5,
  max_depth = c(7)
)
grid
fc<-ncol(train)
library(future)
future::plan("multicore", workers=14)
plan()

set.seed(100)
train.xgb = train(
  x = train[, 1:(fc-1)],
  y = train[, fc],
  trControl = cntrl,
  tuneGrid = grid,
  method = "xgbTree"
)

train.xgb

param <- list(  objective           = 'multi:softprob', 
                booster             = "gbtree",
                eval_metric         = "error",
                eta                 = 0.2, 
                max_depth           = 7, 
                subsample           = 0.5,
                colsample_bytree    = 1,
                gamma               = 0.5,
                min_child_weight = 1,
                num_class=3,
                lambda=0.01
)

x <- as.matrix(train[, 1:(fc-1)])

train.mat <- xgb.DMatrix(data = x, 
                         label = y)
set.seed(300)
# future::plan("multisession", workers = 30) 
plan()
xgb.fit1 <- xgb.train(params = param, data = train.mat, nrounds = 10000)
xgb.fit50<-xgb.fit1

library(InformationValue)
pred <- predict(xgb.fit1, x)
matp_train<-matrix(data=pred, nrow = nrow(train), ncol = 3, byrow = T, dimnames = NULL)
rownames(matp_train)<-rownames(train)
colnames(matp_train)<-c("Normal","TumorEdge","TumorCore")
head(matp_train)
class_train<-apply(matp_train,1,function(x){
  which.max(x)-1
})
head(class_train)
table(y,class_train)

data.testMat <- as.matrix(test[, 1:(fc-1)])
xgb.data.test <- predict(xgb.fit1, data.testMat)
y.test <- as.factor(test[, fc])
matp_test<-matrix(data=xgb.data.test, nrow = nrow(test), ncol = 3, byrow = T, dimnames = NULL)
rownames(matp_test)<-rownames(test)
colnames(matp_test)<-c("Normal","TumorEdge","TumorCore")
head(matp_test)
class_test<-apply(matp_test,1,function(x){
  which.max(x)-1
})

table(y.test,class_test)

impMatrix <- xgb.importance(feature_names = dimnames(x)[[2]], model = xgb.fit1)
impMatrix 

which(impMatrix$Feature%in%c("SPP1","MIF","CD74"))
save(xgb.fit,data,train,test,param,file = "xgboost_noCCL5.Rdata")
write.table(impMatrix,file = "impMatrix.txt",quote = F)

#############catboost########################
devtools::install_github('catboost/catboost', subdir = 'catboost/R-package')
library(caret)
library(catboost)

set.seed(300)

set.seed(3)

#fold_len_multiplier=1.1,
# od_type='Iter',
# od_wait=120,

fit_control <- trainControl(method = "cv",
                            number = 10,
                            classProbs = TRUE)

grid <- expand.grid(depth = c(7,8,9),
                    learning_rate = c(0.2,0.1,0.05),
                    iterations = 1500,
                    l2_leaf_reg = c(0.001,0.01,0.1,1),
                    rsm = 0.95,
                    border_count = 64)
grid <- expand.grid(                  #loss_function='MultiClass',
  #logging_level='Verbose',
  #classes_count=3,
  depth = c(7,8,9),
  learning_rate = c(0.2,0.1),
  iterations = 1500,
  l2_leaf_reg = c(0.01,0.05),
  rsm = 0.95,
  # random_seed=2021,
  # metric_period=50,
  border_count = 64)

future::plan("multiprocess", workers = (detectCores() - 2))
report <- caret::train(train[, 1:(ncol(data)-1)], as.factor(make.names(train[, ncol(data)])),
                       method = catboost.caret,
                       #verbose = TRUE, preProc = NULL,
                       tuneGrid = grid, trControl = fit_control)
importance <- varImp(report, scale = FALSE)
print(importance)

param <- list(  loss_function='MultiClass',
                logging_level='Verbose',
                classes_count=3,
                learning_rate = 0.2, 
                depth           = 7, 
                iterations = 10000,
                border_count = 64,
                rsm = 0.95,
                l2_leaf_reg = 0.05,
                metric_period=50
)


features<-as.data.frame(train[, 1:(ncol(data)-1)])
labels <- train$tissue
train_pool <- catboost.load_pool(data = features, label = labels)
catboost <- catboost.train(train_pool,  NULL,params =param)

real_data <- as.data.frame(test[, 1:(ncol(data)-1)])
colnames(real_data)<-colnames(features)#gsub("\\.","-",colnames(real_data))
labelsy <- test$tissue
real_pool <- catboost.load_pool(data=real_data,label=labelsy)


prediction <- catboost.predict(catboost, real_pool,prediction_type = "Probability")

results<-apply(prediction, 1,function(i){
  which.max(i)-1
})

table(labelsy,results)

im<-catboost.get_feature_importance(catboost, 
                                    pool = NULL, 
                                    type = 'FeatureImportance',
                                    thread_count = -1)
im<-as.data.frame(im)
im$gene<-rownames(im)
im<-im[order(im$V1,decreasing = T),]
im<-im[im$V1!="0",]
head(im)

#catboost::catboost.save_model(catboost$finalModel, "catboost-1")
save(list=ls(),file = "catboost-1.Rdata")
