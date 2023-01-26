######################################################
#install.packages("lightgbm", repos = "https://cran.r-project.org")
#install.packages("Ckmeans.1d.dp")
library(lightgbm)
library(tidyr)
library(tibble)
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(SeuratObject)
library(car) #package to calculate Variance Inflation Factor
library(corrplot) #correlation plots
library(glmnet) #allows ridge regression, LASSO and elastic net
library(reshape2)
library(pROC)
library(rpart) #tissueification and regression trees
library(partykit) #treeplots
library(randomForest) #random forests
library(xgboost) #gradient boosting 
library(caret) #tune hyper-parameters
library(Ckmeans.1d.dp)
rm(list = ls())
setwd("D:/Bioinfrolf/HCC-SC/Myeloid/model")
load("D:/Bioinfrolf/Bioinfrolf/data/HAR-MYE/05.Mono_Macr-tmp.Rdata")
#load("~/Bioinfrolf/Rdata/sceforMatrix.Rdata")
#g1<-rownames(Mono_Macr)
Genes<-read.table("Genes1.txt",header = T)
Genes<-Genes$x

matrix=as.matrix(Mono_Macr@assays$RNA@data)
matrix[1:5,1:5]
g<-rownames(matrix)[match(Genes,rownames(matrix))]%>%na.omit()
#gd<-setdiff(Genes,g)

#write.table(gd,file="gd.txt",quote = F,row.names = F,col.names = F)
matrix<-matrix[g,]

dim(matrix)
results<-Mono_Macr@meta.data$tissue_sub
matrix<-t(matrix)

data<-as.data.frame(matrix)
rownames(data)[1:5]
head(rownames(Mono_Macr@meta.data))
data$tissue<-results
data$tissue<-as.factor(data$tissue)
rm(Mono_Macr)
table(data$tissue)

load("Boruta-sub.Rdata")
library(Boruta)
fNames <- getSelectedAttributes(feature.selection) #withTentative = TRUE
fNames<-gsub("`","",fNames)
fNamesAd<-names(feature.selection$finalDecision[feature.selection$finalDecision!="Rejected"])
fNamesAd<-gsub("`","",fNamesAd)
orignal<-data
data<-data[,fNamesAd]

data$tissue<-results
data$tissue<-as.factor(data$tissue)
rm(matrix,Mono_Macr)
table(data$tissue)

#y <- ifelse(data$tissue == "Tumor", 1, 0)

data$tissue <- ifelse(data$tissue == "Normal", 0, ifelse(data$tissue == "TumorEdge",1,2))
#colnames(data)<-make.names(colnames(data))

set.seed(123) #random number generator
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1, ] #the training data set
test <- data[ind==2, ] #the test data set
str(test) #confirm it worked
table(train$tissue)
table(test$tissue)

ncol(data)
x <- as.matrix(train[, 1:(ncol(data)-1)])
y <- train[, ncol(data)]
################################################
library(lightgbm)
write.csv(as.matrix(train[, 1:(ncol(train)-1)]),file="train.csv")
write.csv(train$tissue,file="train_label.csv")
write.csv(as.matrix(test[, 1:(ncol(test)-1)]),file="test.csv")
write.csv(test$tissue,file="test_label.csv")

# dtrain <- list(data=Matrix(as.matrix(train[, 1:(ncol(train)-1)]),sparse = T),
#               label=train$tissue)
# dtest <- list(data=Matrix(as.matrix(test[, 1:(ncol(test)-1)]),sparse = T),
#               label=test$tissue)
dtrain1 <- lgb.Dataset(dtrain$data, label = train$label)
dtest1 <- lgb.Dataset(dtest$data, label = test$label)
valids <- list(test=dtest1)

param=list(
  boosting_type='gbdt',  
  objective = "multiclass",  
  num_class=3,  
  metric='multi_logloss',  
  num_leaves= 120,  
  min_data_in_leaf=100,  
  learning_rate=0.06,  
  feature_fraction=0.8,  
  bagging_fraction=0.8,  
  bagging_freq=5,  
  lambda_l1= 0.4,  
  lambda_l2= 0.5,  
  min_gain_to_split=0.2,  
  verbose =-1
)

param=list(
  boosting_type='gbdt',
  objective = "multiclass",
  num_class=3,
  #min_data_in_leaf=2,
  feature_pre_filter=T,
  #boost_from_average=T
  #, metric = "l2"
  #, nfold=5
  #, min_data=1
  #, learning_rate=1
  #, early_stopping_rounds=10
)

lgb1 <-lgb.train(
  params = param
  , data = dtrain1
  ,nrounds = 5L
  , valids = valids
)

im<-lgb.importance(lgb1)
pre.lgb=predict(lgb1,bia7)
