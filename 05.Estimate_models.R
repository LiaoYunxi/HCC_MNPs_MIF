library(doParallel)
all_cores <- parallel::detectCores(logical = FALSE)
registerDoParallel(cores = all_cores)

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
library(leaps) #best subsets regression
library(glmnet) #allows ridge regression, LASSO and elastic net
library(reshape2)

library(rpart) #tissueification and regression trees
library(partykit) #treeplots

library(caret) #tune hyper-parameters
library(Ckmeans.1d.dp)

library(xgboost)
library(InformationValue)
library(pROC)

data.testMat <- as.matrix(test[, 1:(fc-1)])
probtest<-predict(xgb.fit,data.testMat)
probtest <- matrix(probtest, ncol=num_class, byrow=TRUE)
t = table(class_test,y.test)
acc = sum(diag(t))/nrow(test) *100
print(paste("模型准确率为：",round(acc,4),'%',sep=''))

rocmat<-matrix(data=0, nrow = nrow(probtest), ncol = 3, 
               byrow = FALSE, dimnames = NULL)
rownames(rocmat)<-rownames(test)
colnames(rocmat)<-c("Normal","TumorEdge","TumorCore")

for(i in 1:length(y.test)){
  if(y.test[i]==0){
    rocmat[i,1]<-1
  }else if(y.test[i]==1){
    rocmat[i,2]<-1
  }else(rocmat[i,3]<-1)
}

rownames(probtest)<-rownames(test)
colnames(probtest)<-c("Normal","TumorEdge","TumorCore")


rocmatt<-as.character(t(rocmat))
mattp<-as.character(t(probtest))%>%as.numeric()
roc_micro<-roc(rocmatt,mattp)
optimalCutoff(rocmatt,mattp)
misClassError(rocmatt,mattp)
confusionMatrix(rocmatt,mattp, threshold = 0.5399555)
precision(rocmatt,mattp, threshold = 0.999)

library(modEvA)
f1_fun = function(pre,y){
  class = sort(unique(y))
  tp=NA
  fp=NA
  fn=NA
  for(i in 1:length(class)){
    tp[i] = sum(pre==class[i] & y==class[i])
    fp[i] = sum(pre==class[i] & y!=class[i])
    fn[i] = sum(pre!=class[i] & y==class[i])
  }
  f1 = 2*tp/(2*tp+fp+fn)
  names(f1) = class
  print(table(pre,y))
  print('-------------f1--------------------')
  print(f1)
  print('--------------mean(f1)-------------------')
  print(mean(f1))
  print(tp)
  print(fp)
  print(fn)
}
f1_fun(class_test,y.test)
R=(346+965+400)/((346+965+400)+(66+151+68))
P=(346+965+400)/((346+965+400)+(41+114+130))
F1=R*P*2/(R+P)

aupr=AUC(obs=rocmatt,pred=mattp,curve = "PR", simplif=TRUE, main = "PR curve",
         curve.col="#E2D200",curve.lwd=5,interval=0.0001)

prob<-c()
for(i in 1:nrow(rocmat)){
  a<-which(rocmat[i,]==1)
  prob[i]=probtest[i,a]
}
head(prob)

prob1<-probtest[,1]
prob2<-probtest[,2]
prob3<-probtest[,3]
rocmat1<-rocmat[,1]
rocmat2<-rocmat[,2]
rocmat3<-rocmat[,3]
roc1<-roc(rocmat1,prob1)
roc2<-roc(rocmat2,prob2)
roc3<-roc(rocmat3,prob3)

optimalCutoff(rocmat1,prob1)
misClassError(rocmat1,prob1)
confusionMatrix(rocmat1,prob1, threshold = 0.5194881)
precision(rocmat1,prob1, threshold = 0.9)
aupr1=AUC(obs=rocmat1,pred=prob1,curve = "PR", simplif=TRUE, main = "PR curve",
          curve.col="#00A08A",curve.lwd=5,interval=0.0001)
optimalCutoff(rocmat2,prob2)
misClassError(rocmat2,prob2)
confusionMatrix(rocmat2,prob2, threshold = 0.6698815)
precision(rocmat2,prob2, threshold = 0.6698815)
aupr1=AUC(obs=rocmat2,pred=prob2,curve = "PR", simplif=TRUE, main = "PR curve",
          curve.col="#F98400",curve.lwd=5,interval=0.0001)
optimalCutoff(rocmat3,prob3)
misClassError(rocmat3,prob3)
confusionMatrix(rocmat3,prob3, threshold = 0.3285257)
precision(rocmat3,prob3, threshold =0.3285257)
aupr1=AUC(obs=rocmat3,pred=prob3,curve = "PR", simplif=TRUE, main = "PR curve",
          curve.col="#FF0000",curve.lwd=5,interval=0.0001)

xgboost_roc <- multiclass.roc(y.test, as.numeric(prob))
#绘制ROC曲线和AUC值
plot.roc(xgboost_roc$rocs[[1]], col="#5BBCD6",add=F)
plot.roc(xgboost_roc$rocs[[2]],col="#B40F20",add=T)
plot.roc(xgboost_roc$rocs[[3]],col="#E2D200",add=T)
auc(xgboost_roc)                         

roc_macro<-data.frame(row.names = rownames(data.testMat),sensitivities=rep(0,nrow(data.testMat)),
                      specificities=rep(0,nrow(data.testMat)),CancerType="macro average")
for(i in 1:nrow(test)){
  roc_macro$sensitivities[i]<-(roc1$sensitivities[i]+roc2$sensitivities[i]+roc3$sensitivities[i])/num_class
  roc_macro$specificities[i]<-(roc1$specificities[i]+roc2$specificities[i]+roc3$specificities[i])/num_class
}
which(is.na(roc_macro$sensitivities))
auc_macro<-sum(as.numeric(roc_macro$sensitivities),na.rm = T)/1996

roc_micro<-data.frame(sensitivities=roc_micro$sensitivities,
                      specificities=roc_micro$specificities,CancerType="micro average")

roc_normal<-data.frame(sensitivities=roc1$sensitivities,
                       specificities=roc1$specificities,CancerType="adjacent normal")
roc_edge<-data.frame(sensitivities=roc2$sensitivities,
                     specificities=roc2$specificities,CancerType="peripheral tumor")
roc_core<-data.frame(sensitivities=roc3$sensitivities,
                     specificities=roc3$specificities,CancerType="core tumor")


col5<-c("macro average"="#5BBCD6","micro average"="#E2D200","adjacent normal"="#00A08A",
        "peripheral tumor"="#F98400","core tumor"="#FF0000")
df<-rbind(roc_macro,roc_micro,roc_normal,roc_edge,roc_core)
df$CancerType<-factor(df$CancerType,levels = c("adjacent normal","peripheral tumor","core tumor",
                                               "macro average","micro average"))
unique(df$CancerType)

jpeg(file="multiROCtest_201_data.jpeg",width =7,height = 5,units = "in", res = 2000)
ggplot(data=df, aes(1-specificities,sensitivities,colour =CancerType))+
  geom_line(size=2,alpha=0.7,key_glyph="smooth")+
  scale_color_manual(values = col5)+
  annotate(geom = "segment", x = 0, y = 0, xend = 1, yend = 1,
           colour="grey",linetype="dashed",size=2)+
  labs(title="ROC Curve", x="1-Specificity (FPR)", y="Sensitivity (TPR)") +
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black'),
        axis.text.x = element_text(hjust = 1,vjust = 1))+
  coord_cartesian(xlim=c(0,1), ylim = c(0,1))#+NoLegend()
dev.off()

ggplot()+annotate(geom = "segment", x = 0, y = 0, xend = 1, yend = 1,
                  colour="grey",linetype="dashed",size=2)+
  geom_line(data=roc_normal, aes(1-specificities,sensitivities),colour ="#00A08A",size=2)+
  geom_line(data=roc_edge, aes(1-specificities,sensitivities),colour ="#F98400",size=2)+
  geom_line(data=roc_core, aes(1-specificities,sensitivities),colour ="#FF0000",size=2)+
  geom_line(data=roc_micro, aes(1-specificities,sensitivities),colour ="#E2D200",size=2) +
  geom_line(data=roc_macro, aes(1-specificities,sensitivities),colour ="#5BBCD6",size=2) +
  labs(title="ROC Curve", x="1-Specificity (FPR)", y="Sensitivity (TPR)") +
  theme(panel.grid = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1), 
        panel.background = element_blank(), 
        plot.title = element_text(size = 15, hjust = 0.5),
        axis.text = element_text(size = 15, color = 'black'), 
        axis.title = element_text(size = 15, color = 'black'),
        axis.text.x = element_text(hjust = 1,vjust = 1))+
  coord_cartesian(xlim=c(0,1), ylim = c(0,1))
