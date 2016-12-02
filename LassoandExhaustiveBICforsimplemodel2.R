setwd("/Users/pihu_yadav/Desktop/Mice_Gene_Expression_Problem-master")
full_data<-read.csv("Kidney_2.csv")
n<-length(full_data$Mapk1)


#Lasso with Cross Validation

library(leaps)

full_data<-subset(full_data, select = -c(Gene.Name) )


x_data<-subset(full_data, select = -c(Mapk1,Gene.Name) )

y_data<-full_data$"Mapk1"

best_subset<-leaps(x_data,y_data,nbest=1,method="adjr2",names=colnames(x_data))

library(glmnet)
x_data<-subset(full_data, select = -c(Mapk1,Gene.Name) )

y_data<-full_data$"Mapk1"
X <- as.matrix(x_data)
Y<-y_data
cvfit <- glmnet::cv.glmnet(X, Y, nfolds=40)
coef(cvfit, s = "lambda.min")
Y_pred<-predict(cvfit, X, s = "lambda.min")
mean((Y_pred-Y)^2)


#Best Subset Selection on entire Dataset

install.packages("leaps")

library(leaps)



regsubsets.out <-

    regsubsets(Mapk1~ Cdc42+Pla2g6+Akt2+Plcg2+Rac2+Rik+Mapkapk2 +Pik3cd+Pla2g5+Sphk2+Map2k1+Pik3r3+Ptk2+Nras+Nos3+Pik3r1+Pik3ca+Ppp3cb+Map2k2+Nfatc4+Mapk13+Rac1+Nfat5       

,

               data = full_data,

               nbest = 1,       # 1 best model for each number of predictors

               nvmax = NULL,    # NULL for no limit on number of variables

               force.in = NULL, force.out = NULL,

               method = "exhaustive")

summary(regsubsets.out)
summary.regsubsets(regsubsets.out)
plot(regsubsets.out)
coef(regsubsets.out,4)

#Calculating Correlations
x_data<-subset(full_data, select = -c(Gene.Name) )
cormat=cor(x_data,x_data)
num0.5=sum(cormat>0.5)
m=-0.5
numless0.5=sum(cormat<m)
# genes that are negatively correlated with a negative value of less than -0.5 are Pik3cd & Rac1 and Plcg2 & Nos3
# genes that are positively correlated with a value of greater than 0.5
# Map2k1 and  Pik3r3 and Nfat5 with Rik
# Rac1 with Pla2g6
# Map2k1 with Pik3r3
# Nfat5 with Pik3r3
# Rak1 with Ptk2

#Best Subset Selection with test set and training set separate
arr=rep(0,40)
for (k in 1:40)
{
set.seed(k)
train=sample(c(1:40),30, replace=FALSE)
train=sort(train)
train=array(train)
#you may also use train_index=sample(c(1:263),131) to create training data
test=full_data[-train,]
regfit.best=regsubsets(Mapk1~ Cdc42+Pla2g6+Akt2+Plcg2+Rac2+Rik+Mapkapk2 +Pik3cd+Pla2g5+Sphk2+Map2k1+Pik3r3+Ptk2+Nras+Nos3+Pik3r1+Pik3ca+Ppp3cb+Map2k2+Nfatc4+Mapk13+Rac1+Nfat5       

,

               data = full_data[train,], nbest = 1,       # 1 best model for each number of predictors

               nvmax = NULL,    # NULL for no limit on number of variables

               force.in = NULL, force.out = NULL,

               method = "exhaustive")
plot(regfit.best)

test.mat=model.matrix(Mapk1~ Cdc42+Pla2g6+Akt2+Plcg2+Rac2+Rik+Mapkapk2 +Pik3cd+Pla2g5+Sphk2+Map2k1+Pik3r3+Ptk2+Nras+Nos3+Pik3r1+Pik3ca+Ppp3cb+Map2k2+Nfatc4+Mapk13+Rac1+Nfat5 ,data=test) # create an X matrix of test data
val.errors=rep(NA,23)
for(i in 1:23){
   coefi=coef(regfit.best,id=i)
   pred=test.mat[,names(coefi)]%*%coefi
   val.errors[i]=mean((test$Mapk1-pred)^2)
}
val.errors
arr[k]=which.min(val.errors)

}
arr

