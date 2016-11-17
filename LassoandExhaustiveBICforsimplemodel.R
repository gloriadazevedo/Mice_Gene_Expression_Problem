setwd("/Users/pihu_yadav/Desktop/Mice_Gene_Expression_Problem-master")
full_data<-read.csv("Kidney_2.csv")

library(glmnet)
x_data<-subset(full_data, select = -c(Mapk1,Gene.Name) )

y_data<-full_data$"Mapk1"
X <- as.matrix(x_data)
Y<-y_data
cvfit <- glmnet::cv.glmnet(X, Y)
coef(cvfit, s = "lambda.1se")
Y_pred<-predict(cvfit, X)
mean((Y_pred-Y)^2)


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
plot(regsubsets.out)
coef(regsubsets.out,4)
