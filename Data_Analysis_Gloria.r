###File for the data analysis and model evaluation##

#Naive way is to implement best subset selection since there are only 24 variables
library(leaps)

#The predictors are everything except for the field Mapk1
#Also can't have categorical data in there
x_data<-subset(full_data, select = -c(Mapk1,Gene.Name) )
y_data<-full_data$"Mapk1"
best_subset<-leaps(x_data,y_data,nbest=1,method="adjr2",names=colnames(x_data))
#Thus we get the best model for each of the sample sizes--then can use linear regression
#to get the coefficients for it.
#Use adjusted R^2 as the criterion since that takes into account model size

#However, in the future, this may not work if we have more than say 40 gene expressions
#so we also implement the forward-backward selection method for the model
library(MASS)
lm_fit<-lm(y_data~.,data=x_data)
step_lm<-stepAIC(lm_fit,direction="both")
#Resulting fit using AIC but still has 12 variables: 
# y_data ~ Cdc42 + Pla2g6 + Akt2 + Plcg2 + Rac2 + Sphk2 + Map2k1 + 
   # Ptk2 + Nos3 + Pik3ca + Ppp3cb + Rac1

#Fit this model to get the coefficients
lm_fit_AIC_result<-lm(y_data ~ Cdc42 + Pla2g6 + Akt2 + Plcg2 + Rac2 + Sphk2 + Map2k1 + 
   Ptk2 + Nos3 + Pik3ca + Ppp3cb + Rac1,data=x_data)
   
#Estimate training error for the model using quadratic loss
y_predict<-predict(lm_fit_AIC_result,x_data)
training_error<-sum((y_predict-y_data)^2)
training_error #Output = 0.1521644

#However we want to minimize the test error for new gene expressions
#which is not the same objective as minimizing the training error
#Since our data set is so small, we opt to use the Leave-One-Out-Cross-Validation 
#for training and testing the data
#LOOCV is a special case of cross validation where the sample size is 1 but K is the number
#of groups so K is actually n=40
library(boot)
#First need to derive a glm model
glm_model<-glm(Mapk1 ~ Cdc42 + Pla2g6 + Akt2 + Plcg2 + Rac2 + Sphk2 + Map2k1 + 
   Ptk2 + Nos3 + Pik3ca + Ppp3cb + Rac1,data=full_data)
loocv_fit<-cv.glm(full_data,glm_model,K=40) #Variable lengths differing for column Cdc42 
#The resulting delta[1] is the prediction error
loocv_fit$delta[1] 
#Output is 0.00842281 which is much lower than the forward-backward stepwise 
#selection using the AIC
