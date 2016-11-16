###File for the data analysis and model evaluation##

#Predefine some commonly used variables
n<-length(full_data$Mapk1)

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
loocv_fit<-cv.glm(full_data,glm_model,K=n) #Variable lengths differing for column Cdc42 
#The resulting delta[1] is the prediction error
loocv_fit$delta[1] 
#Output is 0.00842281 which is much lower than the forward-backward stepwise 
#selection using the AIC

#Using the glmnet procedure to not have to pre-specify a glm model
#implement the k-fold using LOOCV since there is only 40 points
#alpha=0 corresponds to quadratic regularizer
library(glmnet)
cv.out<-cv.glmnet(data.matrix(x_data),y_data,alpha=0,nfold=n) 
bestlambda<-cv.out$lambda.min
#Output 0.1362245

#Refit the data using this best lambda
ridge.model<-glmnet(data.matrix(x_data),y_data,alpha=0,lambda=bestlambda)

#Calculate training error and use sum of squares of differences
y_new<-predict(ridge.model,s=bestlambda,newx=data.matrix(x_data))
training_error<-sum((y_new-y_data)^2)
training_error
#Output 0.2185162

#alpha=1 corresponds to the lasso regularizer
#Lasso has a higher chance of computing a sparse model (less predictors)
cv.out<-cv.glmnet(data.matrix(x_data),y_data,alpha=1,nfold=40) 
bestlambda<-cv.out$lambda.min
bestlambda
#Output is 0.02493913

#Refit the data using this best lambda (with lasso regularization)
ridge.model<-glmnet(data.matrix(x_data),y_data,alpha=1,lambda=bestlambda)

#Calculate training error and use sum of squares of differences
y_new<-predict(ridge.model,s=bestlambda,newx=data.matrix(x_data))
training_error<-sum((y_new-y_data)^2)
training_error
#Output 0.3166494
#Thus we get a higher training error using a lasso regularizer than when 
#we use a quadratic regularizer.  

#Can use cross validation and BIC to figure out optimal model sizes
#Implement LOOCV with steps by hand
test_error<-0
#Pull all data except for "Expression.Name" 
keep_columns<-colnames(full_data)[-1]


#Using LOOCV since there are only 40 data points
for (i in 1:n){
#Define the training and test sets
x_train<-full_data[-i,keep_columns]
x_test<-full_data[i,keep_columns]
y_train<-full_data[-i,"Mapk1"]
y_test<-full_data[i,"Mapk1"]

#Initial model is full
initial_model<-lm(Mapk1~.,data=x_train)
fitted_model<-step(initial_model,direction="both",k=log(n))

#Make a prediction based on the fitted model then calculate the test error
y_predict<-predict(fitted_model,x_test)
test_error<-test_error+(y_predict-y_test)^2
}
test_error
#Output = 1.008477 

#Maybe can implement KNN--but not enough data (discuss!)
