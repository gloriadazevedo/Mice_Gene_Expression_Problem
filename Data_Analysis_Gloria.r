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
mis_classified_vector<-rep(1:15)

#Columns to include (everything except "Gene.Name","Mapk1")
include_col<-colnames(full_data)
include_col<-include_col[-1]
include_col<-include_col[-5]

for (i in 1:length(mis_classified_vector)){
	model_fit<-knn.cv(full_data[include_col],full_data$"Mapk1",i)
	#Use this if the output is 0 or 1
	#mis_classified<-sum(abs(as.numeric(model_fit)-full_data$"Mapk1"))
	#Squared error loss
	mis_classified<-sum((as.numeric(model_fit)-full_data$"Mapk1")^2)
	mis_classified_vector[i]<-mis_classified
}

############################
#Post discussion and figuring out which predictor sets are best#
#Want to calculate coefficients for these linear model  
#Evaluate the fit via R^2 and Adjusted R^2
model_1<-lm(Mapk1~Akt2+Rik+Pik3r3+ Rac1,data=full_data)

summary(model_1)
# Call:
# lm(formula = Mapk1 ~ Akt2 + Rik + Pik3r3 + Rac1, data = full_data)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -0.16348 -0.06566  0.00381  0.05495  0.17092 

# Coefficients:
            # Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.44461    0.15706  -2.831 0.007643 ** 
# Akt2        -0.40548    0.16922  -2.396 0.022047 *  
# Rik          0.22072    0.10178   2.169 0.036994 *  
# Pik3r3       0.24399    0.11350   2.150 0.038565 *  
# Rac1         0.31716    0.08464   3.747 0.000644 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.08281 on 35 degrees of freedom
# Multiple R-squared:  0.6101,    Adjusted R-squared:  0.5656 
# F-statistic: 13.69 on 4 and 35 DF,  p-value: 8.093e-07

model_2<-lm(Mapk1~Rik+Pik3r3+Rac1+Nfat5,data=full_data)

summary(model_2)

# Call:
# lm(formula = Mapk1 ~ Rik + Pik3r3 + Rac1 + Nfat5, data = full_data)

# Residuals:
      # Min        1Q    Median        3Q       Max 
# -0.217531 -0.060483  0.007193  0.051996  0.164377 

# Coefficients:
            # Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  -0.3109     0.1979  -1.571  0.12515   
# Rik           0.1741     0.1143   1.523  0.13665   
# Pik3r3        0.2331     0.1304   1.788  0.08249 . 
# Rac1          0.2570     0.0909   2.827  0.00772 **
# Nfat5         0.0864     0.1397   0.619  0.54019   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.08886 on 35 degrees of freedom
# Multiple R-squared:  0.5511,    Adjusted R-squared:  0.4998 
# F-statistic: 10.74 on 4 and 35 DF,  p-value: 8.703e-06

model_3<-lm(Mapk1~Pik3r3+Rac1,data=full_data)

summary(model_3)

# Call:
# lm(formula = Mapk1 ~ Pik3r3 + Rac1, data = full_data)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -0.19118 -0.05534 -0.00354  0.06817  0.20241 

# Coefficients:
            # Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.05643    0.10239  -0.551 0.584876    
# Pik3r3       0.40410    0.09469   4.268 0.000132 ***
# Rac1         0.27832    0.09016   3.087 0.003821 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.0908 on 37 degrees of freedom
# Multiple R-squared:  0.5045,    Adjusted R-squared:  0.4777 
# F-statistic: 18.83 on 2 and 37 DF,  p-value: 2.284e-06
