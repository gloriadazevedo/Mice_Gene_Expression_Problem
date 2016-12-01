#In order to get a simple model, two approaches have been applied: Using a lasso regularizer and using BIC (which penalizes having a large number of parameters) for exhaustive subset selection, to find the most accurate and simple model
setwd("/Users/pihu_yadav/Desktop/Mice_Gene_Expression_Problem-master")
full_data<-read.csv("Kidney_2.csv")

#Performing Lasso Regression with 10 fold cross-validation
library(glmnet)
x_data<-subset(full_data, select = -c(Mapk1,Gene.Name) )
y_data<-full_data$"Mapk1"
X <- as.matrix(x_data)
Y<-y_data
cvfit <- glmnet::cv.glmnet(X, Y)
coef(cvfit, s = "lambda.1se")
# 24 x 1 sparse Matrix of class "dgCMatrix"
                      # 1
# (Intercept) -0.04389607
# Cdc42        .         
# Pla2g6       .         
# Akt2         .         
# Plcg2        .         
# Rac2         .         
# Rik          0.06664806
# Mapkapk2     .         
# Pik3cd       .         
# Pla2g5       .         
# Sphk2        .         
# Map2k1       .         
# Pik3r3       0.15013675
# Ptk2         .         
# Nras         .         
# Nos3         .         
# Pik3r1       .         
# Pik3ca       .         
# Ppp3cb       .         
# Map2k2       .         
# Nfatc4       .         
# Mapk13       .         
# Rac1         0.08441959
# Nfat5        0.01435892  

#Need to use the lambda.min which corresponds that gives minimum cvm ( mean cross-validated error ).
coef(cvfit, s = "lambda.min")
# 24 x 1 sparse Matrix of class "dgCMatrix"
                  # 1
# (Intercept) -0.1454
# Cdc42        .     
# Pla2g6       .     
# Akt2        -0.0006
# Plcg2        .     
# Rac2         .     
# Rik          0.1115
# Mapkapk2     .     
# Pik3cd      -0.0376
# Pla2g5       .     
# Sphk2        .     
# Map2k1       .     
# Pik3r3       0.1652
# Ptk2         .     
# Nras         .     
# Nos3         .     
# Pik3r1       .     
# Pik3ca       .     
# Ppp3cb       .     
# Map2k2       .     
# Nfatc4       .     
# Mapk13       .     
# Rac1         0.1388
# Nfat5        0.0569
min_lambda_model<-lm(Mapk1 ~ Akt2+Rik+Pik3cd+Pik3r3+Rac1+Nfat5,data=full_data)
summary(min_lambda_model)

cvfit$lambda.1se
[1] 0.03971015
#Cross-validation error is 0.013 as given by cvfit.cvm
#From this value we see that the lasso actually performs a good prediction, as the error is small.

#Using exhaustive best subset selection with BIC as the selection criteria
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
#Subset selection object
#Call: regsubsets.formula(Mapk1 ~ Cdc42 + Pla2g6 + Akt2 + Plcg2 + Rac2 + 
#    Rik + Mapkapk2 + Pik3cd + Pla2g5 + Sphk2 + Map2k1 + Pik3r3 + 
#    Ptk2 + Nras + Nos3 + Pik3r1 + Pik3ca + Ppp3cb + Map2k2 + 
#    Nfatc4 + Mapk13 + Rac1 + Nfat5, data = full_data, nbest = 1, 
#    nvmax = NULL, force.in = NULL, force.out = NULL, method = "exhaustive")
#23 Variables  (and intercept)
#         Forced in Forced out
#Cdc42        FALSE      FALSE
#Pla2g6       FALSE      FALSE
#Akt2         FALSE      FALSE
#Plcg2        FALSE      FALSE
#Rac2         FALSE      FALSE
#Rik          FALSE      FALSE
#Mapkapk2     FALSE      FALSE
#Pik3cd       FALSE      FALSE
#Pla2g5       FALSE      FALSE
#Sphk2        FALSE      FALSE
#Map2k1       FALSE      FALSE
#Pik3r3       FALSE      FALSE
#Ptk2         FALSE      FALSE
#Nras         FALSE      FALSE
#Nos3         FALSE      FALSE
#Pik3r1       FALSE      FALSE
#Pik3ca       FALSE      FALSE
#Ppp3cb       FALSE      FALSE
#Map2k2       FALSE      FALSE
#Nfatc4       FALSE      FALSE
#Mapk13       FALSE      FALSE
#Rac1         FALSE      FALSE
#Nfat5        FALSE      FALSE
#1 subsets of each size up to 23
#Selection Algorithm: exhaustive
#          Cdc42 Pla2g6 Akt2 Plcg2 Rac2 Rik Mapkapk2 Pik3cd Pla2g5 Sphk2 Map2k1
#1  ( 1 )  " "   " "    " "  " "   " "  " " " "      " "    " "    " "   " "   
#2  ( 1 )  " "   " "    " "  " "   " "  " " " "      " "    " "    " "   " "   
#3  ( 1 )  " "   " "    "*"  " "   " "  "*" " "      " "    " "    " "   " "   
#4  ( 1 )  " "   " "    "*"  " "   " "  "*" " "      " "    " "    " "   " "   
#5  ( 1 )  " "   " "    "*"  " "   " "  "*" " "      " "    " "    " "   " "   
#6  ( 1 )  " "   " "    "*"  " "   " "  "*" " "      " "    " "    "*"   " "   
#7  ( 1 )  " "   " "    "*"  "*"   "*"  "*" " "      " "    " "    "*"   " "   
#8  ( 1 )  " "   " "    "*"  "*"   "*"  "*" " "      " "    " "    "*"   " "   
#9  ( 1 )  "*"   " "    "*"  "*"   "*"  " " " "      " "    " "    "*"   " "   
#10  ( 1 ) "*"   " "    "*"  "*"   "*"  " " " "      " "    " "    "*"   " "   
#11  ( 1 ) "*"   " "    "*"  "*"   "*"  " " " "      " "    " "    "*"   " "   
#12  ( 1 ) "*"   "*"    "*"  "*"   "*"  "*" " "      " "    "*"    "*"   "*"   
#13  ( 1 ) "*"   "*"    "*"  "*"   "*"  " " " "      " "    "*"    "*"   "*"   
#14  ( 1 ) "*"   "*"    "*"  "*"   "*"  " " " "      " "    "*"    "*"   "*"   
#15  ( 1 ) "*"   "*"    "*"  "*"   "*"  " " " "      " "    "*"    "*"   "*"   
#16  ( 1 ) "*"   "*"    "*"  "*"   "*"  "*" " "      " "    "*"    "*"   "*"   
#17  ( 1 ) "*"   "*"    "*"  "*"   "*"  "*" " "      " "    "*"    "*"   "*"   
#18  ( 1 ) "*"   "*"    "*"  "*"   "*"  "*" "*"      " "    "*"    "*"   "*"   
#19  ( 1 ) "*"   "*"    "*"  "*"   "*"  "*" "*"      " "    "*"    "*"   "*"   
#20  ( 1 ) "*"   "*"    "*"  "*"   "*"  "*" "*"      " "    "*"    "*"   "*"   
#21  ( 1 ) "*"   "*"    "*"  "*"   "*"  "*" "*"      " "    "*"    "*"   "*"   
#22  ( 1 ) "*"   "*"    "*"  "*"   "*"  "*" "*"      " "    "*"    "*"   "*"   
#23  ( 1 ) "*"   "*"    "*"  "*"   "*"  "*" "*"      "*"    "*"    "*"   "*"   
#          Pik3r3 Ptk2 Nras Nos3 Pik3r1 Pik3ca Ppp3cb Map2k2 Nfatc4 Mapk13 Rac1
#1  ( 1 )  "*"    " "  " "  " "  " "    " "    " "    " "    " "    " "    " " 
#2  ( 1 )  "*"    " "  " "  " "  " "    " "    " "    " "    " "    " "    "*" 
#3  ( 1 )  " "    " "  " "  " "  " "    " "    " "    " "    " "    " "    "*" 
#4  ( 1 )  "*"    " "  " "  " "  " "    " "    " "    " "    " "    " "    "*" 
#5  ( 1 )  "*"    " "  " "  " "  "*"    " "    " "    " "    " "    " "    "*" 
#6  ( 1 )  "*"    " "  " "  " "  "*"    " "    " "    " "    " "    " "    "*" 
#7  ( 1 )  " "    " "  " "  " "  "*"    " "    " "    " "    " "    " "    "*" 
#8  ( 1 )  "*"    " "  " "  " "  "*"    " "    " "    " "    " "    " "    "*" 
#9  ( 1 )  " "    "*"  " "  " "  " "    " "    "*"    " "    "*"    " "    "*" 
#10  ( 1 ) " "    "*"  " "  "*"  " "    " "    "*"    " "    "*"    " "    "*" 
#11  ( 1 ) " "    "*"  " "  "*"  " "    "*"    "*"    " "    "*"    " "    "*" 
#12  ( 1 ) " "    "*"  " "  "*"  " "    " "    " "    " "    " "    " "    "*" 
#13  ( 1 ) " "    "*"  " "  "*"  " "    "*"    "*"    " "    " "    " "    "*" 
#14  ( 1 ) " "    "*"  " "  "*"  " "    "*"    "*"    " "    "*"    " "    "*" 
#15  ( 1 ) "*"    "*"  " "  "*"  " "    "*"    "*"    " "    "*"    " "    "*" 
#16  ( 1 ) " "    "*"  " "  "*"  "*"    "*"    "*"    " "    "*"    " "    "*" 
#17  ( 1 ) "*"    "*"  " "  "*"  "*"    "*"    "*"    " "    "*"    " "    "*" 
#18  ( 1 ) "*"    "*"  " "  "*"  "*"    "*"    "*"    " "    "*"    " "    "*" 
#19  ( 1 ) "*"    "*"  " "  "*"  "*"    "*"    "*"    "*"    "*"    " "    "*" 
#20  ( 1 ) "*"    "*"  " "  "*"  "*"    "*"    "*"    "*"    "*"    "*"    "*" 
#21  ( 1 ) "*"    "*"  "*"  "*"  "*"    "*"    "*"    "*"    "*"    "*"    "*" 
#22  ( 1 ) "*"    "*"  "*"  "*"  "*"    "*"    "*"    "*"    "*"    "*"    "*" 
#23  ( 1 ) "*"    "*"  "*"  "*"  "*"    "*"    "*"    "*"    "*"    "*"    "*" 
#          Nfat5
#1  ( 1 )  " "  
#2  ( 1 )  " "  
#3  ( 1 )  " "  
#4  ( 1 )  " "  
#5  ( 1 )  " "  
#6  ( 1 )  " "  
#7  ( 1 )  " "  
#8  ( 1 )  " "  
#9  ( 1 )  " "  
#10  ( 1 ) " "  
#11  ( 1 ) " "  
#12  ( 1 ) " "  
#13  ( 1 ) " "  
#14  ( 1 ) " "  
#15  ( 1 ) " "  
#16  ( 1 ) " "  
#17  ( 1 ) " "  
#18  ( 1 ) " "  
#19  ( 1 ) " "  
#20  ( 1 ) " "  
#21  ( 1 ) " "  
#22  ( 1 ) "*"  
#23  ( 1 ) "*"  

plot(regsubsets.out)
#From the plot we see that the model with the lowest value of BIC is the model with only 4 variables
coef(regsubsets.out,4)
#(Intercept)        Akt2         Rik      Pik3r3        Rac1 
# -0.4446112  -0.4054803   0.2207164   0.2439920   0.3171645 

#Want to know the specific BIC value for the top two models 
#since the graphic doesn't show more decimal places.

#Top model based on picture
best_1<-lm(Mapk1~Akt2+Rik+Pik3r3+Rac1,data=full_data)
best_1_BIC<-BIC(best_1)
best_1_BIC #Output is -68.98881?!?!

#Second top model based on picture