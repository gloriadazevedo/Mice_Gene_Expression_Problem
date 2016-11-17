#In order to get a simple model, two approaches have been applied: Using a lasso regularizer and using BIC (which penalizes having a large number of parameters) for exhaustive subset selection, to find the most accurate and simple model
setwd("/Users/pihu_yadav/Desktop/Mice_Gene_Expression_Problem-master")
full_data<-read.csv("Kidney_2.csv")

#Performing Lasso Regression
library(glmnet)
x_data<-subset(full_data, select = -c(Mapk1,Gene.Name) )
y_data<-full_data$"Mapk1"
X <- as.matrix(x_data)
Y<-y_data
cvfit <- glmnet::cv.glmnet(X, Y)
coef(cvfit, s = "lambda.1se")
#24 x 1 sparse Matrix of class "dgCMatrix"
                     1
#(Intercept) 0.01017884
#Cdc42       .         
#Pla2g6      .         
#Akt2        .         
#Plcg2       .         
#Rac2        .         
#Rik         0.04474769
#Mapkapk2    .         
#Pik3cd      .         
#Pla2g5      .         
#Sphk2       .         
#Map2k1      .         
#Pik3r3      0.13293097
#Ptk2        .         
#Nras        .         
#Nos3        .         
#Pik3r1      .         
#Pik3ca      .         
#Ppp3cb      .         
#Map2k2      .         
#Nfatc4      .         
#Mapk13      .         
#Rac1        0.04906813
#Nfat5       .         

Y_pred<-predict(cvfit, X)
mean((Y_pred-Y)^2)
#[1] 0.01096602
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
