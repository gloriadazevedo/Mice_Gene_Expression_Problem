###File for the data analysis and model evaluation##

#Naive way is to implement best subset selection since there are only 24 variables
library(leaps)

#The predictors are everything except for the field Mapk1
#Also can't have categorical data in there
x_data<-subset(full_data, select = -c(Mapk1,Gene.Name) )
y_data<-full_data$"Mapk1"
best_subset<-leaps(x_data,y_data,nbest=1,method="adjr2",names=colnames(x_data))