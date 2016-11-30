###Import data script for the Mice Gene Expression project##
##Note that the data is already edited through Microsoft Excel to have the Expressions 
#as the data points and the values of each of the genes for each expression as the fields
#for each of the expressions

#This script should be run in R before the analysis for standard import 
#instead of doing it by hand every time

#Set the working directory for ease of use
setwd("/Users/glori/Documents/GitHub/Mice_Gene_Expression_Problem")

#Actual import: .csv file
full_data<-read.table("Kidney_2.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)

#The data has some really nice characteristics like there are no missing values
#or other weird inconsistencies (so far)
#However, note that our sample size is fairly small for the number of fields (columns) 
#thus we will need to use some resampling techniques such as cross validation and bootstrap
#as well as potentially force the model to be small for lower test error

#Libraries used
library(class)