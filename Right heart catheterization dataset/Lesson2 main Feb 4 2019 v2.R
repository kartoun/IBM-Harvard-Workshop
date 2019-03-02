#Inverse Probability of Treatment Weighting (IPTW)

#Let's set up our R environment
rm(list = ls())
WORK_DIR = "C:/IBM/Harvard_Course_2019/Casualty/Code/"

#Load packages
library(tableone)
library(ipw) #to get inverse probability weights.
library(sandwich) #for robust variance estimation.
library(survey) #to get weighted estimators.

#Read the data:
#http://biostat.mc.vanderbilt.edu/wiki/Main/DataSets
#load(url("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.sav"))
rhc <- read.csv(paste(WORK_DIR, "rhc.csv", sep=""), sep = ','); rhc$X<-NULL

#View the data:
#View(rhc)

#Create a data set with just the variables we are interested at and convert them to a numeric data type:
ARF<-as.numeric(rhc$cat1 == 'ARF')
CHF<-as.numeric(rhc$cat1 == 'CHF')
Cirr<-as.numeric(rhc$cat1 == 'Cirrhosis')
colcan<-as.numeric(rhc$cat1 == 'Colon Cancer')
Coma<-as.numeric(rhc$cat1 == 'Coma')
COPD<-as.numeric(rhc$cat1 == 'COPD')
lungcan<-as.numeric(rhc$cat1 == 'Lung Cancer')
MOSF<-as.numeric(rhc$cat1 == 'MOSF w/Malignancy')
sepsis<-as.numeric(rhc$cat1 == 'MOSF w/Sepsis')
female<-as.numeric(rhc$sex == 'Female')
died<-as.numeric(rhc$death == 'Yes')
age<-rhc$age
treatment<-as.numeric(rhc$swang1 == 'RHC')
meanbp1<-rhc$meanbp1
aps<-rhc$aps1

#Creating the dataset that we will use:
mydata <- cbind(ARF, CHF, Cirr, colcan, Coma, lungcan, MOSF, sepsis,
             age, female, meanbp1, aps, treatment, died)
mydata = data.frame(mydata)
#View(mydata)
#dim(mydata)
#summary(mydata)

#Fitting a propensity score model using logistic regression:
psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+
               sepsis+age+female+meanbp1+aps,
             family = binomial(), data = mydata)

#Show coefficients. Who is more likely to get treated?
summary(psmodel)

#Get propensity scores:
ps = predict(psmodel, mydata, type="response", na.action = na.pass)

ps=data.frame(ps);
mydata <- data.frame(cbind(mydata, ps))

df_treated = mydata[which(mydata$treatment=='1'),]
hist(df_treated$ps, xlim = c(0, 1), ylim = c(0, 400), col = "red")

df_untreated = mydata[which(mydata$treatment=='0'),]
hist(df_untreated$ps, xlim = c(0, 1), ylim = c(0, 400), col = "green")

hist(df_treated$ps, xlim = c(0, 1), ylim = c(0, 400), col=rgb(0,0,1,1/4))
hist(df_untreated$ps, xlim = c(0, 1), ylim = c(0, 400), col=rgb(1,0,0,1/4), add = T)

#Geting into weights

#Create weights:
weight <- ifelse(mydata$treatment == 1, 1 / mydata$ps, 1 / (1 - mydata$ps))

#apply weights to the data:
weighteddata <- svydesign(ids=~1, data = mydata, weights = ~weight)

#Weighted table 1:
xvars <- c("ARF", "CHF", "Cirr", "colcan", "Coma", "lungcan", "MOSF", "sepsis", "age", "female", "meanbp1", "aps")
weightedtable <- svyCreateTableOne(xvars, strata='treatment', data = weighteddata, test=FALSE)

#Show table with SMD:
print(weightedtable, smd = TRUE) #ignore sample sizes, and ignore standard deviations.

#Stratified by treatment
#                     0               1               SMD 
#n                   5760.80         5660.88               
#ARF (mean (sd))        0.43 (0.50)     0.45 (0.50)   0.021
#CHF (mean (sd))        0.08 (0.27)     0.08 (0.27)   0.001
#Cirr (mean (sd))       0.04 (0.19)     0.04 (0.18)   0.017
#colcan (mean (sd))     0.00 (0.04)     0.00 (0.06)   0.039
#Coma (mean (sd))       0.07 (0.26)     0.07 (0.25)   0.025
#lungcan (mean (sd))    0.01 (0.08)     0.01 (0.08)   0.010
#MOSF (mean (sd))       0.07 (0.25)     0.07 (0.26)   0.008
#sepsis (mean (sd))     0.22 (0.41)     0.22 (0.41)   0.004
#age (mean (sd))       61.37 (17.59)   61.52 (15.22)  0.010
#female (mean (sd))     0.44 (0.50)     0.44 (0.50)  <0.001
#meanbp1 (mean (sd))   78.28 (38.20)   78.14 (38.34)  0.004
#aps (mean (sd))       55.05 (20.43)   55.42 (19.84)  0.018

final_model = svyglm(formula = died ~ treatment+ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis+age+female+meanbp1+aps
, weighteddata)
#View(weighteddata$variables)

exp(final_model$coefficients)
summary(final_model)

results = data.frame(summary(final_model)$coefficients)
results$Estimate = round(exp(results$Estimate), 2)

results

