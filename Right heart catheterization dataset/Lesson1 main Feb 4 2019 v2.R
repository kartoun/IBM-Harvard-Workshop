#Propensity score matching in R

#Let's set up our R environment
rm(list = ls())
WORK_DIR = "C:/IBM/Harvard_Course_2019/Casualty/Code/"

#Load packages
library(tableone)
library(MatchIt)

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
pscore=predict(psmodel, mydata, type="response", na.action = na.pass)

pscore=data.frame(pscore);
mydata <- data.frame(cbind(mydata, pscore))

df_treated = mydata[which(mydata$treatment=='1'),]
hist(df_treated$pscore, xlim = c(0, 1), ylim = c(0, 400), col = "red")

df_untreated = mydata[which(mydata$treatment=='0'),]
hist(df_untreated$pscore, xlim = c(0, 1), ylim = c(0, 400), col = "green")

hist(df_treated$pscore, xlim = c(0, 1), ylim = c(0, 400), col=rgb(0,0,1,1/4))
hist(df_untreated$pscore, xlim = c(0, 1), ylim = c(0, 400), col=rgb(1,0,0,1/4), add = T)

#Use matchit for propensity score, nearesr neighbor matching
#This function caclulates the propensity scores and then match them
m.out = matchit(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+
                  sepsis+age+female+meanbp1+aps
                 , data = mydata, method = "nearest")

summary(m.out)

#propensity score plots
plot(m.out, type = "hist")

#Get ids of treated and matched untreated:
match_matrix = m.out$match.matrix

#Example 1:

#Treated:
mydata[5729,]

#Matched untreared:
mydata[1909,]

#Example 2:

#Treated:
mydata[5726,]

#Matched untreared:
mydata[4276,]

############################################

# A different type of matching using a caliper:

library(Matching)

psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+
               sepsis+age+female+meanbp1+aps,
             family = binomial(), data = mydata)

#Get propensity scores:
pscore=predict(psmodel, mydata, type="response", na.action = na.pass)

#do greedy matching on logit(PS)

logit_pscore = log(pscore / (1 - pscore))

psmatch = Match(Tr = mydata$treatment, M = 1, X = logit_pscore, replace = FALSE)
#psmatch = Match(Tr = mydata$treatment, M = 1, X = logit_pscore, replace = FALSE, caliper = .2)

matched <- mydata[unlist(psmatch[c("index.treated", "index.control")]), ]

xvars <- c("ARF", "CHF", "Cirr", "colcan", "Coma", "lungcan", "MOSF", "sepsis", "age", "female", "meanbp1", "aps")

tab1_before_matching = CreateTableOne(vars = xvars, strata = "treatment"
                      , data = mydata, test = FALSE)

print(tab1_before_matching, smd = TRUE)

tab1_after_matching = CreateTableOne(vars = xvars, strata = "treatment"
                             , data = matched, test = FALSE)

print(tab1_after_matching, smd = TRUE)

#Using a caliper:
#Add "caliper = .2" to the line starting with "psmatch" above.
#What is "0.2"? Answer: 0.2 standard deviation units.
#Standard mean deviation threshold of 0.1 is an acceptable threshold.

#### outcome analysis
y_trt <- matched$died[matched$treatment == 1]
dim(data.frame(y_trt))

y_con <- matched$died[matched$treatment == 0]
dim(data.frame(y_con))

#### pairwise difference
diffy <- y_trt - y_con
#### paired t-t
t.test(diffy)

####################################################################

final_model = glm(died ~ treatment+ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis+age+female+meanbp1+aps
                         , matched, family = binomial())

exp(final_model$coefficients)
summary(final_model)

results = data.frame(summary(final_model)$coefficients)
results$Estimate = round(exp(results$Estimate), 2)

results







