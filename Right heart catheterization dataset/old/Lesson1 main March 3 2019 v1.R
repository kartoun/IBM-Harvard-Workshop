#Propensity score matching in R

#Let's set up our R environment
rm(list = ls())
WORK_DIR = "C:/IBM/Harvard_Course_2019/Casualty/Code/"

#Load packages
library(tableone) # The tableone package is an R package that eases the construction of "Table 1", i.e., patient baseline characteristics table commonly found in biomedical research papers.
library(MatchIt) # Selects matched samples of the original treated and control groups with similar covariate distributions.

#Read the data:
#http://biostat.mc.vanderbilt.edu/wiki/Main/DataSets
#load(url("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.sav"))
rhc <- read.csv(paste(WORK_DIR, "rhc_uk.csv", sep=""), sep = ','); rhc$X<-NULL #Note that the 4 coulmns with dates may look like they are integers - use Excel to format these columns to standard dates.
col_names = data.frame(colnames(rhc))

#View the data:
#View(rhc)

#Create a data set with just the variables we are interested at and convert them to a numeric data type:
ARF<-as.numeric(rhc$cat1 == 'ARF') #Acute renal failure
CHF<-as.numeric(rhc$cat1 == 'CHF') #Congestive heart failure
Cirr<-as.numeric(rhc$cat1 == 'Cirrhosis')
colcan<-as.numeric(rhc$cat1 == 'Colon Cancer')
Coma<-as.numeric(rhc$cat1 == 'Coma')
COPD<-as.numeric(rhc$cat1 == 'COPD')
lungcan<-as.numeric(rhc$cat1 == 'Lung Cancer')
MOSF<-as.numeric(rhc$cat1 == 'MOSF w/Malignancy') #Multiple organ system failure
sepsis<-as.numeric(rhc$cat1 == 'MOSF w/Sepsis')
female<-as.numeric(rhc$sex == 'Female')
died<-as.numeric(rhc$death == 'Yes')
age<-rhc$age
treatment<-as.numeric(rhc$swang1 == 'RHC')
meanbp1<-rhc$meanbp1
aps<-rhc$aps1 #APACHE score

print(tab1_before_matching, smd = TRUE)

#Creating the dataset that we will use:
mydata <- cbind(ARF, CHF, Cirr, colcan, Coma, lungcan, MOSF, sepsis,
             age, female, meanbp1, aps, treatment, died)
mydata = data.frame(mydata)
#View(mydata)
#dim(mydata)
#summary(mydata)

#Fitting a propensity score model:
psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+
               sepsis+age+female+meanbp1+aps,
             family = binomial(), data = mydata)

#Show coefficients. Who is more likely to get treated?
summary(psmodel)

pscore<-psmodel$fitted.values # Or use: pscore=predict(psmodel, mydata, type="response", na.action = na.pass)

#Get propensity scores:

pscore=data.frame(pscore);
mydata <- data.frame(cbind(mydata, pscore))

df_treated = mydata[which(mydata$treatment=='1'),]
hist(df_treated$pscore, xlim = c(0, 1), ylim = c(0, 400), col = "red")

df_untreated = mydata[which(mydata$treatment=='0'),]
hist(df_untreated$pscore, xlim = c(0, 1), ylim = c(0, 400), col = "green")

hist(df_treated$pscore, xlim = c(0, 1), ylim = c(0, 400), col=rgb(0,0,1,1/4))
hist(df_untreated$pscore, xlim = c(0, 1), ylim = c(0, 400), col=rgb(1,0,0,1/4), add = T)

#Use matchit for propensity score, nearest neighbor matching
#This function caclulates the propensity scores and then match them
m.out = matchit(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+
                  sepsis+age+female+meanbp1+aps
                 , data = mydata, method = "nearest") # Another method is called "optimal".

summary(m.out)

#propensity score plots
plot(m.out, type = "jitter")
plot(m.out, type = "hist")

#Get ids of treated and matched untreated:
match_matrix = data.frame(m.out$match.matrix)

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
#psmatch = Match(Tr = mydata$treatment, M = 1, X = logit_pscore, replace = FALSE, caliper = 0.2)

#The following is the final data frame that contains cases and their correspoding matched controls.
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
#dim(data.frame(y_trt))

y_con <- matched$died[matched$treatment == 0]
#dim(data.frame(y_con))

#### pairwise difference
diffy <- y_trt - y_con
#### paired t-t
t.test(diffy)

####################################################################









