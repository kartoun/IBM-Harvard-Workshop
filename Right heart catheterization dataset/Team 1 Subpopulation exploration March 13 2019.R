#Propensity score matching in R

#Let's set up our R environment
rm(list = ls()) # Clean all varaibles.
if(!is.null(dev.list())) dev.off() # Clean all plots.
WORK_DIR = "C:/IBM/Harvard_Course_2019/Casualty/Code/"

#Load packages
library(tableone) # The tableone package is an R package that eases the construction of "Table 1", i.e., patient baseline characteristics table commonly found in biomedical research papers.
library(MatchIt) # Selects matched samples of the original treated and control groups with similar covariate distributions.
library(ggplot2) # A library that allows us to plot odds ratios.

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
#colcan<-as.numeric(rhc$cat1 == 'Colon Cancer')
Coma<-as.numeric(rhc$cat1 == 'Coma')
COPD<-as.numeric(rhc$cat1 == 'COPD')
MOSF<-as.numeric(rhc$cat1 == 'MOSF w/Malignancy') #Multiple organ system failure
sepsis<-as.numeric(rhc$cat1 == 'MOSF w/Sepsis')
female<-as.numeric(rhc$sex == 'Female')
died<-as.numeric(rhc$death == 'Yes')
age<-rhc$age
treatment<-as.numeric(rhc$swang1 == 'RHC')
meanbp1<-rhc$meanbp1
aps<-rhc$aps1 #APACHE score
hrt1<-rhc$hrt1
white<-as.numeric(rhc$race == 'white')
black<-as.numeric(rhc$race == 'black')
bili1<-rhc$bili1
Estimated_MELD<- 10 * ((0.957 * log(rhc$crea1, 2.71828)) + (0.378 * log(rhc$bili1, 2.718)) + (1.12 * log(1.499619615, 2.718))) + 6.43

L_Estimated_MELD_Plus = 8.53499496 + 2.06503238 * log10(1 + rhc$bili1) + 2.5967965 * log10(1 + rhc$crea1) - 6.34990436 * log10(1 + rhc$alb1) + 2.99724802 * log10(1 + 1.499619615) + 1.92811726 * log10(1 + rhc$wblc1) + 0.04070442 * rhc$age - 6.47834101 * log10(1 + rhc$sod1)
Estimated_MELD_Plus = exp(L_Estimated_MELD_Plus) / (1 + exp(L_Estimated_MELD_Plus))

liverhx<-as.numeric(rhc$liverhx)

#Creating the dataset that we will use:
mydata <- cbind(ARF, CHF, Cirr, Coma, COPD, MOSF, sepsis,
             age, female, meanbp1, aps, hrt1, white, black, bili1, Estimated_MELD, Estimated_MELD_Plus, liverhx, treatment, died)
mydata = data.frame(mydata)

#Here you can define a sub-population:

#Is RHC harmful for women?
#mydata <- mydata[ which(mydata$female == '1'), ]  

mydata <- mydata[ which(mydata$Estimated_MELD > 20), ]

#mydata <- mydata[ which(mydata$female == '1' 
#                        & mydata$age > 65
#                        & mydata$Cirr == '1'), ]

#What is the prevalence of a variable? We want to remove covariates with
#low prevalence (e.g., less than 1.0%, less than 5.0%)
prev_of_a_variable_in_percent <- 100 * nrow(mydata[ which(mydata$sepsis == '1'), ]) / nrow(mydata)
prev_of_a_variable_in_percent

library(Matching)

psmodel<-glm(treatment~ARF+CHF+Cirr+Coma+COPD+MOSF+
             sepsis+age+meanbp1+aps+hrt1+white+black+bili1+Estimated_MELD+Estimated_MELD_Plus+liverhx,
             family = binomial(), data = mydata)

#Get propensity scores for each patient:
pscore=predict(psmodel, mydata, type="response", na.action = na.pass)

logit_pscore = log(pscore / (1 - pscore))

psmatch = Match(Tr = mydata$treatment, M = 1, X = logit_pscore, replace = FALSE, caliper = 0.2)

#The following is the final data frame that contains cases and their correspoding matched controls.
matched <- mydata[unlist(psmatch[c("index.treated", "index.control")]), ]

xvars <- c("ARF", "CHF", "Cirr", "Coma", "COPD", "MOSF", "sepsis", "age", "meanbp1", "aps", "hrt1", "white", "black", "bili1", "Estimated_MELD", "Estimated_MELD_Plus", "liverhx")

tab1_before_matching = CreateTableOne(vars = xvars, strata = "treatment"
                      , data = mydata, test = FALSE)

print(tab1_before_matching, smd = TRUE)

tab1_after_matching = CreateTableOne(vars = xvars, strata = "treatment"
                             , data = matched, test = FALSE)

print(tab1_after_matching, smd = TRUE)

final_model <- glm(died ~ treatment+ARF+CHF+Cirr+Coma+COPD+MOSF+
               sepsis+age+meanbp1+aps+hrt1+white+black+bili1+Estimated_MELD+Estimated_MELD_Plus+liverhx,
               family = quasibinomial(), data = matched)

summary(final_model)





