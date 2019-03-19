#Propensity score matching in R

#Let's set up our R environment
rm(list = ls()) # Clean all varaibles.
if(!is.null(dev.list())) dev.off() # Clean all plots.
WORK_DIR = "C:/IBM/Harvard_Course_2019/Casualty/Code/"
library(Matching)

rhc <- read.csv(paste(WORK_DIR, "rhc_uk.csv", sep=""), sep = ','); rhc$X<-NULL #Note that the 4 coulmns with dates may look like they are integers - use Excel to format these columns to standard dates.

ARF<-as.numeric(rhc$cat1 == 'ARF') #Acute renal failure
CHF<-as.numeric(rhc$cat1 == 'CHF') #Congestive heart failure
Cirr<-as.numeric(rhc$cat1 == 'Cirrhosis')
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

psmodel<-glm(treatment~ARF+CHF+Cirr+Coma+COPD+MOSF+
             sepsis+age+female+meanbp1+aps+hrt1+white+black+bili1+Estimated_MELD+Estimated_MELD_Plus+liverhx,
             family = binomial(), data = mydata)

#Get propensity scores for each patient:
pscore=predict(psmodel, mydata, type="response", na.action = na.pass)
#pscore

logit_pscore = log(pscore / (1 - pscore))

results_df <- data.frame(OR = double(), CI_Low = double(), CI_High = double())
num_iterations = 5

for (i in 1:num_iterations)
{

psmatch = Match(Tr = mydata$treatment, M = 5, X = logit_pscore, replace = FALSE, caliper = 0.2)

matched <- mydata[unlist(psmatch[c("index.treated", "index.control")]), ]

xvars <- c("ARF", "CHF", "Cirr", "Coma", "COPD", "MOSF", "sepsis", "age", "female", "meanbp1", "aps", "hrt1", "white", "black", "bili1", "Estimated_MELD", "Estimated_MELD_Plus", "liverhx")

final_model <- glm(died ~ treatment+ARF+CHF+Cirr+Coma+COPD+MOSF+
               sepsis+age+female+meanbp1+aps+hrt1+white+black+bili1+Estimated_MELD+Estimated_MELD_Plus+liverhx,
               family = binomial(), data = matched)

#summary(final_model)

results = data.frame(summary(final_model)$coefficients)
results$Estimate = exp(results$Estimate)

conf_int = confint(final_model)
conf_int = exp(conf_int)

df_ORs_Conf_Int = cbind(results$Estimate
                        ,conf_int[,1]
                        ,conf_int[,2])

df_ORs_Conf_Int_sorted = df_ORs_Conf_Int[order(-results$Estimate), ]
df_ORs_Conf_Int_sorted = data.frame(df_ORs_Conf_Int_sorted[!(row.names(df_ORs_Conf_Int_sorted) == '(Intercept)'), ])

current_OR = df_ORs_Conf_Int_sorted[c("treatment"),]$X1
current_CI_Low = df_ORs_Conf_Int_sorted[c("treatment"),]$X2
current_CI_High = df_ORs_Conf_Int_sorted[c("treatment"),]$X3

tmpdf = data.frame(current_OR, current_CI_Low, current_CI_High)
results_df <- rbind(results_df, tmpdf)
print(results_df)

}
  
#Present final results:
names(results_df)<-c("OR", "CI_Low", "CI_High")
mean(results_df$OR)
mean(results_df$CI_Low)
mean(results_df$CI_High)


