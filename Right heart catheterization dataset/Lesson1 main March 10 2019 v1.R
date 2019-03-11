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
MOSF<-as.numeric(rhc$cat1 == 'MOSF w/Malignancy') #Multiple organ system failure
sepsis<-as.numeric(rhc$cat1 == 'MOSF w/Sepsis')
female<-as.numeric(rhc$sex == 'Female')
died<-as.numeric(rhc$death == 'Yes')
age<-rhc$age
treatment<-as.numeric(rhc$swang1 == 'RHC')
meanbp1<-rhc$meanbp1
aps<-rhc$aps1 #APACHE score

hrt1<-rhc$hrt1
#Creating the dataset that we will use:
mydata <- cbind(ARF, CHF, Cirr, colcan, Coma, MOSF, sepsis,
             age, female, meanbp1, aps, hrt1, treatment, died)
mydata = data.frame(mydata)
#View(mydata)
#dim(mydata)
#summary(mydata)

#Fitting a propensity score model:
psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+MOSF+
               sepsis+age+female+meanbp1+aps+hrt1,
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
m.out = matchit(treatment~ARF+CHF+Cirr+colcan+Coma+MOSF+
                  sepsis+age+female+meanbp1+aps+hrt1
                 , data = mydata, method = "nearest") # Another method is called "optimal".

summary(m.out)

#propensity score plots
plot(m.out, type = "jitter")
plot(m.out, type = "hist")

### density plot r ###

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

#mydata <- mydata[ which(mydata$female == '1' 
#                        & mydata$age > 65
#                        & mydata$Cirr == '1'), ]

#mydata <- mydata[ which(mydata$female == '0'), ]

psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+MOSF+
               sepsis+age+female+meanbp1+aps+hrt1,
             family = binomial(), data = mydata)

#Get propensity scores:
pscore=predict(psmodel, mydata, type="response", na.action = na.pass)

#do greedy matching on logit(PS)

logit_pscore = log(pscore / (1 - pscore))

#psmatch = Match(Tr = mydata$treatment, M = 1, X = logit_pscore, replace = FALSE)
psmatch = Match(Tr = mydata$treatment, M = 1, X = logit_pscore, replace = FALSE, caliper = 0.2)

#The following is the final data frame that contains cases and their correspoding matched controls.
matched <- mydata[unlist(psmatch[c("index.treated", "index.control")]), ]

xvars <- c("ARF", "CHF", "Cirr", "colcan", "Coma", "MOSF", "sepsis", "age", "female", "meanbp1", "aps", "hrt1")

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

final_model <- glm(died ~ treatment+ARF+CHF+Cirr+colcan+Coma+MOSF+
               sepsis+age+female+meanbp1+aps+hrt1,
               family = binomial(), data = matched)


#final_model <- glm(died ~ treatment+ARF+CHF+Cirr+colcan+Coma+MOSF+
#                     sepsis+age+female+meanbp1+aps+hrt1,
#                   family = binomial(), data = mydata)

# plot odds ratios
library(ggplot2)

results = data.frame(summary(final_model)$coefficients)
results$Estimate = exp(results$Estimate)

conf_int = confint(final_model)
conf_int = exp(conf_int)

df_ORs_Conf_Int = cbind(results$Estimate
                        ,conf_int[,1]
                        ,conf_int[,2])

df_ORs_Conf_Int_sorted = df_ORs_Conf_Int[order(-results$Estimate), ]

row.names(df_ORs_Conf_Int_sorted)

df_ORs_Conf_Int_sorted = df_ORs_Conf_Int_sorted[!(row.names(df_ORs_Conf_Int_sorted) == '(Intercept)'), ]

boxLabels = rownames(df_ORs_Conf_Int_sorted)

yAxis = length(boxLabels):1

df <- data.frame(
  yAxis = length(boxLabels):1,
  boxOdds = c(df_ORs_Conf_Int_sorted[, 1]),
  boxCILow = c(df_ORs_Conf_Int_sorted[, 2]),
  boxCIHigh = c(df_ORs_Conf_Int_sorted[, 3])
)

p <- ggplot(df, aes(x = boxOdds , y = yAxis))
p + geom_vline(aes(xintercept = 1.0), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = 1.0, height = .2, color = "gray50") +
  geom_point(size = 5.0, color = "black") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = yAxis, labels = boxLabels) +
  scale_x_continuous(breaks = seq(0, 7, 0.5) ) +
  ylab("") +
  xlab("Odds ratio and 95% confidence intervals") +
  theme(axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title = element_text(face = "bold", size = 10))

dev.copy(jpeg, filename = 'Atrius_HTN_odds_ratio_plot.jpg', width = 1440, height = 900);
dev.off()

########################################################################




