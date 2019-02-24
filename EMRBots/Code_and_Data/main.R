#Let's set up our R environment

rm(list = ls())
library(plyr)
library(dplyr)
library(lubridate)

#Change the folder locally:
WORK_DIR = "C:/IBM/Harvard_Course_2019/EMRBots/Code/"

#Let's read in our dataset into a data frame
d.dx <- read.csv(paste(WORK_DIR, "AdmissionsDiagnosesCorePopulatedTable.txt", sep=""), row.names = NULL, sep = '\t', fileEncoding="UTF-8-BOM")
names(d.dx)
#PatientID, AdmissionID, PrimaryDiagnosisCode, PrimaryDiagnosisDescription

d.dx$PatientID_AND_AdmissionID <- paste(d.dx$PatientID, d.dx$AdmissionID, sep = "_")

# Get all patients with malignant neoplasm
d.mp <- d.dx %>%
  filter(grepl("Malignant neoplasm", PrimaryDiagnosisDescription))
head(d.mp)

d.mp_ids <- d.mp %>%
  select(PatientID_AND_AdmissionID) %>% unique()
head(d.mp_ids)

#There is a difference between patients and patient admissions.
cat("Number of admissions with malignant neoplasm",
    d.mp_ids %>% select(PatientID_AND_AdmissionID) %>% unique() %>% nrow(), "\n")
## Number of admissions with malignant neoplasm 4375

cat("Number of patients with malignant neoplasm",
    d.mp %>% select(PatientID) %>% unique() %>% nrow(), "\n")
## Number of patients with malignant neoplasm 3589

d.no_mp_ids <- d.dx %>%
  filter(!(PatientID_AND_AdmissionID %in% unlist(d.mp_ids$PatientID_AND_AdmissionID))) %>%
  select(PatientID_AND_AdmissionID) %>% unique()
cat("Number of admissions without malignant neoplasm",
    d.no_mp_ids %>% select(PatientID_AND_AdmissionID) %>% unique() %>% nrow(), "\n")
## Number of admissions without malignant neoplasm 31768

d.mp_ids$outcome <- 1
d.no_mp_ids$outcome <- 0
d.mp_outcome <- rbind(d.mp_ids, d.no_mp_ids)

write.csv(d.mp_outcome, paste(WORK_DIR, "mp_outcome.csv", sep=""),
          row.names = FALSE)

#Using merge()
d.enc_info <- read.table(paste(WORK_DIR, "AdmissionsCorePopulatedTable.txt", sep=""),
                         sep = "\t", header = TRUE, fileEncoding="UTF-8-BOM")
d.enc_info$PatientID_AND_AdmissionID <- paste(d.enc_info$PatientID, d.enc_info$AdmissionID, sep = "_")
names(d.enc_info)

d.mp_outcome <- read.csv(paste(WORK_DIR, "mp_outcome.csv", sep=""), fileEncoding="UTF-8-BOM")
names(d.mp_outcome)

d.demo <- read.table(paste(WORK_DIR, "PatientCorePopulatedTable.txt", sep=""),
                     sep = "\t", header = TRUE, fileEncoding="UTF-8-BOM")
names(d.demo)

d.cohort <- merge(d.enc_info, d.mp_outcome, by = c("PatientID_AND_AdmissionID"))
d.cohort <- merge(d.cohort, d.demo, by = c("PatientID"))

names(d.cohort)

#Using mutate()
before_mutate = head(d.cohort)
Sys.setenv(tz = "America/Chicago")
d.cohort <- d.cohort %>%
  mutate(AdmissionStartDate = ymd_hms(AdmissionStartDate),
         AdmissionEndDate = ymd_hms(AdmissionEndDate),
         PatientDateOfBirth = ymd_hms(PatientDateOfBirth))

after_mutate = head(d.cohort)

#Calculating age and LOS.
d.cohort <- d.cohort %>%
  mutate(PatientAge = interval(PatientDateOfBirth, AdmissionStartDate) / dyears(1))
d.cohort <- d.cohort %>%
  mutate(LOS = interval(AdmissionStartDate, AdmissionEndDate) / ddays(1))

#Check for consistency: continous ariables
quantile(d.cohort$PatientAge)

#Check for consistency: categorical variables
summary(d.cohort$PatientRace)

#e can start constructing our Table 1.
d.1 <- filter(d.cohort, outcome == 1) #4375 patients
d.0 <- filter(d.cohort, outcome == 0) #36143 patients

cat("Mean age, outcome 1: ", mean(d.1$PatientAge), "\n")
## Mean age, outcome 1: 42.06632

cat("Mean age, outcome 0: ", mean(d.0$PatientAge), "\n")
## Mean age, outcome 0: 41.70317

cat("SD age, outcome 1: ", sd(d.1$PatientAge), "\n")
## SD age, outcome 1: 18.17401

cat("SD age, outcome 1: ", sd(d.0$PatientAge), "\n")
## SD age, outcome 1: 18.04217

#Comparison of continous variables
#Assess differences in age
print(t.test(d.1$PatientAge, d.0$PatientAge))

#Comparision of categorical variables
gender.table <- with(d.cohort, table(outcome, PatientGender))
gender.table

chisq.test(gender.table)

d.most_recent_labs <- read.table(paste(WORK_DIR, "MostRecentLabsPerAdmission.txt", sep=""),
                     sep = "\t", header = TRUE, fileEncoding="UTF-8-BOM")

d.most_recent_labs$PatientID_AND_AdmissionID <- paste(d.most_recent_labs$PatientID, d.most_recent_labs$AdmissionID, sep = "_")

d.cohort <- merge(d.cohort, d.most_recent_labs, by = c("PatientID_AND_AdmissionID"))

names(d.cohort)

d.cohort_for_logistic_regression <- d.cohort[,c("outcome","PatientGender","PatientAge","PatientRace","LOS", "Most.Recent.CBC..ABSOLUTE.NEUTROPHIL", "Most.Recent.CBC..HEMATOCRIT", "Most.Recent.CBC..HEMOGLOBIN", "Most.Recent.CBC..PLATELET.COUNT", "Most.Recent.CBC..RED.BLOOD.CELL.COUNT", "Most.Recent.CBC..WHITE.BLOOD.CELL.COUNT", "Most.Recent.METABOLIC..ALBUMIN","Most.Recent.METABOLIC..BILI.TOTAL","Most.Recent.METABOLIC..CALCIUM","Most.Recent.METABOLIC..POTASSIUM", "Most.Recent.METABOLIC..SODIUM","Most.Recent.URINALYSIS..PH")]

lm_model <- glm(as.factor(outcome) ~
            as.factor(PatientGender)
          + as.numeric(PatientAge)
          + as.factor(PatientRace)
          + as.numeric(LOS)
          + as.numeric(Most.Recent.CBC..ABSOLUTE.NEUTROPHIL)
          + as.numeric(Most.Recent.CBC..HEMATOCRIT)
          + as.numeric(Most.Recent.CBC..HEMOGLOBIN)
          + as.numeric(Most.Recent.CBC..PLATELET.COUNT)
          + as.numeric(Most.Recent.CBC..RED.BLOOD.CELL.COUNT)
          + as.numeric(Most.Recent.CBC..WHITE.BLOOD.CELL.COUNT)
          + as.numeric(Most.Recent.METABOLIC..ALBUMIN)
          + as.numeric(Most.Recent.METABOLIC..BILI.TOTAL)
          + as.numeric(Most.Recent.METABOLIC..CALCIUM)
          + as.numeric(Most.Recent.METABOLIC..POTASSIUM)
          + as.numeric(Most.Recent.METABOLIC..SODIUM)
          + as.numeric(Most.Recent.URINALYSIS..PH)
          , data = d.cohort_for_logistic_regression, family = "binomial")

summary(lm_model)

exp(cbind("Odds ratio" = coef(lm_model), confint.default(lm_model, level = 0.95)))


