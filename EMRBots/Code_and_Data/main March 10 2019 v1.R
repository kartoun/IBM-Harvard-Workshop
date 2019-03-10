#Let's set up our R environment

rm(list = ls()) # Clean all variables
library(plyr) # A data manipulation library.
library(dplyr) # A data manipulation library.
library(lubridate) # A library that allows handling dates.

#Change the folder locally:
WORK_DIR = "C:/IBM/Harvard_Course_2019/EMRBots/Code/"

#Let's read in our dataset into a data frame
d.dx <- read.csv(paste(WORK_DIR, "AdmissionsDiagnosesCorePopulatedTable.txt", sep=""), row.names = NULL, sep = '\t', fileEncoding="UTF-8-BOM")
names(d.dx)
#PatientID, AdmissionID, PrimaryDiagnosisCode, PrimaryDiagnosisDescription

#The following command adds an additional column to the data frame. This addtional column contains a combination of the patient id with an admission id.
d.dx$PatientID_AND_AdmissionID <- paste(d.dx$PatientID, d.dx$AdmissionID, sep = "_")

#Get all patients with malignant neoplasm. Basicall this function filters out all admissions that are not related to "Malignant neoplasm".
d.mp <- d.dx %>% filter(grepl("Malignant neoplasm", PrimaryDiagnosisDescription)) #grepl returns TRUE if a string contains the pattern, otherwise FALSE;

#Present top few records.
head(d.mp, 10)

# Store all unique combinations of patient ids and admission ids.
d.mp_ids <- d.mp %>%
  select(PatientID_AND_AdmissionID) %>% unique()

head(d.mp_ids, 10)

#Count total number of all patients of interest and their assocaited admissions.
#This number equals to the total number of admissions with malignant neoplasm.
cat("Number of admissions with malignant neoplasm",
    d.mp_ids %>% select(PatientID_AND_AdmissionID) %>% unique() %>% nrow(), "\n")

## Number of admissions with malignant neoplasm 4375.

#Count total number of unique patients. It turns out that some patients had 2 or more admissions for malignant neoplasm.
cat("Number of patients with malignant neoplasm",
    d.mp %>% select(PatientID) %>% unique() %>% nrow(), "\n")

## Number of unique patients with malignant neoplasm 3589.

##Let's identify all patients and their associated admissions with no malignant neoplasm.
d.no_mp_ids <- d.dx %>%
  filter(!(PatientID_AND_AdmissionID %in% unlist(d.mp_ids$PatientID_AND_AdmissionID))) %>%
  select(PatientID_AND_AdmissionID) %>% unique()

cat("Number of admissions without malignant neoplasm",
    d.no_mp_ids %>% select(PatientID_AND_AdmissionID) %>% unique() %>% nrow(), "\n")
## Number of admissions without malignant neoplasm 31768

d.mp_ids$outcome <- 1 # This commands adds an outcome column equals to 1.
d.no_mp_ids$outcome <- 0 # This commands adds an outcome column equals to 0.

#This combines the outcome with its combination of patient id and admission id.
d.mp_outcome <- rbind(d.mp_ids, d.no_mp_ids)

#Write the combination patient id, admission id, and outcome, to a file.
write.csv(d.mp_outcome, paste(WORK_DIR, "mp_outcome.csv", sep=""),
          row.names = FALSE)

#Load all admission details including start and end date per admission. 
d.enc_info <- read.table(paste(WORK_DIR, "AdmissionsCorePopulatedTable.txt", sep=""),
                         sep = "\t", header = TRUE, fileEncoding="UTF-8-BOM")

#Add another column to the data frame.
#The new column contains a unique combination of patient id and admission id.
d.enc_info$PatientID_AND_AdmissionID <- paste(d.enc_info$PatientID, d.enc_info$AdmissionID, sep = "_")
names(d.enc_info)

#Load the data frame that contains a uniqe combination of patient id, admission id, and outcome.
d.mp_outcome <- read.csv(paste(WORK_DIR, "mp_outcome.csv", sep=""), fileEncoding="UTF-8-BOM")
names(d.mp_outcome)

#Load demograhpic details for all patients.
d.demo <- read.table(paste(WORK_DIR, "PatientCorePopulatedTable.txt", sep=""),
                     sep = "\t", header = TRUE, fileEncoding="UTF-8-BOM")
names(d.demo)

#Create a cohort that contains all unique combinations of patient id, admission id, addmision start and end dates, and additional demographic details.
d.cohort <- merge(d.enc_info, d.mp_outcome, by = c("PatientID_AND_AdmissionID"))
d.cohort <- merge(d.cohort, d.demo, by = c("PatientID"))

names(d.cohort)
summary(d.cohort)

#Convert dates and times to a more readable format using the command "mutate".
before_mutate = head(d.cohort)
Sys.setenv(tz = "America/Chicago")
d.cohort <- d.cohort %>%
  mutate(AdmissionStartDate = ymd_hms(AdmissionStartDate),
         AdmissionEndDate = ymd_hms(AdmissionEndDate),
         PatientDateOfBirth = ymd_hms(PatientDateOfBirth))

after_mutate = head(d.cohort)

#Calculating age and adding a new column to the data frame.
d.cohort <- d.cohort %>%
  mutate(PatientAge = interval(PatientDateOfBirth, AdmissionStartDate) / dyears(1))

#Calculating LOS and addting a new column to the data frame.
d.cohort <- d.cohort %>%
  mutate(LOS = interval(AdmissionStartDate, AdmissionEndDate) / ddays(1))

#In how many of the admissions the patients were 65 years or older?
nrow(d.cohort[d.cohort$PatientAge >= 65, ])

#How races are distributed?
summary(d.cohort$PatientRace)

#Let's compare the populations (cancer vs. no cancer).
d.1 <- filter(d.cohort, outcome == 1) #4375 admissions.
d.0 <- filter(d.cohort, outcome == 0) #31768 admissions.

cat("Mean age, outcome 1: ", mean(d.1$PatientAge), "\n")
## Mean age, outcome 1: 42.06632

cat("Mean age, outcome 0: ", mean(d.0$PatientAge), "\n")
## Mean age, outcome 0: 41.70317

cat("SD age, outcome 1: ", sd(d.1$PatientAge), "\n")
## SD age, outcome 1: 18.17401

cat("SD age, outcome 1: ", sd(d.0$PatientAge), "\n")
## SD age, outcome 1: 18.04217

#Comparison of continuous variables
#Assess differences in age
print(t.test(d.1$PatientAge, d.0$PatientAge))

#Comparision of categorical variables
gender.table <- with(d.cohort, table(outcome, PatientGender))
gender.table

chisq.test(gender.table)

#Incorporating most recent lab values per admission.
#Note that the table "MostRecentLabsPerAdmission.txt" was created using a seperated program (written in SQL).
#The code for this program is available in the courses' github site.
d.most_recent_labs <- read.table(paste(WORK_DIR, "MostRecentLabsPerAdmission.txt", sep=""),
                     sep = "\t", header = TRUE, fileEncoding="UTF-8-BOM")

d.most_recent_labs$PatientID_AND_AdmissionID <- paste(d.most_recent_labs$PatientID, d.most_recent_labs$AdmissionID, sep = "_")

#Adding the most recent labs for each admission in the main cohort.
d.cohort <- merge(d.cohort, d.most_recent_labs, by = c("PatientID_AND_AdmissionID"))

names(d.cohort)
summary(d.cohort)

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


