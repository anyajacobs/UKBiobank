library(tidyverse)
library(MatchIt)

# Load data if needed
mildmatching_full <- readRDS("data_work/Matching/mildmatching_full.rds")
severematching_full <- readRDS("data_work/Matching/severematching_full.rds")
mildmatching_notownsend <- readRDS("data_work/Matching/mildmatching_notownsend.rds")
severematching_notownsend <- readRDS("data_work/Matching/severematching_notownsend.rds")
# obtain matched data:
mfmatcheddata <- match.data(mildmatching_full) # mf = mild (asthma) fully (matched) cohort with all 5 covariates (matching covariates: age, sex, ethnicity, townsend, centre)
sfmatcheddata <- match.data(severematching_full) # sf = severe (asthma) fully (matched) cohort with all 5 covariates
mrmatcheddata <- match.data(mildmatching_notownsend) # mr = mild (asthma) reduced (matched) cohort with 4 covariates (matching covariates: age, sex, ethnicity, centre)
srmatcheddata <- match.data(severematching_notownsend) # sr = severe (asthma) reduced (matched) cohort with 4 covariates


# ANALYSES ####################################################

## merge hospital data with matched populations ####
# create LONG versions of matched population datasets: #

#subset entire matched population data for variables of interest
# starting with MILD FULLY MATCHED population
mfmatcheddata <- mfmatcheddata %>% 
  dplyr::select(eid,sex,age.group,age.recruit,recruit.date,death.date,
                cancer.baseline.all,MI.baseline,stroke.baseline,
                diabetes.fo.pre,CKD,hypertension,PVD,smoke,ethnicity,
                BMI_cat,COPD.baseline,sleep_apnoea.baseline,mental,
                townsend,qol_score,
                England,Wales,Scotland,centre,dsource,censordate,
                fu_duration,fu_year,year_1,year_2,year_3,year_4,year_5,
                year_6,year_7,year_8,year_9,year_10,year_11,year_12,year_13,
                year_14,year_15,year_16,
                diagnosis_all_med,asthma_mild_ONLY,diagnosis_severe_med) %>% group_by(eid)

nrow(mfmatcheddata) # 119,754 individuals in total (19,959 with mild asthma + 99,795 controls)
# this is the WIDE form of the dataset.

### ~ transform wide to long format #####
mf_long <- mfmatcheddata %>% gather(paste0("year_", 1:16), key=fu_year, value = fu_completion) %>% arrange(eid)

colnames(mf_long)
head(mf_long[,25:31], 20)


## ~ merge with hospital data ###### 
# admissions, hospital days, costs

mf_longhosp <- merge(mf_long, finalhospdata_byICD, by=c("eid", "fu_year"), all.x = T) %>% arrange(eid, fu_year)
head(mf_longhosp, 20)

mf_longhosp <- mf_longhosp %>% arrange(eid, fu_year)

# make NA's 0
mf_longhosp$admi_n[is.na(mf_longhosp$admi_n)] <- 0 
mf_longhosp$hdays_byfuyear[is.na(mf_longhosp$hdays_byfuyear)] <- 0 
mf_longhosp$totcost_byfuyear[is.na(mf_longhosp$totcost_byfuyear)] <- 0 
mf_longhosp$respadmi_byfuyear[is.na(mf_longhosp$respadmi_byfuyear)] <- 0
# so we still have the OG NA's in resp_admi and asthma_admi (and in ICD_diag01) ~ [these aren't YEARLY though]
# but the yearly resp_admi's NA's have been turned to 0, so note they are non-distinguishable from years with no hospital visits at all

## ~ subset to remove years where completion is 0 ######
mf_longhosp <- subset(mf_longhosp, fu_completion>0)

### ~ remove duplicate fu_years - make eid AND fu_year distinct #####
mf_longhosp <- mf_longhosp %>% distinct(eid, fu_year, .keep_all = T)


head(mf_longhosp, 50)

head(mf_longhosp[,c("eid", "diagnosis_all_med", "asthma_mild_ONLY", "diagnosis_severe_med", 
                    "fu_year", "admi_n", "hdays_byfuyear", "totcost_byfuyear", 
                    "ICD_diag01", "respadmi_byfuyear")],100)

# ^ we now have our final dataset (for the MILD FULLY MATCHED population) in long format.

# analysis:
table(mf_longhosp$admi_n) # 1,054,095 years with no hospital admission




#  !! LOAD DATASETS HERE to skip above steps: 
mf_longhosp <- readRDS(file.path(work_data, "mf_longhosp.rds"))
sf_longhosp <- readRDS(file.path(work_data, "sf_longhosp.rds"))
mr_longhosp <- readRDS(file.path(work_data, "mr_longhosp.rds"))
sr_longhosp <- readRDS(file.path(work_data, "sr_longhosp.rds"))



## create binary "died" variable ####
## determine no. of deaths in each follow-up year ##

# make new column of binary death (1 or 0) for yes or no
mf_longhosp$died <- mf_longhosp$death.date
mf_longhosp$died <- as.integer(mf_longhosp$died)
mf_longhosp$died[is.na(mf_longhosp$died)] <- 0
mf_longhosp$died[mf_longhosp$died>0] <- 1
# note this will mean for every person who died, every one of their F-U years
# will have '1' under died column


## count the no. of people who died or were censored in each year of followup:


#y1:
first <- mf_longhosp %>% filter(fu_duration<=365.25) # (we're only looking at eid's who lasted up to 1 year/ had data collection for only 1 yr)
table(first$died) # 719 death years (1's), 32 non-death years (0's) - thus these 32 are lost due to censoring
# we create the column "yr_censored" to reflect this.
# repeat for each year.
first$yr_censored <- 0
first$yr_censored[first$died==0] <- 1

n_distinct(first$eid) # note, only 243 distinct eid's (the 1's are repeated for each person who died)
first <- first %>% distinct(eid, .keep_all = T) # we apply 'distinct eid' from now on
table(first$died) # unique deaths = 218. unique censored = 25.


#y2:
second <- mf_longhosp %>% filter(fu_duration>365.25 & fu_duration<=730.5) %>% distinct(eid, .keep_all = T)
table(second$died) # 326 deaths, 18 censored. 
second$yr_censored[second$died==0] <- 1
second$yr_censored[is.na(second$yr_censored)] <- 0

#y3:
third <- mf_longhosp %>% filter(fu_duration>730.5 & fu_duration<=1095.75) %>% distinct(eid, .keep_all = T)
table(third$died) # 424 died in their 3rd year of follow-up. 
# 39 lost due to reasons other than death (i.e. this is their last yr due to losstofollowup or censored) 
third$yr_censored[third$died==0] <- 1
third$yr_censored[is.na(third$yr_censored)] <- 0

# ...

# y15:
fifteen <- distincteid %>% filter(fu_duration>5113.5 & fu_duration<=5478.75) 
table(fifteen$deathyesorno) # total 11 participants left. all 11 died in this yr (15) of follow-up.


# summarise longest year of followup for each id:

mf_longhospDISTINCT <- mf_longhosp %>% distinct(eid, .keep_all = T)

table(mf_longhosp$fu_year, useNA = "ifany") 
# 146,145 people (85% of our study population) have 10-years follow-up 
#' ^ (no. this is not unique ID)


ten <- mf_longhosp %>% filter(fu_duration>3652.5 & fu_duration<=4017.75) %>% distinct(eid, .keep_all = T)
table(ten$died) # 940 died in their 10th year of follow-up, 35,602 are censored in this year.
# 39 lost due to reasons other than death (i.e. this is their last yr due to losstofollowup or censored) 
ten$yr_censored[ten$died==0] <- 1
ten$yr_censored[is.na(ten$yr_censored)] <- 0





### ~ calculate mean duration of follow-up #####

## First, subset the separate populations
# subset the mild asthma population
submf_mildasthma <- mf_longhosp[mf_longhosp$diagnosis_all_med == 1, ]
# subset for the control population
submf_mildcontrols <- mf_longhosp[mf_longhosp$diagnosis_all_med == 0, ] 

## Calculate mean duration of follow-up for each group
summary(submf_mildasthma$fu_duration/365.25) # 10.72 years
summary(submf_mildcontrols$fu_duration/365.25) # 10.73 years
sd(submf_mildasthma$fu_duration/365.25) # (1.36)
sd(submf_mildcontrols$fu_duration/365.25) # (1.34)

table(submf_mildasthma$fu_year)
table(submf_mildcontrols$fu_year)



## ~ exclude fu_years>10 (& further filtering by fu_completion #####

## Exclude years of follow-up beyond 10 and use this for all following analyses.

mf_longhosp10y <- mf_longhosp %>% filter(fu_year != "year_11")
mf_longhosp10y <- mf_longhosp10y %>% filter(fu_year != "year_12")
mf_longhosp10y <- mf_longhosp10y %>% filter(fu_year != "year_13")
mf_longhosp10y <- mf_longhosp10y %>% filter(fu_year != "year_14")
mf_longhosp10y <- mf_longhosp10y %>% filter(fu_year != "year_15")
mf_longhosp10y <- mf_longhosp10y %>% filter(fu_year != "year_16")

nrow(mf_longhosp)
nrow(mf_longhosp10y)
n_distinct(mf_longhosp10y$eid)

# remove partial years of followup for those who are alive [i.e. the year they were censored] &
# keep all years for those that died, including partial years [i.e. year of death]
# (years of 0 completion were already removed earlier).
mf_longhosp10y <- mf_longhosp10y %>% filter((fu_completion==1 & died==0) | (fu_completion<=1 & died==1))

summary(mf_longhosp10y$fu_completion) 

head(mf_longhosp10y)
# ^ this is the final version of dataset to be used for modelling.



### Annual hospital data (descriptive) #######################################################


# Re-subset with mf_longhosp10y:
# subset the mild asthma population
submf_mildasthma10 <- mf_longhosp10y[mf_longhosp10y$asthma_mild_ONLY== 1, ]
# subset for the control population
submf_mildcontrols10 <- mf_longhosp10y[mf_longhosp10y$asthma_mild_ONLY == 0, ] 

library(tableone)
cont_variables <- c("admi_n", "hdays_byfuyear", "totcost_byfuyear", "respadmi_byfuyear")
cat_variables <- c("resp_admi", "asthma_admi", "ICD_diag01")


# Calculate annual hospital data with new subsets (<= 10-yr followup)
table1_mf_asthma <- CreateContTable(vars = cont_variables, data = submf_mildasthma10)
t1mf_asthmaicd <- CreateCatTable(vars = cat_variables, data = submf_mildasthma10)

table1_mf_controls <- CreateContTable(vars = cont_variables, data = submf_mildcontrols10)
t1mf_controlsicd <- CreateCatTable(vars = cat_variables, data = submf_mildcontrols10)

# for mild asthmatics...
table1_mf_asthma
t1mf_asthmaicd # 6.4% of all admissions are respiratory, 1.3% of all admissions are due to asthma
table(submf_mildasthma10$ICD_diag01, useNA = "ifany") 

# for non-asthmatic controls...
table1_mf_controls
t1mf_controlsicd # 2.7% of all admissions are respiratory, 0.1% are asthma


# confirm significance:
library(MASS)
# compare between mild asthma + non-asthma controls:
t.test(mf_longhosp10y$asthma_mild_ONLY, mf_longhosp10y$admi_n)
t.test(mf_longhosp10y$asthma_mild_ONLY, mf_longhosp10y$hdays_byfuyear)
t.test(mf_longhosp10y$asthma_mild_ONLY, mf_longhosp10y$totcost_byfuyear)
t.test(mf_longhosp10y$asthma_mild_ONLY, mf_longhosp10y$resp_admi)
t.test(mf_longhosp10y$asthma_mild_ONLY, mf_longhosp10y$asthma_admi)
# ALL are highly significant (p<0.000000000000001)

pairwise.t.test(mf_longhosp10y$asthma_mild_ONLY, mf_longhosp10y$ICD_diag01, p.adjust.method = "none")

## Hospital outcomes BY FU year:

## calculate mean costs by FU year:  [below was repeated for no. admissions (n_admi) + days (hdays_byfuyear)]
# mean costs by FU year for asthma cohort (mild):
aggregate(submf_mildasthma10$totcost_byfuyear,list(submf_mildasthma10$fu_year),mean)
aggregate(submf_mildasthma10$totcost_byfuyear,list(submf_mildasthma10$fu_year),sd)
# mean costs by FU year for control cohort:
aggregate(submf_mildcontrols10$totcost_byfuyear,list(submf_mildcontrols10$fu_year),mean)
aggregate(submf_mildcontrols10$totcost_byfuyear,list(submf_mildcontrols10$fu_year),sd)


## calculate mean hospital outcomes by ICD10 admission chapter: [below was repeated for no. admissions + days]

# mean costs by FU year for asthma cohort (mild):
aggregate(submf_mildasthma10$totcost_byfuyear,list(submf_mildasthma10$ICD_diag01),mean)
aggregate(submf_mildasthma10$totcost_byfuyear,list(submf_mildasthma10$ICD_diag01),sd)
# mean costs by FU year for control cohort:
aggregate(submf_mildcontrols10$totcost_byfuyear,list(submf_mildcontrols10$ICD_diag01),mean)
aggregate(submf_mildcontrols10$totcost_byfuyear,list(submf_mildcontrols10$ICD_diag01),sd)


## proportion of years of followup with hospital costs [example below = for the asthma [mild] cohort in mf dataset]
propfcosts <- submf_mildasthma10 %>% distinct(eid, fu_year, .keep_all = T)

# make new binary column (1 or 0) for yes or no for costs
propfcosts$costyesorno <- propfcosts$totcost_byfuyear
propfcosts$costyesorno[is.na(propfcosts$costyesorno)] <- 0
propfcosts$costyesorno[propfcosts$costyesorno>0] <- 1

#
propfcosts$costyesorno 
aggregate(propfcosts$costyesorno,list(propfcosts$fu_year),mean)
head(propfcosts[,c("eid", "fu_year", "totcost_byfuyear", "costyesorno")], 100)

table(propfcosts$costyesorno) # 23.8% (44,948 FU years have costs, out of 188,619 FU years) - mild asthma
# for the controls, 18.7% (176,401 FU years have costs, out of 944,469 years) 
nrow(propfcosts)



# REPEATS FOR OTHER DATASETS # -------------------------------------------------

### ~~~~ FULL SEVERE MATCHING ~~~~ ###########

#subset entire matched population data for variables of interest
sfmatcheddata <- sfmatcheddata %>% 
  dplyr::select(eid,sex,age.group,age.recruit,recruit.date,death.date,
                cancer.baseline.all,MI.baseline,stroke.baseline,
                diabetes.fo.pre,CKD,hypertension,PVD,smoke,ethnicity,
                BMI_cat,COPD.baseline,sleep_apnoea.baseline,mental,
                townsend,qol_score,
                England,Wales,Scotland,centre,dsource,censordate,
                fu_duration,fu_year,year_1,year_2,year_3,year_4,year_5,
                year_6,year_7,year_8,year_9,year_10,year_11,year_12,year_13,
                year_14,year_15,year_16,
                diagnosis_all_med,asthma_mild_ONLY,diagnosis_severe_med) %>% 
  group_by(eid)

nrow(sfmatcheddata) # 30,432 individuals in total (5,072 with mod-severe asthma + 25,360 controls)

########## ~ transform wide to long format ##
sf_long <- sfmatcheddata %>% gather(paste0("year_", 1:16), key=fu_year, value = fu_completion) %>% arrange(eid)

colnames(sf_long)
head(sf_long[,25:31], 20)


########## ~ merge with hospital data ##
# admissions, hospital days, costs

sf_longhosp <- merge(sf_long, finalhospdata_byICD, by=c("eid", "fu_year"), all.x = T) %>% arrange(eid, fu_year)
head(sf_longhosp, 20)

sf_longhosp <- sf_longhosp %>% arrange(eid, fu_year)

# make NA's 0
sf_longhosp$admi_n[is.na(sf_longhosp$admi_n)] <- 0 
sf_longhosp$hdays_byfuyear[is.na(sf_longhosp$hdays_byfuyear)] <- 0 
sf_longhosp$totcost_byfuyear[is.na(sf_longhosp$totcost_byfuyear)] <- 0 
sf_longhosp$respadmi_byfuyear[is.na(sf_longhosp$respadmi_byfuyear)] <- 0

########## ~ subset to remove years where completion is 0 ##
sf_longhosp <- subset(sf_longhosp, fu_completion>0)


############# ~ remove duplicate fu_years - make eid AND fu_year distinct #
sf_longhosp <- sf_longhosp %>% distinct(eid, fu_year, .keep_all = T)

head(sf_longhosp, 50)

head(sf_longhosp[,c("eid", "fu_year", "diagnosis_all_med", "asthma_mild_ONLY", "diagnosis_severe_med", 
                    "admi_n", "hdays_byfuyear", "totcost_byfuyear", "respadmi_byfuyear")],100)

# ^ we now have our final dataset (for the SEVERE FULLY MATCHED population) in long format.

## create binary "died" variable ##

# make new column of binary death (1 or 0) for yes or no
sf_longhosp$died <- sf_longhosp$death.date
sf_longhosp$died <- as.integer(sf_longhosp$died)
sf_longhosp$died[is.na(sf_longhosp$died)] <- 0
sf_longhosp$died[sf_longhosp$died>0] <- 1
# note this will mean for every person who died, every one of their F-U years
# will have '1' under died column


## count the no. of people who died or were censored in each year of followup:


#y1:
first <- sf_longhosp %>% filter(fu_duration<=365.25) # (we're only looking at eid's who lasted up to 1 year/ had data collection for only 1 yr)
table(first$died) # 719 death years (1's), 32 non-death years (0's) - thus these 32 are lost due to censoring
# we create the column "yr_censored" to reflect this.
# repeat for each year.
first$yr_censored <- 0
first$yr_censored[first$died==0] <- 1

n_distinct(first$eid) # note, only 243 distinct eid's (the 1's are repeated for each person who died)
first <- first %>% distinct(eid, .keep_all = T) # we apply 'distinct eid' from now on
table(first$died) # unique deaths = 218. unique censored = 25.


#y2:
second <- sf_longhosp %>% filter(fu_duration>365.25 & fu_duration<=730.5) %>% distinct(eid, .keep_all = T)
table(second$died) # 326 deaths, 18 censored. 
second$yr_censored[second$died==0] <- 1
second$yr_censored[is.na(second$yr_censored)] <- 0

#y3:
third <- sf_longhosp %>% filter(fu_duration>730.5 & fu_duration<=1095.75) %>% distinct(eid, .keep_all = T)
table(third$died) # 424 died in their 3rd year of follow-up. 
# 39 lost due to reasons other than death (i.e. this is their last yr due to losstofollowup or censored) 
third$yr_censored[third$died==0] <- 1
third$yr_censored[is.na(third$yr_censored)] <- 0

# ...

# y15:
fifteen <- distincteid %>% filter(fu_duration>5113.5 & fu_duration<=5478.75) 
table(fifteen$deathyesorno) # total 11 participants left. all 11 died in this yr (15) of follow-up.


# summarise longest year of followup for each id:

sf_longhospDISTINCT <- sf_longhosp %>% distinct(eid, .keep_all = T)

table(sf_longhosp$fu_year, useNA = "ifany") 
# 146,145 people (85% of our study population) have 10-years follow-up 
#' ^ (no. this is not unique ID)


ten <- sf_longhosp %>% filter(fu_duration>3652.5 & fu_duration<=4017.75) %>% distinct(eid, .keep_all = T)
table(ten$died) # 940 died in their 10th year of follow-up, 35,602 are censored in this year.
# 39 lost due to reasons other than death (i.e. this is their last yr due to losstofollowup or censored) 
ten$yr_censored[ten$died==0] <- 1
ten$yr_censored[is.na(ten$yr_censored)] <- 0





## ~ calculate mean duration of follow-up ##

# subset mod-severe asthma + control populations
subSF_asthma <- sf_longhosp[sf_longhosp$diagnosis_severe_med == 1, ]
subSF_controls <- sf_longhosp[sf_longhosp$diagnosis_severe_med == 0, ]

## Calculate mean duration of follow-up for each group
summary(subSF_asthma$fu_duration/365.25) # 10.58 years
summary(subSF_controls$fu_duration/365.25) # 10.73 years
sd(subSF_asthma$fu_duration/365.25) # (1.58)
sd(subSF_controls$fu_duration/365.25) # (1.38)


## ~ exclude fu_years>10 (& further filtering by fu_completion ##

## Exclude years of follow-up beyond 10 and use this for all following analyses.

sf_longhosp10y <- sf_longhosp %>% filter(fu_year != "year_11")
sf_longhosp10y <- sf_longhosp10y %>% filter(fu_year != "year_12")
sf_longhosp10y <- sf_longhosp10y %>% filter(fu_year != "year_13")
sf_longhosp10y <- sf_longhosp10y %>% filter(fu_year != "year_14")
sf_longhosp10y <- sf_longhosp10y %>% filter(fu_year != "year_15")
sf_longhosp10y <- sf_longhosp10y %>% filter(fu_year != "year_16")

nrow(sf_longhosp)
nrow(sf_longhosp10y) 
n_distinct(sf_longhosp10y$eid)

# remove partial years of followup for those who are alive [i.e. the year they were censored] &
# keep all years for those that died, including partial years [i.e. year of death]
# (years of 0 completion were already removed earlier).
sf_longhosp10y <- sf_longhosp10y %>% filter((fu_completion==1 & died==0) | (fu_completion<=1 & died==1))

summary(sf_longhosp10y$fu_completion) 

head(sf_longhosp10y)



### Annual hospital data (descriptive) for SF dataset ##

# Re-subset with sf_longhosp10y:
# subset the mild asthma population
subSF_asthma10 <- sf_longhosp10y[sf_longhosp10y$diagnosis_severe_med== 1, ]
# subset for the control population
subSF_controls10 <- sf_longhosp10y[sf_longhosp10y$diagnosis_severe_med == 0, ] 

library(tableone)
cont_variables <- c("admi_n", "hdays_byfuyear", "totcost_byfuyear")
cat_variables <- c("resp_admi", "asthma_admi", "ICD_diag01")


# Re-calculate annual hospital data with new subsets (<= 10-yr followup)
table1_sf_asthma <- CreateContTable(vars = cont_variables, data = subSF_asthma10)
table1_sf_asthma
table1_sf_controls <- CreateContTable(vars = cont_variables, data = subSF_controls10)
table1_sf_controls

## admissions by ICD chapter
table1_sf_asthmaICD <- CreateCatTable(vars = cat_variables, data = subSF_asthma10)
table1_sf_controlsICD <- CreateCatTable(vars = cat_variables, data = subSF_controls10)
# for severe asthmatics...
table1_sf_asthmaICD # 12.5% of all admissions are respiratory etc
table(subSF_asthma10$asthma_admi, useNA = "ifany") # 2.1% of all admissions are due to asthma

# for non-asthmatic controls...
table1_sf_controlsICD # 2.8% of all admissions are respiratory, 0.1% asthma


# confirm significance:
# compare between severe asthma + non-asthma controls:
t.test(sf_longhosp10y$diagnosis_severe_med, sf_longhosp10y$admi_n)
t.test(sf_longhosp10y$diagnosis_severe_med, sf_longhosp10y$hdays_byfuyear)
t.test(sf_longhosp10y$diagnosis_severe_med, sf_longhosp10y$totcost_byfuyear)
t.test(sf_longhosp10y$diagnosis_severe_med, sf_longhosp10y$resp_admi)
t.test(sf_longhosp10y$diagnosis_severe_med, sf_longhosp10y$asthma_admi)
# ALL are highly significant (p<2.2e-16)


## Costs by FU year:

## calculate mean costs/ admissions/ days by FU year:  

# mean costs by FU year for asthma cohort (severe): [below was repeated for no. admissions + days]
aggregate(subSF_asthma10$totcost_byfuyear,list(subSF_asthma10$fu_year),mean)
aggregate(subSF_asthma10$totcost_byfuyear,list(subSF_asthma10$fu_year),sd)

# mean costs by FU year for control cohort:
aggregate(subSF_controls10$totcost_byfuyear,list(subSF_controls10$fu_year),mean)
aggregate(subSF_controls10$totcost_byfuyear,list(subSF_controls10$fu_year),sd)


## calculate mean hospital outcomes by ICD10 admission chapter:

# mean costs by FU year for asthma cohort (mild): [below was repeated for no. admissions + days]
aggregate(subSF_asthma10$totcost_byfuyear,list(subSF_asthma10$ICD_diag01),mean)
aggregate(subSF_asthma10$totcost_byfuyear,list(subSF_asthma10$ICD_diag01),sd)
# mean costs by FU year for control cohort:
aggregate(subSF_controls10$totcost_byfuyear,list(subSF_controls10$ICD_diag01),mean)
aggregate(subSF_controls10$totcost_byfuyear,list(subSF_controls10$ICD_diag01),sd)



## proportion of years of followup with hospital costs [example below = for the asthma [severe] cohort in sf dataset]

# make new binary column (1 or 0) for yes or no for costs
subSF_asthma10$costyesorno <- subSF_asthma10$totcost_byfuyear
subSF_asthma10$costyesorno[is.na(subSF_asthma10$costyesorno)] <- 0
subSF_asthma10$costyesorno[subSF_asthma10$costyesorno>0] <- 1

table(subSF_asthma10$costyesorno) # 32.5%. (15,257 FU years have costs, out of 46,980 years)
table(subSF_controls10$costyesorno) # 20.7% (49,671 years have costs, out of 239,500)

summary(as.numeric(subSF_asthma10$fu_year))

head(sf_longhosp10y)




### ~~~~ REDUCED MILD MATCHING ~~~~ (no townsend) ###########

#subset entire matched population data for variables of interest

mrmatcheddata <- mrmatcheddata %>% 
  dplyr::select(eid,sex,age.group,age.recruit,recruit.date,death.date,
                cancer.baseline.all,MI.baseline,stroke.baseline,
                diabetes.fo.pre,CKD,hypertension,PVD,smoke,ethnicity,
                BMI_cat,COPD.baseline,sleep_apnoea.baseline,mental,
                townsend,qol_score,
                England,Wales,Scotland,centre,dsource,censordate,
                fu_duration,fu_year,year_1,year_2,year_3,year_4,year_5,
                year_6,year_7,year_8,year_9,year_10,year_11,year_12,year_13,
                year_14,year_15,year_16,
                diagnosis_all_med,asthma_mild_ONLY,diagnosis_severe_med) %>% 
  group_by(eid)

nrow(mrmatcheddata) # 119,754 individuals in total

########## ~ transform wide to long format ##
mr_long <- mrmatcheddata %>% gather(paste0("year_", 1:16), key=fu_year, value = fu_completion) %>% arrange(eid)

########## ~ merge with hospital data ##
# admissions, hospital days, costs
mr_longhosp <- merge(mr_long, finalhospdata_byICD, by=c("eid", "fu_year"), all.x = T) %>% arrange(eid, fu_year)
head(mr_longhosp, 20)

mr_longhosp <- mr_longhosp %>% arrange(eid, fu_year)

# make NA's 0
mr_longhosp$admi_n[is.na(mr_longhosp$admi_n)] <- 0 
mr_longhosp$hdays_byfuyear[is.na(mr_longhosp$hdays_byfuyear)] <- 0 
mr_longhosp$totcost_byfuyear[is.na(mr_longhosp$totcost_byfuyear)] <- 0 
mr_longhosp$respadmi_byfuyear[is.na(mr_longhosp$respadmi_byfuyear)] <- 0
mr_longhosp$ICD1_admibyfuyear[is.na(mr_longhosp$ICD1_admibyfuyear)] <- 0
mr_longhosp$ICD2_admibyfuyear[is.na(mr_longhosp$ICD2_admibyfuyear)] <- 0
mr_longhosp$ICD3_admibyfuyear[is.na(mr_longhosp$ICD3_admibyfuyear)] <- 0
mr_longhosp$ICD4_admibyfuyear[is.na(mr_longhosp$ICD4_admibyfuyear)] <- 0
mr_longhosp$ICD5_admibyfuyear[is.na(mr_longhosp$ICD5_admibyfuyear)] <- 0
mr_longhosp$ICD6_admibyfuyear[is.na(mr_longhosp$ICD6_admibyfuyear)] <- 0
mr_longhosp$ICD7_admibyfuyear[is.na(mr_longhosp$ICD7_admibyfuyear)] <- 0
mr_longhosp$ICD8_admibyfuyear[is.na(mr_longhosp$ICD8_admibyfuyear)] <- 0
mr_longhosp$ICD9_admibyfuyear[is.na(mr_longhosp$ICD9_admibyfuyear)] <- 0
mr_longhosp$ICD10_admibyfuyear[is.na(mr_longhosp$ICD10_admibyfuyear)] <- 0
mr_longhosp$ICD11_admibyfuyear[is.na(mr_longhosp$ICD11_admibyfuyear)] <- 0
mr_longhosp$ICD12_admibyfuyear[is.na(mr_longhosp$ICD12_admibyfuyear)] <- 0
mr_longhosp$ICD13_admibyfuyear[is.na(mr_longhosp$ICD13_admibyfuyear)] <- 0
mr_longhosp$ICD14_admibyfuyear[is.na(mr_longhosp$ICD14_admibyfuyear)] <- 0
mr_longhosp$ICD15_admibyfuyear[is.na(mr_longhosp$ICD15_admibyfuyear)] <- 0
mr_longhosp$ICD16_admibyfuyear[is.na(mr_longhosp$ICD16_admibyfuyear)] <- 0
mr_longhosp$ICD17_admibyfuyear[is.na(mr_longhosp$ICD17_admibyfuyear)] <- 0
mr_longhosp$ICD18_admibyfuyear[is.na(mr_longhosp$ICD18_admibyfuyear)] <- 0
mr_longhosp$ICD19_admibyfuyear[is.na(mr_longhosp$ICD19_admibyfuyear)] <- 0
mr_longhosp$ICD20_admibyfuyear[is.na(mr_longhosp$ICD20_admibyfuyear)] <- 0
mr_longhosp$ICD21_admibyfuyear[is.na(mr_longhosp$ICD21_admibyfuyear)] <- 0
mr_longhosp$ICD22_admibyfuyear[is.na(mr_longhosp$ICD22_admibyfuyear)] <- 0

########## ~ subset to remove years where completion is 0 ##
mr_longhosp <- subset(mr_longhosp, fu_completion>0)

############# ~ remove duplicate fu_years - make eid AND fu_year distinct #
mr_longhosp <- mr_longhosp %>% distinct(eid, fu_year, .keep_all = T)
head(mr_longhosp, 50)

head(mr_longhosp[,c("eid", "fu_year", "diagnosis_all_med", "asthma_mild_ONLY", "diagnosis_severe_med", 
                    "admi_n", "hdays_byfuyear", "totcost_byfuyear")],100)

# ^ final dataset for the MILD REDUCED MATCHED population in long format.


## create binary "died" variable ##
# make new column of binary death (1 or 0) for yes or no
mr_longhosp$died <- mr_longhosp$death.date
mr_longhosp$died <- as.integer(mr_longhosp$died)
mr_longhosp$died[is.na(mr_longhosp$died)] <- 0
mr_longhosp$died[mr_longhosp$died>0] <- 1
# note this will mean for every person who died, every one of their F-U years
# will have '1' under died column


## count the no. of people who died or were censored in each year of followup:


#y1:
first <- mr_longhosp %>% filter(fu_duration<=365.25) # (we're only looking at eid's who lasted up to 1 year/ had data collection for only 1 yr)
table(first$died) # 719 death years (1's), 32 non-death years (0's) - thus these 32 are lost due to censoring
# we create the column "yr_censored" to reflect this.
# repeat for each year.
first$yr_censored <- 0
first$yr_censored[first$died==0] <- 1

n_distinct(first$eid) # note, only 243 distinct eid's (the 1's are repeated for each person who died)
first <- first %>% distinct(eid, .keep_all = T) # we apply 'distinct eid' from now on
table(first$died) # unique deaths = 218. unique censored = 25.


#y2:
second <- mr_longhosp %>% filter(fu_duration>365.25 & fu_duration<=730.5) %>% distinct(eid, .keep_all = T)
table(second$died) # 326 deaths, 18 censored. 
second$yr_censored[second$died==0] <- 1
second$yr_censored[is.na(second$yr_censored)] <- 0

#y3:
third <- mr_longhosp %>% filter(fu_duration>730.5 & fu_duration<=1095.75) %>% distinct(eid, .keep_all = T)
table(third$died) # 424 died in their 3rd year of follow-up. 
# 39 lost due to reasons other than death (i.e. this is their last yr due to losstofollowup or censored) 
third$yr_censored[third$died==0] <- 1
third$yr_censored[is.na(third$yr_censored)] <- 0

# ...

# y15:
fifteen <- distincteid %>% filter(fu_duration>5113.5 & fu_duration<=5478.75) 
table(fifteen$deathyesorno) # total 11 participants left. all 11 died in this yr (15) of follow-up.


# summarise longest year of followup for each id:

mr_longhospDISTINCT <- mr_longhosp %>% distinct(eid, .keep_all = T)

table(mr_longhosp$fu_year, useNA = "ifany") 
# 146,145 people (85% of our study population) have 10-years follow-up 
#' ^ (no. this is not unique ID)

ten <- mr_longhosp %>% filter(fu_duration>3652.5 & fu_duration<=4017.75) %>% distinct(eid, .keep_all = T)
table(ten$died) # 940 died in their 10th year of follow-up, 35,602 are censored in this year.
# 39 lost due to reasons other than death (i.e. this is their last yr due to losstofollowup or censored) 
ten$yr_censored[ten$died==0] <- 1
ten$yr_censored[is.na(ten$yr_censored)] <- 0



## ~ calculate mean duration of follow-up ##

# subset mod-severe asthma + control populations
submr_asthma <- mr_longhosp[mr_longhosp$diagnosis_all_med == 1, ]
submr_controls <- mr_longhosp[mr_longhosp$diagnosis_all_med == 0, ]

## Calculate mean duration of follow-up for each group
summary(submr_asthma$fu_duration/365.25) # 10.72 years (THE SAME AS FULL MILD ASTHMA)
summary(submr_controls$fu_duration/365.25) # 10.73 years
sd(submr_asthma$fu_duration/365.25) # (1.36)
sd(submr_controls$fu_duration/365.25) # (1.33)


###### mr_longhosp10y ########

## ~ exclude fu_years>10 (& further filtering by fu_completion ##

## Exclude years of follow-up beyond 10 and use this for all following analyses.

mr_longhosp10y <- mr_longhosp %>% filter(fu_year != "year_11")
mr_longhosp10y <- mr_longhosp10y %>% filter(fu_year != "year_12")
mr_longhosp10y <- mr_longhosp10y %>% filter(fu_year != "year_13")
mr_longhosp10y <- mr_longhosp10y %>% filter(fu_year != "year_14")
mr_longhosp10y <- mr_longhosp10y %>% filter(fu_year != "year_15")
mr_longhosp10y <- mr_longhosp10y %>% filter(fu_year != "year_16")

# remove partial years of followup for those who are alive [i.e. the year they were censored] &
# keep all years for those that died, including partial years [i.e. year of death]
# (years of 0 completion were already removed earlier).
mr_longhosp10y <- mr_longhosp10y %>% filter((fu_completion==1 & died==0) | (fu_completion<=1 & died==1))


### Annual hospital data (descriptive) for mr dataset ##

## annual hospital data ("admi_n", "hdays_byfuyear", "totcost_byfuyear")
## overall AND by FU_year, and admissions by ICD chapter data, prop. years of FU...
## are exactly the same as for the full mild matching 

## HOWEVER, data for the controls have changed! (slightly different control cohort due to not matching for townsend)

install.packages("tableone")
library(tableone)
cont_variables <- c("admi_n", "hdays_byfuyear", "totcost_byfuyear", "respadmi_byfuyear")
cat_variables <- c("resp_admi", "asthma_admi", "ICD_diag01")

# subset the mild asthma population with mr_longhosp10y:
subMR_mildasthma10 <- mr_longhosp10y[mr_longhosp10y$asthma_mild_ONLY== 1, ]
# Subset control population:
subMR_controls10 <- mr_longhosp10y[mr_longhosp10y$asthma_mild_ONLY == 0, ] 

## Annual hospital data (<= 10-yr followup)
# # Asthma population
table1_mr_asthma <- CreateContTable(vars = cont_variables, data = subMR_mildasthma10)
table1_mr_asthma 
# admissions by ICD chapter
table1_mr_asthmaICD <- CreateCatTable(vars = cat_variables, data = subMR_mildasthma10)
table1_mr_asthmaICD # resp admi = 6.3%, asthma admi = 1.3%
# # Control population
table1_mr_controls <- CreateContTable(vars = cont_variables, data = subMR_controls10)
table1_mr_controls # mean annual data has slightly reduced compared to the fully mild matched controls
# admissions by ICD chapter
table1_mr_controlsICD <- CreateCatTable(vars = cat_variables, data = subMR_controls10)
table1_mr_controlsICD # resp admi = 2.6%, asthma admi = 0.1%


# Prop. of fu_years w/ hospital costs
# make new binary column (1 or 0) for yes or no for costs
subMR_controls10$costyesorno <- subMR_controls10$totcost_byfuyear
subMR_controls10$costyesorno[is.na(subMR_controls10$costyesorno)] <- 0
subMR_controls10$costyesorno[subMR_controls10$costyesorno>0] <- 1

table(subMR_controls10$costyesorno) # 18.4% (174,055 FU years have costs -
# slightly different to control cohort of FULL matched (18.7%)


## Costs by FU year:
## calculate mean costs by FU year:  [below was repeated for no. admissions + days]
# mean costs by FU year for control cohort:
aggregate(subMR_mildasthma10$totcost_byfuyear,list(subMR_mildasthma10$fu_year),mean)
aggregate(subMR_mildasthma10$totcost_byfuyear,list(subMR_mildasthma10$fu_year),sd)


## Mean hospital data by ICD-10 chapter:

# mean admissions by FU year for asthma cohort (mild):
# divide by 10 to get the mean ANNUAL data.
(aggregate(subMR_mildasthma10$admi_n,list(subMR_mildasthma10$ICD_diag01),mean))/10
(aggregate(subMR_mildasthma10$admi_n,list(subMR_mildasthma10$ICD_diag01),sd))/10
# mean admissions by FU year for control cohort:
(aggregate(subMR_controls10$admi_n,list(subMR_controls10$ICD_diag01),mean))/10
(aggregate(subMR_controls10$admi_n,list(subMR_controls10$ICD_diag01),sd))/10
# mean days in hospital by FU year for asthma cohort (mild):
(aggregate(subMR_mildasthma10$hdays_byfuyear,list(subMR_mildasthma10$ICD_diag01),mean))/10
(aggregate(subMR_mildasthma10$hdays_byfuyear,list(subMR_mildasthma10$ICD_diag01),sd))/10
# mean days in hospital by FU year for control cohort:
(aggregate(subMR_controls10$hdays_byfuyear,list(subMR_controls10$ICD_diag01),mean))/10
(aggregate(subMR_controls10$hdays_byfuyear,list(subMR_controls10$ICD_diag01),sd))/10
# mean costs by FU year for asthma cohort (mild):
(aggregate(subMR_mildasthma10$totcost_byfuyear,list(subMR_mildasthma10$ICD_diag01),mean))/10
(aggregate(subMR_mildasthma10$totcost_byfuyear,list(subMR_mildasthma10$ICD_diag01),sd))/10
# mean costs by FU year for control cohort:
(aggregate(subMR_controls10$totcost_byfuyear,list(subMR_controls10$ICD_diag01),mean))/10
(aggregate(subMR_controls10$totcost_byfuyear,list(subMR_controls10$ICD_diag01),sd))/10


## proportion of years of followup with hospital costs [example below = for the asthma [mild] cohort in mf dataset]
propfcosts <- subMR_mildasthma10 %>% distinct(eid, fu_year, .keep_all = T)

# make new binary column (1 or 0) for yes or no for costs
propfcosts$costyesorno <- propfcosts$totcost_byfuyear
propfcosts$costyesorno[is.na(propfcosts$costyesorno)] <- 0
propfcosts$costyesorno[propfcosts$costyesorno>0] <- 1

#
propfcosts$costyesorno 
aggregate(propfcosts$costyesorno,list(propfcosts$fu_year),mean)
head(propfcosts[,c("eid", "fu_year", "totcost_byfuyear", "costyesorno")], 100)

table(propfcosts$costyesorno) # 23.8% (44,948 FU years have costs, out of 188,619 FU years) - mild asthma
# for the controls, 18.7% (176,401 FU years have costs, out of 944,469 years) 
nrow(propfcosts)




### ~~~~ REDUCED SEVERE MATCHING ~~~~ (no townsend) ###########

#subset entire matched population data for variables of interest
srmatcheddata <- srmatcheddata %>% 
  dplyr::select(eid,sex,age.group,age.recruit,recruit.date,death.date,
                cancer.baseline.all,MI.baseline,stroke.baseline,
                diabetes.fo.pre,CKD,hypertension,PVD,smoke,ethnicity,
                BMI_cat,COPD.baseline,sleep_apnoea.baseline,mental,
                townsend,qol_score,
                England,Wales,Scotland,centre,dsource,censordate,
                fu_duration,fu_year,year_1,year_2,year_3,year_4,year_5,
                year_6,year_7,year_8,year_9,year_10,year_11,year_12,year_13,
                year_14,year_15,year_16,
                diagnosis_all_med,asthma_mild_ONLY,diagnosis_severe_med) %>% 
  group_by(eid)

nrow(srmatcheddata) # 30,432 individuals in total (5,072 with mod-severe asthma + 25,360 controls)

########## ~ transform wide to long format ##
sr_long <- srmatcheddata %>% gather(paste0("year_", 1:16), key=fu_year, value = fu_completion) %>% arrange(eid)

########## ~ merge with hospital data ##
# admissions, hospital days, costs

sr_longhosp <- merge(sr_long, finalhospdata_byICD, by=c("eid", "fu_year"), all.x = T) %>% arrange(eid, fu_year)

sr_longhosp <- sr_longhosp %>% arrange(eid, fu_year)

# make NA's 0
sr_longhosp$admi_n[is.na(sr_longhosp$admi_n)] <- 0 
sr_longhosp$hdays_byfuyear[is.na(sr_longhosp$hdays_byfuyear)] <- 0 
sr_longhosp$totcost_byfuyear[is.na(sr_longhosp$totcost_byfuyear)] <- 0 
sr_longhosp$respadmi_byfuyear[is.na(sr_longhosp$respadmi_byfuyear)] <- 0
sr_longhosp$ICD1_admibyfuyear[is.na(sr_longhosp$ICD1_admibyfuyear)] <- 0
sr_longhosp$ICD2_admibyfuyear[is.na(sr_longhosp$ICD2_admibyfuyear)] <- 0
sr_longhosp$ICD3_admibyfuyear[is.na(sr_longhosp$ICD3_admibyfuyear)] <- 0
sr_longhosp$ICD4_admibyfuyear[is.na(sr_longhosp$ICD4_admibyfuyear)] <- 0
sr_longhosp$ICD5_admibyfuyear[is.na(sr_longhosp$ICD5_admibyfuyear)] <- 0
sr_longhosp$ICD6_admibyfuyear[is.na(sr_longhosp$ICD6_admibyfuyear)] <- 0
sr_longhosp$ICD7_admibyfuyear[is.na(sr_longhosp$ICD7_admibyfuyear)] <- 0
sr_longhosp$ICD8_admibyfuyear[is.na(sr_longhosp$ICD8_admibyfuyear)] <- 0
sr_longhosp$ICD9_admibyfuyear[is.na(sr_longhosp$ICD9_admibyfuyear)] <- 0
sr_longhosp$ICD10_admibyfuyear[is.na(sr_longhosp$ICD10_admibyfuyear)] <- 0
sr_longhosp$ICD11_admibyfuyear[is.na(sr_longhosp$ICD11_admibyfuyear)] <- 0
sr_longhosp$ICD12_admibyfuyear[is.na(sr_longhosp$ICD12_admibyfuyear)] <- 0
sr_longhosp$ICD13_admibyfuyear[is.na(sr_longhosp$ICD13_admibyfuyear)] <- 0
sr_longhosp$ICD14_admibyfuyear[is.na(sr_longhosp$ICD14_admibyfuyear)] <- 0
sr_longhosp$ICD15_admibyfuyear[is.na(sr_longhosp$ICD15_admibyfuyear)] <- 0
sr_longhosp$ICD16_admibyfuyear[is.na(sr_longhosp$ICD16_admibyfuyear)] <- 0
sr_longhosp$ICD17_admibyfuyear[is.na(sr_longhosp$ICD17_admibyfuyear)] <- 0
sr_longhosp$ICD18_admibyfuyear[is.na(sr_longhosp$ICD18_admibyfuyear)] <- 0
sr_longhosp$ICD19_admibyfuyear[is.na(sr_longhosp$ICD19_admibyfuyear)] <- 0
sr_longhosp$ICD20_admibyfuyear[is.na(sr_longhosp$ICD20_admibyfuyear)] <- 0
sr_longhosp$ICD21_admibyfuyear[is.na(sr_longhosp$ICD21_admibyfuyear)] <- 0
sr_longhosp$ICD22_admibyfuyear[is.na(sr_longhosp$ICD22_admibyfuyear)] <- 0

########## ~ subset to remove years where completion is 0 ##
sr_longhosp <- subset(sr_longhosp, fu_completion>0)

############# ~ remove duplicate fu_years - make eid AND fu_year distinct #
sr_longhosp <- sr_longhosp %>% distinct(eid, fu_year, .keep_all = T)

head(sr_longhosp, 50)

head(sr_longhosp[,c("eid", "fu_year", "diagnosis_all_med", "asthma_mild_ONLY", "diagnosis_severe_med", 
                    "admi_n", "hdays_byfuyear", "totcost_byfuyear", "respadmi_byfuyear")],100)

# ^ we now have our final dataset for the SEVERE REDUCED MATCHED population in long format.



## create binary "died" variable ##

# make new column of binary death (1 or 0) for yes or no
sr_longhosp$died <- sr_longhosp$death.date
sr_longhosp$died <- as.integer(sr_longhosp$died)
sr_longhosp$died[is.na(sr_longhosp$died)] <- 0
sr_longhosp$died[sr_longhosp$died>0] <- 1
# note this will mean for every person who died, every one of their F-U years
# will have '1' under died column


## count the no. of people who died or were censored in each year of followup:

#y1:
first <- sr_longhosp %>% filter(fu_duration<=365.25) # (we're only looking at eid's who lasted up to 1 year/ had data collection for only 1 yr)
table(first$died) # 719 death years (1's), 32 non-death years (0's) - thus these 32 are lost due to censoring
# we create the column "yr_censored" to reflect this.
# repeat for each year.
first$yr_censored <- 0
first$yr_censored[first$died==0] <- 1

n_distinct(first$eid) # note, only 243 distinct eid's (the 1's are repeated for each person who died)
first <- first %>% distinct(eid, .keep_all = T) # we apply 'distinct eid' from now on
table(first$died) # unique deaths = 218. unique censored = 25.

#y2:
second <- sr_longhosp %>% filter(fu_duration>365.25 & fu_duration<=730.5) %>% distinct(eid, .keep_all = T)
table(second$died) # 326 deaths, 18 censored. 
second$yr_censored[second$died==0] <- 1
second$yr_censored[is.na(second$yr_censored)] <- 0

#y3:
third <- sr_longhosp %>% filter(fu_duration>730.5 & fu_duration<=1095.75) %>% distinct(eid, .keep_all = T)
table(third$died) # 424 died in their 3rd year of follow-up. 
# 39 lost due to reasons other than death (i.e. this is their last yr due to losstofollowup or censored) 
third$yr_censored[third$died==0] <- 1
third$yr_censored[is.na(third$yr_censored)] <- 0

# ...

# y15:
fifteen <- distincteid %>% filter(fu_duration>5113.5 & fu_duration<=5478.75) 
table(fifteen$deathyesorno) # total 11 participants left. all 11 died in this yr (15) of follow-up.


# summarise longest year of followup for each id:

sr_longhospDISTINCT <- sr_longhosp %>% distinct(eid, .keep_all = T)

table(sr_longhosp$fu_year, useNA = "ifany") 
# 146,145 people (85% of our study population) have 10-years follow-up 
#' ^ (no. this is not unique ID)


ten <- sr_longhosp %>% filter(fu_duration>3652.5 & fu_duration<=4017.75) %>% distinct(eid, .keep_all = T)
table(ten$died) # 940 died in their 10th year of follow-up, 35,602 are censored in this year.
# 39 lost due to reasons other than death (i.e. this is their last yr due to losstofollowup or censored) 
ten$yr_censored[ten$died==0] <- 1
ten$yr_censored[is.na(ten$yr_censored)] <- 0



## ~ calculate mean duration of follow-up ##

# subset mod-severe asthma + control populations
subsr_asthma <- sr_longhosp[sr_longhosp$diagnosis_severe_med == 1, ]
subsr_controls <- sr_longhosp[sr_longhosp$diagnosis_severe_med == 0, ]

## Calculate mean duration of follow-up for each group
summary(subsr_asthma$fu_duration/365.25) # 10.58 years
summary(subsr_controls$fu_duration/365.25) # 10.73 years
sd(subsr_asthma$fu_duration/365.25) # (1.58)
sd(subsr_controls$fu_duration/365.25) # (1.37)



###### sr_longhosp10y ########

## ~ exclude fu_years>10 (& further filtering by fu_completion ##

## Exclude years of follow-up beyond 10 and use this for all following analyses.

sr_longhosp10y <- sr_longhosp %>% filter(fu_year != "year_11")
sr_longhosp10y <- sr_longhosp10y %>% filter(fu_year != "year_12")
sr_longhosp10y <- sr_longhosp10y %>% filter(fu_year != "year_13")
sr_longhosp10y <- sr_longhosp10y %>% filter(fu_year != "year_14")
sr_longhosp10y <- sr_longhosp10y %>% filter(fu_year != "year_15")
sr_longhosp10y <- sr_longhosp10y %>% filter(fu_year != "year_16")

nrow(sr_longhosp)
nrow(sr_longhosp10y) 
n_distinct(sr_longhosp10y$eid)

# remove partial years of followup for those who are alive [i.e. the year they were censored] &
# keep all years for those that died, including partial years [i.e. year of death]
# (years of 0 completion were already removed earlier).
sr_longhosp10y <- sr_longhosp10y %>% filter((fu_completion==1 & died==0) | (fu_completion<=1 & died==1))

summary(sr_longhosp10y$fu_completion) 

head(sr_longhosp10y)



### Annual hospital data (descriptive) for SR dataset ##

## annual hospital data ("admi_n", "hdays_byfuyear", "totcost_byfuyear")
## overall AND by FU_year, and admissions by ICD chapter data, prop. years of FU...
## of the reduced SEVERE cohort, are exactly the same as for the FULL severe matching 

## admissions by ICD chapter
# for severe asthmatics...
table(subsr_asthma10$resp_admi, useNA = "ifany") # 12.5% of all admissions are respiratory etc etc...


## HOWEVER, data for the controls have changed! (slightly different control cohort due to not matching for townsend)

cont_variables <- c("admi_n", "hdays_byfuyear", "totcost_byfuyear")
cat_variables <- c("resp_admi", "asthma_admi", "ICD_diag01")

# Re-subset asthma population with sr_longhosp10y:
subSR_asthma10 <- sr_longhosp10y[sr_longhosp10y$diagnosis_severe_med == 1, ] 
# Re-subset control population with sr_longhosp10y:
subSR_controls10 <- sr_longhosp10y[sr_longhosp10y$diagnosis_severe_med == 0, ] 

# Re-calculate annual hospital data (<= 10-yr followup)
table1_sr_asthma <- CreateContTable(vars = cont_variables, data = subSR_asthma10)
table1_sr_asthma # mean annual data has slightly reduced compared to the fully severe matched controls
table1_sr_asthmaICD <- CreateCatTable(vars = cat_variables, data = subSR_asthma10)
table1_sr_asthmaICD # resp admi = 2.7%, asthma admi = 0.1%
# for controls:
table1_sr_controls <- CreateContTable(vars = cont_variables, data = subSR_controls10)
table1_sr_controls # mean annual data has slightly reduced compared to the fully severe matched controls
table1_sr_controlsICD <- CreateCatTable(vars = cat_variables, data = subSR_controls10)
table1_sr_controlsICD # resp admi = 2.7%, asthma admi = 0.1%


# prop. of fu_years w/ hospital costs
# make new binary column (1 or 0) for yes or no for costs
subSR_asthma10$costyesorno <- subSR_asthma10$totcost_byfuyear
subSR_asthma10$costyesorno[is.na(subSR_asthma10$costyesorno)] <- 0
subSR_asthma10$costyesorno[subSR_asthma10$costyesorno>0] <- 1
table(subSR_asthma10$costyesorno) # 32.5% (the same as FULL severe cohort)
subSR_controls10$costyesorno <- subSR_controls10$totcost_byfuyear
subSR_controls10$costyesorno[is.na(subSR_controls10$costyesorno)] <- 0
subSR_controls10$costyesorno[subSR_controls10$costyesorno>0] <- 1
table(subSR_controls10$costyesorno) # 20.3% (48,703 FU years have costs -
# slightly different to control cohort of FULL matched (20.9%)


## Costs by FU year:

## calculate mean costs by FU year:  [below was repeated for no. admissions + days]
# mean costs by FU year for asthma cohort:
aggregate(subSR_asthma10$totcost_byfuyear,list(subSR_asthma10$fu_year),mean)
aggregate(subSR_asthma10$totcost_byfuyear,list(subSR_asthma10$fu_year),sd)
# mean costs by FU year for control cohort:
aggregate(subSR_controls10$totcost_byfuyear,list(subSR_controls10$fu_year),mean)
aggregate(subSR_controls10$totcost_byfuyear,list(subSR_controls10$fu_year),sd)


## Mean hospital data by ICD-10 chapter:

# mean admissions by FU year for asthma cohort (severe):
# divide by 10 to get the mean ANNUAL data.
(aggregate(subSR_asthma10$admi_n,list(subSR_asthma10$ICD_diag01),mean))/10
(aggregate(subSR_asthma10$admi_n,list(subSR_asthma10$ICD_diag01),sd))/10
# mean admissions by FU year for control cohort:
(aggregate(subSR_controls10$admi_n,list(subSR_controls10$ICD_diag01),mean))/10
(aggregate(subSR_controls10$admi_n,list(subsR_controls10$ICD_diag01),sd))/10
# mean days in hospital by FU year for asthma cohort (severe):
(aggregate(subSR_asthma10$hdays_byfuyear,list(subSR_asthma10$ICD_diag01),mean))/10
(aggregate(subSR_asthma10$hdays_byfuyear,list(subSR_asthma10$ICD_diag01),sd))/10
# mean days in hospital by FU year for control cohort:
(aggregate(subSR_controls10$hdays_byfuyear,list(subSR_controls10$ICD_diag01),mean))/10
(aggregate(subSR_controls10$hdays_byfuyear,list(subSR_controls10$ICD_diag01),sd))/10
# mean costs by FU year for asthma cohort (severe):
(aggregate(subSR_asthma10$totcost_byfuyear,list(subSR_asthma10$ICD_diag01),mean))/10
(aggregate(subSR_asthma10$totcost_byfuyear,list(subSR_asthma10$ICD_diag01),sd))/10
# mean costs by FU year for control cohort:
(aggregate(subSR_controls10$totcost_byfuyear,list(subSR_controls10$ICD_diag01),mean))/10
(aggregate(subSR_controls10$totcost_byfuyear,list(subSR_controls10$ICD_diag01),sd))/10


# save ----
saveRDS(mr_longhosp, file = file.path(work_data, "mr_longhosp.rds"), compress = F)
mf_longhosp <- readRDS(file.path(work_data, "mf_longhosp.rds"))
sf_longhosp <- readRDS(file.path(work_data, "sf_longhosp.rds"))
mr_longhosp <- readRDS(file.path(work_data, "mr_longhosp.rds"))
sr_longhosp <- readRDS(file.path(work_data, "sr_longhosp.rds"))
# mf = mild (asthma) fully (matched) with all 5 covariates, plus hospital data
# sf = severe (asthma) fully (matched) with all 5 covariates, plus hospital data
# mr = mild (asthma) reduced (matched) with 4 covariates, plus hospital data
# sr = severe (asthma) reduced (matched) with 4 covariates, plus hospital data
# duplicate fu_years have been removed (using distinct(eid, fu_year)
# fu_years with 0 completion have been removed
# note, these datasets need to be converted to "mf_longhosp10y" etc before being used in models.
