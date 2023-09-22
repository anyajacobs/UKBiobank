library(tidyverse)

UKB_a <- readRDS(file.path(work_data, "UKB_a.rds.rds"))

# MATCHING ################################################

install.packages("MatchIt")
library(MatchIt)
?method_nearest

# Remove NA'S from matching covariates (matchit package does not allow missing values in the covariates)

summary(UKB_a) # There are only NAs in the ethnicity and townsend covariates
# sum(is.na(UKB_a$ethnicity)) = 2,776 ethnicity NA's in entire UKB cohort.
# sum(is.na(UKB_a$townsend)) = 623 townsend NA's (or 620 left once all ethnicity NA's have been dropped)
# table(UKB_a$diagnosis_all_med, UKB_a$ethnicity, useNA = "ifany") # 130 asthmatics have N/A ethnicity
# table(UKB_a$diagnosis_all_med, UKB_a$townsend, useNA = "ifany") # 33 asthmatics have N/A townsend
UKB_a_drop_ethnicitytownsendNAs <- UKB_a %>% drop_na(ethnicity)
UKB_a_drop_ethnicitytownsendNAs <- UKB_a_drop_ethnicitytownsendNAs %>% drop_na(townsend)

# check:
nrow(UKB_a_drop_ethnicitytownsendNAs) #rows have dropped to 499,108 (from 502,504)
table(UKB_a_drop_ethnicitytownsendNAs$diagnosis_all_med) 
#  FINAL asthma population for matching is 25,031 (down from 25,194)
#  final mild = 19,959. final mod-severe = 5,072 (after excluding townsend + ethnicity NA's)

library(Rcpp)
library("MatchIt")

# factor relevant variables (except for centre; use as.character)

UKB_a_drop_ethnicitytownsendNAs$centre <- as.character(UKB_a_drop_ethnicitytownsendNAs$centre)
UKB_a_drop_ethnicitytownsendNAs$townsend <- as.factor(UKB_a_drop_ethnicitytownsendNAs$townsend) 
UKB_a_drop_ethnicitytownsendNAs$sex <- as.factor(UKB_a_drop_ethnicitytownsendNAs$sex)
UKB_a_drop_ethnicitytownsendNAs$ethnicity <- as.factor(UKB_a_drop_ethnicitytownsendNAs$ethnicity)


# Match the asthma population separately by severity status
# 5:1 nearest neighbour matching without replacement. 
  
### Full matching (using all 5 OG covariates)  ####

## Full mild - match mild asthmatic population to controls:

# first remove ALL mod_severe asthma people from dataset, to prevent 
# matching of mild asthmatics to controls that are severe asthmatics (who happen to have similar age, sex etc)

subsetnosevere_UKB_a_NAsdropped <- UKB_a_drop_ethnicitytownsendNAs[UKB_a_drop_ethnicitytownsendNAs$diagnosis_severe_med == 0, ]

mildmatching_full <- matchit(asthma_mild_ONLY ~ age.recruit + sex + ethnicity + townsend + centre, 
                              data = subsetnosevere_UKB_a_NAsdropped, method = "nearest", distance = "mahalanobis", 
                              ratio = 5)

mildmatching_full # 5:1 nearest neighbour matching without replacement. target estimand: ATT
#saveRDS(mildmatching_full, file = file.path(work_data, "mildmatching_full.rds"))
mildmatching_full <- readRDS("data_work/Matching/mildmatching_full.rds")

summary(mildmatching_full) # balance assessment
# Focus on the 'percent balance improvement' table: 
# Note that variance ratios are not computed for binary covariates because they are a function
# of the prevalence in each group, which is captured in the mean difference and eCDF statistics

# A better way of assessing balance is to use the cobalt package:
install.packages("cobalt")
library(cobalt)
bal.tab(mildmatching_full, un = TRUE, binary = "std") 
# un = TRUE - means the balance statistics before matching are also displayed.
# Standard mean differences (SMD) after matching = Diff.Adj, before matching = Diff.un
# SMD's after matching are very close to zero = very well balanced.

# obtain matched data:
mfmatcheddata <- match.data(mildmatching_full)

##' IF FOLLOWING ERROR OCCURS: THEN ADD ',DATA = '
##' Error: A valid dataset could not be found. Please supply an argument to 'data' containing the original dataset used in the matching.
##' then use: mfmatcheddata <- match.data(mildmatching_full, data = subsetnosevere_UKB_a_NAsdropped)
##'   (i.e., supply the original data source after the matching object)

head(mfmatcheddata)
summary(mfmatcheddata)
# obtain matched data for control + asthma groups separately if needed:
TREATEDmildmatching_full <- match.data(mildmatching_full, group = "treat")
CONTROLmildmatching_full <- match.data(mildmatching_full, group = "control")
#summary(TREATEDmildmatching_full)
#summary(CONTROLmildmatching_full)
nrow(TREATEDmildmatching_full) # = 19,959
nrow(CONTROLmildmatching_full) # = 99,795


## Full severe matching ###

# Now subset to remove all eid's used in the previous matching (remove all mild's and their matched controls)

subsetsnomild_UKB_a_NAsdropped <-
  UKB_a_drop_ethnicitytownsendNAs[!UKB_a_drop_ethnicitytownsendNAs$eid 
                                             %in% mfmatcheddata$eid, ]

nrow(subsetsnomild_UKB_a_NAsdropped) # rows = 379,354

# COMPARE EID'S TO ENSURE NO OVERLAP:
overlap <- mfmatcheddata[mfmatcheddata$eid %in% subsetsnomild_UKB_a_NAsdropped$eid, ]
# ^ 0 rows = no overlap.

severematching_full <- matchit(diagnosis_severe_med ~ age.recruit + sex + ethnicity + townsend + centre, 
                                data = subsetsnomild_UKB_a_NAsdropped, method = "nearest", distance = "mahalanobis", 
                                ratio = 5) 

# saveRDS(severematching_full, file = file.path(work_data, "severematching_full.rds")) 
severematching_full <- readRDS("data_work/Matching/severematching_full.rds")

bal.tab(severematching_full, un = TRUE, binary = "std") 

# obtain matched data:
sfmatcheddata <- match.data(severematching_full)
# obtain matched data for control + asthma groups separately if needed:
TREATEDseverematching_full <- match.data(severematching_full, group = "treat")
CONTROLseverematching_full <- match.data(severematching_full, group = "control")
nrow(TREATEDseverematching_full) # = 5,072
nrow(CONTROLseverematching_full) # = 25,360

# check there are no mild asthmatics: table(sfmatcheddata$asthma_mild_ONLY) = 0


### Reduced/ simplified matching (remove townsend as a covariate):   ######

# mild:
mildmatching_notownsend <- matchit(asthma_mild_ONLY ~ age.recruit + sex + ethnicity + centre, 
                                   data = subsetnosevere_UKB_a_NAsdropped, method = "nearest", distance = "mahalanobis", 
                                   ratio = 5) 

#severe:
# again subset to remove all eid's used in mild reduced matching first
subsetsnoMR_UKB_a_NAsdropped <-
  UKB_a_drop_ethnicitytownsendNAs[!UKB_a_drop_ethnicitytownsendNAs$eid 
                                             %in% mrmatcheddata$eid, ]
nrow(subsetsnoMR_UKB_a_NAsdropped) # rows = 379,354

severematching_notownsend <- matchit(diagnosis_severe_med ~ age.recruit + sex + ethnicity + centre, 
                                      data = subsetsnoMR_UKB_a_NAsdropped, method = "nearest", distance = "mahalanobis", 
                                      ratio = 5) 

View(mildmatching_notownsend)

mildmatching_notownsend <- readRDS("data_work/Matching/mildmatching_notownsend.rds")
# ^ saved in 'Matching' folder
# obtain matched data:
mrmatcheddata <- match.data(mildmatching_notownsend)
# obtain matched data for control + asthma groups separately if needed:
TREATEDmildmatching_notownsend <- match.data(mildmatching_notownsend, group = "treat")
CONTROLmildmatching_notownsend <- match.data(mildmatching_notownsend, group = "control")
nrow(TREATEDmildmatching_notownsend) # = 19,959
nrow(CONTROLmildmatching_notownsend) # = 99,795


severematching_notownsend <- readRDS("data_work/Matching/severematching_notownsend.rds")
# obtain matched data:
srmatcheddata <- match.data(severematching_notownsend)
# obtain matched data for control + asthma groups separately if needed:
TREATEDseverematching_notownsend <- match.data(severematching_notownsend, group = "treat")
CONTROLseverematching_notownsend <- match.data(severematching_notownsend, group = "control")
nrow(TREATEDseverematching_notownsend) # = 5,072
nrow(CONTROLseverematching_notownsend) # = 25,360


# ptm <- proc.time()
# print(Sys.time())


# # Delete the only one case with nearly all information missing (id=3462265)
# & delete participants that have left UKB:
# droplist <- c(3462265, 1330343, 2376313, 3075648, 3123970, 3385704, 3387728, 
#              3809403, 4236833, 4602053, 5063590, 5158442)
# '%ni%' <- Negate('%in%')
# UKB_a_imp <- UKB_a_imp[UKB_a_imp$id %ni% droplist, ]



# BASELINE CHARACTERISTICS ########################################

install.packages("tableone")
library(tableone)

head(UKB_a) # check variables
colnames(UKB_a) 

### Table 1 - full matchings ######

#### full mild matching ######

#subset entire matched population data for variables of interest    (tip, use "dplyr::" if needed, to recognise %>%)
mfmatcheddata <- mfmatcheddata %>% 
  dplyr::select(eid,sex,age.group,age.recruit,recruit.date,cancer.baseline.all,MI.baseline,
                stroke.baseline,diabetes.fo.pre,CKD,hypertension,
                PVD,smoke,ethnicity,BMI_cat,COPD.baseline,
                sleep_apnoea.baseline,mental,townsend,
                MO, SC, UA, PD, AD, EQ5D, qol_score,
                England,Wales,Scotland,centre,
                diagnosis_all_med,asthma_mild_ONLY,diagnosis_severe_med,death.date) %>% group_by(eid)

nrow(mfmatcheddata) # 119,754 individuals in total (19,959 with asthma + 99,795 controls)
table(mfmatcheddata$asthma_mild_ONLY) # no. = 19,959
table(mfmatcheddata$diagnosis_severe_med) #check there are no severe ppl that have been used as controls

# currently in wide format. 
head(mfmatcheddata,20)

# merge with hospital costs
#wide_costs <- merge(mfmatcheddata, sumtot_cost, by=c("eid", "fu_year"), all.x = T) %>% arrange(eid, fu_year)
#head(wide_costs, 20)
# merge with admissions 
#finalwide <- merge(wide_costs, admissions_and_days, by=c("eid", "fu_year"), all.x = T)
# make NA's 0
#finalwide$admi_n[is.na(finalwide$admi_n)] <- 0 
#finalwide$admi_tot_hospdays[is.na(finalwide$admi_tot_hospdays)] <- 0 
#finalwide$totcost_byfuyear[is.na(finalwide$totcost_byfuyear)] <- 0 


cat_variables <- c("sex", "age.group", "ethnicity", "COPD.baseline", 
                   "diabetes.fo.pre", "hypertension", "MI.baseline",
                   "stroke.baseline", "cancer.baseline.all", "sleep_apnoea.baseline", 
                   "CKD", "PVD", "smoke", "BMI_cat", "mental", "townsend",
                   "MO", "SC", "UA", "PD", "AD", "EQ5D",
                   "England", "Wales", "Scotland", "centre",
                   "diagnosis_all_med", "asthma_mild_ONLY", "diagnosis_severe_med")

cont_variables <- c("admi_n", "admi_tot_hospdays", 
                    "totcost_byfuyear", "respadmi_byfuyear", 
                    "age.recruit", "qol_score")

# subset the mild asthma population
submf_asthmamild <- mfmatcheddata[mfmatcheddata$asthma_mild_ONLY == 1, ]
nrow(submf_asthmamild) # check: no. of people s/ id's = 19,959
head(submf_asthmamild)

# summarise variables for mild population

table1_mild_cat <- CreateCatTable(vars = cat_variables, data = submf_asthmamild)
table1_mild_cont <- CreateContTable(vars = cont_variables, data = submf_asthmamild)
summary(table1_mild_cat)
table1_mild_cat
table1_mild_cont

# summarise variables for the control(mild) population
submf_control <- mfmatcheddata[mfmatcheddata$asthma_mild_ONLY == 0, ]
nrow(submf_control) # 99,795 controls
n_distinct(submf_control$eid) # 

table1_controlmild_cat <- CreateCatTable(vars = cat_variables, data = submf_control)
table1_controlmild_cont <- CreateContTable(vars = cont_variables, data = submf_control)
table1_controlmild_cat
table1_controlmild_cont





#### full severe matching ######

# summarise variables for moderate-severe population

#subset entire matched population data for variables of interest    (tip, use "dplyr::" if needed, to recognise %>%)
sfmatcheddata <- sfmatcheddata %>% 
  dplyr::select(eid,sex,age.group,age.recruit,recruit.date,cancer.baseline.all,MI.baseline,
                stroke.baseline,diabetes.fo.pre,CKD,hypertension,
                PVD,smoke,ethnicity,BMI_cat,COPD.baseline,
                sleep_apnoea.baseline,mental,townsend,
                MO, SC, UA, PD, AD, EQ5D, qol_score,
                England,Wales,Scotland,centre,
                diagnosis_all_med,asthma_mild_ONLY,diagnosis_severe_med,death.date) %>% group_by(eid)

subsf_asthmasevere <- sfmatcheddata[sfmatcheddata$diagnosis_severe_med == 1, ]
nrow(subsf_asthmasevere) # check: no. of people s/ id's = 5,072

table1_severe_cat <- CreateCatTable(vars = cat_variables, data = subsf_asthmasevere)
table1_severe_cont <- CreateContTable(vars = cont_variables, data = subsf_asthmasevere)
table1_severe_cat
table1_severe_cont


# summarise variables for the control(severe) population
subsf_controls <- sfmatcheddata[sfmatcheddata$diagnosis_severe_med == 0, ] 
nrow(subsf_controls) # 25,360 controls

table1_controlsevere_cat <- CreateCatTable(vars = cat_variables, data = subsf_controls)
table1_controlsevere_cont <- CreateContTable(vars = cont_variables, data = subsf_controls)
table1_controlsevere_cat
table1_controlsevere_cont




### Table 1 - reduced matchings ######

#### reduced mild matching ######

#subset entire matched population data for variables of interest    (tip, use "dplyr::" if needed, to recognise %>%)
mrmatcheddata <- mrmatcheddata %>% 
  dplyr::select(eid,sex,age.group,age.recruit,recruit.date,cancer.baseline.all,MI.baseline,
                stroke.baseline,diabetes.fo.pre,CKD,hypertension,
                PVD,smoke,ethnicity,BMI_cat,COPD.baseline,
                sleep_apnoea.baseline,mental,townsend,
                MO, SC, UA, PD, AD, EQ5D, qol_score,
                England,Wales,Scotland,centre,
                diagnosis_all_med,asthma_mild_ONLY,diagnosis_severe_med,death.date) %>% group_by(eid)


nrow(mrmatcheddata) # 119,754 individuals in total (n with asthma + n controls)
table(mrmatcheddata$asthma_mild_ONLY) # no. = 
table(mrmatcheddata$diagnosis_severe_med) #check there are no severe ppl that have been used as controls


# subset the mild asthma population (reduced)
submr_asthmamild <- mrmatcheddata[mrmatcheddata$asthma_mild_ONLY == 1, ]
nrow(submr_asthmamild) # check: no. of people s/ id's = 19,959
head(submr_asthmamild)
# summarise variables for mild population (reduced)
table1_mrasthma_cat <- CreateCatTable(vars = cat_variables, data = submr_asthmamild)
table1_mrasthma_cont <- CreateContTable(vars = cont_variables, data = submr_asthmamild)
table1_mrasthma_cat
table1_mrasthma_cont 

# summarise variables for the controls (for mild reduced) population
submr_controls <- mrmatcheddata[mrmatcheddata$diagnosis_all_med == 0, ] 
nrow(submr_controls) # 99,795 controls
n_distinct(submr_controls) # 

table1_mrcontrols_cat <- CreateCatTable(vars = cat_variables, data = submr_controls)
table1_mrcontrols_cont <- CreateContTable(vars = cont_variables, data = submr_controls)
table1_mrcontrols_cat
table1_mrcontrols_cont



#### reduced severe matching ######

#subset entire matched population data for variables of interest    (tip, use "dplyr::" if needed, to recognise %>%)
srmatcheddata <- srmatcheddata %>% 
  dplyr::select(eid,sex,age.group,age.recruit,recruit.date,cancer.baseline.all,MI.baseline,
                stroke.baseline,diabetes.fo.pre,CKD,hypertension,
                PVD,smoke,ethnicity,BMI_cat,COPD.baseline,
                sleep_apnoea.baseline,mental,townsend,
                MO, SC, UA, PD, AD, EQ5D, qol_score,
                England,Wales,Scotland,centre,
                diagnosis_all_med,asthma_mild_ONLY,diagnosis_severe_med,death.date) %>% group_by(eid)


nrow(srmatcheddata) # 30,432 individuals in total (n with asthma + n controls)
table(srmatcheddata$diagnosis_severe_med) # no. = 5,072
table(srmatcheddata$asthma_mild_ONLY) #check there are no mild ppl that have been used as controls


# subset the severe asthma population (reduced)
subsetsr_asthma <- srmatcheddata[srmatcheddata$diagnosis_severe_med == 1, ]
# summarise variables for severe population (reduced)
table1_srasthma_cat <- CreateCatTable(vars = cat_variables, data = subsetsr_asthma)
table1_srasthma_cont <- CreateContTable(vars = cont_variables, data = subsetsr_asthma)
table1_srasthma_cat
table1_srasthma_cont 

# summarise variables for the controls (for severe reduced) population
subsetsr_controls <- srmatcheddata[srmatcheddata$diagnosis_severe_med == 0, ] 

table1_srcontrols_cat <- CreateCatTable(vars = cat_variables, data = subsetsr_controls)
table1_srcontrols_cont <- CreateContTable(vars = cont_variables, data = subsetsr_controls)
table1_srcontrols_cat
table1_srcontrols_cont




### Statistical testing ###########################

library(MASS)

## for full matchings:

#' compare mild asthma to their controls:
chisq.test(mfmatcheddata$asthma_mild_ONLY, mfmatcheddata$COPD.baseline, correct=FALSE)
chisq.test(mfmatcheddata$asthma_mild_ONLY, mfmatcheddata$diabetes.fo.pre, correct=FALSE)
chisq.test(mfmatcheddata$asthma_mild_ONLY, mfmatcheddata$hypertension, correct=FALSE)
chisq.test(mfmatcheddata$asthma_mild_ONLY, mfmatcheddata$MI.baseline, correct=FALSE)
chisq.test(mfmatcheddata$asthma_mild_ONLY, mfmatcheddata$stroke.baseline, correct=FALSE)
chisq.test(mfmatcheddata$asthma_mild_ONLY, mfmatcheddata$cancer.baseline.all, correct=FALSE)
chisq.test(mfmatcheddata$asthma_mild_ONLY, mfmatcheddata$sleep_apnoea.baseline, correct=FALSE)
chisq.test(mfmatcheddata$asthma_mild_ONLY, mfmatcheddata$CKD, correct=FALSE)
chisq.test(mfmatcheddata$asthma_mild_ONLY, mfmatcheddata$PVD, correct=FALSE)
chisq.test(mfmatcheddata$asthma_mild_ONLY, mfmatcheddata$mental, correct=FALSE)
chisq.test(mfmatcheddata$asthma_mild_ONLY, mfmatcheddata$qol_score, correct=FALSE)

# for categorical variables:
chisq.test(table(mfmatcheddata$asthma_mild_ONLY, mfmatcheddata$smoke))
chisq.test(table(mfmatcheddata$asthma_mild_ONLY, mfmatcheddata$BMI_cat))


#' compare moderate-severe asthma to their controls
chisq.test(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$COPD.baseline, correct=FALSE)
chisq.test(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$diabetes.fo.pre, correct=FALSE)
chisq.test(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$hypertension, correct=FALSE)
chisq.test(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$MI.baseline, correct=FALSE)
chisq.test(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$stroke.baseline, correct=FALSE)
chisq.test(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$cancer.baseline.all, correct=FALSE)
chisq.test(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$sleep_apnoea.baseline, correct=FALSE)
chisq.test(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$CKD, correct=FALSE)
chisq.test(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$PVD, correct=FALSE)
chisq.test(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$mental, correct=FALSE)
chisq.test(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$qol_score, correct=FALSE)


## We see that the p-value is 1 for variables that were matched for:
#chisq.test(table(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$age.recruit))
#chisq.test(table(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$sex))
#chisq.test(table(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$ethnicity))
#chisq.test(table(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$centre))
#chisq.test(table(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$townsend)) 

# TABLE for categorical variables:
chisq.test(table(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$smoke))
chisq.test(table(sfmatcheddata$diagnosis_severe_med, sfmatcheddata$BMI_cat))






##' for reduced matchings:

#' compare mild asthma to their controls:
chisq.test(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$COPD.baseline, correct=FALSE)
chisq.test(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$diabetes.fo.pre, correct=FALSE)
chisq.test(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$hypertension, correct=FALSE)
chisq.test(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$MI.baseline, correct=FALSE)
chisq.test(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$stroke.baseline, correct=FALSE)
chisq.test(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$cancer.baseline.all, correct=FALSE)
chisq.test(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$sleep_apnoea.baseline, correct=FALSE)
chisq.test(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$CKD, correct=FALSE)
chisq.test(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$PVD, correct=FALSE)
chisq.test(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$mental, correct=FALSE)
chisq.test(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$qol_score, correct=FALSE)

# for categorical variables:
chisq.test(table(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$smoke))
chisq.test(table(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$BMI_cat))
chisq.test(table(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$EQ5D.MO))
chisq.test(table(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$EQ5D.SC))
chisq.test(table(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$EQ5D.UA))
chisq.test(table(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$EQ5D.PD))
chisq.test(table(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$EQ5D.AD))

# Additional statistical test for missing matching variable townsend:
chisq.test(table(mrmatcheddata$asthma_mild_ONLY, mrmatcheddata$townsend)) 



#' compare mod-severe asthma to their controls:
chisq.test(srmatcheddata$diagnosis_severe_med, srmatcheddata$COPD.baseline, correct=FALSE)
chisq.test(srmatcheddata$diagnosis_severe_med, srmatcheddata$diabetes.fo.pre, correct=FALSE)
chisq.test(srmatcheddata$diagnosis_severe_med, srmatcheddata$hypertension, correct=FALSE)
chisq.test(srmatcheddata$diagnosis_severe_med, srmatcheddata$MI.baseline, correct=FALSE)
chisq.test(srmatcheddata$diagnosis_severe_med, srmatcheddata$stroke.baseline, correct=FALSE)
chisq.test(srmatcheddata$diagnosis_severe_med, srmatcheddata$cancer.baseline.all, correct=FALSE)
chisq.test(srmatcheddata$diagnosis_severe_med, srmatcheddata$sleep_apnoea.baseline, correct=FALSE)
chisq.test(srmatcheddata$diagnosis_severe_med, srmatcheddata$CKD, correct=FALSE)
chisq.test(srmatcheddata$diagnosis_severe_med, srmatcheddata$PVD, correct=FALSE)
chisq.test(srmatcheddata$diagnosis_severe_med, srmatcheddata$mental, correct=FALSE)
chisq.test(srmatcheddata$diagnosis_severe_med, srmatcheddata$qol_score, correct=FALSE)

# for categorical variables:
chisq.test(table(srmatcheddata$diagnosis_severe_med, srmatcheddata$smoke))
chisq.test(table(srmatcheddata$diagnosis_severe_med, srmatcheddata$BMI_cat))
chisq.test(table(srmatcheddata$diagnosis_severe_med, srmatcheddata$EQ5D.MO))
chisq.test(table(srmatcheddata$diagnosis_severe_med, srmatcheddata$EQ5D.SC))
chisq.test(table(srmatcheddata$diagnosis_severe_med, srmatcheddata$EQ5D.UA))
chisq.test(table(srmatcheddata$diagnosis_severe_med, srmatcheddata$EQ5D.PD))
chisq.test(table(srmatcheddata$diagnosis_severe_med, srmatcheddata$EQ5D.AD))
# Additional statistical test for missing matching variable townsend:
chisq.test(table(srmatcheddata$diagnosis_severe_med, srmatcheddata$townsend)) 


# saved data ----
mildmatching_full <- readRDS("data_work/Matching/mildmatching_full.rds")
severematching_full <- readRDS("data_work/Matching/severematching_full.rds")
mildmatching_notownsend <- readRDS("data_work/Matching/mildmatching_notownsend.rds")
severematching_notownsend <- readRDS("data_work/Matching/severematching_notownsend.rds")