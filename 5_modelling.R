library(tidyverse)

# Load data
mf_longhosp <- readRDS(file.path(work_data, "mf_longhosp.rds")) # mf = mild (asthma) fully (matched) cohort with all 5 covariates, plus hospital data
sf_longhosp <- readRDS(file.path(work_data, "sf_longhosp.rds")) # sf = severe (asthma) fully (matched) cohort with all 5 covariates, plus hospital data
mr_longhosp <- readRDS(file.path(work_data, "mr_longhosp.rds")) # mr = mild (asthma) reduced (matched) cohort with 4 covariates, plus hospital data
sr_longhosp <- readRDS(file.path(work_data, "sr_longhosp.rds")) # sr = severe (asthma) reduced (matched) cohort with 4 covariates, plus hospital data
# note: duplicate fu_years have already been removed (using distinct(eid, fu_year), and fu_years with 0 completion have been removed

# Convert each to 10y versions: 
# ~ exclude fu_years>10 & further filtering by fu_completion

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

summary(mr_longhosp10y$fu_completion) 
# ^ repeat above for sr_longhosp


## Modelling ############################# 

### REDUCED MATCHED MODELS ########################################### 
# "reduced" meaning asthma patients were matched to controls based on 4 variables (age, sex, ethnicity, centre/location), excluding townsend
# these were the final models used for analyses. The fully matched models (below) were sensitivity analyses.

# First factor key variables that will be used in the adjusted models (minimally, intermediate & fully adjusted).

mr_longhosp10y$ethnicity <- factor(mr_longhosp10y$ethnicity, ordered = FALSE,
                                   levels = c("White", "Black", "South Asian", "Others"))
mr_longhosp10y$smoke <- factor(mr_longhosp10y$smoke, ordered = FALSE,
                               levels = c("Never", "Current", "Previous"))
mr_longhosp10y$BMI_cat <- factor(mr_longhosp10y$BMI_cat, ordered = FALSE,
                                 levels = c("18.5-25", "<18.5", "25-30", 
                                            "30-35", "35-40", "40+"))
mr_longhosp10y$townsend <- factor(mr_longhosp10y$townsend, ordered = FALSE,
                                  levels = c(1, 2, 3, 4, 5))
mr_longhosp10y$centre <- factor(mr_longhosp10y$centre, ordered = FALSE,
                                levels = c(11010, 11002, 10003, 11001, 11003, 11004, 
                                           11005, 11006, 11007, 11008, 11009, 11011, 
                                           11012, 11013, 11014, 
                                           11016, 11017, 11018, 11020, 11021,
                                           11022, 11023)) 

####~ Mild vs no asthma ###########################################

####### ~~~ MAIN models #######################################
# Main models have NO interaction terms

## 1 - unadjusted models:

library(MASS)
install.packages("texreg")
library(texreg)

# Different models were tested - negative binomial (from MASS package) outperformed
# all other models for admissions.
# zero-inflated negative binomial models were built using pscl and boot packages

# NON- ADJUSTED #
# Admissions
# negative binomial model:
admi1_mr_x <- glm.nb(admi_n ~ asthma_mild_ONLY + fu_year, data = mr_longhosp10y)
# zero-inflated negative binomial model (ZINB)
admi2_mr_x <- zeroinfl(admi_n ~ asthma_mild_ONLY + fu_year, data = mr_longhosp10y,
                       dist = "negbin")
screenreg(list(admi1_mr_x, admi2_mr_x)) # log likelihoods are the same but AIC is lower for the nb model which is better.

# - DAYS - 
# negative binomial
days1_mr_x <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY + fu_year, data = mr_longhosp10y)
# to test zero-inflated models we need to round days to integers:
mr_longhosp10y$rounded_days <- round(mr_longhosp10y$hdays_byfuyear, digits = 0)
# zero-inflated negative binomial:
days2_mr_x <- zeroinfl(rounded_days ~ asthma_mild_ONLY + fu_year, data = mr_longhosp10y,
                       dist = "negbin")
# zero-inflated poisson (ZIP): 
days3_mr_x <- zeroinfl(rounded_days ~ asthma_mild_ONLY + fu_year, data = mr_longhosp10y,
                       dist = "poisson")
screenreg(list(days1_mr_x, days2_mr_x, days3_mr_x)) 
# ZIP model is the worst (smallest LogLik)
# NB model is the best.

# exponentiated coefficients for ZINB model:
expcoef <- exp(coef(days2_mr_x))
expcoef2 <- matrix(expcoef, ncol = 2)
rownames(expcoef2) <- names(coef(expcoef))
colnames(expcoef) <- c("Count_model", "Zero_inflation_model")
expcoef
# clustered confidence intervals:
exp((days2_mr_x__cis <- coefci(days2_mr_x, vcov = vcovCL, cluster = ~eid)))


# - COSTS - #
# GLM with gaussian family (linear regression) and identity link:
library(miceadds)
library(sandwich)
costs1_mr_x <- glm.cluster(formula = totcost_byfuyear ~ asthma_mild_ONLY + fu_year,
                           cluster ="eid", 
                           family = gaussian(link = "identity"),
                           data = mr_longhosp10y)
# negative binomial
costs2_mr_x <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY + fu_year, data = mr_longhosp10y)
# zero-inflated negative binomial (with rounded costs variable to make all integers)
mr_longhosp10y$rounded_costs <- round(mr_longhosp10y$totcost_byfuyear, digits = 0)
costs3_mr_x <- zeroinfl(rounded_costs ~ asthma_mild_ONLY + fu_year, 
                        data = mr_longhosp10y, dist = "negbin")
# zero-inflated poisson: 
costs4_mr_x <- zeroinfl(rounded_costs ~ asthma_mild_ONLY + fu_year, data = mr_longhosp10y,
                        dist = "poisson")

screenreg(list(costs1_mr_x, costs2_mr_x, costs3_mr_x, costs4_mr_x))
# standard NB (negative binomial) has the largest log likelihood, 
# beating linear regression, zero-inflated poisson (ZIP), and zero-inflated negative binomial (ZINB)



## 2 - minimally adjusted models:
# include the matched covariates (remember, no townsend)

# negative binomial models for admissions, days, and costs:

# FIRST centre age:
mr_longhosp10y$age.centred <- scale(mr_longhosp10y$age.recruit, center = TRUE, scale = FALSE)
# ^ age is centred at 56-years.

admi1_mr_adj_x <- glm.nb(admi_n ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)

days1_mr_adj_x <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
days2_mr_adj_x <- zeroinfl(rounded_days ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y,
                           dist = "negbin")
costs1_mr_adj_x <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)

screenreg(list(admi1_mr_adj_x, days1_mr_adj_x, costs1_mr_adj_x))

(expcoef <- exp(coef(days2_mr_adj_x)))


## 3 - maximally adjusted models:
# include the matched covariates AND townsend, BMI, smoking, and comorbidities

# negative binomial model for admissions
admi1_mr_maxadj_x <- glm.nb(admi_n ~ asthma_mild_ONLY + fu_year + age.centred
                            + sex + ethnicity + townsend + centre
                            + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                            + hypertension + MI.baseline + stroke.baseline 
                            + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                            + PVD + mental, 
                            data = mr_longhosp10y)

# negative binomial model for days
days1_mr_maxadj_x <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY + fu_year + age.centred
                            + sex + ethnicity + townsend + centre
                            + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                            + hypertension + MI.baseline + stroke.baseline 
                            + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                            + PVD + mental, 
                            data = mr_longhosp10y)
# ZINB test
days2_mr_maxadj_x <- zeroinfl(rounded_days ~ asthma_mild_ONLY + fu_year + age.centred
                              + sex + ethnicity + townsend + centre
                              + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                              + hypertension + MI.baseline + stroke.baseline 
                              + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                              + PVD + mental, 
                              data = mr_longhosp10y, dist = "negbin")


# negative binomial model for costs
costs1_mr_maxadj_x <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY + fu_year + age.centred
                             + sex + ethnicity + townsend + centre
                             + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                             + hypertension + MI.baseline + stroke.baseline 
                             + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                             + PVD + mental, 
                             data = mr_longhosp10y)


## 4 - INTERMEDIATE adjusted models: (only for the "reduced" matched populations) #

# adjust for all 4 matching covariates AND townsend ~ by comparing these models to 
# minimally adjusted models we can estimate the contribution of deprivation

admi1_mr_interadj_x <- glm.nb(admi_n ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
# compare all 4 models of admissions: [all neg. binomial]
screenreg(list(admi1_mr_x, admi1_mr_adj_x, admi1_mr_interadj_x, admi1_mr_maxadj_x))

days1_mr_interadj_x <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)

costs1_mr_interadj_x <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)

screenreg(list(admi1_mr_interadj_x, days1_mr_interadj_x, costs1_mr_interadj_x))


##### clustered standard errors: ###

## calculate clustered standard errors + confidence intervals:

install.packages("lmtest")
library(lmtest)
library(sandwich)
library(boot)

# for REDUCED matched models
# for mild asthma vs no-asthma models:
# COEF TEST (standard errors are clustered by eid) to obtain TRUE significance level:
round(admi1_mr_x_ coef <- coeftest(admi1_mr_x, vcov. = vcovCL, cluster = ~eid), digits = 2)
(admi1_mr_adj_x_coef <- coeftest(admi1_mr_adj_x, vcov. = vcovCL, cluster = ~eid))
(admi1_mr_interadj_x_coef <- coeftest(admi1_mr_interadj_x, vcov. = vcovCL, cluster = ~eid))
(admi1_mr_maxadj_x_coef <- coeftest(admi1_mr_maxadj_x, vcov. = vcovCL, cluster = ~eid))

(days1_mr_x_coef <- coeftest(days1_mr_x, vcov. = vcovCL, cluster = ~eid))
(days1_mr_adj_x_coef <- coeftest(days1_mr_adj_x, vcov. = vcovCL, cluster = ~eid))
(days1_mr_interadj_x_coef <- coeftest(days1_mr_interadj_x, vcov. = vcovCL, cluster = ~eid))
(days1_mr_maxadj_x_coef <- coeftest(days1_mr_maxadj_x, vcov. = vcovCL, cluster = ~eid))

(costs2_mr_x_coef <- coeftest(costs2_mr_x, vcov. = vcovCL, cluster = ~eid))
(costs1_mr_adj_x_coef <- coeftest(costs1_mr_adj_x, vcov. = vcovCL, cluster = ~eid))
(costs1_mr_interadj_x_coef <- coeftest(costs1_mr_interadj_x, vcov. = vcovCL, cluster = ~eid))
(costs1_mr_maxadj_x_coef <- coeftest(costs1_mr_maxadj_x, vcov. = vcovCL, cluster = ~eid))

# EXPONENTIATED COEFFICIENTS:  (this converts coefficients to Incident Rate Ratios)
as.matrix(round(exp(coef(admi1_mr_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(admi1_mr_adj_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(admi1_mr_interadj_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(admi1_mr_maxadj_x)), digits = 2), ncol = 2, rownames = TRUE)

as.matrix(round(exp(coef(days1_mr_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_mr_adj_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_mr_interadj_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_mr_maxadj_x)), digits = 2), ncol = 2, rownames = TRUE)

as.matrix(round(exp(coef(costs2_mr_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_mr_adj_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_mr_interadj_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_mr_maxadj_x)), digits = 2), ncol = 2, rownames = TRUE)

# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((admi1_mr_x_cis <- coefci(admi1_mr_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((admi1_mr_adj_x_cis <- coefci(admi1_mr_adj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((admi1_mr_interadj_x_cis <- coefci(admi1_mr_interadj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((admi1_mr_maxadj_x_cis <- coefci(admi1_mr_maxadj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)

round(exp((days1_mr_x_cis <- coefci(days1_mr_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days1_mr_adj_x_cis <- coefci(days1_mr_adj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days1_mr_interadj_x_cis <- coefci(days1_mr_interadj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days1_mr_maxadj_x_cis <- coefci(days1_mr_maxadj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)

round(exp((costs2_mr_x_cis <- coefci(costs2_mr_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((costs1_mr_adj_x_cis <- coefci(costs1_mr_adj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((costs1_mr_interadj_x_cis <- coefci(costs1_mr_interadj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((costs1_mr_maxadj_x_cis <- coefci(costs1_mr_maxadj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)




####### ~~~ check fu_year interaction #######################################

# Compare models with and without interaction for each level of adjustment,
# then compare USING LIKELIHOOD RATIO TESTS

# admissions first

admi1_mr <- glm.nb(admi_n ~ asthma_mild_ONLY * fu_year, data = mr_longhosp10y)
lrtest(admi1_mr_x, admi1_mr)
# 1st model (with interaction) has a larger log likelihood than without interaction (-829609 vs -829615)
# however, the p value is not significant 
# this means both models fit the data equally well

admi1_mr_adj <- glm.nb(admi_n ~ asthma_mild_ONLY * fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
lrtest(admi1_mr_adj_x, admi1_mr_adj)
# again, though 1st model does have larger LogLik,
# p>0.05 ~ the p value isn't significant
# thus the model WITHOUT interaction is better

# [INTER-ADJUSTED]:
admi1_mr_interadj <- glm.nb(admi_n ~ asthma_mild_ONLY * fu_year + age.centred
                            + sex + ethnicity + townsend + centre, 
                            data = mr_longhosp10y)
lrtest(admi1_mr_interadj_x, admi1_mr_interadj)
# p value is not significant

# maximally adjusted:
admi1_mr_maxadj <- glm.nb(admi_n ~ asthma_mild_ONLY * fu_year + age.centred
                          + sex + ethnicity + townsend + centre
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental, 
                          data = mr_longhosp10y)
lrtest(admi1_mr_maxadj_x, admi1_mr_maxadj)
# p value is NOT significant 


# now length of stay:
days1_mr <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * fu_year, data = mr_longhosp10y)
lrtest(days1_mr_x, days1_mr)
# p value is NOT significant

days1_mr_adj <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * fu_year + age.centred
                       + sex + ethnicity + centre, 
                       data = mr_longhosp10y)
lrtest(days1_mr_adj_x, days1_mr_adj)
#!!  p value IS significant (p<0.001) ***
# thus model WITH interaction is better
days1_mr_maxadj <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * fu_year + age.centred
                          + sex + ethnicity + townsend + centre 
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental,
                          data = mr_longhosp10y)
lrtest(days1_mr_maxadj_x, days1_mr_maxadj)
# p value IS significant (p<0.001) ***
# thus the model WITH interaction is better.

# [INTER-ADJUSTED]:
days1_mr_interadj <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * fu_year + age.centred
                            + sex + ethnicity + townsend + centre, 
                            data = mr_longhosp10y)
lrtest(days1_mr_interadj_x, days1_mr_interadj)
# p value IS significant (p<0.01) **
# thus the model WITH interaction is better.



## now COSTS: ##

costs1_mr <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * fu_year, data = mr_longhosp10y)
lrtest(costs2_mr_x, costs1_mr)
# p value is NOT significant
# thus the model WITHOUT interaction is better
costs1_mr_adj <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * fu_year + age.centred
                        + sex + ethnicity + centre, 
                        data = mr_longhosp10y)
lrtest(costs1_mr_adj_x, costs1_mr_adj)
# p value is NOT significant
# thus the model WITHOUT interaction is better
costs1_mr_maxadj <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * fu_year + age.centred
                           + sex + ethnicity + townsend + centre
                           + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                           + hypertension + MI.baseline + stroke.baseline 
                           + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                           + PVD + mental,
                           data = mr_longhosp10y)
lrtest(costs1_mr_maxadj_x, costs1_mr_maxadj)
# p value is NOT significant
# thus the model WITHOUT interaction is better
# [INTER-ADJUSTED]:
costs1_mr_interadj <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * fu_year + age.centred
                             + sex + ethnicity + townsend + centre, 
                             data = mr_longhosp10y)
lrtest(costs1_mr_interadj_x, costs1_mr_interadj)
# p value is NOT significant
# thus the model WITHOUT interaction is better


##### clustered standard errors: ###

## calculate clustered standard errors + confidence intervals:

install.packages("lmtest")
library(lmtest)
library(sandwich)
library(boot)

# mild asthma vs no-asthma models:
# COEF TEST (clustered by eid) to obtain TRUE significance level:
round((admi1_mr_coef <- coeftest(admi1_mr, vcov. = vcovCL, cluster = ~eid)), digits = 2)
(admi1_mr_adj_coef <- coeftest(admi1_mr_adj, vcov. = vcovCL, cluster = ~eid))
(admi1_mr_interadj_coef <- coeftest(admi1_mr_interadj, vcov. = vcovCL, cluster = ~eid))
(admi1_mr_maxadj_coef <- coeftest(admi1_mr_maxadj, vcov. = vcovCL, cluster = ~eid))

(days1_mr_coef <- coeftest(days1_mr, vcov. = vcovCL, cluster = ~eid))
(days1_mr_adj_coef <- coeftest(days1_mr_adj, vcov. = vcovCL, cluster = ~eid))
(days1_mr_interadj_coef <- coeftest(days1_mr_interadj, vcov. = vcovCL, cluster = ~eid))
(days1_mr_maxadj_coef <- coeftest(days1_mr_maxadj, vcov. = vcovCL, cluster = ~eid))

(costs2_mr_coef <- coeftest(costs2_mr, vcov. = vcovCL, cluster = ~eid))
(costs1_mr_adj_coef <- coeftest(costs1_mr_adj, vcov. = vcovCL, cluster = ~eid))
(costs1_mr_interadj_coef <- coeftest(costs1_mr_interadj, vcov. = vcovCL, cluster = ~eid))
(costs1_mr_maxadj_coef <- coeftest(costs1_mr_maxadj, vcov. = vcovCL, cluster = ~eid))

# EXPONENTIATED COEFFICIENTS:  (this converts coefficients to Incident Rate Ratios)
as.matrix(round(exp(coef(admi1_mr)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(admi1_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(admi1_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(admi1_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)

as.matrix(round(exp(coef(days1_mr)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)

as.matrix(round(exp(coef(costs1_mr)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)

# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((admi1_mr_cis <- coefci(admi1_mr, vcov = vcovCL, cluster = ~eid))), digits = 2)  




############ ~~~ check continuous fu_year interaction ###########################################

# first create new continuous variable for follow-up year (for models with only a single interaction term)
mr_longhosp10y <- mr_longhosp10y %>% 
  mutate(new_fu_year = case_when(mr_longhosp10y$fu_year == "year_1" ~ 1,
                                 mr_longhosp10y$fu_year == "year_2" ~ 2,
                                 mr_longhosp10y$fu_year == "year_3" ~ 3,
                                 mr_longhosp10y$fu_year == "year_4" ~ 4,
                                 mr_longhosp10y$fu_year == "year_5" ~ 5,
                                 mr_longhosp10y$fu_year == "year_6" ~ 6,
                                 mr_longhosp10y$fu_year == "year_7" ~ 7,
                                 mr_longhosp10y$fu_year == "year_8" ~ 8,
                                 mr_longhosp10y$fu_year == "year_9" ~ 9,
                                 mr_longhosp10y$fu_year == "year_10" ~ 10))

## Model follow-up year as a single continuous variable for all outcomes 
## using the 'new_fu_year' variable:
# admissions:

admi_mr_n <- glm.nb(admi_n ~ asthma_mild_ONLY * new_fu_year, data = mr_longhosp10y)

admi_mr_adj_n <- glm.nb(admi_n ~ asthma_mild_ONLY * new_fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)

admi_mr_interadj_n <- glm.nb(admi_n ~ asthma_mild_ONLY * new_fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)

admi_mr_maxadj_n <- glm.nb(admi_n ~ asthma_mild_ONLY * new_fu_year + age.centred
                           + sex + ethnicity + townsend + centre
                           + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                           + hypertension + MI.baseline + stroke.baseline 
                           + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                           + PVD + mental, 
                           data = mr_longhosp10y)

screenreg(list(admi_mr_n, admi_mr_adj_n, admi_mr_interadj_n, admi_mr_maxadj_n))


# days: #
days_mr_n <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * new_fu_year, data = mr_longhosp10y)

days_mr_adj_n <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * 
                          new_fu_year + age.centred + sex 
                        + ethnicity + centre, data = mr_longhosp10y)

days_mr_interadj_n <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * new_fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)

days_mr_maxadj_n <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * new_fu_year + age.centred
                           + sex + ethnicity + townsend + centre
                           + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                           + hypertension + MI.baseline + stroke.baseline 
                           + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                           + PVD + mental, 
                           data = mr_longhosp10y)

screenreg(list(days_mr_adj_n, days_mr_interadj_n, days_mr_maxadj_n))


# costs: #
costs_mr_n <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * new_fu_year, data = mr_longhosp10y)

costs_mr_adj_n <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * 
                           new_fu_year + age.centred + sex 
                         + ethnicity + centre, data = mr_longhosp10y)

costs_mr_interadj_n <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * 
                                new_fu_year + age.centred + sex 
                              + ethnicity + centre + townsend, data = mr_longhosp10y)

costs_mr_maxadj_n <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * new_fu_year + age.centred
                            + sex + ethnicity + townsend + centre
                            + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                            + hypertension + MI.baseline + stroke.baseline 
                            + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                            + PVD + mental, 
                            data = mr_longhosp10y)

screenreg(list(costs_mr_n, costs_mr_adj_n, costs_mr_interadj_n, costs_mr_maxadj_n))



# CLUSTERED STANDARD ERRORS:

library(lmtest)
library(sandwich)
library(boot)
(admi_mr_n_coef <- coeftest(admi_mr_n, vcov. = vcovCL, cluster = ~eid))
(admi_mr_adj_n_coef <- coeftest(admi_mr_adj_n, vcov. = vcovCL, cluster = ~eid))
(admi_mr_interadj_n_coef <- coeftest(admi_mr_interadj_n, vcov. = vcovCL, cluster = ~eid))
(admi_mr_maxadj_n_coef <- coeftest(admi_mr_maxadj_n, vcov. = vcovCL, cluster = ~eid))

(days_mr_n_coef <- coeftest(days_mr_n, vcov. = vcovCL, cluster = ~eid))
(days_mr_adj_n_coef <- coeftest(days_mr_adj_n, vcov. = vcovCL, cluster = ~eid))
(days_mr_interadj_n_coef <- coeftest(days_mr_interadj_n, vcov. = vcovCL, cluster = ~eid))
(days_mr_maxadj_n_coef <- coeftest(days_mr_maxadj_n, vcov. = vcovCL, cluster = ~eid))

(costs_mr_n_coef <- coeftest(costs_mr_n, vcov. = vcovCL, cluster = ~eid))
(costs_mr_adj_n_coef <- coeftest(costs_mr_adj_n, vcov. = vcovCL, cluster = ~eid))
(costs_mr_interadj_n_coef <- coeftest(costs_mr_interadj_n, vcov. = vcovCL, cluster = ~eid))
(costs_mr_maxadj_n_coef <- coeftest(costs_mr_maxadj_n, vcov. = vcovCL, cluster = ~eid))


#Likelihood ratio tests comparing categorical fu_year & continuous new_fu_year:
lrtest(admi1_mr, admi_mr_n)    # (p<0.001) *** = categorical fu_year = better
lrtest(admi1_mr_adj, admi_mr_adj_n) # (p<0.001) *** = categorical fu_year = better
lrtest(admi1_mr_interadj, admi_mr_interadj_n) # (p<0.001) *** = categorical fu_year = better
lrtest(admi1_mr_maxadj, admi_mr_maxadj_n) # (p<0.001) *** = categorical fu_year = better

lrtest(days1_mr, days_mr_n)    # (p<0.01) ** = categorical fu_year = better
lrtest(days1_mr_adj, days_mr_adj_n) # (p<0.001) *** = categorical fu_year = better
lrtest(days1_mr_interadj, days_mr_interadj_n) # (p<0.001) *** = categorical fu_year = better
lrtest(days1_mr_maxadj, days_mr_maxadj_n) # (p<0.001) *** = categorical fu_year = better

lrtest(costs2_mr, costs_mr_n)    # NO significance. Tiny difference in LogLik scores.
lrtest(costs1_mr_adj, costs_mr_adj_n) # (p<0.001) *** = categorical fu_year = better
lrtest(costs1_mr_interadj, costs_mr_interadj_n) # NO significance.
lrtest(costs1_mr_maxadj, costs_mr_maxadj_n) # NO significance.

# EXPONENTIATED COEFFICIENTS:  (this converts coefficients to Incident Rate Ratios)
as.matrix(round(exp(coef(admi_mr_adj_n)), digits = 2), ncol = 2, rownames = TRUE)
# ^ repeat for all models

# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((admi_mr_adj_n_cis <- coefci(admi_mr_adj_n, vcov = vcovCL, cluster = ~eid))), digits = 2)
# ^ repeat for all models


####### ~~~ check deprivation interaction #######################################

#  - admissions - 
a_mr_dep <- glm.nb(admi_n ~ asthma_mild_ONLY * townsend, data = mr_longhosp10y)
a_mr_adj_dep <- glm.nb(admi_n ~ asthma_mild_ONLY * townsend + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
a_mr_interadj_dep <- glm.nb(admi_n ~ asthma_mild_ONLY * townsend + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
a_mr_maxadj_dep <- glm.nb(admi_n ~ asthma_mild_ONLY * townsend + age.centred
                          + sex + ethnicity + townsend + centre
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental, data = mr_longhosp10y)

# lrtest comparing with and without townsend (categorical) interaction:
# i.e., comparing interaction with FU_year OR townsend
lrtest(admi1_mr, a_mr_dep)
# The model WITHOUT townsend interaction has a larger LogLik than with deprivation interaction.
# p value IS significant (p<0.001) ***
# thus the model WITHOUT townsend interaction is better.
lrtest(admi1_mr_adj, a_mr_adj_dep) # p value IS significant (p<0.001) ***, WITHOUT townsend interaction is better.
lrtest(admi1_mr_interadj, a_mr_interadj_dep) # significant *** = WITHOUT townsend is better
lrtest(admi1_mr_maxadj, a_mr_maxadj_dep) # significant *** = WITHOUT townsend is better


#  - days -
h_mr_dep <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * townsend, data = mr_longhosp10y)
h_mr_adj_dep <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * townsend + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
h_mr_interadj_dep <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * townsend + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
h_mr_maxadj_dep <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * townsend + age.centred
                          + sex + ethnicity + townsend + centre
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental, 
                          data = mr_longhosp10y)

# lrtest comparing models with and without townsend (categorical) interaction:
lrtest(days1_mr, h_mr_dep) 
lrtest(days1_mr_adj, h_mr_adj_dep) 
lrtest(days1_mr_maxadj, h_mr_maxadj_dep) # significant *** = WITHOUT townsend is better


#  - costs -
c_mr_dep <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * townsend, data = mr_longhosp10y)
c_mr_adj_dep <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * townsend + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
c_mr_interadj_dep <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * townsend + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
c_mr_maxadj_dep <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * townsend + age.centred
                          + sex + ethnicity + townsend + centre
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental, 
                          data = mr_longhosp10y)
# lrtest comparing models with and without townsend (categorical) interaction:
lrtest(costs2_mr, c_mr_dep) # significant *** = WITHOUT townsend is better
lrtest(costs1_mr_adj, c_mr_adj_dep) # significant *** = WITHOUT townsend is better
lrtest(costs1_mr_maxadj, c_mr_maxadj_dep) # significant *** = WITHOUT townsend is better


# -- cLUSTERED STANDARD ERRORS --:
(a_mr_dep_coef <- coeftest(a_mr_dep, vcov. = vcovCL, cluster = ~eid))
(a_mr_adj_dep_coef <- coeftest(a_mr_adj_dep, vcov. = vcovCL, cluster = ~eid))
(a_mr_interadj_dep_coef <- coeftest(a_mr_interadj_dep, vcov. = vcovCL, cluster = ~eid))
(a_mr_maxadj_dep_coef <- coeftest(a_mr_maxadj_dep, vcov. = vcovCL, cluster = ~eid))

(h_mr_dep_coef <- coeftest(h_mr_dep, vcov. = vcovCL, cluster = ~eid))
(h_mr_adj_dep_coef <- coeftest(h_mr_adj_dep, vcov. = vcovCL, cluster = ~eid))
(h_mr_interadj_dep_coef <- coeftest(h_mr_interadj_dep, vcov. = vcovCL, cluster = ~eid))
(h_mr_maxadj_dep_coef <- coeftest(h_mr_maxadj_dep, vcov. = vcovCL, cluster = ~eid))

(c_mr_dep_coef <- coeftest(c_mr_dep, vcov. = vcovCL, cluster = ~eid))
(c_mr_adj_dep_coef <- coeftest(c_mr_adj_dep, vcov. = vcovCL, cluster = ~eid))
(c_mr_interadj_dep_coef <- coeftest(c_mr_interadj_dep, vcov. = vcovCL, cluster = ~eid))
(c_mr_maxadj_dep_coef <- coeftest(c_mr_maxadj_dep, vcov. = vcovCL, cluster = ~eid))





######## ~~~ check continuous deprivation interaction #######################################

# convert townsend from factored into a count variable:

mr_longhosp10y$new_townsend <- as.numeric(mr_longhosp10y$townsend)

# test all model specifications except for inter-adjusted (because we already include townsend as interaction term)

# - admissions - 
a_mr_depcont <- glm.nb(admi_n ~ asthma_mild_ONLY * new_townsend, data = mr_longhosp10y)
a_mr_adj_depcont <- glm.nb(admi_n ~ asthma_mild_ONLY * new_townsend + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
a_mr_interadj_depcont <- glm.nb(admi_n ~ asthma_mild_ONLY * new_townsend + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
a_mr_maxadj_depcont <- glm.nb(admi_n ~ asthma_mild_ONLY * new_townsend + age.centred
                              + sex + ethnicity + townsend + centre
                              + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                              + hypertension + MI.baseline + stroke.baseline 
                              + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                              + PVD + mental, data = mr_longhosp10y)
# lrtest comparing with and without townsend (numeric/ continuous) interaction:
lrtest(admi1_mr, a_mr_depcont) # The model WITHOUT townsend interaction has a larger LogLik than with deprivation interaction.
# p value IS significant (p<0.001) ***, thus the model WITHOUT townsend interaction is better.
lrtest(admi1_mr_adj, a_mr_adj_depcont) # p value IS significant (p<0.001) ***, WITHOUT townsend interaction is better.
lrtest(admi1_mr_maxadj, a_mr_maxadj_depcont) # significant *** = WITHOUT townsend is better


#screenreg to get log likelihoods:
screenreg(list(a_mr_depcont, a_mr_adj_depcont,a_mr_maxadj_depcont))

# CLUSTERED STANDARD ERRORS:
round((a_mr_depcont_coef <- coeftest(a_mr_depcont, vcov. = vcovCL, cluster = ~eid)), digits = 2)
(a_mr_adj_depcont_coef <- coeftest(a_mr_adj_depcont, vcov. = vcovCL, cluster = ~eid))
(a_mr_interadj_depcont_coef <- coeftest(a_mr_interadj_depcont, vcov. = vcovCL, cluster = ~eid))
(a_mr_maxadj_depcont_coef <- coeftest(a_mr_maxadj_depcont, vcov. = vcovCL, cluster = ~eid))
# EXPONENTIATED COEFFICIENTS:
as.matrix(round(exp(coef(a_mr_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(a_mr_adj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(a_mr_interadj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(a_mr_maxadj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((a_mr_depcont_cis <- coefci(a_mr_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((a_mr_adj_depcont_cis <- coefci(a_mr_adj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((a_mr_interadj_depcont_cis <- coefci(a_mr_interadj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((a_mr_maxadj_depcont_cis <- coefci(a_mr_maxadj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)


# - days - 
d_mr_depcont <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * new_townsend, data = mr_longhosp10y)
d_mr_adj_depcont <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * new_townsend + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
d_mr_interadj_depcont <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * new_townsend + age.centred + sex + ethnicity + centre + new_townsend, data = mr_longhosp10y)
d_mr_maxadj_depcont <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * new_townsend + age.centred
                              + sex + ethnicity + townsend + centre
                              + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                              + hypertension + MI.baseline + stroke.baseline 
                              + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                              + PVD + mental, data = mr_longhosp10y)
#screenreg to get log likelihoods:
screenreg(list(d_mr_depcont, d_mr_adj_depcont,d_mr_maxadj_depcont))
# lrtest comparing with and without townsend (numeric/ continuous) interaction:
lrtest(days1_mr, d_mr_depcont) 
lrtest(days1_mr_adj, d_mr_adj_depcont) 
lrtest(days1_mr_maxadj, d_mr_maxadj_depcont) # significant *** = WITHOUT townsend is better

# CLUSTERED STANDARD ERRORS:
(d_mr_depcont_coef <- coeftest(d_mr_depcont, vcov. = vcovCL, cluster = ~eid))
(d_mr_adj_depcont_coef <- coeftest(d_mr_adj_depcont, vcov. = vcovCL, cluster = ~eid))
(d_mr_interadj_depcont_coef <- coeftest(d_mr_interadj_depcont, vcov. = vcovCL, cluster = ~eid))
(d_mr_maxadj_depcont_coef <- coeftest(d_mr_maxadj_depcont, vcov. = vcovCL, cluster = ~eid))
# EXPONENTIATED COEFFICIENTS:
as.matrix(round(exp(coef(d_mr_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(d_mr_adj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(d_mr_interadj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(d_mr_maxadj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((d_mr_depcont_cis <- coefci(d_mr_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((d_mr_adj_depcont_cis <- coefci(d_mr_adj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((d_mr_interadj_depcont_cis <- coefci(d_mr_interadj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((d_mr_maxadj_depcont_cis <- coefci(d_mr_maxadj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)


# - costs - 
c_mr_depcont <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * new_townsend, data = mr_longhosp10y)
c_mr_adj_depcont <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * new_townsend + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
c_mr_interadj_depcont <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * new_townsend + age.centred + sex + ethnicity + centre + new_townsend, data = mr_longhosp10y)
c_mr_maxadj_depcont <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * new_townsend + age.centred
                              + sex + ethnicity + townsend + centre
                              + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                              + hypertension + MI.baseline + stroke.baseline 
                              + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                              + PVD + mental, data = mr_longhosp10y)
#screenreg to get log likelihoods:
screenreg(list(c_mr_depcont, c_mr_adj_depcont,c_mr_maxadj_depcont))
# lrtest comparing with and without townsend (numeric/ continuous) interaction:
lrtest(costs1_mr, c_mr_depcont) # 
lrtest(costs1_mr_adj, c_mr_adj_depcont) # 
lrtest(costs1_mr_maxadj, c_mr_maxadj_depcont) # 

# CLUSTERED STANDARD ERRORS:
(c_mr_depcont_coef <- coeftest(c_mr_depcont, vcov. = vcovCL, cluster = ~eid))
(c_mr_adj_depcont_coef <- coeftest(c_mr_adj_depcont, vcov. = vcovCL, cluster = ~eid))
(c_mr_interadj_depcont_coef <- coeftest(c_mr_interadj_depcont, vcov. = vcovCL, cluster = ~eid))
(c_mr_maxadj_depcont_coef <- coeftest(c_mr_maxadj_depcont, vcov. = vcovCL, cluster = ~eid))
# EXPONENTIATED COEFFICIENTS:
as.matrix(round(exp(coef(c_mr_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(c_mr_adj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(c_mr_interadj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(c_mr_maxadj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((c_mr_depcont_cis <- coefci(c_mr_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((c_mr_adj_depcont_cis <- coefci(c_mr_adj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((c_mr_interadj_depcont_cis <- coefci(c_mr_interadj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((c_mr_maxadj_depcont_cis <- coefci(c_mr_maxadj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)





############## ~~~ by RESPIRATORY admissions ###########################################

# (comparing NB and ZINB models showed that, similar to overall admissions,
# NB is the better model ~ log likelihoods were identical but NB'S AIC was lower)

respadmi_mr <- glm.nb(respadmi_byfuyear ~ asthma_mild_ONLY * fu_year, data = mr_longhosp10y)

respadmi_mr_adj <- glm.nb(respadmi_byfuyear ~ asthma_mild_ONLY * fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)

respadmi_mr_interadj <- glm.nb(respadmi_byfuyear ~ asthma_mild_ONLY * fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)

respadmi_mr_maxadj <- glm.nb(respadmi_byfuyear ~ asthma_mild_ONLY * fu_year + age.centred
                             + sex + ethnicity + townsend + centre
                             + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                             + hypertension + MI.baseline + stroke.baseline 
                             + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                             + PVD + mental, 
                             data = mr_longhosp10y)

screenreg(list(respadmi_mr, respadmi_mr_adj, respadmi_mr_interadj, respadmi_mr_maxadj))


install.packages("lmtest")
library(lmtest)
library(sandwich)
library(boot)

# use coeftest to get true significance levels (because clustering by eid's updates the SE's)
round((respadmi_mr_coef <- coeftest(respadmi_mr, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((respadmi_mr_adj_coef <- coeftest(respadmi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((respadmi_mr_adj_interadj_coef <- coeftest(respadmi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((respadmi_mr_adj_maxadj_coef <- coeftest(respadmi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
# exponentiated coefficients:
as.matrix(round(exp(coef(respadmi_mr)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(respadmi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(respadmi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(respadmi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
# ^ round coefficients to 2 digits (avoids the e+01 notation etc with long no.'s)

# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((respadmi_mr_cis <- coefci(respadmi_mr, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((respadmi_mr_adj_cis <- coefci(respadmi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((respadmi_mr_interadj_cis <- coefci(respadmi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((respadmi_mr_maxadj_cis <- coefci(respadmi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)


############## ~~~ by other ICD admissions ###########################################

#mr_longhosp10y$ICD_diag01 <- factor(mr_longhosp10y$ICD_diag01, ordered = FALSE,
#                                   levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 21, 22))

#a_mr_icd <- glm.nb(admi_n ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + ICD_diag01, data = mr_longhosp10y)
#screenreg(a_mr_icd)


# All miminally adjusted models: 
ICD1admi_mr_adj <- glm.nb(ICD1_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD2admi_mr_adj <- glm.nb(ICD2_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD3admi_mr_adj <- glm.nb(ICD3_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD4admi_mr_adj <- glm.nb(ICD4_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD5admi_mr_adj <- glm.nb(ICD5_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD6admi_mr_adj <- glm.nb(ICD6_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD7admi_mr_adj <- glm.nb(ICD7_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD8admi_mr_adj <- glm.nb(ICD8_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD9admi_mr_adj <- glm.nb(ICD9_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD10admi_mr_adj <- glm.nb(ICD10_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD11admi_mr_adj <- glm.nb(ICD11_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD12admi_mr_adj <- glm.nb(ICD12_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD13admi_mr_adj <- glm.nb(ICD13_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD14admi_mr_adj <- glm.nb(ICD14_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD15admi_mr_adj <- glm.nb(ICD15_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD16admi_mr_adj <- glm.nb(ICD16_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD17admi_mr_adj <- glm.nb(ICD17_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD18admi_mr_adj <- glm.nb(ICD18_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD19admi_mr_adj <- glm.nb(ICD19_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD20admi_mr_adj <- glm.nb(ICD20_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD21admi_mr_adj <- glm.nb(ICD21_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)
ICD22admi_mr_adj <- glm.nb(ICD22_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre, data = mr_longhosp10y)

# use coeftest to get true significance levels (because clustering by eid's updates the SE's)
# NB: no admissions recorded associated with ICD 16, 20 or 22. not enough admi's for 15 --
round((ICD1admi_mr_adj_coef <- coeftest(ICD1admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD2admi_mr_adj_coef <- coeftest(ICD2admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD3admi_mr_adj_coef <- coeftest(ICD3admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD4admi_mr_adj_coef <- coeftest(ICD4admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD5admi_mr_adj_coef <- coeftest(ICD5admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD6admi_mr_adj_coef <- coeftest(ICD6admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD7admi_mr_adj_coef <- coeftest(ICD7admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD8admi_mr_adj_coef <- coeftest(ICD8admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD9admi_mr_adj_coef <- coeftest(ICD9admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD10admi_mr_adj_coef <- coeftest(ICD10admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD11admi_mr_adj_coef <- coeftest(ICD11admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD12dmi_mr_adj_coef <- coeftest(ICD12admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD13admi_mr_adj_coef <- coeftest(ICD13admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD14admi_mr_adj_coef <- coeftest(ICD14admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD15admi_mr_adj_coef <- coeftest(ICD15admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD16admi_mr_adj_coef <- coeftest(ICD16admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD17admi_mr_adj_coef <- coeftest(ICD17admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD18admi_mr_adj_coef <- coeftest(ICD18admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD19admi_mr_adj_coef <- coeftest(ICD19admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD20admi_mr_adj_coef <- coeftest(ICD20admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD21admi_mr_adj_coef <- coeftest(ICD21admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD22admi_mr_adj_coef <- coeftest(ICD22admi_mr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
# exponentiated coefficients:
as.matrix(round(exp(coef(ICD1admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD2admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD3admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD4admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD5admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD6admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD7admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD8admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD9admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD10admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD11admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD12admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD13admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD14admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD15admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD16admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD17admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD18admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD19admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD20admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD21admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD22admi_mr_adj)), digits = 2), ncol = 2, rownames = TRUE)
# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((ICD1admi_mr_adj_cis <- coefci(ICD1admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD2admi_mr_adj_cis <- coefci(ICD2admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD3admi_mr_adj_cis <- coefci(ICD3admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD4admi_mr_adj_cis <- coefci(ICD4admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD5admi_mr_adj_cis <- coefci(ICD5admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD6admi_mr_adj_cis <- coefci(ICD6admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD7admi_mr_adj_cis <- coefci(ICD7admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD8admi_mr_adj_cis <- coefci(ICD8admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD9admi_mr_adj_cis <- coefci(ICD9admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD10admi_mr_adj_cis <- coefci(ICD10admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD11admi_mr_adj_cis <- coefci(ICD11admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD12admi_mr_adj_cis <- coefci(ICD12admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD13admi_mr_adj_cis <- coefci(ICD13admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD14admi_mr_adj_cis <- coefci(ICD14admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD15admi_mr_adj_cis <- coefci(ICD15admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD16admi_mr_adj_cis <- coefci(ICD16admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD17admi_mr_adj_cis <- coefci(ICD17admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD18admi_mr_adj_cis <- coefci(ICD18admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD19admi_mr_adj_cis <- coefci(ICD19admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD20admi_mr_adj_cis <- coefci(ICD20admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD21admi_mr_adj_cis <- coefci(ICD21admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD22admi_mr_adj_cis <- coefci(ICD22admi_mr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)

# All intermediate adjusted models: 
#  models not run for ICD 15, 16, 20 & 22 bc either no or not enough admi's --
ICD1admi_mr_interadj <- glm.nb(ICD1_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD2admi_mr_interadj <- glm.nb(ICD2_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD3admi_mr_interadj <- glm.nb(ICD3_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD4admi_mr_interadj <- glm.nb(ICD4_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD5admi_mr_interadj <- glm.nb(ICD5_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD6admi_mr_interadj <- glm.nb(ICD6_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD7admi_mr_interadj <- glm.nb(ICD7_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD8admi_mr_interadj <- glm.nb(ICD8_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD9admi_mr_interadj <- glm.nb(ICD9_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD10admi_mr_interadj <- glm.nb(ICD10_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD11admi_mr_interadj <- glm.nb(ICD11_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD12admi_mr_interadj <- glm.nb(ICD12_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD13admi_mr_interadj <- glm.nb(ICD13_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD14admi_mr_interadj <- glm.nb(ICD14_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD17admi_mr_interadj <- glm.nb(ICD17_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD18admi_mr_interadj <- glm.nb(ICD18_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD19admi_mr_interadj <- glm.nb(ICD19_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
ICD21admi_mr_interadj <- glm.nb(ICD21_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + centre + townsend, data = mr_longhosp10y)
# use coeftest to get true significance levels (because clustering by eid's updates the SE's)
round((ICD1admi_mr_interadj_coef <- coeftest(ICD1admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD2admi_mr_interadj_coef <- coeftest(ICD2admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD3admi_mr_interadj_coef <- coeftest(ICD3admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD4admi_mr_interadj_coef <- coeftest(ICD4admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD5admi_mr_interadj_coef <- coeftest(ICD5admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD6admi_mr_interadj_coef <- coeftest(ICD6admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD7admi_mr_interadj_coef <- coeftest(ICD7admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD8admi_mr_interadj_coef <- coeftest(ICD8admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD9admi_mr_interadj_coef <- coeftest(ICD9admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD10admi_mr_interadj_coef <- coeftest(ICD10admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD11admi_mr_interadj_coef <- coeftest(ICD11admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD12dmi_mr_interadj_coef <- coeftest(ICD12admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD13admi_mr_interadj_coef <- coeftest(ICD13admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD14admi_mr_interadj_coef <- coeftest(ICD14admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD17admi_mr_interadj_coef <- coeftest(ICD17admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD18admi_mr_interadj_coef <- coeftest(ICD18admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD19admi_mr_interadj_coef <- coeftest(ICD19admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD21admi_mr_interadj_coef <- coeftest(ICD21admi_mr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
# exponentiated coefficients:
as.matrix(round(exp(coef(ICD1admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD2admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD3admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD4admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD5admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD6admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD7admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD8admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD9admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD10admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD11admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD12admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD13admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD14admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD17admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD18admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD19admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD21admi_mr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((ICD1admi_mr_interadj_cis <- coefci(ICD1admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD2admi_mr_interadj_cis <- coefci(ICD2admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD3admi_mr_interadj_cis <- coefci(ICD3admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD4admi_mr_interadj_cis <- coefci(ICD4admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD5admi_mr_interadj_cis <- coefci(ICD5admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD6admi_mr_interadj_cis <- coefci(ICD6admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD7admi_mr_interadj_cis <- coefci(ICD7admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD8admi_mr_interadj_cis <- coefci(ICD8admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD9admi_mr_interadj_cis <- coefci(ICD9admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD10admi_mr_interadj_cis <- coefci(ICD10admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD11admi_mr_interadj_cis <- coefci(ICD11admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD12admi_mr_interadj_cis <- coefci(ICD12admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD13admi_mr_interadj_cis <- coefci(ICD13admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD14admi_mr_interadj_cis <- coefci(ICD14admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD17admi_mr_interadj_cis <- coefci(ICD17admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD18admi_mr_interadj_cis <- coefci(ICD18admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD19admi_mr_interadj_cis <- coefci(ICD19admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD21admi_mr_interadj_cis <- coefci(ICD21admi_mr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)

# Maximally adjusted models: 
#  models not run for ICD 15, 16, 20 & 22 bc either no or not enough admi's --
ICD1admi_mr_maxadj <- glm.nb(ICD1_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD2admi_mr_maxadj <- glm.nb(ICD2_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD3admi_mr_maxadj <- glm.nb(ICD3_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD4admi_mr_maxadj <- glm.nb(ICD4_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD5admi_mr_maxadj <- glm.nb(ICD5_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD6admi_mr_maxadj <- glm.nb(ICD6_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD7admi_mr_maxadj <- glm.nb(ICD7_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD8admi_mr_maxadj <- glm.nb(ICD8_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD9admi_mr_maxadj <- glm.nb(ICD9_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD10admi_mr_maxadj <- glm.nb(ICD10_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD11admi_mr_maxadj <- glm.nb(ICD11_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD12admi_mr_maxadj <- glm.nb(ICD12_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD13admi_mr_maxadj <- glm.nb(ICD13_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD14admi_mr_maxadj <- glm.nb(ICD14_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD17admi_mr_maxadj <- glm.nb(ICD17_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y)
ICD18admi_mr_maxadj <- glm.nb(ICD18_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
ICD19admi_mr_maxadj <- glm.nb(ICD19_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y)
ICD21admi_mr_maxadj <- glm.nb(ICD21_admibyfuyear ~ asthma_mild_ONLY + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = mr_longhosp10y) 
# use coeftest to get true significance levels 
round((ICD1admi_mr_maxadj_coef <- coeftest(ICD1admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD2admi_mr_maxadj_coef <- coeftest(ICD2admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD3admi_mr_maxadj_coef <- coeftest(ICD3admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD4admi_mr_maxadj_coef <- coeftest(ICD4admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD5admi_mr_maxadj_coef <- coeftest(ICD5admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD6admi_mr_maxadj_coef <- coeftest(ICD6admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD7admi_mr_maxadj_coef <- coeftest(ICD7admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD8admi_mr_maxadj_coef <- coeftest(ICD8admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD9admi_mr_maxadj_coef <- coeftest(ICD9admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD10admi_mr_maxadj_coef <- coeftest(ICD10admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD11admi_mr_maxadj_coef <- coeftest(ICD11admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD12dmi_mr_maxadj_coef <- coeftest(ICD12admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD13admi_mr_maxadj_coef <- coeftest(ICD13admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD14admi_mr_maxadj_coef <- coeftest(ICD14admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD17admi_mr_maxadj_coef <- coeftest(ICD17admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD18admi_mr_maxadj_coef <- coeftest(ICD18admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD19admi_mr_maxadj_coef <- coeftest(ICD19admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD21admi_mr_maxadj_coef <- coeftest(ICD21admi_mr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
# exponentiated coefficients:
as.matrix(round(exp(coef(ICD1admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD2admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD3admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD4admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD5admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD6admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD7admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD8admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD9admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD10admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD11admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD12admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD13admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD14admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD17admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD18admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD19admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD21admi_mr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((ICD1admi_mr_maxadj_cis <- coefci(ICD1admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD2admi_mr_maxadj_cis <- coefci(ICD2admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD3admi_mr_maxadj_cis <- coefci(ICD3admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD4admi_mr_maxadj_cis <- coefci(ICD4admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD5admi_mr_maxadj_cis <- coefci(ICD5admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD6admi_mr_maxadj_cis <- coefci(ICD6admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD7admi_mr_maxadj_cis <- coefci(ICD7admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD8admi_mr_maxadj_cis <- coefci(ICD8admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD9admi_mr_maxadj_cis <- coefci(ICD9admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD10admi_mr_maxadj_cis <- coefci(ICD10admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD11admi_mr_maxadj_cis <- coefci(ICD11admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD12admi_mr_maxadj_cis <- coefci(ICD12admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD13admi_mr_maxadj_cis <- coefci(ICD13admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD14admi_mr_maxadj_cis <- coefci(ICD14admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD17admi_mr_maxadj_cis <- coefci(ICD17admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD18admi_mr_maxadj_cis <- coefci(ICD18admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD19admi_mr_maxadj_cis <- coefci(ICD19admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD21admi_mr_maxadj_cis <- coefci(ICD21admi_mr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)




###~ Severe vs no asthma ###########################################

# (factor variables first)
sr_longhosp10y$ethnicity <- factor(sr_longhosp10y$ethnicity, ordered = FALSE, levels = c("White", "Black", "South Asian", "Others"))
sr_longhosp10y$smoke <- factor(sr_longhosp10y$smoke, ordered = FALSE, levels = c("Never", "Current", "Previous"))
sr_longhosp10y$BMI_cat <- factor(sr_longhosp10y$BMI_cat, ordered = FALSE, levels = c("18.5-25", "<18.5", "25-30", "30-35", "35-40", "40+"))
sr_longhosp10y$townsend <- factor(sr_longhosp10y$townsend, ordered = FALSE, levels = c(1, 2, 3, 4, 5))
sr_longhosp10y$centre <- factor(sr_longhosp10y$centre, ordered = FALSE, levels = c(11010, 11002, 10003, 11001, 11003, 11004, 11005, 11006, 11007, 11008, 11009, 11011, 11012, 11013, 11014, 11016, 11017, 11018, 11020, 11021, 11022, 11023))
# centre age:
sr_longhosp10y$age.centred <- scale(sr_longhosp10y$age.recruit, center = TRUE, scale = FALSE)

####### ~~~ MAIN MODELS #######################################
# (without fu_year interaction) #

# Admissions first:

admi1_sr_x <- glm.nb(admi_n ~ diagnosis_severe_med + fu_year, data = sr_longhosp10y) # NON-adjusted

admi1_sr_adj_x <- glm.nb(admi_n ~ diagnosis_severe_med + fu_year + age.centred + # MINIMALLY-adjusted
                           sex + ethnicity + centre, 
                         data = sr_longhosp10y)

admi1_sr_maxadj_x <- glm.nb(admi_n ~ diagnosis_severe_med + fu_year + age.centred # MAXIMALLY-adjusted
                            + sex + ethnicity + townsend + centre
                            + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                            + hypertension + MI.baseline + stroke.baseline 
                            + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                            + PVD + mental, 
                            data = sr_longhosp10y)

# [INTER-ADJUSTED]:
# adjust for all 4 matching covariates AND townsend ~ by comparing these models to 
# minimally adjusted models we can estimate the contribution of deprivation
admi1_sr_interadj_x <- glm.nb(admi_n ~ diagnosis_severe_med + fu_year + age.centred # INTERMEDIATE-adjusted
                              + sex + ethnicity + townsend + centre, 
                              data = sr_longhosp10y)

## Now length of stay: ##

days1_sr_x <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med + fu_year, data = sr_longhosp10y)

days1_sr_adj_x <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med + fu_year + age.centred
                         + sex + ethnicity + centre, 
                         data = sr_longhosp10y)

days1_sr_maxadj_x <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med + fu_year + age.centred
                            + sex + ethnicity + townsend + centre 
                            + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                            + hypertension + MI.baseline + stroke.baseline 
                            + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                            + PVD + mental,
                            data = sr_longhosp10y)
# [INTER-ADJUSTED]:
days1_sr_interadj_x <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med + fu_year + age.centred
                              + sex + ethnicity + townsend + centre, 
                              data = sr_longhosp10y)

## now COSTS: ##

costs1_sr_x <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med + fu_year, data = sr_longhosp10y)

costs1_sr_adj_x <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med + fu_year + age.centred
                          + sex + ethnicity + centre, 
                          data = sr_longhosp10y)

costs1_sr_maxadj_x <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med + fu_year + age.centred
                             + sex + ethnicity + townsend + centre
                             + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                             + hypertension + MI.baseline + stroke.baseline 
                             + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                             + PVD + mental,
                             data = sr_longhosp10y)

# [INTER-ADJUSTED]:
costs1_sr_interadj_x <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med + fu_year + age.centred
                               + sex + ethnicity + townsend + centre, 
                               data = sr_longhosp10y)

# ---- CLUSTERED - STANDARD - ERRORS ---- ! #
# ---- CLUSTERED - CONFIDENCE - INTERVALS ---- ! #

# MOD-SEVERE vs no-asthma models:
library(lmtest)
library(sandwich)
library(boot)

# COEF TEST (standard errors are clustered by eid) to obtain TRUE significance level:
(admi1_sr_x_coef <- coeftest(admi1_sr_x, vcov. = vcovCL, cluster = ~eid))
(admi1_sr_adj_x_coef <- coeftest(admi1_sr_adj_x, vcov. = vcovCL, cluster = ~eid))
(admi1_sr_interadj_x_coef <- coeftest(admi1_sr_interadj_x, vcov. = vcovCL, cluster = ~eid))
(admi1_sr_maxadj_x_coef <- coeftest(admi1_sr_maxadj_x, vcov. = vcovCL, cluster = ~eid))

(days1_sr_x_coef <- coeftest(days1_sr_x, vcov. = vcovCL, cluster = ~eid))
(days1_sr_adj_x_coef <- coeftest(days1_sr_adj_x, vcov. = vcovCL, cluster = ~eid))
(days1_sr_interadj_x_coef <- coeftest(days1_sr_interadj_x, vcov. = vcovCL, cluster = ~eid))
(days1_sr_maxadj_x_coef <- coeftest(days1_sr_maxadj_x, vcov. = vcovCL, cluster = ~eid))

(costs2_sr_x_coef <- coeftest(costs2_sr_x, vcov. = vcovCL, cluster = ~eid))
(costs1_sr_adj_x_coef <- coeftest(costs1_sr_adj_x, vcov. = vcovCL, cluster = ~eid))
(costs1_sr_interadj_x_coef <- coeftest(costs1_sr_interadj_x, vcov. = vcovCL, cluster = ~eid))
(costs1_sr_maxadj_x_coef <- coeftest(costs1_sr_maxadj_x, vcov. = vcovCL, cluster = ~eid))

# EXPONENTIATED COEFFICIENTS:  (this converts coefficients to Incident Rate Ratios)
as.matrix(round(exp(coef(admi1_sr_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(admi1_sr_adj_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(admi1_sr_interadj_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(admi1_sr_maxadj_x)), digits = 2), ncol = 2, rownames = TRUE)

as.matrix(round(exp(coef(days1_sr_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_sr_adj_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_sr_interadj_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_sr_maxadj_x)), digits = 2), ncol = 2, rownames = TRUE)

as.matrix(round(exp(coef(costs2_sr_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_sr_adj_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_sr_interadj_x)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_sr_maxadj_x)), digits = 2), ncol = 2, rownames = TRUE)

# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((admi1_sr_x_cis <- coefci(admi1_sr_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((admi1_sr_adj_x_cis <- coefci(admi1_sr_adj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((admi1_sr_interadj_x_cis <- coefci(admi1_sr_interadj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((admi1_sr_maxadj_x_cis <- coefci(admi1_sr_maxadj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)

round(exp((days1_sr_x_cis <- coefci(days1_sr_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days1_sr_adj_x_cis <- coefci(days1_sr_adj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days1_sr_interadj_x_cis <- coefci(days1_sr_interadj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days1_sr_maxadj_x_cis <- coefci(days1_sr_maxadj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)

round(exp((costs2_sr_x_cis <- coefci(costs2_sr_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((costs1_sr_adj_x_cis <- coefci(costs1_sr_adj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((costs1_sr_interadj_x_cis <- coefci(costs1_sr_interadj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((costs1_sr_maxadj_x_cis <- coefci(costs1_sr_maxadj_x, vcov = vcovCL, cluster = ~eid))), digits = 2)



####### ~~~ check fu_year interaction #######################################

# Compare models with and without interaction for each level of adjustment,
# then compare USING LIKELIHOOD RATIO TESTS

## 1 - unadjusted models:

# admissions
admi1_sr <- glm.nb(admi_n ~ diagnosis_severe_med * fu_year, data = sr_longhosp10y)
screenreg(admi1_sr)
# days
days1_sr <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * fu_year, data = sr_longhosp10y)
screenreg(m1days_sr)
# costs
costs1_sr <- glm.cluster(formula = totcost_byfuyear ~ diagnosis_severe_med * fu_year,
                         cluster ="eid", 
                         family = gaussian(link = "identity"),
                         data = sr_longhosp10y)

costs2_sr <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * fu_year, data = sr_longhosp10y)
screenreg(list(costs1_sr, costs2_sr))
# again, negative binomial outperforms
lrtest(costs2_sr, costs1_sr_x)
# p value is NOT significant
# thus the model WITHOUT interaction is better 


## 2 - minimally adjusted models:
# include the matched covariates (NO TOWNSEND)

# negative binomial models for admissions, days, and costs:

admi1_sr_adj <- glm.nb(admi_n ~ diagnosis_severe_med * fu_year + age.centred
                       + sex + ethnicity + centre, 
                       data = sr_longhosp10y)
screenreg(admi1_sr_adj)


days1_sr_adj <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * fu_year + age.centred
                       + sex + ethnicity + centre, 
                       data = sr_longhosp10y)
screenreg(days1_sr_adj)


costs1_sr_adj <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * fu_year + age.centred
                        + sex + ethnicity + centre, 
                        data = sr_longhosp10y)
screenreg(costs1_sr_adj)
screenreg(list(admi1_sr_adj, days1_sr_adj, costs1_sr_adj))


## 3 - maximally adjusted models:
# include the matched covariates AND townsend, BMI, smoking, and comorbidities

# negative binomial model for admissions
admi1_sr_maxadj <- glm.nb(admi_n ~ diagnosis_severe_med * fu_year + age.centred
                          + sex + ethnicity + townsend + centre
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental, 
                          data = sr_longhosp10y)

# compare all 3 models of admissions: [all neg. binomial]
screenreg(list(admi1_sr, admi1_sr_adj, admi1_sr_maxadj))


# negative binomial model for days
days1_sr_maxadj <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * fu_year + age.centred
                          + sex + ethnicity + townsend + centre
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental, 
                          data = sr_longhosp10y)

# compare all 3 models of admissions: [all neg. binomial]
screenreg(list(days1_sr, days1_sr_adj, days1_sr_maxadj))


# negative binomial model for costs
costs1_sr_maxadj <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * fu_year + age.centred
                           + sex + ethnicity + townsend + centre
                           + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                           + hypertension + MI.baseline + stroke.baseline 
                           + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                           + PVD + mental, 
                           data = sr_longhosp10y)

# compare all 3 models of admissions: [all neg. binomial]
screenreg(list(costs1_sr, costs1_sr_adj, costs1_sr_maxadj))


## 4 - INTERMEDIATE adjusted models:

admi1_sr_interadj <- glm.nb(admi_n ~ diagnosis_severe_med * fu_year + age.centred
                            + sex + ethnicity + centre + townsend, 
                            data = sr_longhosp10y)
screenreg(admi1_sr_interadj)


days1_sr_interadj <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * fu_year + age.centred
                            + sex + ethnicity + centre + townsend, 
                            data = sr_longhosp10y)
screenreg(days1_sr_interadj)


costs1_sr_interadj <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * fu_year + age.centred
                             + sex + ethnicity + centre + townsend, 
                             data = sr_longhosp10y)
screenreg(costs1_sr_interadj)



####### ~~~ check deprivation interaction #######################################

a_sr_dep <- glm.nb(admi_n ~ diagnosis_severe_med * townsend, data = sr_longhosp10y)
a_sr_adj_dep <- glm.nb(admi_n ~ diagnosis_severe_med * townsend + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
a_sr_interadj_dep <- glm.nb(admi_n ~ diagnosis_severe_med * townsend + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
a_sr_maxadj_dep <- glm.nb(admi_n ~ diagnosis_severe_med * townsend + age.centred
                          + sex + ethnicity + townsend + centre
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental, 
                          data = sr_longhosp10y)
# lrtest comparing with and without townsend (categorical) interaction:
lrtest(admi1_sr, a_sr_dep) # significant *** = WITHOUT townsend is better
lrtest(admi1_sr_adj, a_sr_adj_dep) # significant *** = WITHOUT townsend is better
lrtest(admi1_sr_maxadj, a_sr_maxadj_dep) # significant *** = WITHOUT townsend is better


h_sr_dep <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * townsend, data = sr_longhosp10y)
h_sr_adj_dep <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * townsend + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
h_sr_interadj_dep <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * townsend + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
h_sr_maxadj_dep <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * townsend + age.centred
                          + sex + ethnicity + townsend + centre
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental, 
                          data = sr_longhosp10y)
# lrtest comparing with and without townsend (categorical) interaction:
lrtest(days1_sr, h_sr_dep) # *** significant
lrtest(days1_sr_adj, h_sr_adj_dep) # *** significant
lrtest(days1_sr_maxadj, h_sr_maxadj_dep) # significant *** = WITHOUT townsend is better

c_sr_dep <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * townsend, data = sr_longhosp10y)
c_sr_adj_dep <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * townsend + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
c_sr_interadj_dep <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * townsend + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
c_sr_maxadj_dep <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * townsend + age.centred
                          + sex + ethnicity + townsend + centre
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental, 
                          data = sr_longhosp10y)
# lrtest comparing with and without townsend (categorical) interaction:
lrtest(costs2_sr, c_sr_dep) # significant *** = WITHOUT townsend is better
lrtest(costs1_sr_adj, c_sr_adj_dep) # significant *** = WITHOUT townsend is better
lrtest(costs1_sr_maxadj, c_sr_maxadj_dep) # significant *** = WITHOUT townsend is better


# cLUSTERED STANDARD ERRORS:
(a_sr_dep_coef <- coeftest(a_sr_dep, vcov. = vcovCL, cluster = ~eid))
(a_sr_adj_dep_coef <- coeftest(a_sr_adj_dep, vcov. = vcovCL, cluster = ~eid))
(a_sr_interadj_dep_coef <- coeftest(a_sr_interadj_dep, vcov. = vcovCL, cluster = ~eid))
(a_sr_maxadj_dep_coef <- coeftest(a_sr_maxadj_dep, vcov. = vcovCL, cluster = ~eid))

(h_sr_dep_coef <- coeftest(h_sr_dep, vcov. = vcovCL, cluster = ~eid))
(h_sr_adj_dep_coef <- coeftest(h_sr_adj_dep, vcov. = vcovCL, cluster = ~eid))
(h_sr_interadj_dep_coef <- coeftest(h_sr_interadj_dep, vcov. = vcovCL, cluster = ~eid))
(h_sr_maxadj_dep_coef <- coeftest(h_sr_maxadj_dep, vcov. = vcovCL, cluster = ~eid))

(c_sr_dep_coef <- coeftest(c_sr_dep, vcov. = vcovCL, cluster = ~eid))
(c_sr_adj_dep_coef <- coeftest(c_sr_adj_dep, vcov. = vcovCL, cluster = ~eid))
(c_sr_interadj_dep_coef <- coeftest(c_sr_interadj_dep, vcov. = vcovCL, cluster = ~eid))
(c_sr_maxadj_dep_coef <- coeftest(c_sr_maxadj_dep, vcov. = vcovCL, cluster = ~eid))





######## ~~~ continuous deprivation interaction #######################################

# convert townsend from factored into numeric variable:

sr_longhosp10y$new_townsend <- as.numeric(sr_longhosp10y$townsend)

# test all model specifications except for inter-adjusted (because we already include townsend as interaction term)

# - admissions - 
a_sr_depcont <- glm.nb(admi_n ~ diagnosis_severe_med * new_townsend, data = sr_longhosp10y)
a_sr_adj_depcont <- glm.nb(admi_n ~ diagnosis_severe_med * new_townsend + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
a_sr_interadj_depcont <- glm.nb(admi_n ~ diagnosis_severe_med * new_townsend + age.centred + sex + ethnicity + centre + new_townsend, data = sr_longhosp10y)
a_sr_maxadj_depcont <- glm.nb(admi_n ~ diagnosis_severe_med * new_townsend + age.centred
                              + sex + ethnicity + townsend + centre
                              + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                              + hypertension + MI.baseline + stroke.baseline 
                              + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                              + PVD + mental, data = sr_longhosp10y)
#screenreg to get log likelihoods:
screenreg(list(a_sr_depcont, a_sr_adj_depcont,a_sr_maxadj_depcont))
# lrtest comparing with and without townsend (numeric/ continuous) interaction:
lrtest(admi1_sr, a_sr_depcont) # *** = WITHOUT townsend interaction is better.
lrtest(admi1_sr_adj, a_sr_adj_depcont) # *** = WITHOUT townsend interaction is better.
lrtest(admi1_sr_maxadj, a_sr_maxadj_depcont) # significant *** = WITHOUT townsend is better

# cLUSTERED STANDARD ERRORS:
(a_sr_depcont_coef <- coeftest(a_sr_depcont, vcov. = vcovCL, cluster = ~eid))
(a_sr_adj_depcont_coef <- coeftest(a_sr_adj_depcont, vcov. = vcovCL, cluster = ~eid))
(a_sr_interadj_depcont_coef <- coeftest(a_sr_interadj_depcont, vcov. = vcovCL, cluster = ~eid))
(a_sr_maxadj_depcont_coef <- coeftest(a_sr_maxadj_depcont, vcov. = vcovCL, cluster = ~eid))
# EXPONENTIATED COEFFICIENTS:
as.matrix(round(exp(coef(a_sr_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(a_sr_adj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(a_sr_interadj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(a_sr_maxadj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((a_sr_depcont_cis <- coefci(a_sr_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((a_sr_adj_depcont_cis <- coefci(a_sr_adj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((a_sr_interadj_depcont_cis <- coefci(a_sr_interadj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((a_sr_maxadj_depcont_cis <- coefci(a_sr_maxadj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)



# - days - 
d_sr_depcont <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * new_townsend, data = sr_longhosp10y)
d_sr_adj_depcont <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * new_townsend + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
d_sr_interadj_depcont <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * new_townsend + age.centred + sex + ethnicity + centre + new_townsend, data = sr_longhosp10y)
d_sr_maxadj_depcont <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * new_townsend + age.centred
                              + sex + ethnicity + townsend + centre
                              + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                              + hypertension + MI.baseline + stroke.baseline 
                              + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                              + PVD + mental, data = sr_longhosp10y)
#screenreg to get log likelihoods:
screenreg(list(d_sr_depcont, d_sr_adj_depcont,d_sr_maxadj_depcont))
# lrtest comparing with and without townsend (numeric/ continuous) interaction:
lrtest(days1_sr, d_sr_depcont) 
lrtest(days1_sr_adj, d_sr_adj_depcont) 
lrtest(days1_sr_maxadj, d_sr_maxadj_depcont) # significant *** = WITHOUT townsend is better

# cLUSTERED STANDARD ERRORS:
(d_sr_depcont_coef <- coeftest(d_sr_depcont, vcov. = vcovCL, cluster = ~eid))
(d_sr_adj_depcont_coef <- coeftest(d_sr_adj_depcont, vcov. = vcovCL, cluster = ~eid))
(d_sr_interadj_depcont_coef <- coeftest(d_sr_interadj_depcont, vcov. = vcovCL, cluster = ~eid))
(d_sr_maxadj_depcont_coef <- coeftest(d_sr_maxadj_depcont, vcov. = vcovCL, cluster = ~eid))
# EXPONENTIATED COEFFICIENTS:
as.matrix(round(exp(coef(d_sr_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(d_sr_adj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(d_sr_interadj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(d_sr_maxadj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((d_sr_depcont_cis <- coefci(d_sr_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((d_sr_adj_depcont_cis <- coefci(d_sr_adj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((d_sr_interadj_depcont_cis <- coefci(d_sr_interadj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((d_sr_maxadj_depcont_cis <- coefci(d_sr_maxadj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)



# - costs - 
c_sr_depcont <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * new_townsend, data = sr_longhosp10y)
c_sr_adj_depcont <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * new_townsend + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
c_sr_interadj_depcont <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * new_townsend + age.centred + sex + ethnicity + centre + new_townsend, data = sr_longhosp10y)
c_sr_maxadj_depcont <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * new_townsend + age.centred
                              + sex + ethnicity + townsend + centre
                              + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                              + hypertension + MI.baseline + stroke.baseline 
                              + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                              + PVD + mental, data = sr_longhosp10y)
#screenreg to get log likelihoods:
screenreg(list(c_sr_depcont, c_sr_adj_depcont,c_sr_maxadj_depcont))
# lrtest comparing with and without townsend (numeric/ continuous) interaction:
lrtest(costs2_sr, c_sr_depcont) # significant *** = WITHOUT townsend is better
lrtest(costs1_sr_adj, c_sr_adj_depcont) # significant *** = WITHOUT townsend is better
lrtest(costs1_sr_maxadj, c_sr_maxadj_depcont) # significant *** = WITHOUT townsend is better

# cLUSTERED STANDARD ERRORS:
(c_sr_depcont_coef <- coeftest(c_sr_depcont, vcov. = vcovCL, cluster = ~eid))
(c_sr_adj_depcont_coef <- coeftest(c_sr_adj_depcont, vcov. = vcovCL, cluster = ~eid))
(c_sr_interadj_depcont_coef <- coeftest(c_sr_interadj_depcont, vcov. = vcovCL, cluster = ~eid))
(c_sr_maxadj_depcont_coef <- coeftest(c_sr_maxadj_depcont, vcov. = vcovCL, cluster = ~eid))
# EXPONENTIATED COEFFICIENTS:
as.matrix(round(exp(coef(c_sr_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(c_sr_adj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(c_sr_interadj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(c_sr_maxadj_depcont)), digits = 2), ncol = 2, rownames = TRUE)
# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((c_sr_depcont_cis <- coefci(c_sr_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((c_sr_adj_depcont_cis <- coefci(c_sr_adj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((c_sr_interadj_depcont_cis <- coefci(c_sr_interadj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((c_sr_maxadj_depcont_cis <- coefci(c_sr_maxadj_depcont, vcov = vcovCL, cluster = ~eid))), digits = 2)



############ ~~~ by continuous fu_year ###########################################

# first create new continuous variable for follow-up year
sr_longhosp10y <- sr_longhosp10y %>% 
  mutate(new_fu_year = case_when(sr_longhosp10y$fu_year == "year_1" ~ 1,
                                 sr_longhosp10y$fu_year == "year_2" ~ 2,
                                 sr_longhosp10y$fu_year == "year_3" ~ 3,
                                 sr_longhosp10y$fu_year == "year_4" ~ 4,
                                 sr_longhosp10y$fu_year == "year_5" ~ 5,
                                 sr_longhosp10y$fu_year == "year_6" ~ 6,
                                 sr_longhosp10y$fu_year == "year_7" ~ 7,
                                 sr_longhosp10y$fu_year == "year_8" ~ 8,
                                 sr_longhosp10y$fu_year == "year_9" ~ 9,
                                 sr_longhosp10y$fu_year == "year_10" ~ 10))

## Model follow-up year as a single continuous variable for all outcomes 
## using the 'new_fu_year' variable:
# admissions:

admi_sr_n <- glm.nb(admi_n ~ diagnosis_severe_med * new_fu_year, data = sr_longhosp10y)

admi_sr_adj_n <- glm.nb(admi_n ~ diagnosis_severe_med * new_fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)

admi_sr_interadj_n <- glm.nb(admi_n ~ diagnosis_severe_med * new_fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)

admi_sr_maxadj_n <- glm.nb(admi_n ~ diagnosis_severe_med * new_fu_year + age.centred
                           + sex + ethnicity + townsend + centre
                           + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                           + hypertension + MI.baseline + stroke.baseline 
                           + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                           + PVD + mental, 
                           data = sr_longhosp10y)

screenreg(list(admi_sr_n, admi_sr_adj_n, admi_sr_interadj_n, admi_sr_maxadj_n))


# days: #
days_sr_n <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * new_fu_year, data = sr_longhosp10y)

days_sr_adj_n <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * 
                          new_fu_year + age.centred + sex 
                        + ethnicity + centre, data = sr_longhosp10y)

days_sr_interadj_n <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * new_fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)

days_sr_maxadj_n <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * new_fu_year + age.centred
                           + sex + ethnicity + townsend + centre
                           + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                           + hypertension + MI.baseline + stroke.baseline 
                           + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                           + PVD + mental, 
                           data = sr_longhosp10y)

screenreg(list(days_sr_n, days_sr_adj_n, days_sr_interadj_n, days_sr_maxadj_n))


# costs: #
costs_sr_n <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * new_fu_year, data = sr_longhosp10y)

costs_sr_adj_n <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * 
                           new_fu_year + age.centred + sex 
                         + ethnicity + centre, data = sr_longhosp10y)

costs_sr_interadj_n <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * new_fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)


costs_sr_maxadj_n <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * new_fu_year + age.centred
                            + sex + ethnicity + townsend + centre
                            + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                            + hypertension + MI.baseline + stroke.baseline 
                            + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                            + PVD + mental, 
                            data = sr_longhosp10y)

screenreg(list(costs_sr_n, costs_sr_adj_n, costs_sr_interadj_n, costs_sr_maxadj_n))

# cLUSTERED STANDARD ERRORS:
(admi_sr_n_coef <- coeftest(admi_sr_n, vcov. = vcovCL, cluster = ~eid))
(admi_sr_adj_n_coef <- coeftest(admi_sr_adj_n, vcov. = vcovCL, cluster = ~eid))
(admi_sr_interadj_n_coef <- coeftest(admi_sr_interadj_n, vcov. = vcovCL, cluster = ~eid))
(admi_sr_maxadj_n_coef <- coeftest(admi_sr_maxadj_n, vcov. = vcovCL, cluster = ~eid))

(days_sr_n_coef <- coeftest(days_sr_n, vcov. = vcovCL, cluster = ~eid))
(days_sr_adj_n_coef <- coeftest(days_sr_adj_n, vcov. = vcovCL, cluster = ~eid))
(days_sr_interadj_n_coef <- coeftest(days_sr_interadj_n, vcov. = vcovCL, cluster = ~eid))
(days_sr_maxadj_n_coef <- coeftest(days_sr_maxadj_n, vcov. = vcovCL, cluster = ~eid))

(costs_sr_n_coef <- coeftest(costs_sr_n, vcov. = vcovCL, cluster = ~eid))
(costs_sr_adj_n_coef <- coeftest(costs_sr_adj_n, vcov. = vcovCL, cluster = ~eid))
(costs_sr_interadj_n_coef <- coeftest(costs_sr_interadj_n, vcov. = vcovCL, cluster = ~eid))
(costs_sr_maxadj_n_coef <- coeftest(costs_sr_maxadj_n, vcov. = vcovCL, cluster = ~eid))



#Likelihood ratio tests comparing continuous new_fu_year interaction with NO fu_year interaction!!:
lrtest(admi1_sr_x, admi_sr_n)    # NON-SIGNIFICANT = WITHOUT interaction is better
lrtest(admi1_sr_adj_x, admi_sr_adj_n) # NON-SIGNIFICANT = WITHOUT interaction is better
lrtest(admi1_sr_interadj_x, admi_sr_interadj_n) # NON-SIGNIFICANT = WITHOUT interaction is better
lrtest(admi1_sr_maxadj_x, admi_sr_maxadj_n) # (p<0.01) ** 

lrtest(days1_sr_x, days_sr_n) # SIGNIFICANT (p<0.001) *** non-interaction model has larger LogLik = WITHOUT interaction is better
lrtest(days1_sr_adj_x, days_sr_adj_n) # ""  ""
lrtest(days1_sr_interadj_x, days_sr_interadj_n) # ""  ""
lrtest(days1_sr_maxadj_x, days_sr_maxadj_n) # ""  ""

lrtest(costs1_sr_x, costs_sr_n)    # NO significance.
lrtest(costs1_sr_adj_x, costs_sr_adj_n) # NO significance.
lrtest(costs1_sr_interadj_x, costs_sr_interadj_n) # NO significance.
lrtest(costs1_sr_maxadj_x, costs_sr_maxadj_n) # NO significance.


#Likelihood ratio tests comparing categorical fu_year & continuous new_fu_year:
lrtest(admi1_sr, admi_sr_n)    # (p<0.01) ** = categorical fu_year = better
lrtest(admi1_sr_adj, admi_sr_adj_n) # (p<0.01) ** = categorical fu_year = better
lrtest(admi1_sr_interadj, admi_sr_interadj_n) # (p<0.01) ** = categorical fu_year = better
lrtest(admi1_sr_maxadj, admi_sr_maxadj_n) # (p<0.01) ** = categorical fu_year = better

lrtest(days1_sr, days_sr_n)    # (p<0.001) *** = categorical fu_year = better
lrtest(days1_sr_adj, days_sr_adj_n) # (p<0.001) *** = categorical fu_year = better
lrtest(days1_sr_interadj, days_sr_interadj_n) # (p<0.001) *** = categorical fu_year = better
lrtest(days1_sr_maxadj, days_sr_maxadj_n) # (p<0.001) *** = categorical fu_year = better

lrtest(costs2_sr, costs_sr_n)    # NO significance. Tiny difference in LogLik scores.
lrtest(costs1_sr_adj, costs_sr_adj_n) # NO significance.
lrtest(costs1_sr_interadj, costs_sr_interadj_n) # NO significance.
lrtest(costs1_sr_maxadj, costs_sr_maxadj_n) # NO significance.




############# ~~~ by RESPIRATORY admissions ###########################################

respadmi_sr <- glm.nb(respadmi_byfuyear ~ diagnosis_severe_med * fu_year, data = sr_longhosp10y)

respadmi_sr_adj <- glm.nb(respadmi_byfuyear ~ diagnosis_severe_med * fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)

respadmi_sr_interadj <- glm.nb(respadmi_byfuyear ~ diagnosis_severe_med * fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)

respadmi_sr_maxadj <- glm.nb(respadmi_byfuyear ~ diagnosis_severe_med * fu_year + age.centred
                             + sex + ethnicity + townsend + centre
                             + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                             + hypertension + MI.baseline + stroke.baseline 
                             + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                             + PVD + mental, 
                             data = sr_longhosp10y)

# screenreg to get Log likelihoods:
screenreg(list(respadmi_sr, respadmi_sr_adj, respadmi_sr_interadj, respadmi_sr_maxadj))

install.packages("lmtest")
library(lmtest)
library(sandwich)
library(boot)

# use coeftest to get true significance levels (because clustering by eid's updates the SE's)
round((respadmi_sr_coef <- coeftest(respadmi_sr, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((respadmi_sr_adj_coef <- coeftest(respadmi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((respadmi_sr_adj_interadj_coef <- coeftest(respadmi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((respadmi_sr_adj_maxadj_coef <- coeftest(respadmi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
# exponentiated coefficients:
as.matrix(round(exp(coef(respadmi_sr)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(respadmi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(respadmi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(respadmi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((respadmi_sr_cis <- coefci(respadmi_sr, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((respadmi_sr_adj_cis <- coefci(respadmi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((respadmi_sr_interadj_cis <- coefci(respadmi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((respadmi_sr_maxadj_cis <- coefci(respadmi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)

############## ~~~ by other ICD admissions ###########################################

# All miminally adjusted models: 
ICD1admi_sr_adj <- glm.nb(ICD1_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD2admi_sr_adj <- glm.nb(ICD2_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD3admi_sr_adj <- glm.nb(ICD3_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD4admi_sr_adj <- glm.nb(ICD4_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD5admi_sr_adj <- glm.nb(ICD5_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD6admi_sr_adj <- glm.nb(ICD6_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD7admi_sr_adj <- glm.nb(ICD7_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD8admi_sr_adj <- glm.nb(ICD8_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD9admi_sr_adj <- glm.nb(ICD9_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD10admi_sr_adj <- glm.nb(ICD10_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD11admi_sr_adj <- glm.nb(ICD11_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD12admi_sr_adj <- glm.nb(ICD12_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD13admi_sr_adj <- glm.nb(ICD13_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD14admi_sr_adj <- glm.nb(ICD14_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD15admi_sr_adj <- glm.nb(ICD15_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD16admi_sr_adj <- glm.nb(ICD16_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD17admi_sr_adj <- glm.nb(ICD17_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD18admi_sr_adj <- glm.nb(ICD18_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD19admi_sr_adj <- glm.nb(ICD19_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD20admi_sr_adj <- glm.nb(ICD20_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
ICD21admi_sr_adj <- glm.nb(ICD21_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre, data = sr_longhosp10y)
# use coeftest to get true significance levels
round((ICD1admi_sr_adj_coef <- coeftest(ICD1admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD2admi_sr_adj_coef <- coeftest(ICD2admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD3admi_sr_adj_coef <- coeftest(ICD3admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD4admi_sr_adj_coef <- coeftest(ICD4admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD5admi_sr_adj_coef <- coeftest(ICD5admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD6admi_sr_adj_coef <- coeftest(ICD6admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD7admi_sr_adj_coef <- coeftest(ICD7admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD8admi_sr_adj_coef <- coeftest(ICD8admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD9admi_sr_adj_coef <- coeftest(ICD9admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD10admi_sr_adj_coef <- coeftest(ICD10admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD11admi_sr_adj_coef <- coeftest(ICD11admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD12dmi_sr_adj_coef <- coeftest(ICD12admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD13admi_sr_adj_coef <- coeftest(ICD13admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD14admi_sr_adj_coef <- coeftest(ICD14admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD15admi_sr_adj_coef <- coeftest(ICD15admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD16admi_sr_adj_coef <- coeftest(ICD16admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD17admi_sr_adj_coef <- coeftest(ICD17admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD18admi_sr_adj_coef <- coeftest(ICD18admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD19admi_sr_adj_coef <- coeftest(ICD19admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD20admi_sr_adj_coef <- coeftest(ICD20admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD21admi_sr_adj_coef <- coeftest(ICD21admi_sr_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
# exponentiated coefficients:
as.matrix(round(exp(coef(ICD1admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD2admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD3admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD4admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD5admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD6admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD7admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD8admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD9admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD10admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD11admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD12admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD13admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD14admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD15admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD16admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD17admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD18admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD19admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD20admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD21admi_sr_adj)), digits = 2), ncol = 2, rownames = TRUE)
# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((ICD1admi_sr_adj_cis <- coefci(ICD1admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD2admi_sr_adj_cis <- coefci(ICD2admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD3admi_sr_adj_cis <- coefci(ICD3admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD4admi_sr_adj_cis <- coefci(ICD4admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD5admi_sr_adj_cis <- coefci(ICD5admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD6admi_sr_adj_cis <- coefci(ICD6admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD7admi_sr_adj_cis <- coefci(ICD7admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD8admi_sr_adj_cis <- coefci(ICD8admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD9admi_sr_adj_cis <- coefci(ICD9admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD10admi_sr_adj_cis <- coefci(ICD10admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD11admi_sr_adj_cis <- coefci(ICD11admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD12admi_sr_adj_cis <- coefci(ICD12admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD13admi_sr_adj_cis <- coefci(ICD13admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD14admi_sr_adj_cis <- coefci(ICD14admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD15admi_sr_adj_cis <- coefci(ICD15admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD16admi_sr_adj_cis <- coefci(ICD16admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD17admi_sr_adj_cis <- coefci(ICD17admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD18admi_sr_adj_cis <- coefci(ICD18admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD19admi_sr_adj_cis <- coefci(ICD19admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD20admi_sr_adj_cis <- coefci(ICD20admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD21admi_sr_adj_cis <- coefci(ICD21admi_sr_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)

# All intermediate adjusted models: 
#  models not run for ICD 15, 16, 20 & 22 bc either no or not enough admi's --
ICD1admi_sr_interadj <- glm.nb(ICD1_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD2admi_sr_interadj <- glm.nb(ICD2_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD3admi_sr_interadj <- glm.nb(ICD3_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD4admi_sr_interadj <- glm.nb(ICD4_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD5admi_sr_interadj <- glm.nb(ICD5_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD6admi_sr_interadj <- glm.nb(ICD6_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD7admi_sr_interadj <- glm.nb(ICD7_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD8admi_sr_interadj <- glm.nb(ICD8_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD9admi_sr_interadj <- glm.nb(ICD9_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD10admi_sr_interadj <- glm.nb(ICD10_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD11admi_sr_interadj <- glm.nb(ICD11_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD12admi_sr_interadj <- glm.nb(ICD12_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD13admi_sr_interadj <- glm.nb(ICD13_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD14admi_sr_interadj <- glm.nb(ICD14_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD17admi_sr_interadj <- glm.nb(ICD17_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD18admi_sr_interadj <- glm.nb(ICD18_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD19admi_sr_interadj <- glm.nb(ICD19_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
ICD21admi_sr_interadj <- glm.nb(ICD21_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + centre + townsend, data = sr_longhosp10y)
# use coeftest to get true significance levels (because clustering by eid's updates the SE's)
round((ICD1admi_sr_interadj_coef <- coeftest(ICD1admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD2admi_sr_interadj_coef <- coeftest(ICD2admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD3admi_sr_interadj_coef <- coeftest(ICD3admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD4admi_sr_interadj_coef <- coeftest(ICD4admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD5admi_sr_interadj_coef <- coeftest(ICD5admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD6admi_sr_interadj_coef <- coeftest(ICD6admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD7admi_sr_interadj_coef <- coeftest(ICD7admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD8admi_sr_interadj_coef <- coeftest(ICD8admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD9admi_sr_interadj_coef <- coeftest(ICD9admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD10admi_sr_interadj_coef <- coeftest(ICD10admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD11admi_sr_interadj_coef <- coeftest(ICD11admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD12dmi_sr_interadj_coef <- coeftest(ICD12admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD13admi_sr_interadj_coef <- coeftest(ICD13admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD14admi_sr_interadj_coef <- coeftest(ICD14admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD17admi_sr_interadj_coef <- coeftest(ICD17admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD18admi_sr_interadj_coef <- coeftest(ICD18admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD19admi_sr_interadj_coef <- coeftest(ICD19admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD21admi_sr_interadj_coef <- coeftest(ICD21admi_sr_interadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
# exponentiated coefficients:
as.matrix(round(exp(coef(ICD1admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD2admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD3admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD4admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD5admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD6admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD7admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD8admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD9admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD10admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD11admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD12admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD13admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD14admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD17admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD18admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD19admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD21admi_sr_interadj)), digits = 2), ncol = 2, rownames = TRUE)
# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((ICD1admi_sr_interadj_cis <- coefci(ICD1admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD2admi_sr_interadj_cis <- coefci(ICD2admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD3admi_sr_interadj_cis <- coefci(ICD3admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD4admi_sr_interadj_cis <- coefci(ICD4admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD5admi_sr_interadj_cis <- coefci(ICD5admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD6admi_sr_interadj_cis <- coefci(ICD6admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD7admi_sr_interadj_cis <- coefci(ICD7admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD8admi_sr_interadj_cis <- coefci(ICD8admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD9admi_sr_interadj_cis <- coefci(ICD9admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD10admi_sr_interadj_cis <- coefci(ICD10admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD11admi_sr_interadj_cis <- coefci(ICD11admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD12admi_sr_interadj_cis <- coefci(ICD12admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD13admi_sr_interadj_cis <- coefci(ICD13admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD14admi_sr_interadj_cis <- coefci(ICD14admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD17admi_sr_interadj_cis <- coefci(ICD17admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD18admi_sr_interadj_cis <- coefci(ICD18admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD19admi_sr_interadj_cis <- coefci(ICD19admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD21admi_sr_interadj_cis <- coefci(ICD21admi_sr_interadj, vcov = vcovCL, cluster = ~eid))), digits = 2)

# Maximally adjusted models: 
#  models not run for ICD 15, 16, 20 & 22 bc either no or not enough admi's --
ICD1admi_sr_maxadj <- glm.nb(ICD1_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD2admi_sr_maxadj <- glm.nb(ICD2_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD3admi_sr_maxadj <- glm.nb(ICD3_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD4admi_sr_maxadj <- glm.nb(ICD4_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD5admi_sr_maxadj <- glm.nb(ICD5_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD6admi_sr_maxadj <- glm.nb(ICD6_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD7admi_sr_maxadj <- glm.nb(ICD7_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD8admi_sr_maxadj <- glm.nb(ICD8_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD9admi_sr_maxadj <- glm.nb(ICD9_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD10admi_sr_maxadj <- glm.nb(ICD10_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD11admi_sr_maxadj <- glm.nb(ICD11_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD12admi_sr_maxadj <- glm.nb(ICD12_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD13admi_sr_maxadj <- glm.nb(ICD13_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD14admi_sr_maxadj <- glm.nb(ICD14_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD17admi_sr_maxadj <- glm.nb(ICD17_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y)
ICD18admi_sr_maxadj <- glm.nb(ICD18_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
ICD19admi_sr_maxadj <- glm.nb(ICD19_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y)
ICD21admi_sr_maxadj <- glm.nb(ICD21_admibyfuyear ~ diagnosis_severe_med + fu_year + age.centred + sex + ethnicity + townsend + centre + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre + hypertension + MI.baseline + stroke.baseline + cancer.baseline.all + sleep_apnoea.baseline + CKD + PVD + mental, data = sr_longhosp10y) 
# use coeftest to get true significance levels 
round((ICD1admi_sr_maxadj_coef <- coeftest(ICD1admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD2admi_sr_maxadj_coef <- coeftest(ICD2admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD3admi_sr_maxadj_coef <- coeftest(ICD3admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD4admi_sr_maxadj_coef <- coeftest(ICD4admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD5admi_sr_maxadj_coef <- coeftest(ICD5admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD6admi_sr_maxadj_coef <- coeftest(ICD6admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD7admi_sr_maxadj_coef <- coeftest(ICD7admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD8admi_sr_maxadj_coef <- coeftest(ICD8admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD9admi_sr_maxadj_coef <- coeftest(ICD9admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD10admi_sr_maxadj_coef <- coeftest(ICD10admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD11admi_sr_maxadj_coef <- coeftest(ICD11admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD12dmi_sr_maxadj_coef <- coeftest(ICD12admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD13admi_sr_maxadj_coef <- coeftest(ICD13admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD14admi_sr_maxadj_coef <- coeftest(ICD14admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD17admi_sr_maxadj_coef <- coeftest(ICD17admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD18admi_sr_maxadj_coef <- coeftest(ICD18admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD19admi_sr_maxadj_coef <- coeftest(ICD19admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((ICD21admi_sr_maxadj_coef <- coeftest(ICD21admi_sr_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
# exponentiated coefficients:
as.matrix(round(exp(coef(ICD1admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD2admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD3admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD4admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD5admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD6admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD7admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD8admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD9admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD10admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD11admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD12admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD13admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD14admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD17admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD18admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD19admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(ICD21admi_sr_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
# EXPONENTIATED CLUSTERED CONFIDENCE INTERVALS:
round(exp((ICD1admi_sr_maxadj_cis <- coefci(ICD1admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD2admi_sr_maxadj_cis <- coefci(ICD2admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD3admi_sr_maxadj_cis <- coefci(ICD3admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD4admi_sr_maxadj_cis <- coefci(ICD4admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD5admi_sr_maxadj_cis <- coefci(ICD5admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD6admi_sr_maxadj_cis <- coefci(ICD6admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD7admi_sr_maxadj_cis <- coefci(ICD7admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD8admi_sr_maxadj_cis <- coefci(ICD8admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD9admi_sr_maxadj_cis <- coefci(ICD9admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD10admi_sr_maxadj_cis <- coefci(ICD10admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD11admi_sr_maxadj_cis <- coefci(ICD11admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD12admi_sr_maxadj_cis <- coefci(ICD12admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD13admi_sr_maxadj_cis <- coefci(ICD13admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD14admi_sr_maxadj_cis <- coefci(ICD14admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD17admi_sr_maxadj_cis <- coefci(ICD17admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD18admi_sr_maxadj_cis <- coefci(ICD18admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD19admi_sr_maxadj_cis <- coefci(ICD19admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((ICD21admi_sr_maxadj_cis <- coefci(ICD21admi_sr_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)




### FULLY MATCHED MODELS ###########################################

####~ Mild vs no asthma ###########################################

# FACTOR VARIABLES FIRST!
mf_longhosp10y$ethnicity <- factor(mf_longhosp10y$ethnicity, ordered = FALSE,
                                   levels = c("White", "Black", "South Asian", "Others"))

mf_longhosp10y$smoke <- factor(mf_longhosp10y$smoke, ordered = FALSE,
                               levels = c("Never", "Current", "Previous"))

mf_longhosp10y$BMI_cat <- factor(mf_longhosp10y$BMI_cat, ordered = FALSE,
                                 levels = c("18.5-25", "<18.5", "25-30", 
                                            "30-35", "35-40", "40+"))

mf_longhosp10y$townsend <- factor(mf_longhosp10y$townsend, ordered = FALSE,
                                  levels = c(1, 2, 3, 4, 5))

mf_longhosp10y$centre <- factor(mf_longhosp10y$centre, ordered = FALSE,
                                levels = c(11010, 11002, 10003, 11001, 11003, 11004, 
                                           11005, 11006, 11007, 11008, 11009, 
                                           11011, 11012, 11013, 11014, 
                                           11016, 11017, 11018, 11020, 11021,
                                           11022, 11023))

# create new continuous variable for follow-up year (for models with only a single interaction term)
mf_longhosp10y <- mf_longhosp10y %>% 
  mutate(new_fu_year = case_when(mf_longhosp10y$fu_year == "year_1" ~ 1,
                                 mf_longhosp10y$fu_year == "year_2" ~ 2,
                                 mf_longhosp10y$fu_year == "year_3" ~ 3,
                                 mf_longhosp10y$fu_year == "year_4" ~ 4,
                                 mf_longhosp10y$fu_year == "year_5" ~ 5,
                                 mf_longhosp10y$fu_year == "year_6" ~ 6,
                                 mf_longhosp10y$fu_year == "year_7" ~ 7,
                                 mf_longhosp10y$fu_year == "year_8" ~ 8,
                                 mf_longhosp10y$fu_year == "year_9" ~ 9,
                                 mf_longhosp10y$fu_year == "year_10" ~ 10))

## 1 - unadjusted models:

# use glm.nb function from MASS package to estimate a negative binomial regression

library(MASS)
library(texreg)

# no. admissions:
admi1_mf <- glm.nb(admi_n ~ asthma_mild_ONLY * fu_year, data = mf_longhosp10y)

summary(admi1_mf)
screenreg(admi1_mf)

# now try zero-inflated negative binomial model
# first load necessary packages
install.packages("pscl")
install.packages("boot")
library(pscl)
library(boot)

admi2_mf <- zeroinfl(admi_n ~ asthma_mild_ONLY + fu_year, data = mf_longhosp10y,
                     dist = "negbin")
# compare both models:
screenreg(list(admi1_mf, admi2_mf))
# negative binomial is better than zero-inflated


## length of stay in hospital - negative binomial:
days1_mf <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * fu_year, data = mf_longhosp10y)
summary(days1_mf)
screenreg(days1_mf)


# costs:
# GLM with gamma distribution and log link:

costs1_mf <- glm(formula = cost2 ~ asthma_mild_ONLY + fu_year,
                 family = Gamma(link = "log"),
                 data = mf_longhosp10y)
summary(costs1_mf)
# Error in eval(family$initialize) : 
# non-positive values not allowed for the 'Gamma' family

mf_longhosp10y$cost2 <- mf_longhosp10y$totcost_byfuyear + 0.1
mf_longhosp10y$totcost_byfuyear
mf_longhosp10y$cost2
summary(mf_longhosp10y$cost2)
# repeat glm with above -> still doesn't work

# GLM - gaussian family (linear regression) and identity link:
# use glm.cluster to get clustered robust standard errors [available from miceadds package]
library(miceadds)
install.packages("sandwich")
library(sandwich)

costs2_mf <- glm.cluster(formula = totcost_byfuyear ~ asthma_mild_ONLY * fu_year,
                         cluster ="eid", 
                         family = gaussian(link = "identity"),
                         data = mf_longhosp10y)
summary(costs2_mf)

# negative binomial regression model for costs
costs3_mf <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * fu_year, data = mf_longhosp10y)

# compare linear regression and negative binomial for costs
screenreg(list(costs2_mf, costs3_mf))
# negative binomial model has a larger log likelihood, so fits better




## 2 - minimally adjusted models:

# include 4 out of 5 of the matched covariates (age, sex, ethnicity, location)
# NO DEPRIVATION!!

# negative binomial model for admissions
admi1_mf_adj <- glm.nb(admi_n ~ asthma_mild_ONLY * fu_year + age.centred
                       + sex + ethnicity + centre, 
                       data = mf_longhosp10y)

# compare non-adjusted and minimally adjusted models of admissions: [both neg. binomial]
screenreg(list(admi1_mf, admi1_mf_adj))

## negative binomial length of stay in hospital:
days1_mf_adj <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * fu_year + age.centred
                       + sex + ethnicity + centre, 
                       data = mf_longhosp10y)
summary(days1_mf_adj)

# to test zero-inflated models we need to round days to integers:
mf_longhosp10y$rounded_days <- round(mf_longhosp10y$hdays_byfuyear, digits = 0)
days3_mf_adj <- zeroinfl(rounded_days ~ asthma_mild_ONLY * fu_year 
                         + age.centred + sex + ethnicity + townsend + centre, 
                         data = mf_longhosp10y, 
                         dist = "negbin")
# compare both models:
screenreg(list(days1_mf_adj, days3_mf_adj))


# negative binomial regression model for costs
costs1_mf_adj <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * fu_year + age.centred
                        + sex + ethnicity + centre, 
                        data = mf_longhosp10y)

screenreg(costs1_mf_adj)
# compare non-adjusted and minimally adjusted models of hospital days: [both neg. binomial]
screenreg(list(costs3_mf, costs1_mf_adj))




## 3 - maximally adjusted models:

# include the matched covariates (age, sex, ethnicity, location, deprivation by townsend score)
# AND include BMI, smoking, and comorbidities


# negative binomial model for admissions
admi1_mf_maxadj <- glm.nb(admi_n ~ asthma_mild_ONLY * fu_year + age.centred
                          + sex + ethnicity + townsend + centre
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental, 
                          data = mf_longhosp10y)

# compare all 3 models of admissions: [all neg. binomial]
screenreg(list(admi1_mf, admi1_mf_adj, admi1_mf_maxadj))


## negative binomial length of stay in hospital:
days1_mf_maxadj <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * fu_year + age.centred
                          + sex + ethnicity + townsend + centre 
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental,
                          data = mf_longhosp10y)
summary(days1_mf_maxadj)
# compare all 3 models of hospital days: [all neg. binomial]
screenreg(list(days1_mf, days1_mf_adj, days1_mf_maxadj))



# negative binomial regression model for costs
costs1_mf_maxadj <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * fu_year + age.centred
                           + sex + ethnicity + townsend + centre
                           + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                           + hypertension + MI.baseline + stroke.baseline 
                           + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                           + PVD + mental,
                           data = mf_longhosp10y)

screenreg(costs1_mf_maxadj)

# compare all 3 models of hospital days: [all neg. binomial]
screenreg(list(costs3_mf, costs1_mf_adj, costs1_mf_maxadj))



install.packages("lmtest")
library(lmtest)
library(sandwich)
library(boot)

# use coeftest to get true significance levels (because clustering by eid's updates the SE's)
round((admi1_mf_coef <- coeftest(admi1_mf, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((admi1_mf_adj_coef <- coeftest(admi1_mf_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((admi1_mf_maxadj_coef <- coeftest(admi1_mf_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((days1_mf_coef <- coeftest(days1_mf, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((days1_mf_adj_coef <- coeftest(days1_mf_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((days3_mf_adj_coef <- coeftest(days3_mf_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((days1_mf_maxadj_coef <- coeftest(days1_mf_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((costs1_mf_coef <- coeftest(costs1_mf, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((costs1_mf_adj_coef <- coeftest(costs1_mf_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((costs1_mf_maxadj_coef <- coeftest(costs1_mf_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
# exponentiated coefficients:
as.matrix(round(exp(coef(admi1_mf)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(admi1_mf_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(admi1_mf_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_mf)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_mf_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days3_mf_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_mf_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_mf)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_mf_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_mf_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
# exponentiated clustered confidence intervals:
round(exp((admi1_mf_cis <- coefci(admi1_mf, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((admi1_mf_adj_cis <- coefci(admi1_mf_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((admi1_mf_maxadj_cis <- coefci(admi1_mf_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days1_mf_cis <- coefci(days1_mf, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days1_mf_adj_cis <- coefci(days1_mf_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days3_mf_adj_cis <- coefci(days3_mf_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days1_mf_maxadj_cis <- coefci(days1_mf_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((costs1_mf_cis <- coefci(costs1_mf, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((costs1_mf_adj_cis <- coefci(costs1_mf_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((costs1_mf_maxadj_cis <- coefci(costs1_mf_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)


### Likelihood ratio tests ###
## Repeat models WITHOUT interaction terms. 
# PERFORM LIKELIHOOD RATIO TESTS:

####### ~~~ without fu_year interaction #######################################
# Compare models with and without interaction for each level of adjustment:

library(lmtest)

# admissions first:

noint_admi_mf <- glm.nb(admi_n ~ asthma_mild_ONLY + fu_year, data = mf_longhosp10y)

(noint_admi_mf_coef <- coeftest(noint_admi_mf, vcov. = vcovCL, cluster = ~eid))
lrtest(m1admi_coef, m1admi2_coef)
lrtest(m1admi_mf, noint_admi_mf)
# p value is significant @ 0.001 
# thus the model WITH interaction is better

m1admi_adj_2 <- glm.nb(admi_n ~ asthma_mild_ONLY + fu_year + age.centred + 
                         sex + ethnicity + townsend + centre, 
                       data = finaldataset10y)
lrtest(m1admi_adj, m1admi_adj_2)
# p value is significant @ 0.001
# thus the model WITH interaction is better

m1admi_maxadj_2 <- glm.nb(admi_n ~ asthma_mild_ONLY + fu_year + age.centred
                          + sex + ethnicity + townsend + centre
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental, 
                          data = finaldataset10y)
lrtest(m1admi_maxadj, m1admi_maxadj_2)
# p value is NOT significant for 2nd model,
# thus the model WITHOUT interaction is better



# now length of stay:

m1days_2 <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY + fu_year, data = finaldataset10y)
lrtest(m1days, m1days_2)
# p value is NOT significant
# thus the model WITHOUT interaction is better
m1days_adj_2 <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY + fu_year + age.centred
                       + sex + ethnicity + townsend + centre, 
                       data = finaldataset10y)
lrtest(m1days_adj, m1days_adj_2)
# p value is NOT significant
# thus the model WITHOUT interaction is better
m1days_maxadj_2 <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY + fu_year + age.centred
                          + sex + ethnicity + townsend + centre 
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental,
                          data = finaldataset10y)
lrtest(m1days_maxadj, m1days_maxadj_2)
# p value is NOT significant
# thus the model WITHOUT interaction is better



# now costs:

m3costs_2 <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY + fu_year, data = finaldataset10y)
lrtest(m3costs, m3costs_2)
# p value is NOT significant
# thus the model WITHOUT interaction is better
m1costs_adj_2 <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY + fu_year + age.centred
                        + sex + ethnicity + townsend + centre, 
                        data = finaldataset10y)
lrtest(m1costs_adj, m1costs_adj_2)
# p value is NOT significant
# thus the model WITHOUT interaction is better
m1costs_maxadj_2 <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY + fu_year + age.centred
                           + sex + ethnicity + townsend + centre
                           + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                           + hypertension + MI.baseline + stroke.baseline 
                           + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                           + PVD + mental,
                           data = finaldataset10y)
lrtest(m1costs_maxadj, m1costs_maxadj_2)
# p value is NOT significant
# thus the model WITHOUT interaction is better






############ ~~~ by continuous fu_year ###########################################

## Model follow-up year as a single continuous variable for all outcomes 
## using the 'new_fu_year' variable:
# admissions:
admimf_single <- glm.nb(admi_n ~ asthma_mild_ONLY * new_fu_year, data = mf_longhosp10y)

admimf_adj_single <- glm.nb(admi_n ~ asthma_mild_ONLY * 
                              new_fu_year + age.centred + sex 
                            + ethnicity + townsend + centre, data = mf_longhosp10y)

admimf_maxadj_single <- glm.nb(admi_n ~ asthma_mild_ONLY * new_fu_year + age.centred
                               + sex + ethnicity + townsend + centre
                               + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                               + hypertension + MI.baseline + stroke.baseline 
                               + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                               + PVD + mental, 
                               data = mf_longhosp10y)

screenreg(list(admimf_single, admimf_adj_single, admimf_maxadj_single))

# clustered SE's:
(a_mr_dep_coef <- coeftest(a_mr_dep, vcov. = vcovCL, cluster = ~eid))
(a_mr_adj_dep_coef <- coeftest(a_mr_adj_dep, vcov. = vcovCL, cluster = ~eid))
(a_mr_interad_dep_coef <- coeftest(a_mr_interad_dep, vcov. = vcovCL, cluster = ~eid))
(a_mr_maxadj_dep_coef <- coeftest(a_mr_maxadj_dep, vcov. = vcovCL, cluster = ~eid))


# days:
daysmf_single <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * new_fu_year, data = mf_longhosp10y)

daysmf_adj_single <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * 
                              new_fu_year + age.centred + sex 
                            + ethnicity + townsend + centre, data = mf_longhosp10y)

daysmf_maxadj_single <- glm.nb(hdays_byfuyear ~ asthma_mild_ONLY * new_fu_year + age.centred
                               + sex + ethnicity + townsend + centre
                               + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                               + hypertension + MI.baseline + stroke.baseline 
                               + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                               + PVD + mental, 
                               data = mf_longhosp10y)

screenreg(list(daysmf_single, daysmf_adj_single, daysmf_maxadj_single))

# costs:
costsmf_single <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * new_fu_year, data = mf_longhosp10y)

costsmf_adj_single <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * 
                               new_fu_year + age.centred + sex 
                             + ethnicity + townsend + centre, data = mf_longhosp10y)

costsmf_maxadj_single <- glm.nb(totcost_byfuyear ~ asthma_mild_ONLY * new_fu_year + age.centred
                                + sex + ethnicity + townsend + centre
                                + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                                + hypertension + MI.baseline + stroke.baseline 
                                + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                                + PVD + mental, 
                                data = mf_longhosp10y)

screenreg(list(costsmf_single, costsmf_adj_single, costsmf_maxadj_single))




###~ Severe vs no asthma ###########################################


# FACTOR VARIABLES FIRST!
sf_longhosp10y$ethnicity <- factor(sf_longhosp10y$ethnicity, ordered = FALSE,
                                   levels = c("White", "Black", "South Asian", "Others"))
sf_longhosp10y$smoke <- factor(sf_longhosp10y$smoke, ordered = FALSE,
                               levels = c("Never", "Current", "Previous"))
sf_longhosp10y$BMI_cat <- factor(sf_longhosp10y$BMI_cat, ordered = FALSE,
                                 levels = c("18.5-25", "<18.5", "25-30", 
                                            "30-35", "35-40", "40+"))
sf_longhosp10y$townsend <- factor(sf_longhosp10y$townsend, ordered = FALSE,
                                  levels = c(1, 2, 3, 4, 5))
sf_longhosp10y$centre <- factor(sf_longhosp10y$centre, ordered = FALSE,
                                levels = c(11010, 11002, 10003, 11001, 11003, 11004, 
                                           11005, 11006, 11007, 11008, 11009, 
                                           11011, 11012, 11013, 11014, 
                                           11016, 11017, 11018, 11020, 11021,
                                           11022, 11023))

## 1 - unadjusted models:

library(MASS)
library(texreg)

# admissions
admi1_sf <- glm.nb(admi_n ~ diagnosis_severe_med * fu_year, data = sf_longhosp10y)
screenreg(admi1_sf)
# days
days1_sf <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * fu_year, data = sf_longhosp10y)
screenreg(days1_sf)
# costs
costs1_sf <- glm.cluster(formula = totcost_byfuyear ~ diagnosis_severe_med * fu_year,
                         cluster ="eid", 
                         family = gaussian(link = "identity"),
                         data = sf_longhosp10y)

costs2_sf <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * fu_year, data = sf_longhosp10y)
screenreg(list(costs1_sf, costs2_sf))
# again, negative binomial has much larger log likelihood


## 2 - minimally adjusted models:

# include 4 out of 5 of the matched covariates (age, sex, ethnicity, location)
# NO DEPRIVATION!!

# negative binomial models for admissions, days, and costs:

admi1_sf_adj <- glm.nb(admi_n ~ diagnosis_severe_med * fu_year + age.centred
                       + sex + ethnicity + centre, 
                       data = sf_longhosp10y)
screenreg(admi1_sf_adj)


days1_sf_adj <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * fu_year + age.centred
                       + sex + ethnicity + centre, 
                       data = sf_longhosp10y)
screenreg(days1_sf_adj)


costs1_sf_adj <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * fu_year + age.centred
                        + sex + ethnicity + centre, 
                        data = sf_longhosp10y)
screenreg(costs1_sf_adj)





## 3 - maximally adjusted models:

# include the matched covariates (age, sex, ethnicity, location, deprivation by townsend score)
# AND include BMI, smoking, and comorbidities

# negative binomial model for admissions
admi1_sf_maxadj <- glm.nb(admi_n ~ diagnosis_severe_med * fu_year + age.centred
                          + sex + ethnicity + townsend + centre
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental, 
                          data = sf_longhosp10y)

# compare all 3 models of admissions: [all neg. binomial]
screenreg(list(admi1_sf, admi1_sf_adj, admi1_sf_maxadj))


## negative binomial length of stay in hospital:
days1_sf_maxadj <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * fu_year + age.centred
                          + sex + ethnicity + townsend + centre 
                          + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                          + hypertension + MI.baseline + stroke.baseline 
                          + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                          + PVD + mental,
                          data = sf_longhosp10y)
summary(days1_sf_maxadj)
# compare all 3 models of hospital days: [all neg. binomial]
screenreg(list(days1_sf, days1_sf_adj, days1_sf_maxadj))



# negative binomial regression model for costs
costs1_sf_maxadj <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * fu_year + age.centred
                           + sex + ethnicity + townsend + centre
                           + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                           + hypertension + MI.baseline + stroke.baseline 
                           + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                           + PVD + mental,
                           data = sf_longhosp10y)

screenreg(costs1_sf_maxadj)

# compare all 3 models of hospital days: [all neg. binomial]
screenreg(list(costs3_sf, costs1_sf_adj, costs1_sf_maxadj))


# use coeftest to get true significance levels (because clustering by eid's updates the SE's)
round((admi1_sf_coef <- coeftest(admi1_sf, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((admi1_sf_adj_coef <- coeftest(admi1_sf_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((admi1_sf_maxadj_coef <- coeftest(admi1_sf_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((days1_sf_coef <- coeftest(days1_sf, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((days1_sf_adj_coef <- coeftest(days1_sf_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((days3_sf_adj_coef <- coeftest(days3_sf_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((days1_sf_maxadj_coef <- coeftest(days1_sf_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((costs1_sf_coef <- coeftest(costs1_sf, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((costs1_sf_adj_coef <- coeftest(costs1_sf_adj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
round((costs1_sf_maxadj_coef <- coeftest(costs1_sf_maxadj, vcov. = vcovCL, cluster = ~eid)), digits = 2)
# exponentiated coefficients:
as.matrix(round(exp(coef(admi1_sf)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(admi1_sf_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(admi1_sf_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_sf)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_sf_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days3_sf_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(days1_sf_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_sf)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_sf_adj)), digits = 2), ncol = 2, rownames = TRUE)
as.matrix(round(exp(coef(costs1_sf_maxadj)), digits = 2), ncol = 2, rownames = TRUE)
# exponentiated clustered confidence intervals:
round(exp((admi1_sf_cis <- coefci(admi1_sf, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((admi1_sf_adj_cis <- coefci(admi1_sf_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((admi1_sf_maxadj_cis <- coefci(admi1_sf_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days1_sf_cis <- coefci(days1_sf, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days1_sf_adj_cis <- coefci(days1_sf_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days3_sf_adj_cis <- coefci(days3_sf_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((days1_sf_maxadj_cis <- coefci(days1_sf_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((costs1_sf_cis <- coefci(costs1_sf, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((costs1_sf_adj_cis <- coefci(costs1_sf_adj, vcov = vcovCL, cluster = ~eid))), digits = 2)
round(exp((costs1_sf_maxadj_cis <- coefci(costs1_sf_maxadj, vcov = vcovCL, cluster = ~eid))), digits = 2)




####### ~~~ without fu_year interaction #######################################
# Compare models with and without interaction for each level of adjustment,
# then compare with LIKELIHOOD RATIO TESTS

library(lmtest)

# admissions

admi2_severity <- glm.nb(admi_n ~ diagnosis_severe_med + fu_year, data = sf_longhosp10y)
(admi2_severity_coef <- coeftest(admi2_severity, vcov. = vcovCL, cluster = ~eid))
lrtest(admi1_severity_coef, admi2_severity_coef)
lrtest(admi1_severity, admi2_severity)
# p value is significant @ 0.001
# thus the model WITH interaction is better

admi2_severity_adj <- glm.nb(admi_n ~ diagnosis_severe_med + fu_year + age.centred
                             + sex + ethnicity + townsend + centre, 
                             data = sf_longhosp10y)
lrtest(admi1_severity_adj, admi2_severity_adj)
# p value is significant @ 0.01
# thus the model WITH interaction is better

admi2_severity_maxadj <- glm.nb(admi_n ~ diagnosis_severe_med + fu_year + age.centred
                                + sex + ethnicity + townsend + centre
                                + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                                + hypertension + MI.baseline + stroke.baseline 
                                + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                                + PVD + mental, 
                                data = sf_longhosp10y)
lrtest(admi1_severity_maxadj, admi2_severity_maxadj)
# p value is not significant
# thus the model WITHOUT interaction is better


## length of stay

days2_severity <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med + fu_year, data = sf_longhosp10y)
lrtest(days1_severity, days2_severity)
# p value is significant @ 0.001
# thus the model WITH interaction is better

days2_severity_adj <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med + fu_year + age.centred
                             + sex + ethnicity + townsend + centre, 
                             data = sf_longhosp10y)
lrtest(days1_severity_adj, days2_severity_adj)
# p value is significant @ 0.001
# thus model WITH interaction is better

days2_severity_maxadj <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med + fu_year + age.centred
                                + sex + ethnicity + townsend + centre
                                + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                                + hypertension + MI.baseline + stroke.baseline 
                                + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                                + PVD + mental, 
                                data = sf_longhosp10y)
lrtest(days1_severity_maxadj, days2_severity_maxadj)
# # p value is significant @ 0.01
# thus the model WITH interaction is better


## costs
costs3_severity <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med + fu_year, data = sf_longhosp10y)
lrtest(costs2_severity, costs3_severity)
# p value = non-significant
# thus model WITHOUT interaction is better

costs2_severity_adj <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med + fu_year + age.centred
                              + sex + ethnicity + townsend + centre, 
                              data = sf_longhosp10y)
lrtest(costs1_severity_adj, costs2_severity_adj)
# p value = non-significant
# thus model WITHOUT interaction is better

costs2_severity_maxadj <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med + fu_year + age.centred
                                 + sex + ethnicity + townsend + centre
                                 + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                                 + hypertension + MI.baseline + stroke.baseline 
                                 + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                                 + PVD + mental, 
                                 data = sf_longhosp10y)
lrtest(costs1_severity_maxadj, costs2_severity_maxadj)
# # p value = non-significant. 
# thus the model WITHOUT interaction is better



############ ~~~ by continuous fu_year ###########################################

## Model follow-up year as a single continuous variable for all outcomes 
## using the 'new_fu_year' variable:
# admissions:
admi_severity_single <- glm.nb(admi_n ~ diagnosis_severe_med * new_fu_year, data = sf_longhosp10y)

admi_severity_adj_single <- glm.nb(admi_n ~ diagnosis_severe_med * new_fu_year +
                                     age.centred + sex 
                                   + ethnicity + townsend + centre, data = sf_longhosp10y)

admi_severity_maxadj_single <- glm.nb(admi_n ~ diagnosis_severe_med * new_fu_year + age.centred
                                      + sex + ethnicity + townsend + centre
                                      + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                                      + hypertension + MI.baseline + stroke.baseline 
                                      + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                                      + PVD + mental, 
                                      data = sf_longhosp10y)

screenreg(list(admi_severity_single, admi_severity_adj_single, admi_severity_maxadj_single))

# days:
days_severity_single <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * new_fu_year, data = sf_longhosp10y)

days_severity_adj_single <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * 
                                     new_fu_year + age.centred + sex 
                                   + ethnicity + townsend + centre, data = sf_longhosp10y)

days_severity_maxadj_single <- glm.nb(hdays_byfuyear ~ diagnosis_severe_med * new_fu_year + age.centred
                                      + sex + ethnicity + townsend + centre
                                      + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                                      + hypertension + MI.baseline + stroke.baseline 
                                      + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                                      + PVD + mental, 
                                      data = sf_longhosp10y)

screenreg(list(days_severity_single, days_severity_adj_single, days_severity_maxadj_single))

# costs:
costs_severity_single <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * new_fu_year, data = sf_longhosp10y)

costs_severity_adj_single <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * 
                                      new_fu_year + age.centred + sex 
                                    + ethnicity + townsend + centre, data = sf_longhosp10y)

costs_severity_maxadj_single <- glm.nb(totcost_byfuyear ~ diagnosis_severe_med * new_fu_year + age.centred
                                       + sex + ethnicity + townsend + centre
                                       + smoke + BMI_cat + COPD.baseline + diabetes.fo.pre
                                       + hypertension + MI.baseline + stroke.baseline 
                                       + cancer.baseline.all + sleep_apnoea.baseline + CKD 
                                       + PVD + mental, 
                                       data = sf_longhosp10y)

screenreg(list(costs_severity_single, costs_severity_adj_single, costs_severity_maxadj_single))
