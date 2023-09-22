# UKBiobank
Hospital analyses of asthma patients in the UK Biobank cohort

Please progress through the code in the following order of scripts:
1_define_variables.R 
2_matching.R
3_hospital_variables.R
4_descriptive_analyses.R
5_modelling.R

• 1_define_variables.R: this script loads and cleans data and defines key variables for the population (age, sex, MI, stroke, 
  cancer, diabetes, hypertension, PVD, smoking status, ethnicity, BMI, CKD, COPD, sleep apnoea, deprivation (Townsend), mental health illness, location (UKB assessment centre)),
  censoring date and follow_up completion. The asthma population is then defined.

• 2_matching.R: this script documents the matching algorithm for matching mild and moderate-severe asthma patients separately to controls using R’s “MatchIt” package, giving 
  the final study population. Baseline characteristics are then defined for the asthma cohorts and control cohorts separately, using the "tableone" package.

• 3_hospital_variables.R: this script defines the following hospital variables: hospital admissions by ICD-10 chapter (including respiratory), annual no. of admissions, annual 
  no. of days spent in hospital and annual hospital costs

• 4_descriptive_analyses.R: matched data are merged with hospital data and duration of follow_up is limited to 10-years before performing descriptive analyses of the following:   
  i) mean hospital costs, no. admissions and no. hospital days by follow_up year and ii) mean hospital costs, admissions and days by ICD-10 chapter 

• 5_modelling.R: Negative binomial regression models of hospital admissions, days spent in hospital and hospital costs for people with mild asthma compared to matched controls, 
  and separately for people with moderate-severe asthma compared to matched controls. Three levels of covariate adjustments were pre-specified: (1) minimally adjusted models with 
  adjustments for the matching factors (age, sex, ethnicity, and location) and calendar year; (2) intermediately-adjusted models with adjustments for quintile of socioeconomic 
  deprivation (as measured by the Townsend deprivation index) in addition to the matching factors and calendar year; and (3) fully adjusted models with adjustments for smoking 
  status, BMI category and comorbidities in addition to socioeconomic deprivation, matching factors and calendar year.  
  Additional analyses assess the impact of asthma on different types (ICD-10) of hospital admission and whether associations between asthma severity and hospital outcomes 
  differed with duration of follow-up or across quintiles of socioeconomic deprivation. 


