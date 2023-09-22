library(tidyverse)

# Load hospital data ################################

# open files:
ukb_annual_cost_hosp <- readRDS("UKB costs processed 2022-10/UKB hosp admission costs/data/ukb_annual_cost_hosp.rds")
summary(ukb_annual_cost_hosp) # only 5 columns. Organised by follow-up YEAR per eid
head(ukb_annual_cost_hosp, 20)
# This dataset contains the annual hospital cost in 2018/2019 value for 
# all the calendar years since participants' entry until 2016/03/13, lost 
# of follow-up or death, whichever earlier. 
# -	year = calendar year since participant's entry 
# -	yearnotobs = proportion of year not being observed. It is 0 before the last year of follow-up, and usually higher than 0 at the last year of follow-up due to censoring (except the year of death, in which yearnotobs is assumed to be 0, implying full year of observation at the year of death). proportion of year not being observed. It is 0 before the last year of follow-up, and usually higher than 0 at the last year of follow-up due to censoring (except the year of death, in which yearnotobs is assumed to be 0, implying full year of observation at the year of death).
# -	cost = annual hospital cost

hesin2_cost_hrg <- readRDS("UKB costs processed 2022-10/UKB hosp admission costs/data/hesin2_cost_hrg.rds")
summary(hesin2_cost_hrg) # 14 columns - contains "ins_index", epiorder, hrg, tot, core, c.type, dur etc. 
# organised by ins_index per eid 

hesin2_prehrg <- readRDS("UKB costs processed 2022-10/UKB hosp admission costs/data/hesin2_prehrg.rds")
colnames(hesin2_prehrg) # lotsa columns - contains DIAGS, OPERS, ins_index, admidate, disdate, epistart, epiend
# also organised by ins_index per eid 
# there are 14 DIAG_nn columns which can be populated per episode
# the field DIAG_01 contains the PRIMARY diagnosis for each admission
# (subsequent DIAG_nn's represent secondary/subsidiary diagnoses)

## merge datasets ######

# ^ combine datasets to get necessary variables
# ensure merging to keep a structure where data is organised by eid and ins_index  
ukbhospcost2022 <- merge(hesin2_prehrg[,c("eid", "ins_index", "admidate", "disdate", 
                                          "epistart", "epiend", "DIAG01")], hesin2_cost_hrg[,c(
                                            "eid", "ins_index", "tot")], by=c("eid", "ins_index"), all.x = T)
head(ukbhospcost2022)
# this is essentially a mini dataset with only the variables of interest
# we can merge it back later by id or id+instance

# use dsource variable from hesin to get data sources 
hesin <- read.delim(file.path(raw_data, "hesin.txt"))

# merge the 'dsource' variable into UKBhospcost2022, by id and instance
# maintain the structure of UKBhospcost2022
ukbhospcost2022 <- merge(ukbhospcost2022, hesin[,c("eid", "ins_index", "dsource")], by=c("eid", "ins_index"), all.x = T)
head(ukbhospcost2022)

## merge recruit date & other columns from main dataset into ukbhospcost2022 dataset
ukbhospcost2022 <- merge(ukbhospcost2022, UKB_a_2ndversion[, c("eid", "recruit.date",
                                                               "lossfudate", "death.date", 
                                                               "censordate", "fu_duration", 
                                                               "fu_year")], by = "eid", all.x = T)
# summarise admission and recruit dates and difference
summary(as.numeric(ukbhospcost2022$admidate-ukbhospcost2022$recruit.date)/365.25)
 

### admissions by ICD-10 chapter #########

# Categorise by ICD-10 chapter of primary diagnosis in first episode of admission. 
# Chapter X = diseases of the respiratory system, codes J00 - J99
# Asthma = codes J45 & J46
# make new column of binary (1 or 0) for "respiratory" yes or no as primary diagnosis of admission
ukbhospcost2022$resp_admi <- 0
ukbhospcost2022$resp_admi[grepl("^J", ukbhospcost2022$DIAG01)] <- 1
table(ukbhospcost2022$resp_admi)

# make binary column yes/ no for asthma as primary diagnosis
ukbhospcost2022$asthma_admi <- 0
ukbhospcost2022$asthma_admi[grepl("^J45|^J46", ukbhospcost2022$DIAG01)] <- 1
table(ukbhospcost2022$asthma_admi)


# Code categorical column containing ICD10 chapters 1-22

# 1 = infectious/ parasitic diseases
# 2 = neoplasms
# 3 = blood diseases
# 4 = endocrine/ nutritional/ metabolic diseases
# 5 = mental/ behavioural disorders
# 6 = nervous system diseases
# 7 = eye diseases
# 8 = ear diseases
# 9 = circulatory diseases
# 10 = respiratory diseases
# 11 = digestive diseases
# 12 = skin diseases
# 13 = musculoskeletal/ connective diseases
# 14 = genitourinary diseases
# 15 = pregnancy/ childbirth/ the puerperium
# 16 = perinatal period conditions 
# 17 = congenital malformations/ deformations/ chromosomal abnormalities
# 18 = symptoms/ signs/ abnormal clinical & lab findings, unclassified
# 19 = injury/ poisoning/ certain other consequences of external causes
# 20 = external causes of morbidity & mortality
# 21 = factors influencing health status and contact with health service
# 22 = codes for special purposes 

# categorical column of overall ICD chapter 
ukbhospcost2022$ICD_diag01 <- 0
ukbhospcost2022$ICD_diag01[grepl("^A|^B", ukbhospcost2022$DIAG01)] <- 1
ukbhospcost2022$ICD_diag01[grepl("^C|^D0|^D1|^D2|^D3|^D4", ukbhospcost2022$DIAG01)] <- 2
ukbhospcost2022$ICD_diag01[grepl("^D5|^D6|^D7|^D8", ukbhospcost2022$DIAG01)] <- 3
ukbhospcost2022$ICD_diag01[grepl("^E", ukbhospcost2022$DIAG01)] <- 4
ukbhospcost2022$ICD_diag01[grepl("^F", ukbhospcost2022$DIAG01)] <- 5
ukbhospcost2022$ICD_diag01[grepl("^G", ukbhospcost2022$DIAG01)] <- 6
ukbhospcost2022$ICD_diag01[grepl("^H0|^H1|^H2|^H3|^H4|^H5", ukbhospcost2022$DIAG01)] <- 7
ukbhospcost2022$ICD_diag01[grepl("^H6|^H7|^H8|^H9", ukbhospcost2022$DIAG01)] <- 8
ukbhospcost2022$ICD_diag01[grepl("^I", ukbhospcost2022$DIAG01)] <- 9
ukbhospcost2022$ICD_diag01[grepl("^J", ukbhospcost2022$DIAG01)] <- 10
ukbhospcost2022$ICD_diag01[grepl("^K", ukbhospcost2022$DIAG01)] <- 11
ukbhospcost2022$ICD_diag01[grepl("^L", ukbhospcost2022$DIAG01)] <- 12
ukbhospcost2022$ICD_diag01[grepl("^M", ukbhospcost2022$DIAG01)] <- 13
ukbhospcost2022$ICD_diag01[grepl("^N", ukbhospcost2022$DIAG01)] <- 14
ukbhospcost2022$ICD_diag01[grepl("^O", ukbhospcost2022$DIAG01)] <- 15
ukbhospcost2022$ICD_diag01[grepl("^P", ukbhospcost2022$DIAG01)] <- 16
ukbhospcost2022$ICD_diag01[grepl("^Q", ukbhospcost2022$DIAG01)] <- 17
ukbhospcost2022$ICD_diag01[grepl("^R", ukbhospcost2022$DIAG01)] <- 18
ukbhospcost2022$ICD_diag01[grepl("^S|^T", ukbhospcost2022$DIAG01)] <- 19
ukbhospcost2022$ICD_diag01[grepl("^V|^Y", ukbhospcost2022$DIAG01)] <- 20
ukbhospcost2022$ICD_diag01[grepl("^Z", ukbhospcost2022$DIAG01)] <- 21
ukbhospcost2022$ICD_diag01[grepl("^U", ukbhospcost2022$DIAG01)] <- 22

table(ukbhospcost2022$ICD_diag01) # 3,636 are empty (0) i.e., they do not have a reason/ primary diagnosis assigned
# 3636 / 2314850 = 0.2% of admissions are undiagnosed. 

# make binary column yes/ no for blood diseases (ICD chapter 3) as primary diagnosis
ukbhospcost2022$ICD3_admi <- 0
ukbhospcost2022$ICD3_admi[grepl("^D5|^D6|^D7|^D8", ukbhospcost2022$DIAG01)] <- 3
table(ukbhospcost2022$ICD3_admi)



#### impute admissions #######

# Now count admission in particular years
# the same method as previously

# (if we want to transform to wide format:)
ukbhospcost2022_wide <- ukbhospcost2022 %>% group_by(eid) %>% distinct(eid, .keep_all = T) 

summary(ukbhospcost2022$admidate); summary(ukbhospcost2022$epistart)
# admidate has many missing values, while epistart has no missing

# implement Boby's algorithm for missing admidate

# copy a new dataset for imputation and order by eid, epistart, epiend
# generate a count of admission for each id (generate an "admcount" column)
ukbhospcost2022_imp <- ukbhospcost2022 %>% 
  arrange(eid, epistart, epiend) %>% 
  select(eid, ins_index, epistart, epiend) %>% 
  group_by(eid) %>% 
  mutate(admcount=n()) 

# ^ here mutate() creates new column admcount filled by function =n(),
# which counts the no. of rows per group, in this case "eid" (we see it's 
# preceded by group_by(eid %>%))


# create an order index for episode with individual (generate an epi_index column)
# so epi_index shows numbers each episode for each individual
ukbhospcost2022_imp <- ukbhospcost2022_imp %>% group_by(eid) %>% mutate(epi_index=1:n())

# create our own admidate2 and disdate2
ukbhospcost2022_imp$admidate2 <- as.Date(NA)

ukbhospcost2022_imp$disdate2 <- as.Date(NA)

# for epi_index==1, impute admidate2 & disdate2 with epistart and epiend (discharge date)
ukbhospcost2022_imp$admidate2[ukbhospcost2022_imp$epi_index==1] <- 
  as.Date(ukbhospcost2022_imp$epistart[ukbhospcost2022_imp$epi_index==1])

ukbhospcost2022_imp$disdate2[ukbhospcost2022_imp$epi_index==1] <- 
  as.Date(ukbhospcost2022_imp$epiend[ukbhospcost2022_imp$epi_index==1])

# for id's with only 1 admission, the work has been complete


# now focus on id's with more than one admission
# split the dataset into two
# imp_1 keeps those with only 1 admission. Their admi + disdates don't need to be imputed. 
# only those with multiple admisions (imp_m) need their admissions (beyond the 1st admission)
# to be imputed, to ensure no overlap of epistart/ end dates etc. 
# work on ukbhospcost2022_imp_m
ukbhospcost2022_imp_1 <- ukbhospcost2022_imp %>% filter(admcount==1)

ukbhospcost2022_imp_m <- ukbhospcost2022_imp %>% filter(admcount>1)

###' HERE, WE SKIP THIS STEP AND USE RUNGUO'S CODE TO IMPUTE THOSE WITH 
###' ADMISSIONS > 1 (in "output" file. This code is in admi_count script)
output <- readRDS("hcost_imp_m.rds")

# combine the imputed dataset and ukbhospcost2022_imp_1 to replace ukbhospcost2022_imp
ukbhospcost2022_imp <- rbind(ukbhospcost2022_imp_1, output[,1:8])

#saveRDS(ukbhospcost2022_imp, file = file.path(work_data, "ukbhospcost2022_imp.rds"), compress = F)
ukbhospcost2022_imp <- readRDS("data_work/ukbhospcost2022_imp.rds")


# merge admidate2 and disdate2 from ukbhospcost2022_imp into ukbhospcost2022
# by eid, ins_index

ukbhospcost2022 <- subset(ukbhospcost2022, select = -c(admidate2, disdate2)) # remove old draft admidate 2 / disdate2

ukbhospcost2022 <- merge(ukbhospcost2022, ukbhospcost2022_imp[, c("eid", "ins_index", "admidate2", "disdate2")], by=c("eid", "ins_index"), all.x = T)
head(ukbhospcost2022)

# and  merge with long format:
ukbhospcost2022_uni_long_imp <- merge(ukbhospcost2022_uni_long, ukbhospcost2022_imp[, c("eid", "ins_index", "admidate2", "disdate2")], by=c("eid", "ins_index"), all.x = T)

# from then on use the imputed admidate2 and disdate2 
# leave the original admidate and disdate as references



##### Year of follow-up ###########
# (with imputed admi dates)

# first generate "fu_days"
ukbhospcost2_uni_long_imp$fu_days = as.numeric(ukbhospcost2_uni_long_imp$admidate2-ukbhospcost2_uni_long_imp$recruit.date)

ukbhospcost2_uni_long_imp <- ukbhospcost2_uni_long_imp %>% 
  mutate(fu_year = 
           case_when(fu_days %/% (365.25) == -3 ~ "year_3_ago",
                     fu_days %/% (365.25) == -2 ~ "year_2_ago",
                     fu_days %/% (365.25) == -1 ~ "year_1_ago",
                     fu_days %/% (365.25) == 0 ~ "year_1",
                     fu_days %/% (365.25) == 1 ~ "year_2",
                     fu_days %/% (365.25) == 2 ~ "year_3",
                     fu_days %/% (365.25) == 3 ~ "year_4",
                     fu_days %/% (365.25) == 4 ~ "year_5",
                     fu_days %/% (365.25) == 5 ~ "year_6",
                     fu_days %/% (365.25) == 6 ~ "year_7",
                     fu_days %/% (365.25) == 7 ~ "year_8",
                     fu_days %/% (365.25) == 8 ~ "year_9",
                     fu_days %/% (365.25) == 9 ~ "year_10",
                     fu_days %/% (365.25) == 10 ~ "year_11",
                     fu_days %/% (365.25) == 11 ~ "year_12",
                     fu_days %/% (365.25) == 12 ~ "year_13",
                     fu_days %/% (365.25) == 13 ~ "year_14",
                     fu_days %/% (365.25) == 14 ~ "year_15"
           ))

# check the counts:
table(ukbhospcost2_mini_long_imp$fu_year, useNA = "ifany")
# e.g. 411,120 episodes at 7 yrs hospital follow-up
# 1,520 episodes at 14 yrs follow-up

#saveRDS(ukbhospcost2_uni_long_imp, file = file.path(work_data, "ukbhospcost2_uni_long_imp.rds"), compress = F)


# repeat for ukbhospcost2022:

ukbhospcost2022$fu_days = as.numeric(ukbhospcost2022$admidate2-ukbhospcost2022$recruit.date)

ukbhospcost2022 <- ukbhospcost2022 %>%  
  mutate(fu_year = 
           case_when(fu_days %/% (365.25) == -3 ~ "year_3_ago",
                     fu_days %/% (365.25) == -2 ~ "year_2_ago",
                     fu_days %/% (365.25) == -1 ~ "year_1_ago",
                     fu_days %/% (365.25) == 0 ~ "year_1",
                     fu_days %/% (365.25) == 1 ~ "year_2",
                     fu_days %/% (365.25) == 2 ~ "year_3",
                     fu_days %/% (365.25) == 3 ~ "year_4",
                     fu_days %/% (365.25) == 4 ~ "year_5",
                     fu_days %/% (365.25) == 5 ~ "year_6",
                     fu_days %/% (365.25) == 6 ~ "year_7",
                     fu_days %/% (365.25) == 7 ~ "year_8",
                     fu_days %/% (365.25) == 8 ~ "year_9",
                     fu_days %/% (365.25) == 9 ~ "year_10",
                     fu_days %/% (365.25) == 10 ~ "year_11",
                     fu_days %/% (365.25) == 11 ~ "year_12",
                     fu_days %/% (365.25) == 12 ~ "year_13",
                     fu_days %/% (365.25) == 13 ~ "year_14",
                     fu_days %/% (365.25) == 14 ~ "year_15"
           ))



#saveRDS(ukbhospcost2022, file = file.path(work_data, "ukbhospcost2022.rds"), compress = F)
#saveRDS(ukbhospcost2_uni_long_imp, file = file.path(work_data, "ukbhospcost2_uni_long_imp.rds"), compress = F)


!!# Load this before each session !!
ukbhospcost2022 <- readRDS("data_work/ukbhospcost2022.rds")
ukbhospcost2_uni_long_imp <- readRDS("data_work/ukbhospcost2_uni_long_imp.rds") 



###### long-form dataset + fu_completion #############

# start to change the data structure
# create a large table where each id has the same no. of rows indicating follow-up years such that we can fill information into the table

# the current data sets have the structure of single eid that could correspond to multiple ins_index 
# however, the lossfudate, death.date, censordate and fu_duration should be unique to each id
# generate a new dataset with unique ids
# in case censordate differ by ins_index within one individual, so duration differ
# keep the longest duration only
## i.e., convert to wide-form dataset (1 row per eid): 
ukbhospcost2022_uni <- ukbhospcost2022 %>% group_by(eid) %>% filter(fu_duration==max(fu_duration)) %>% distinct(eid, .keep_all = T)
summary(ukbhospcost2022_uni)
head(ukbhospcost2022_uni,20) # a wide-form dataset
## now we have a dataset with unique id and follow-up duration
# (by doing this we lose the ins_index structure/ data)


table(ukbhospcost2022_uni$fu_year)

# since the longest fu year is 15.05, generate 16 year columns

for (i in 1:16) {
  
  # first create one column with all 0
  ukbhospcost2022_uni[[paste0("year_", i)]] <- 0
  
  # change the year with full follow up to be 1
  ukbhospcost2022_uni[[paste0("year_", i)]][ukbhospcost2022_uni$fu_duration / 365.25 >= i] <- 1
  
  # change the year with partial follow up to a proportion of the year
  ukbhospcost2022_uni[[paste0("year_", i)]][ukbhospcost2022_uni$fu_duration / 365.25 < i & ukbhospcost2022_uni$fu_duration / 365.25 > i-1] <- 
    (ukbhospcost2022_uni$fu_duration[ukbhospcost2022_uni$fu_duration / 365.25 < i & ukbhospcost2022_uni$fu_duration / 365.25 > i-1] %% 365.25)/365.25 # %%: modulus ex. 366%%365.25=0.75
}

head(ukbhospcost2022_uni, 20) # a wide-form dataset
colnames(ukbhospcost2022_uni) # we can see all the years separately 


### generate "fu_completion" 
# transform wide to long dataset
ukbhospcost2022_uni_long <- ukbhospcost2022_uni %>% gather(paste0("year_", 1:16), key=fu_year, 
                                                           value = fu_completion) %>% arrange(eid) %>% mutate(fu_duration=NULL)
head(ukbhospcost2022_uni_long, 32)
colnames(ukbhospcost2022_uni_long)
ukbhospcost2022_uni_long[10:20,c("eid","fu_year","fu_completion")]

## ^ now we have a long-form dataset
# each individual has 16 rows from year 1 to 16
# fu_completion = 1 means the year is followed; 0 nonfollowed; decimal = partially followed
# current dataset is in the form of a large table where each id has the same no. of rows indicating follow-up years


################## 1) annual costs ############

# sum admission costs - i.e. annual costs by FU year
annualhospcost <- aggregate(x = ukbhospcost2022$tot,
                            by = list(ukbhospcost2022$eid, ukbhospcost2022$fu_year),
                            FUN = sum)

# rename variables
annualhospcost <- annualhospcost %>% rename(eid=Group.1, fu_year=Group.2, totcost_byfuyear=x) %>% arrange(eid, fu_year)
head(annualhospcost)

# make eid and fu_year distinct:
annualhospcost <- annualhospcost %>% distinct(eid, fu_year, .keep_all = T)

# annualhospcost is FINAL COST DATASET
#saveRDS(annualhospcost, file = file.path(work_data, "annualhospcost.rds"), compress = F)

# calculate mean hospital costs by FU year:
aggregate(annualhospcost$totcost_byfuyear,list(annualhospcost$fu_year),mean)

n_distinct(annualhospcost$eid) # there are hospital records for 373,269 people in biobank


# final annual costs:

finalhospdata <- merge(ukbhospcost2022, annualhospcost, by=c("eid", "fu_year"), all.x = T)



#################### 2) no. overall admissions + days in hospital ############

# before keep distinct admidate2, we need to sort 
# ensure we keep the last record (the longest of duration) of records with the same admidate2 for each individual
# sort or not will not influence the count of admission

# first sort by eid and admidate2 in ascending order, but disdate2 in descending order
ukbhospcost2022_sort <- ukbhospcost2022 %>% arrange(eid, admidate2, desc(disdate2))
head(ukbhospcost2022_sort)


# generate a new dataset but only keep distinct admidate 
# i.e. remove overlapping episodes with the same admidate
# as a result of sorting above, actually keeps admidate2 with the latest disdate2
# create new column admi_duration representing the no. of days in hosp
# also use this opportunity to maintain resp_admi & asthma_admi columns such that they can be merged later
admi_count2 <- ukbhospcost2022_sort %>% 
  group_by(eid) %>% 
  distinct(admidate2, .keep_all = T) %>% 
  mutate(admi_duration = as.numeric(disdate2 - admidate2)) %>% 
  select(eid, fu_year, admi_duration, admidate2, disdate2, resp_admi, asthma_admi, ICD_diag01)

summary(admi_count2$admi_duration)
colnames(admi_count2)

# replace 0's (i.e. where someone has been discharged on the same date as admitted) with 0.5

admi_count2$admi_duration[admi_count2$admi_duration==0] <- 0.5
head(admi_count2, 20)


# generate a variable indicating count the number of admission per fu year (admi_n)
admi_count2 <- admi_count2 %>% 
  group_by(eid, fu_year) %>% 
  mutate(admi_n = n()) %>% 
  select(eid, fu_year, admi_n, admi_duration, resp_admi, asthma_admi, ICD_diag01)


# generate a variable indicating count the number of admission per fu year (admi_n)

head(admi_count2$admi_n)


# sum admission duration - i.e. no. days in hospital by FU year
hospdaysby_fuyear2 <- aggregate(x = admi_count2$admi_duration,
                                by = list(admi_count2$eid, admi_count2$fu_year),
                                FUN = sum)

# rename variables
hospdaysby_fuyear2 <- hospdaysby_fuyear2 %>% rename(eid=Group.1, fu_year=Group.2, hdays_byfuyear=x)

# merge the 3 variables together in one dataset:
admi_count2 <- merge(admi_count2, hospdaysby_fuyear2, by=c("eid", "fu_year"), all.x = T)
head(admi_count2,10)


# get mean per FU year
aggregate(hospdaysby_fuyear2$hdays_byfuyear,list(hospdaysby_fuyear2$fu_year),mean)


# code the NA admin_n as 0
admi_count2$admi_n[is.na(admi_count2$admi_n)] <- 0

summary(admi_count2$admi_n)
# ^ we see that the max no. admissions for a person is 210


# transform from long to wide 
admi_count_wide <- admi_count2 %>% 
  group_by(eid, fu_year) %>% 
  mutate(admi_index = 1:n()) %>% 
  spread(admi_index, admi_duration)

head(admi_count_wide,10)
summary(admi_count_wide$hdays_byfuyear)

# rename the columns
names(admi_count_wide)[4:213] <- paste0("admi_", colnames(admi_count_wide[, 4:213]))

head(admi_count_wide,10)
colnames(admi_count_wide)
# now we have had a dataset with id, follow up year, number of admission and the duration of each admission


admi_count_wide2 <- admi_count_wide %>% select("eid", "fu_year", "admi_n", "admi_hdays_byfuyear", resp_admi, asthma_admi)
head(admi_count_wide2)



# Final hospital data only:  (showing only eid's that have experienced hosp admission(s))
finalhospdata <- merge(admi_count2, annualhospcost, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD <- merge(admi_count2, annualhospcost, by=c("eid", "fu_year"), all.x = T)
head(finalhospdata_byICD,20)


#################### 3) no. of respiratory admissions/ year ####################
respadmi_n <- aggregate(x = finalhospdata_byICD$resp_admi, 
                        by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                        FUN = sum)

respadmi_n <- respadmi_n %>% rename(eid=Group.1, fu_year=Group.2, respadmi_byfuyear=x)

finalhospdata_byICD <- merge(finalhospdata_byICD, respadmi_n, by=c("eid", "fu_year"), all.x = T)

finalhospdata_byICD$respadmi_byfuyear[is.na(finalhospdata_byICD$respadmi_byfuyear)] <- 0

summary(finalhospdata_byICD$respadmi_byfuyear) # max no. of RESPIRATORY admissions in any 1 year = 60


##################### 4) no. of other ICD admissions/ year ####################

# ICD chapter 1 admissions
# first make binary column yes/ no for ICD chapter 1 as primary diagnosis
finalhospdata_byICD$ICD1_admi <- 0
finalhospdata_byICD$ICD1_admi[finalhospdata_byICD$ICD_diag01==1] <- 1
table(finalhospdata_byICD$ICD1_admi)

ICD1admi_n <- aggregate(x = finalhospdata_byICD$ICD1_admi, 
                        by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                        FUN = sum)
ICD1admi_n <- ICD1admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD1_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD1admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD1_admibyfuyear[is.na(finalhospdata_byICD$ICD1_admibyfuyear)] <- 0

# ICD chapter 2 admissions
finalhospdata_byICD$ICD2_admi <- 0
finalhospdata_byICD$ICD2_admi[finalhospdata_byICD$ICD_diag01==2] <- 1
table(finalhospdata_byICD$ICD2_admi)

ICD2admi_n <- aggregate(x = finalhospdata_byICD$ICD2_admi, 
                        by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                        FUN = sum)
ICD2admi_n <- ICD2admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD2_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD2admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD2_admibyfuyear[is.na(finalhospdata_byICD$ICD2_admibyfuyear)] <- 0
# ICD chapter 3 admissions (blood diseases)
finalhospdata_byICD$ICD3_admi <- 0
finalhospdata_byICD$ICD3_admi[finalhospdata_byICD$ICD_diag01==3] <- 1
table(finalhospdata_byICD$ICD3_admi)

ICD3admi_n <- aggregate(x = finalhospdata_byICD$ICD3_admi, 
                        by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                        FUN = sum)
ICD3admi_n <- ICD3admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD3_admibyfuyear=x)

finalhospdata_byICD <- merge(finalhospdata_byICD, ICD3admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD3_admibyfuyear[is.na(finalhospdata_byICD$ICD3_admibyfuyear)] <- 0
# ICD chapter 4 admissions
finalhospdata_byICD$ICD4_admi <- 0
finalhospdata_byICD$ICD4_admi[finalhospdata_byICD$ICD_diag01==4] <- 1
table(finalhospdata_byICD$ICD4_admi)
ICD4admi_n <- aggregate(x = finalhospdata_byICD$ICD4_admi, 
                        by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                        FUN = sum)
ICD4admi_n <- ICD4admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD4_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD4admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD4_admibyfuyear[is.na(finalhospdata_byICD$ICD4_admibyfuyear)] <- 0
# ICD chapter 5 admissions
finalhospdata_byICD$ICD5_admi <- 0
finalhospdata_byICD$ICD5_admi[finalhospdata_byICD$ICD_diag01==5] <- 1
table(finalhospdata_byICD$ICD5_admi)
ICD5admi_n <- aggregate(x = finalhospdata_byICD$ICD5_admi, 
                        by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                        FUN = sum)
ICD5admi_n <- ICD5admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD5_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD5admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD5_admibyfuyear[is.na(finalhospdata_byICD$ICD5_admibyfuyear)] <- 0
# ICD chapter 6 admissions
finalhospdata_byICD$ICD6_admi <- 0
finalhospdata_byICD$ICD6_admi[finalhospdata_byICD$ICD_diag01==6] <- 1
table(finalhospdata_byICD$ICD6_admi)
ICD6admi_n <- aggregate(x = finalhospdata_byICD$ICD6_admi, 
                        by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                        FUN = sum)
ICD6admi_n <- ICD6admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD6_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD6admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD6_admibyfuyear[is.na(finalhospdata_byICD$ICD6_admibyfuyear)] <- 0
# ICD chapter 7 admissions
finalhospdata_byICD$ICD7_admi <- 0
finalhospdata_byICD$ICD7_admi[finalhospdata_byICD$ICD_diag01==7] <- 1
table(finalhospdata_byICD$ICD7_admi)
ICD7admi_n <- aggregate(x = finalhospdata_byICD$ICD7_admi, 
                        by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                        FUN = sum)
ICD7admi_n <- ICD7admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD7_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD7admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD7_admibyfuyear[is.na(finalhospdata_byICD$ICD7_admibyfuyear)] <- 0
# ICD chapter 8 admissions
finalhospdata_byICD$ICD8_admi <- 0
finalhospdata_byICD$ICD8_admi[finalhospdata_byICD$ICD_diag01==8] <- 1
table(finalhospdata_byICD$ICD8_admi)
ICD8admi_n <- aggregate(x = finalhospdata_byICD$ICD8_admi, 
                        by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                        FUN = sum)
ICD8admi_n <- ICD8admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD8_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD8admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD8_admibyfuyear[is.na(finalhospdata_byICD$ICD8_admibyfuyear)] <- 0
# ICD chapter 9 admissions
finalhospdata_byICD$ICD9_admi <- 0
finalhospdata_byICD$ICD9_admi[finalhospdata_byICD$ICD_diag01==9] <- 1
table(finalhospdata_byICD$ICD9_admi)
ICD9admi_n <- aggregate(x = finalhospdata_byICD$ICD9_admi, 
                        by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                        FUN = sum)
ICD9admi_n <- ICD9admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD9_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD9admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD9_admibyfuyear[is.na(finalhospdata_byICD$ICD9_admibyfuyear)] <- 0
# ICD chapter 10 admissions (RESPIRATORY)
finalhospdata_byICD$ICD10_admi <- 0
finalhospdata_byICD$ICD10_admi[finalhospdata_byICD$ICD_diag01==10] <- 1
table(finalhospdata_byICD$ICD10_admi)
ICD10admi_n <- aggregate(x = finalhospdata_byICD$ICD10_admi, 
                         by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                         FUN = sum)
ICD10admi_n <- ICD10admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD10_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD10admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD10_admibyfuyear[is.na(finalhospdata_byICD$ICD10_admibyfuyear)] <- 0
# ICD chapter 11 admissions
finalhospdata_byICD$ICD11_admi <- 0
finalhospdata_byICD$ICD11_admi[finalhospdata_byICD$ICD_diag01==11] <- 1
table(finalhospdata_byICD$ICD11_admi)
ICD11admi_n <- aggregate(x = finalhospdata_byICD$ICD11_admi, 
                         by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                         FUN = sum)
ICD11admi_n <- ICD11admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD11_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD11admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD11_admibyfuyear[is.na(finalhospdata_byICD$ICD11_admibyfuyear)] <- 0
# ICD chapter 12 admissions
finalhospdata_byICD$ICD12_admi <- 0
finalhospdata_byICD$ICD12_admi[finalhospdata_byICD$ICD_diag01==12] <- 1
table(finalhospdata_byICD$ICD12_admi)
ICD12admi_n <- aggregate(x = finalhospdata_byICD$ICD12_admi, 
                         by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                         FUN = sum)
ICD12admi_n <- ICD12admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD12_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD12admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD12_admibyfuyear[is.na(finalhospdata_byICD$ICD12_admibyfuyear)] <- 0
# ICD chapter 13 admissions
finalhospdata_byICD$ICD13_admi <- 0
finalhospdata_byICD$ICD13_admi[finalhospdata_byICD$ICD_diag01==13] <- 1
table(finalhospdata_byICD$ICD13_admi)
ICD13admi_n <- aggregate(x = finalhospdata_byICD$ICD13_admi, 
                         by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                         FUN = sum)
ICD13admi_n <- ICD13admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD13_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD13admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD13_admibyfuyear[is.na(finalhospdata_byICD$ICD13_admibyfuyear)] <- 0
# ICD chapter 14 admissions
finalhospdata_byICD$ICD14_admi <- 0
finalhospdata_byICD$ICD14_admi[finalhospdata_byICD$ICD_diag01==14] <- 1
table(finalhospdata_byICD$ICD14_admi)
ICD14admi_n <- aggregate(x = finalhospdata_byICD$ICD14_admi, 
                         by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                         FUN = sum)
ICD14admi_n <- ICD14admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD14_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD14admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD14_admibyfuyear[is.na(finalhospdata_byICD$ICD14_admibyfuyear)] <- 0
# ICD chapter 15 admissions
finalhospdata_byICD$ICD15_admi <- 0
finalhospdata_byICD$ICD15_admi[finalhospdata_byICD$ICD_diag01==15] <- 1
table(finalhospdata_byICD$ICD15_admi)
ICD15admi_n <- aggregate(x = finalhospdata_byICD$ICD15_admi, 
                         by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                         FUN = sum)
ICD15admi_n <- ICD15admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD15_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD15admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD15_admibyfuyear[is.na(finalhospdata_byICD$ICD15_admibyfuyear)] <- 0
# ICD chapter 16 admissions
finalhospdata_byICD$ICD16_admi <- 0
finalhospdata_byICD$ICD16_admi[finalhospdata_byICD$ICD_diag01==16] <- 1
table(finalhospdata_byICD$ICD16_admi)
ICD16admi_n <- aggregate(x = finalhospdata_byICD$ICD16_admi, 
                         by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                         FUN = sum)
ICD16admi_n <- ICD16admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD16_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD16admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD16_admibyfuyear[is.na(finalhospdata_byICD$ICD16_admibyfuyear)] <- 0
# ICD chapter 17 admissions
finalhospdata_byICD$ICD17_admi <- 0
finalhospdata_byICD$ICD17_admi[finalhospdata_byICD$ICD_diag01==17] <- 1
table(finalhospdata_byICD$ICD17_admi)
ICD17admi_n <- aggregate(x = finalhospdata_byICD$ICD17_admi, 
                         by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                         FUN = sum)
ICD17admi_n <- ICD17admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD17_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD17admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD17_admibyfuyear[is.na(finalhospdata_byICD$ICD17_admibyfuyear)] <- 0
# ICD chapter 18 admissions
finalhospdata_byICD$ICD18_admi <- 0
finalhospdata_byICD$ICD18_admi[finalhospdata_byICD$ICD_diag01==18] <- 1
table(finalhospdata_byICD$ICD18_admi)
ICD18admi_n <- aggregate(x = finalhospdata_byICD$ICD18_admi, 
                         by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                         FUN = sum)
ICD18admi_n <- ICD18admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD18_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD18admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD18_admibyfuyear[is.na(finalhospdata_byICD$ICD18_admibyfuyear)] <- 0
# ICD chapter 19 admissions
finalhospdata_byICD$ICD19_admi <- 0
finalhospdata_byICD$ICD19_admi[finalhospdata_byICD$ICD_diag01==19] <- 1
table(finalhospdata_byICD$ICD19_admi)
ICD19admi_n <- aggregate(x = finalhospdata_byICD$ICD19_admi, 
                         by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                         FUN = sum)
ICD19admi_n <- ICD19admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD19_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD19admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD19_admibyfuyear[is.na(finalhospdata_byICD$ICD19_admibyfuyear)] <- 0
# ICD chapter 20 admissions
finalhospdata_byICD$ICD20_admi <- 0
finalhospdata_byICD$ICD20_admi[finalhospdata_byICD$ICD_diag01==20] <- 1
table(finalhospdata_byICD$ICD20_admi)
ICD20admi_n <- aggregate(x = finalhospdata_byICD$ICD20_admi, 
                         by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                         FUN = sum)
ICD20admi_n <- ICD20admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD20_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD20admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD20_admibyfuyear[is.na(finalhospdata_byICD$ICD20_admibyfuyear)] <- 0
# ICD chapter 21 admissions
finalhospdata_byICD$ICD21_admi <- 0
finalhospdata_byICD$ICD21_admi[finalhospdata_byICD$ICD_diag01==21] <- 1
table(finalhospdata_byICD$ICD21_admi)
ICD21admi_n <- aggregate(x = finalhospdata_byICD$ICD21_admi, 
                         by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                         FUN = sum)
ICD21admi_n <- ICD21admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD21_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD21admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD21_admibyfuyear[is.na(finalhospdata_byICD$ICD21_admibyfuyear)] <- 0
# ICD chapter 22 admissions
finalhospdata_byICD$ICD22_admi <- 0
finalhospdata_byICD$ICD22_admi[finalhospdata_byICD$ICD_diag01==22] <- 1
table(finalhospdata_byICD$ICD22_admi)
ICD22admi_n <- aggregate(x = finalhospdata_byICD$ICD22_admi, 
                         by = list(finalhospdata_byICD$eid, finalhospdata_byICD$fu_year), 
                         FUN = sum)
ICD22admi_n <- ICD22admi_n %>% rename(eid=Group.1, fu_year=Group.2, ICD22_admibyfuyear=x)
finalhospdata_byICD <- merge(finalhospdata_byICD, ICD22admi_n, by=c("eid", "fu_year"), all.x = T)
finalhospdata_byICD$ICD22_admibyfuyear[is.na(finalhospdata_byICD$ICD22_admibyfuyear)] <- 0



# merge the count data into the long-form hospital data with each with 16 year rows
# the merge key is eid and fu_year  

ukbhospcost2022_longadmi <- merge(ukbhospcost2022_uni_long, admi_count_wide2, by=c("eid", "fu_year"), all.x = T)
head(ukbhospcost2022_longadmi)


# sort by id and fu year

ukbhospcost2022_longadmi$fu_year <- factor(ukbhospcost2022_longadmi$fu_year, levels = paste0("year_", 1:16))

ukbhospcost2022_longadmi <- ukbhospcost2022_longadmi %>% arrange(eid, fu_year)
head(ukbhospcost2022_longadmi)


# code the admin_n NA's as 0 again
ukbhospcost2022_longadmi$admi_n[is.na(ukbhospcost2022_longadmi$admi_n)] <- 0


# see the final dataset example
head(ukbhospcost2022_longadmi, 20)

ukbhospcost2022_longadmi$admi_duration (# null because we removed admi_duration earlier)
colnames(ukbhospcost2022_longadmi)
  

## subset to remove fu years that are 0.
ukbhospcost2022_longadmi_f <- subset(ukbhospcost2022_longadmi, fu_completion>0)
head(ukbhospcost2022_longadmi_f, 20)
head(ukbhospcost2022_longadmi_f[,c("eid", "fu_year", "fu_completion", "admi_n", "admi_hdays_byfuyear")], 33)


# calculate mean no. admissions by FU year (currently for only those who have hospital records):
aggregate(ukbhospcost2022_longadmi_f$admi_n,list(ukbhospcost2022_longadmi_f$fu_year),mean)



# save ----

ukbhospcost2 <- readRDS(file.path(work_data, "ukbhospcost2.rds"))
# ^ secondary care data with ALL variables

ukbhospcost2022 <- readRDS(file.path(work_data, "ukbhospcost2022.rds"))
# ^ smaller version of secondary care data 
# (now with imputed admidates + disdates = admidate2 + disdate2)
# "fu_days" column is generated from (imputed) admidate2, updating "fu_year"
# ICD-10 respiratory chapters are defined

ukbhospcost2022_uni_long <- readRDS(file.path(work_data, "ukbhospcost2022_uni_long.rds"))
# the above (ukbhospcost2022) converted to wide before being converted to long again
## 1) converted to wide by filtering by max fu_duration (& distinct eid)
## giving a dataset with unique id and follow-up duration
## 16 "year_" columns are added using the (max) fu_duration ~ (still in wide-format)
## 2) converted back to long again, by generating "fu_completion" column
## long-form dataset - each individual has 16 ROWS years 1 - 16
## fu_duration=NULL so this column no longer exists

#saveRDS(finalhospdata, file = file.path(work_data, "finalhospdata.rds"), compress = F)
finalhospdata <- readRDS(file.path(work_data, "finalhospdata.rds"))
finalhospdata_byICD <- readRDS(file.path(work_data, "finalhospdata_byICD.rds"))
# no. admissions, no. days in hospital & costs annually 
# note in these datasets, if there are multiple admissions in a year, then the year is repeated.
