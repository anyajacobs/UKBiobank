install.packages("tidyverse")
library(tidyverse)

load(file.path(work_data, "bd.Rdata"))

# Load data ---------------------------------------------------------------

# read inpatient bulk data  
# diagnosis table
hesin_diag <- read.delim(file.path(raw_data, "hesin_diag.txt"))
# use hesin to incorporate dates
hesin <- read.delim(file.path(raw_data, "hesin.txt"))
# operation table
hesin_oper <- read.delim(file.path(raw_data, "hesin_oper.txt"))
# death tables
death <- read.delim(file.path(raw_data, "death.txt"))
death_cause <- read.delim(file.path(raw_data, "death_cause.txt")) 

## Create UKB_a dataframe ------
# use the main dataset to create a new dataset for recoding variables
# keep the original dataset - "bd" intact
# start with the basic information
UKB_a <- data.frame(id = bd$f.eid, 
                    birth.year = bd$f.34.0.0, 
                    sex = bd$f.31.0.0,
                    age.recruit = bd$f.21022.0.0, 
                    recruit.date = as.Date(bd$f.53.0.0))
# age group 
UKB_a$age.group <- ifelse(UKB_a$age.recruit<45, "<45", 
                          ifelse(UKB_a$age.recruit %in% 45:49, "45-49", 
                                 ifelse(UKB_a$age.recruit %in% 50:54, "50-54",
                                        ifelse(UKB_a$age.recruit %in% 55:59, "55-59",
                                               ifelse(UKB_a$age.recruit %in% 60:64, "60-64",
                                                      ifelse(UKB_a$age.recruit>=65, "65+",NA))))))

# Bulk data ---------------------------------------------------------------

# the variables such as MI.hes, stroke.hes and CRV.hes derived from bulk HES data 
# are basically the same as those from main dataset
# However, bulk data are updated more frequently than main dataset
# I'd recommend to use those from bulk HES data

# code empty cases of icd 10 and icd 9 as NA
hesin_diag$diag_icd10[hesin_diag$diag_icd10==""] <- NA
hesin_diag$diag_icd9[hesin_diag$diag_icd9==""] <- NA

# transform hesin_diag from long to wide by arr_index
# only keep three id vars. and icd10 and icd9
# converting data from long to wide 
hesin_diag <- hesin_diag %>% 
  select("eid", "ins_index", "arr_index", "diag_icd10", "diag_icd9")%>% 
  group_by(eid, ins_index) %>% 
  gather("diag_icd10", "diag_icd9", key = icd, value = code) %>% 
  unite(icd_array, icd, arr_index) %>% 
  spread(icd_array, code)

# remove columns with all NAs
hesin_diag <- hesin_diag[, colSums(is.na(hesin_diag))<nrow(hesin_diag)]

# generate diagnosis date, using episode start date as the primary proxy
# admission date as the secondary proxy.
hesin$date <- hesin$epistart
hesin$date[is.na(hesin$date)] <- hesin[is.na(hesin$date), "admidate"]
# convert integer date to date format
# hesin <- transform(hesin, date=as.Date(as.character(date), "%Y%m%d"))
# RW: 2021-03-23: the 2020 updated hesin use different date format
hesin <- transform(hesin, date=as.Date(as.character(date), "%d/%m/%Y"))

hesin <- hesin %>% select("eid", "ins_index", "date")

# merge hesin into hesin_diag
hesin_diag <- merge(hesin_diag, hesin, by=c("eid", "ins_index"), all.x = T)

# refer to UKB_a for recruitment date
hesin_diag <- merge(hesin_diag, UKB_a[, c("id","recruit.date")], by.x = "eid", 
                    by.y = "id", all.x = T)

# only keep instance where date >= recruitment date
hesin_diag <- hesin_diag[hesin_diag$date>=hesin_diag$recruit.date,]

# again, remove columns with all NAs
hesin_diag <- hesin_diag[, colSums(is.na(hesin_diag))<nrow(hesin_diag)]
# all icd-9 codes are effectively removed because they all happened before recruitment

### code HES event ----
# MI first
# follow-up MI, without old MI
hesin_diag$MI.hes <- 0
# MI codes include: all I21X, all I22X, all I23X, I241in ICD-10
for (i in 0:19) { 
  text2 <- paste0("diag_icd10_",i)
  hesin_diag$MI.hes[grepl("^I21|^I22|^I23|I241", hesin_diag[[text2]])] <- 1 # exclude I252, old MI
}

# only keep those with MI=1
# generate a new instance index, because the old one has been discrete due to filtering
# convert from long to wide
hesin_MI <- hesin_diag %>% filter(MI.hes==1) %>% 
  select("eid", "ins_index", "date") %>% group_by(eid) %>% 
  mutate(ins = 1:n(), ins_index=NULL) %>% spread(ins, date)

# finally select the earliest date among the all
hesin_MI$earliest <- apply(hesin_MI[, -1], 1, function(x) min(x, na.rm = T))
hesin_MI <- hesin_MI[, c(1, 57, 2:56)]

# output
hesin_MI <- hesin_MI %>% select("eid", "earliest") %>% mutate(MI.hes=1) %>% 
  rename(MI.hes.date = earliest)
hesin_MI$MI.hes.date <- as.Date(hesin_MI$MI.hes.date)

UKB_a <- merge(UKB_a, hesin_MI, by.x = "id", by.y = "eid", all.x = T)
UKB_a$MI.hes[is.na(UKB_a$MI.hes)] <- 0
# table(UKB_a$MI.hes, UKB_a$MI.inpatient.post)

# stroke
hesin_diag$stroke.hes <- 0
for (i in 0:19) { 
  text2 <- paste0("diag_icd10_",i)
  hesin_diag$stroke.hes[grepl("^I60|^I61|^I62|^I63|I64", hesin_diag[[text2]])] <- 1
}

hesin_stroke <- hesin_diag %>% filter(stroke.hes==1) %>% 
  select("eid", "ins_index", "date") %>% group_by(eid) %>% 
  mutate(ins = 1:n(), ins_index=NULL) %>% spread(ins, date)

# finally select the earliest date among the all
hesin_stroke$stroke.hes.date <- apply(hesin_stroke[, -1], 1, function(x) min(x, na.rm = T))
hesin_stroke$stroke.hes.date <- as.Date(hesin_stroke$stroke.hes.date)
# output
hesin_stroke <- hesin_stroke %>% select("eid", "stroke.hes.date") %>% mutate(stroke.hes=1)

UKB_a <- merge(UKB_a, hesin_stroke, by.x = "id", by.y = "eid", all.x = T)
UKB_a$stroke.hes[is.na(UKB_a$stroke.hes)] <- 0
# table(UKB_a$stroke.hes, UKB_a$stroke.inpatient.post)

# cancer
hesin_diag$cancer.hes <- 0
for (i in 0:19) { 
  text2 <- paste0("diag_icd10_",i)
  hesin_diag$cancer.hes[grepl("^C", hesin_diag[[text2]]) &
                          !grepl("^C44", hesin_diag[[text2]])] <- 1
}

hesin_cancer <- hesin_diag %>% filter(cancer.hes==1) %>% 
  select("eid", "ins_index", "date") %>% group_by(eid) %>% 
  mutate(ins = 1:n(), ins_index=NULL) %>% spread(ins, date)

# finally select the earliest date among the all
hesin_cancer$cancer.hes.date <- apply(hesin_cancer[, -1], 1, function(x) min(x, na.rm = T))
hesin_cancer$cancer.hes.date <- as.Date(hesin_cancer$cancer.hes.date)
# output
hesin_cancer <- hesin_cancer %>% select("eid", "cancer.hes.date") %>% mutate(cancer.hes=1)

UKB_a <- merge(UKB_a, hesin_cancer, by.x = "id", by.y = "eid", all.x = T)
UKB_a$cancer.hes[is.na(UKB_a$cancer.hes)] <- 0


### apply the same method to operation codes

# code empty cases of oper4 as NA
# oper3 has been coded as NA 
hesin_oper$oper4[hesin_oper$oper4==""] <- NA

hesin_oper <- hesin_oper %>% 
  select("eid", "ins_index", "arr_index", "oper4", "oper3")%>% 
  group_by(eid, ins_index) %>% 
  gather("oper4", "oper3", key = icd, value = code) %>% 
  unite(icd_array, icd, arr_index) %>% spread(icd_array, code)

# remove columns with all NAs
hesin_oper <- hesin_oper[, colSums(is.na(hesin_oper))<nrow(hesin_oper)]

# merge hesin 
hesin_oper <- merge(hesin_oper, hesin, by=c("eid", "ins_index"), all.x = T)

# refer to UKB_a for recruitment date
cruit_date <- UKB_a[, c("id","recruit.date")] 
hesin_oper <- merge(hesin_oper, cruit_date, by.x = "eid", by.y = "id", all.x = T)
rm(cruit_date)

# only keep instance where date >= recruitment date
hesin_oper <- hesin_oper[hesin_oper$date>=hesin_oper$recruit.date,]

# again, remove columns with all NAs
hesin_oper <- hesin_oper[, colSums(is.na(hesin_oper))<nrow(hesin_oper)]
# all opcs 3 codes are effectively removed because they all happened before recruitment

# scan the code row by row to check if there is our interested ones

# CRV
hesin_oper$CRV.hes <- 0
for (i in 0:23) { 
  text2 <- paste0("oper4_",i)
  hesin_oper$CRV.hes[grepl("^K49|^K75|^K76|^K40|^K41|^K42|^K43|^K44|^K45|^K46|K501|K504", 
                           hesin_oper[[text2]])] <- 1 
}

hesin_CRV <- hesin_oper %>% filter(CRV.hes==1) %>% 
  select("eid", "ins_index", "date") %>% group_by(eid) %>% 
  mutate(ins = 1:n(), ins_index=NULL) %>% spread(ins, date)

# finally select the earliest date among the all
hesin_CRV$CRV.hes.date <- apply(hesin_CRV[, -1], 1, function(x) min(x, na.rm = T))
hesin_CRV$CRV.hes.date <- as.Date(hesin_CRV$CRV.hes.date)
# output
hesin_CRV <- hesin_CRV %>% select("eid", "CRV.hes.date") %>% mutate(CRV.hes=1)

UKB_a <- merge(UKB_a, hesin_CRV, by.x = "id", by.y = "eid", all.x = T)
UKB_a$CRV.hes[is.na(UKB_a$CRV.hes)] <- 0

# remove hesin data
rm(hesin_diag, hesin, hesin_oper, hesin_MI, hesin_stroke, hesin_CRV, hesin_cancer)


### death based on bulk data --------------------------------------------------

# directly use bulk data
# bulk deaths include all main dataset deaths and have a few extra

# convert integer date to date format
death <- transform(death, date=as.Date(as.character(date_of_death), "%d/%m/%Y"))
death <- death %>% select(eid, ins_index, date)

death_cause_prim <- death_cause %>% filter(level==1)
death_cause_sec <- death_cause %>% filter(level==2)

# transform from long to wide
death_cause_prim <- death_cause_prim %>% 
  select("eid","ins_index","arr_index","cause_icd10") %>% 
  group_by(eid, ins_index) %>% 
  gather("cause_icd10", key = icd, value = code) %>% 
  unite(icd_array, icd, arr_index) %>% 
  spread(icd_array, code) 

death_cause_sec <- death_cause_sec %>% 
  select("eid","ins_index","arr_index","cause_icd10") %>% 
  group_by(eid, ins_index) %>% 
  gather("cause_icd10", key = icd, value = code) %>% 
  unite(icd_array, icd, arr_index) %>% 
  spread(icd_array, code) 

# merge
death <- merge(death, death_cause_prim, by=c("eid", "ins_index"), all.x = T)
colnames(death)[4] <- "ICD10_prim"
death <- merge(death, death_cause_sec, by=c("eid", "ins_index"), all.x = T)
colnames(death) <- sub("cause_icd10", "ICD10_sec", colnames(death))

# vascular death
# all ICD-10 I category: circulatory
# all ICD-10 R category: unclassfied elsewhere
death$VD <- 0
death$VD[grepl("^I", death$ICD10_prim) | 
           grepl("^R", death$ICD10_prim)] <- 1
# on the condition of CVD as a secondary cause 
# Y832: anastomosis, bypass or graft
# Y835: amputation of limb
# W19: unspecified fall
for (i in 1:14) {
  text1 <- paste0("ICD10_sec_",i)
  death$VD[grepl("^I", death[[text1]]) & (death$ICD10_prim %in% 
                                            c("Y832", "Y835") | grepl("^W19", death$ICD10_prim))] <- 1
}

death_VD <- death %>% filter(VD==1) %>% 
  select(eid, date) %>% 
  distinct(eid, .keep_all = T) %>% 
  mutate(death.vascular=1) %>% 
  rename(death.vascular.date=date)

# merge to UKB_a
UKB_a <- merge(UKB_a, death_VD, by.x = "id", by.y = "eid", all.x = T)
UKB_a$death.vascular[is.na(UKB_a$death.vascular)] <- 0

# all death
death_all <- death %>% select(eid, date) %>% 
  distinct(eid, .keep_all = T) %>% 
  mutate(death.allcause=1) %>% 
  rename(death.date=date)

UKB_a <- merge(UKB_a, death_all, by.x = "id", by.y = "eid", all.x = T)
UKB_a$death.allcause[is.na(UKB_a$death.allcause)] <- 0

# non-vascular death
UKB_a$death.nonvascular <- ifelse(UKB_a$death.allcause==1 & 
                                    UKB_a$death.vascular==0, 1, 0)

UKB_a$death.nonvascular.date <- UKB_a$death.date
UKB_a$death.nonvascular.date[UKB_a$death.vascular==1] <- NA



# -------------------------------------------------------------------------------------- # 
# Main dataset #####
# -------------------------------------------------------------------------------------- # 


### MI ----------------------------------------------------------------------

# baseline MI (myocardial infarction) use UKB algorithm 

# combination of verbal interview and inpatient record before recruitment
UKB_a$MI.baseline <- ifelse(!is.na(bd$f.42001.0.0) & 
                              bd$f.42000.0.0<UKB_a$recruit.date, 1, 0)

#### incident MI
# first
# MI death
# death with a primary or secondary reason as MI
# bd$f.40001 has two instance although
# cases in the second instance is actually included in the first instance

# # primary reason as MI
# UKB_a$MI.death <- 0 # RW 2020-10-12
# for (i in 0:1) {
#   text2 <- paste0("f.40001.",i, ".", 0)
#   UKB_a$MI.death[grepl("^I21|^I22|^I23|I241", bd[[text2]])] <- 1
# }
# 
# # secondary reason as MI
# for (i in 0:1) {
#   for (j in 0: 13) {
#     text2 <- paste0("f.40002.",i, ".", j)
#     UKB_a$MI.death[grepl("^I21|^I22|^I23|I241", bd[[text2]])] <- 1
#   }
# }
# 
# UKB_a$MI.death.date <- as.Date(bd$f.40000.0.0)
# UKB_a$MI.death.date[UKB_a$MI.death==0] <- NA

# use bulk death data instead
death$MI.death <- 0
# primary reason
death$MI.death[grepl("^I21|^I22|^I23|I241", death$ICD10_prim)] <- 1
# secondary reason
for (i in 1:14) {
  text1 <- paste0("ICD10_sec_",i)
  death$MI.death[grepl("^I21|^I22|^I23|I241", death[[text1]])] <- 1
}

MI_death <- death %>% filter(MI.death==1) %>% 
  select(eid, date, MI.death) %>% 
  distinct(eid, .keep_all = T) %>% 
  rename(MI.death.date=date)

# merge to UKB_a
UKB_a <- merge(UKB_a, MI_death, by.x = "id", by.y = "eid", all.x = T)
UKB_a$MI.death[is.na(UKB_a$MI.death)] <- 0


# second
# first occurrence date, including primary care
# check I21,22,23 only, since not all I24 belongs to MI
# I21: 131298
# I22: 131300
# I23: 131302
# because it lack some code, just used for complement only

temp <- data.frame(I21=bd$f.131298.0.0, I22=bd$f.131300.0.0, I23=bd$f.131302.0.0)
# generate MI ever
UKB_a$MI.fo <- 0
UKB_a$MI.fo[rowSums(!is.na(temp)) > 0] <- 1
# given one person can have more than one dates for MI, we should choose the first one
# before this, check "1901-01-01", "1902-02-02", "1903-03-03" and "2037-07-07". recode them as NA
# summary(temp$I21)
# summary(temp$I22)
# summary(temp$I23)
# answer is NO

earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(UKB_a), but keeps the same row names as names
# join using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 

UKB_a$MI.fo.date <- temp$date
UKB_a$MI.fo.date <- as.Date(UKB_a$MI.fo.date)

# generate pre-recruitment and post-recruitment
UKB_a$MI.fo.pre <- ifelse(is.na(UKB_a$MI.fo.date), 0,
                          ifelse(UKB_a$MI.fo.date<UKB_a$recruit.date, 1, 0))
UKB_a$MI.fo.pre.date <- UKB_a$MI.fo.date
UKB_a$MI.fo.pre.date[UKB_a$MI.fo.pre==0] <- NA

UKB_a$MI.fo.post <- ifelse(is.na(UKB_a$MI.fo.date), 0,
                           ifelse(UKB_a$MI.fo.date>=UKB_a$recruit.date, 1, 0))
UKB_a$MI.fo.post.date <- UKB_a$MI.fo.date
UKB_a$MI.fo.post.date[UKB_a$MI.fo.post==0] <- NA


# third
# combine MI death, inpatient MI record and first occurrence post
# date choose whichever the earliest
# data from HES inpatient records contain all inpatient records in the main dataset 
# we'd use HES records instead of that from main dataset inpatient records.  

#UKB_a$MI.all.post <- UKB_a$MI.inpatient.post
UKB_a$MI.all.post <- UKB_a$MI.hes
UKB_a$MI.all.post[UKB_a$MI.death==1] <- 1
UKB_a$MI.all.post[UKB_a$MI.fo.post==1] <- 1

# pick up the earliest date among the three
earliest <- apply(UKB_a[UKB_a$MI.hes==1 | UKB_a$MI.death==1 | UKB_a$MI.fo.post==1, 
                        c("MI.hes.date", "MI.death.date", "MI.fo.post.date")], 1, function(x) min(x, na.rm = T))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))

temp <- as.data.frame(UKB_a$MI.hes.date)
temp$key <- as.numeric(rownames(temp))
temp <- merge(temp, earliest, all = T) 
# check if the inpatient.date is always the earliest
# sum(temp$`UKB_a$MI.inpatient.date`!=temp$date, na.rm = T)

temp$date <- as.Date(temp$date)

# whichever the earliest
UKB_a$MI.all.date <- temp$date

# finally, also complement MI baseline
UKB_a$MI.baseline[UKB_a$MI.fo.pre==1] <- 1



### stroke ------------------------------------------------------------------

# first inpatient record
# We finally use bulk HES data, but don't mute the follows, because baseline stroke 
# refers to the information.
# baseline MI basically follows the UKB algorithm so do not need inpatient data

UKB_a$stroke.inpatient <- 0
# MI codes include: all CXXX in ICD-10
for (i in 0:212) { # f.41270 has 213 arrays
  text2 <- paste0("f.41270.0.",i)
  UKB_a$stroke.inpatient[grepl("^I60|^I61|^I62|^I63|I64 ", bd[[text2]])] <- 1
}
# ICD-9 codes all happened years before recruitment.
# they are included here only to supplement identifying baseline stroke
# so they are not considered in identifying dates
for (i in 0:46) { # f.41271 has 47 arrays
  text2 <- paste0("f.41271.0.",i)
  UKB_a$stroke.inpatient[bd[[text2]]=="4309" |
                           bd[[text2]]=="4319" |
                           bd[[text2]]=="4349" |
                           bd[[text2]]=="4369" |
                           bd[[text2]]=="4320" |
                           bd[[text2]]=="4321"] <- 1
  
}
# diagnosis date
# ICD-10
# create a new date frame to collect all dates linked to stroke diagnoses
# one individual may have more than one diagnoses and linked dates
temp <- as.data.frame(UKB_a$stroke.inpatient)
# above is an easy way to create a date frame with the same structure with UKB_a
for (i in 0:212) { # f.41280 has 213 arrays, corresponding to 41270
  text1 <- paste0("stroke.inpatient.date.", i)
  text2 <- paste0("f.41270.0.",i)
  text3 <- paste0("f.41280.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^I60|^I61|^I62|^I63|I64", bd[[text2]]) &
                  !is.na(bd[[text2]])] <-
    bd[[text3]][grepl("^I60|^I61|^I62|^I63|I64", bd[[text2]]) &
                  !is.na(bd[[text2]])]
}

# if there are multiple date for stroke ICD, use the earliest one after recruitment date
# get rid of date before recruitment date and first column
for (i in 0:212){
  text1 <- paste0("stroke.inpatient.date.", i)
  temp[[text1]][temp[[text1]] < UKB_a$recruit.date] <- NA
}
temp$`UKB_a$stroke.inpatient` <- NULL
# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# select the earliest date
# the following command only apply to rows with at least one non-missing values
earliest <- apply(temp[rowSums(is.na(temp)) <
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(UKB_a), but keeps the same row
# names as names
# join to temp using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T)
# give the value to UKB_a
UKB_a$stroke.inpatient.date <- temp$date
UKB_a$stroke.inpatient.date <- as.Date(UKB_a$stroke.inpatient.date)

# stroke.inpatient.date is only post-recruitment

# dates linked to ICD-9 codes are no later than 1996-03-30
# don't work on ICD-9 dates

UKB_a$stroke.inpatient.post <- UKB_a$stroke.inpatient
UKB_a$stroke.inpatient.post[is.na(UKB_a$stroke.inpatient.date)] <- 0

# baseline stroke
# base on f.20002, self-reported non-cancer illness
text1 <- "stroke.baseline"
UKB_a[[text1]] <- 0
for (j in 0:33) { # gather all records from multiple arrays
  text2 <- paste0("f.20002.0.", j)
  UKB_a[[text1]][bd[[text2]]== 1081 |
                   bd[[text2]]== 1583 |
                   bd[[text2]]== 1086 |
                   bd[[text2]]== 1491 |
                   bd[[text2]]== 1083] <- 1
}

# and inpatient record before recruitment
# is.na(UKB_a$stroke.inpatient.date) means the date precede recruit date
UKB_a$stroke.baseline[UKB_a$stroke.inpatient==1 & is.na(UKB_a$stroke.inpatient.date)] <- 1

# second
# stroke death
# death with a primary or secondary reason as stroke
# bd$f.40001 has two instance although, 
# cases in the second instance is actually included in the first instance
# primary reason as stroke
# UKB_a$stroke.death <- 0 # RW 2020-10-12
# for (i in 0:1) {
#   text2 <- paste0("f.40001.",i, ".", 0)
#   UKB_a$stroke.death[grepl("^I60|^I61|^I62|^I63|I64", bd[[text2]])] <- 1
# }
# 
# # secondary reason as stroke
# for (i in 0:1) {
#   for (j in 0: 13) {
#     text2 <- paste0("f.40002.",i, ".", j)
#     UKB_a$stroke.death[grepl("^I60|^I61|^I62|^I63|I64", bd[[text2]])] <- 1
#   }
# }
# 
# UKB_a$stroke.death.date <- bd$f.40000.0.0
# UKB_a$stroke.death.date[UKB_a$stroke.death==0] <- NA

# use bulk death data instead
death$stroke.death <- 0
# primary reason
death$stroke.death[grepl("^I60|^I61|^I62|^I63|I64", death$ICD10_prim)] <- 1
# secondary reason
for (i in 1:14) {
  text1 <- paste0("ICD10_sec_",i)
  death$stroke.death[grepl("^I60|^I61|^I62|^I63|I64", death[[text1]])] <- 1
}

stroke_death <- death %>% filter(stroke.death==1) %>% 
  select(eid, date, stroke.death) %>% 
  distinct(eid, .keep_all = T) %>% 
  rename(stroke.death.date=date)

# merge to UKB_a
UKB_a <- merge(UKB_a, stroke_death, by.x = "id", by.y = "eid", all.x = T)
UKB_a$stroke.death[is.na(UKB_a$stroke.death)] <- 0


# third
# first occurrence date, including primary care
# check I21,22,23 only, since not all I24 belongs to MI
# I60: 131360
# I61: 131362
# I62: 131364
# I63: 131366
# I64: 131368

temp <- data.frame(I60=bd$f.131360.0.0, I61=bd$f.131362.0.0, I62=bd$f.131364.0.0, 
                   I63=bd$f.131366.0.0, I64=bd$f.131368.0.0)
# generate MI ever
UKB_a$stroke.fo <- 0
UKB_a$stroke.fo[rowSums(!is.na(temp)) > 0] <- 1
# given one person can have more than one dates for stroke, we should choose the first one
# before this, check "1901-01-01", "1902-02-02", "1903-03-03" and "2037-07-07". recode them as NA
# summary(temp$I60)
# summary(temp$I61)
# summary(temp$I62)
# summary(temp$I63)
# summary(temp$I64)
# answer is NO

earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(UKB_a), but keeps the same row names as names
# join using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 

UKB_a$stroke.fo.date <- temp$date
UKB_a$stroke.fo.date <- as.Date(UKB_a$stroke.fo.date)

# generate pre-recruitment and post-recruitment
UKB_a$stroke.fo.pre <- ifelse(is.na(UKB_a$stroke.fo.date), 0,
                              ifelse(UKB_a$stroke.fo.date<UKB_a$recruit.date, 1, 0))
UKB_a$stroke.fo.pre.date <- UKB_a$stroke.fo.date
UKB_a$stroke.fo.pre.date[UKB_a$stroke.fo.pre==0] <- NA

UKB_a$stroke.fo.post <- ifelse(is.na(UKB_a$stroke.fo.date), 0,
                               ifelse(UKB_a$stroke.fo.date>=UKB_a$recruit.date, 1, 0))
UKB_a$stroke.fo.post.date <- UKB_a$stroke.fo.date
UKB_a$stroke.fo.post.date[UKB_a$stroke.fo.post==0] <- NA

# forth
# combine stroke death, inpatient stroke record and first occurrence post
# date choose whichever the earliest

# UKB_a$stroke.all.post <- UKB_a$stroke.inpatient.post
UKB_a$stroke.all.post <- UKB_a$stroke.hes
UKB_a$stroke.all.post[UKB_a$stroke.death==1] <- 1
UKB_a$stroke.all.post[UKB_a$stroke.fo.post==1] <- 1

# pick up the earliest date among the three
earliest <- apply(UKB_a[UKB_a$stroke.hes==1 | UKB_a$stroke.death==1 | UKB_a$stroke.fo.post==1, 
                        c("stroke.hes.date", "stroke.death.date", "stroke.fo.post.date")], 1, function(x) min(x, na.rm = T))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))

temp <- as.data.frame(UKB_a$stroke.hes.date)
temp$key <- as.numeric(rownames(temp))
temp <- merge(temp, earliest, all = T) 

temp$date <- as.Date(temp$date)

# whichever the earliest
UKB_a$stroke.all.date <- temp$date

# finally also complement stroke baseline
UKB_a$stroke.baseline[UKB_a$stroke.fo.pre==1] <- 1


### CRV ---------------------------------------------------------------------

# only use hospital diagnoses
# no baseline CRV

UKB_a$CRV.inpatient <- 0
# MI codes include: all K49X, K501, all K75X, all K76X, all K40X
# all K41X, all K42X, all K43X, all K44X, all K45X, all K46X in OPCS-4
for (i in 0:116) { # f.41272 has 117 arrays
  text2 <- paste0("f.41272.0.",i)
  UKB_a$CRV.inpatient[grepl("^K49|^K75|^K76|^K40|^K41|^K42|^K43|^K44|^K45|^K46|K501|K504", bd[[text2]])] <- 1
}

# CRV date
temp <- as.data.frame(UKB_a$CRV.inpatient)
# OPCS-4 only
# dates linked to OPCS-3 codes are no later than 1989, so don't use them
for (i in 0:116) { # f.41282 has 117 arrays, corresponding to 41272
  text1 <- paste0("CRV.inpatient.date.", i)
  text2 <- paste0("f.41272.0.",i)
  text3 <- paste0("f.41282.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^K49|^K75|^K76|^K40|^K41|^K42|^K43|^K44|^K45|^K46|K501|K504", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("^K49|^K75|^K76|^K40|^K41|^K42|^K43|^K44|^K45|^K46|K501|K504", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}

for (i in 0:116){
  text1 <- paste0("CRV.inpatient.date.", i)
  temp[[text1]][temp[[text1]] < UKB_a$recruit.date] <- NA
}
temp$`UKB_a$CRV.inpatient` <- NULL
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 
# give the value to UKB_a
UKB_a$CRV.inpatient.date <- temp$date
UKB_a$CRV.inpatient.date <- as.Date(UKB_a$CRV.inpatient.date)

# post-recruitment 
UKB_a$CRV.inpatient.post <- UKB_a$CRV.inpatient
UKB_a$CRV.inpatient.post[is.na(UKB_a$CRV.inpatient.date)] <- 0

# CRV.hes include all CRV.inpatient.post
UKB_a$CRV.all.post <- UKB_a$CRV.hes

# so just pick up the earliest date among the two
earliest <- apply(UKB_a[UKB_a$CRV.inpatient.post==1 | UKB_a$CRV.hes, 
                        c("CRV.inpatient.date", "CRV.hes.date")], 1, function(x) min(x, na.rm = T))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))

temp <- as.data.frame(UKB_a$CRV.hes.date)
temp$key <- as.numeric(rownames(temp))
temp <- merge(temp, earliest, all = T) 

temp$date <- as.Date(temp$date)

UKB_a$CRV.all.date <- temp$date

# check 
# identical(UKB_a$CRV.all.date, UKB_a$CRV.hes.date)
# the two are the same, so we actually only need to use CRV.hes data


### cancer ------------------------------------------------------------------

# first
# hospital diagnoses

UKB_a$cancer.inpatient <- 0
# MI codes include: all CXXX in ICD-10
for (i in 0:212) { # f.41270 has 213 arrays
  text2 <- paste0("f.41270.0.",i)
  # all malignant neoplasms but C44 family of non-melanoma skin cancers
  UKB_a$cancer.inpatient[grepl("^C", bd[[text2]]) &
                           !grepl("^C44", bd[[text2]])] <- 1
}
# 140X-208X in ICD-9
# ICD-9 codes all happened years before recruitment. 
# they are included here only to supplement identifying baseline cancer 
# so they are not considered in identifying dates
for (i in 0:46) { # f.41271 has 47 arrays
  text2 <- paste0("f.41271.0.",i)
  # all malignant neoplasms but 173 family of non-melanoma skin cancers
  UKB_a$cancer.inpatient[(grepl("^14", bd[[text2]]) | 
                            grepl("^15", bd[[text2]]) |
                            grepl("^16", bd[[text2]]) |
                            grepl("^17", bd[[text2]]) |
                            grepl("^18", bd[[text2]]) |
                            grepl("^19", bd[[text2]]) |
                            grepl("^20", bd[[text2]])) &
                           !grepl("^173", bd[[text2]])] <- 1
}
# cancer diagnosis date
# ICD-10
# create a new date frame to collect all dates linked to cancer diagnoses
# one individual may have more than one diagnoses and linked dates
temp <- as.data.frame(UKB_a$cancer.inpatient)
# above is an easy way to create a date frame with the same structure with UKB_a
for (i in 0:212) { # f.41280 has 213 arrays, corresponding to 41270
  text1 <- paste0("cancer.inpatient.date.", i)
  text2 <- paste0("f.41270.0.",i)
  text3 <- paste0("f.41280.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}
# if there are multiple date for cancer ICD, use the earliest one after recruitment date
# get rid of date before recruitment date and first column
for (i in 0:212){
  text1 <- paste0("cancer.inpatient.date.", i)
  temp[[text1]][temp[[text1]] < UKB_a$recruit.date] <- NA
}
temp$`UKB_a$cancer.inpatient` <- NULL
# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# select the earliest date
# the following command only apply to rows with at least one non-missing values
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(UKB_a), but keeps the same row
# join to temp using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 
# give the value to UKB_a
UKB_a$cancer.inpatient.post.date <- temp$date
UKB_a$cancer.inpatient.post.date <- as.Date(UKB_a$cancer.inpatient.post.date)

UKB_a$cancer.inpatient.post <- UKB_a$cancer.inpatient
UKB_a$cancer.inpatient.post[is.na(UKB_a$cancer.inpatient.post.date)] <- 0


# inpatient records date before recruitment

# ICD-10
temp <- as.data.frame(UKB_a$cancer.inpatient)
# above is an easy way to create a date frame with the same structure with UKB_a
for (i in 0:212) { # f.41280 has 213 arrays, corresponding to 41270
  text1 <- paste0("cancer.inpatient.date.", i)
  text2 <- paste0("f.41270.0.",i)
  text3 <- paste0("f.41280.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}

for (i in 0:46) { # f.41271 has 47 arrays
  text1 <- paste0("cancer.inpatient.date.", i+213)
  text2 <- paste0("f.41271.0.",i)
  text3 <- paste0("f.41281.0.",i)
  # all malignant neoplasms but 173 family of non-melanoma skin cancers
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][(grepl("^14", bd[[text2]]) | 
                   grepl("^15", bd[[text2]]) |
                   grepl("^16", bd[[text2]]) |
                   grepl("^17", bd[[text2]]) |
                   grepl("^18", bd[[text2]]) |
                   grepl("^19", bd[[text2]]) |
                   grepl("^20", bd[[text2]])) &
                  !grepl("^173", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][(grepl("^14", bd[[text2]]) | 
                   grepl("^15", bd[[text2]]) |
                   grepl("^16", bd[[text2]]) |
                   grepl("^17", bd[[text2]]) |
                   grepl("^18", bd[[text2]]) |
                   grepl("^19", bd[[text2]]) |
                   grepl("^20", bd[[text2]])) &
                  !grepl("^173", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}

# if there are multiple date for cancer ICD, use the earliest one before recruitment date
# get rid of date after recruitment date and first column
for (i in 0:259){
  text1 <- paste0("cancer.inpatient.date.", i)
  temp[[text1]][temp[[text1]] >= UKB_a$recruit.date] <- NA
}
temp$`UKB_a$cancer.inpatient` <- NULL

# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# select the earliest date
# the following command only apply to rows with at least one non-missing values
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(UKB_a), but keeps the same row
# names as names
# join to temp using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 
# give the value to UKB_a
UKB_a$cancer.inpatient.pre.date <- temp$date
UKB_a$cancer.inpatient.pre.date <- as.Date(UKB_a$cancer.inpatient.pre.date)

# there is a date, there is a baseline cancer
UKB_a$cancer.inpatient.pre <- UKB_a$cancer.inpatient
UKB_a$cancer.inpatient.pre[is.na(UKB_a$cancer.inpatient.pre.date)] <- 0

# second
# cancer registry

UKB_a$cancer.registry <- 0
# MI codes include: all CXXX in ICD-10
for (i in 0:16) { # f.40006: type of cancer ICD-10 has 17 instances 
  text2 <- paste0("f.40006.",i,".0")
  # all malignant neoplasms but C44 family of non-melanoma skin cancers
  UKB_a$cancer.registry[grepl("^C", bd[[text2]]) &
                          !grepl("^C44", bd[[text2]])] <- 1
}
# 140X-208X in ICD-9
# ICD-9 codes all happened years before recruitment. 
# they are included here only to supplement identifying baseline cancer 
# so they are not considered in identifying dates
for (i in 0:14) { # f.40013: type of cancer ICD-9 has 15 instances
  text2 <- paste0("f.40013.",i, ".0")
  # all malignant neoplasms but 173 family of non-melanoma skin cancers
  UKB_a$cancer.registry[(grepl("^14", bd[[text2]]) | 
                           grepl("^15", bd[[text2]]) |
                           grepl("^16", bd[[text2]]) |
                           grepl("^17", bd[[text2]]) |
                           grepl("^18", bd[[text2]]) |
                           grepl("^19", bd[[text2]]) |
                           grepl("^20", bd[[text2]])) &
                          !grepl("^173", bd[[text2]])] <- 1
}
# cancer diagnosis date
# ICD-10
# create a new date frame to collect all dates linked to cancer diagnoses
# one individual may have more than one diagnoses and linked dates
temp <- as.data.frame(UKB_a$cancer.registry)
# above is an easy way to create a date frame with the same structure with UKB_a
for (i in 0:16) { # f.40005: date of cancer diagnosis has 213 instances
  text1 <- paste0("cancer.registry.date.", i)
  text2 <- paste0("f.40006.",i,".0")
  text3 <- paste0("f.40005.",i,".0")
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}
# if there are multiple date for cancer ICD, use the earliest one after recruitment date
# get rid of date before recruitment date and first column
for (i in 0:16){
  text1 <- paste0("cancer.registry.date.", i)
  temp[[text1]][temp[[text1]] < UKB_a$recruit.date] <- NA
}
temp$`UKB_a$cancer.registry` <- NULL
# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# select the earliest date
# the following command only apply to rows with at least one non-missing values
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(UKB_a), but keeps the same row
# names as names
# join to temp using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 
# give the value to UKB_a
UKB_a$cancer.registry.post.date <- temp$date
UKB_a$cancer.registry.post.date <- as.Date(UKB_a$cancer.registry.post.date)

UKB_a$cancer.registry.post <- UKB_a$cancer.registry
UKB_a$cancer.registry.post[is.na(UKB_a$cancer.registry.post.date)] <- 0


# cancer register date before recruitment

# ICD-10

temp <- as.data.frame(UKB_a$cancer.registry)
# above is an easy way to create a date frame with the same structure with UKB_a
for (i in 0:16) { # f.40005: date of cancer diagnosis has 17 instances
  text1 <- paste0("cancer.registry.date.", i)
  text2 <- paste0("f.40006.",i,".0")
  text3 <- paste0("f.40005.",i,".0")
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}
# ICD-9
for (i in 0:14) { # f.40013 has 15 arrays
  text1 <- paste0("cancer.registry.date.", i+17)
  text2 <- paste0("f.40013.",i,".0")
  text3 <- paste0("f.40005.",i,".0")
  # all malignant neoplasms but 173 family of non-melanoma skin cancers
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][(grepl("^14", bd[[text2]]) | 
                   grepl("^15", bd[[text2]]) |
                   grepl("^16", bd[[text2]]) |
                   grepl("^17", bd[[text2]]) |
                   grepl("^18", bd[[text2]]) |
                   grepl("^19", bd[[text2]]) |
                   grepl("^20", bd[[text2]])) &
                  !grepl("^173", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][(grepl("^14", bd[[text2]]) | 
                   grepl("^15", bd[[text2]]) |
                   grepl("^16", bd[[text2]]) |
                   grepl("^17", bd[[text2]]) |
                   grepl("^18", bd[[text2]]) |
                   grepl("^19", bd[[text2]]) |
                   grepl("^20", bd[[text2]])) &
                  !grepl("^173", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}
# if there are multiple date for cancer ICD, use the earliest one after recruitment date
# get rid of date after recruitment date and first column
for (i in 0:31){
  text1 <- paste0("cancer.registry.date.", i)
  temp[[text1]][temp[[text1]] >= UKB_a$recruit.date] <- NA
}
temp$`UKB_a$cancer.registry` <- NULL

# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# select the earliest date
# the following command only apply to rows with at least one non-missing values
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(UKB_a), but keeps the same row
# names as names
# join to temp using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 
# give the value to UKB_a
UKB_a$cancer.registry.pre.date <- temp$date
UKB_a$cancer.registry.pre.date <- as.Date(UKB_a$cancer.registry.pre.date)

# cancer registry before recruitment
UKB_a$cancer.registry.pre <- UKB_a$cancer.registry
UKB_a$cancer.registry.pre[is.na(UKB_a$cancer.registry.pre.date)] <- 0


# third
# baseline cancer
# base on f.20001, self-reported cancer
# exclude non-melanoma skin cancers

text1 <- "cancer.baseline"
UKB_a[[text1]] <- 0
for (j in 0:5) {
  text2 <- paste0("f.20001.0.", j)
  UKB_a[[text1]][!is.na(bd[[text2]]) & 
                   bd[[text2]] != 1060 &
                   bd[[text2]] != 1061 &
                   bd[[text2]] != 1062 &
                   bd[[text2]] != 1073] <- 1
}

# cancer baseline date 
temp <- as.data.frame(UKB_a$cancer.baseline)
for (j in 0:5) {
  text1 <- paste0("cancer.baseline.date.", j)
  text2 <- paste0("f.20001.0.", j)
  text3 <- paste0("f.20006.0.", j)
  temp[[text1]] <- NA
  temp[[text1]][!is.na(bd[[text2]]) & 
                  bd[[text2]] != 1060 &
                  bd[[text2]] != 1061 &
                  bd[[text2]] != 1062 &
                  bd[[text2]] != 1073] <- 
    bd[[text3]][!is.na(bd[[text2]]) & 
                  bd[[text2]] != 1060 &
                  bd[[text2]] != 1061 &
                  bd[[text2]] != 1062 &
                  bd[[text2]] != 1073]
}
# select the earliest
temp$`UKB_a$cancer.baseline` <- NULL
# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# select the earliest date
# the following command only apply to rows with at least one non-missing values
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(UKB_a), but keeps the same row
# names as names
# join to temp using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 
# give the value to UKB_a
UKB_a$cancer.baseline.date <- temp$date
# recode date=-1 NA
UKB_a$cancer.baseline.date[UKB_a$cancer.baseline.date==-1] <- NA

# RW 2020-12-17
# because baseline date is only year, assume it is the mid day of the year
# the computer interestingly assume the date as today, though year is different
# so we minus the difference between the current date XXXX-XX-XX and the XXXX-01-01
# then add 182 to the mid of the year
dif <- as.numeric(Sys.Date() - as.Date(paste0(substr(as.character(Sys.Date()), 1, 4), "-01-01")))
UKB_a$cancer.baseline.date <- as.Date(as.character(UKB_a$cancer.baseline.date), format = "%Y")-dif+182

# a problem is that some baseline cancer dates exceed the recruitment date (n =374) 
# we have to assume them 1 day before the recruitment date
UKB_a$cancer.baseline.date[!is.na(UKB_a$cancer.baseline.date) & UKB_a$cancer.baseline.date > UKB_a$recruit.date] <- 
  UKB_a$recruit.date[!is.na(UKB_a$cancer.baseline.date) & UKB_a$cancer.baseline.date > UKB_a$recruit.date]-1

# and inpatient record before recruitment
# is.na(UKB_a$cancer.inpatient.post.date) means the date precede recruit date
UKB_a$cancer.baseline.all <- UKB_a$cancer.baseline
# inpatient
UKB_a$cancer.baseline.all[UKB_a$cancer.inpatient.pre==1] <- 1
# cancer register
UKB_a$cancer.baseline.all[UKB_a$cancer.registry.pre==1] <- 1

# fourth
# cancer death
# death with a primary or secondary reason as cancer
# bd$f.40001 has two instance although, 
# cases in the second instance is actually included in the first instance
# primary reason as cancer
# UKB_a$cancer.death <- 0 # RW 2020-10-12
# for (i in 0:1) {
#   text2 <- paste0("f.40001.",i, ".", 0)
#   UKB_a$cancer.death[grepl("^C", bd[[text2]]) &
#                      !grepl("^C44", bd[[text2]])] <- 1
# }
# # secondary reason as cancer
# for (i in 0:1) {
#   for (j in 0: 13) {
#     text2 <- paste0("f.40002.",i, ".", j)
#     UKB_a$cancer.death[grepl("^C", bd[[text2]]) &
#                        !grepl("^C44", bd[[text2]])] <- 1
#     
#   }
# }
# 
# UKB_a$cancer.death.date <- as.Date(bd$f.40000.0.0)
# UKB_a$cancer.death.date[UKB_a$cancer.death==0] <- NA

# use bulk death data

death$cancer.death <- 0
# primary reason as cancer
death$cancer.death[grepl("^C", death$ICD10_prim) & 
                     !grepl("^C44", death$ICD10_prim)] <- 1
# secondary reason as cancer
for (i in 1:14) {
  text1 <- paste0("ICD10_sec_",i)
  death$cancer.death[grepl("^C", death[[text1]]) & !grepl("^C44", death[[text1]])] <- 1
}

cancer_death <- death %>% filter(cancer.death==1) %>% 
  select(eid, date, cancer.death) %>% 
  distinct(eid, .keep_all = T) %>% 
  rename(cancer.death.date=date)

# merge to UKB_a
UKB_a <- merge(UKB_a, cancer_death, by.x = "id", by.y = "eid", all.x = T)
UKB_a$cancer.death[is.na(UKB_a$cancer.death)] <- 0


# fifth
# combine cancer death, inpatient cancer record and first occurrence post
# date choose whichever the earliest
# as MI and stroke, cancer.hes replaces cancer.inpatient.post

# UKB_a$cancer.all.post <- UKB_a$cancer.inpatient.post
UKB_a$cancer.all.post <- UKB_a$cancer.hes
UKB_a$cancer.all.post[UKB_a$cancer.death==1] <- 1
UKB_a$cancer.all.post[UKB_a$cancer.registry.post==1] <- 1

# pick up the earliest date among the three
earliest <- apply(UKB_a[UKB_a$cancer.hes==1 | UKB_a$cancer.death==1 | 
                          UKB_a$cancer.registry.post==1, c("cancer.hes.date", "cancer.death.date", 
                                                           "cancer.registry.post.date")], 1, function(x) min(x, na.rm = T))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))

temp <- as.data.frame(UKB_a$cancer.hes.date)
temp$key <- as.numeric(rownames(temp))
temp <- merge(temp, earliest, all = T) 

temp$date <- as.Date(temp$date)

# whichever the earliest
UKB_a$cancer.all.date <- temp$date

# baseline cancer diagnosis date
# interview, inpatient, register
# missing values exist
# pick up the earliest date among the three
earliest <- apply(UKB_a[!is.na(UKB_a$cancer.inpatient.pre.date) | 
                          !is.na(UKB_a$cancer.baseline.date) | !is.na(UKB_a$cancer.registry.pre.date),
                        c("cancer.inpatient.pre.date", "cancer.baseline.date", 
                          "cancer.registry.pre.date")], 1, function(x) min(x, na.rm = T))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))

temp <- as.data.frame(UKB_a$cancer.inpatient.post.date)
temp$key <- as.numeric(rownames(temp))
temp <- merge(temp, earliest, all = T) 

# whichever the earliest
UKB_a$cancer.baseline.all.date <- as.Date(temp$date)
# there are a few missing dates 
# because there are 77 missing dates for baseline interview cancer data.
# assume them to be the median of the cancer.baseline.date - "2002-07-01"
UKB_a$cancer.baseline.all.date[is.na(UKB_a$cancer.baseline.all.date) & UKB_a$cancer.baseline.all==1] <- 
  as.Date(as.numeric(summary(UKB_a$cancer.baseline.date)["Median"]), origin= "1970-01-01")

# incident cancer that has no baseline cancer

UKB_a$cancer.incident.only <- UKB_a$cancer.all.post
UKB_a$cancer.incident.only[UKB_a$cancer.baseline.all==1] <- 0
UKB_a$cancer.incident.only.date <- UKB_a$cancer.all.date
UKB_a$cancer.incident.only.date[UKB_a$cancer.baseline.all==1] <- NA


### diabetes ----------------------------------------------------------------

###
# incorporate diabetes medication info from primary care data 
# there are four variables that record prescriptions:
# read 2 code, BNF code, dmd code and drug name
# All have considerable missing values
# The best combination is read 2 + drug name, with missingness for only five ids, for whom other code variables are also missing

# insulin prescription
# read2
# f1... short-acting insulin preparations
# f2... medium/long-acting insulins
# fw... short with intermediate-acting insulins
gp_scr <- read.delim(file.path(raw_data, "gp_scripts.txt"))
read2 <- c("^f1", "^f2", "^fw")
rd2ptn <- paste(read2, collapse = "|")

# drug name
drug_name <- c("insulin","hypurin","neusulin","quicksol","velosulin","actrapid",
               "humulin","novopen","penject","humaject","pur-in","autopen","BD pen",
               "BD ultra pen","insuman","diapen","exubera","humalog","novorapid",
               "apidra","rapitard","penmix","lentard","neulente","tempulin",
               "monotard","semitard","ultratard","insulatard","monophane","neuphane",
               "initard","mixtard","protaphane","actraphane","protaphane","isophane",
               "lantus","toujeo","abasaglar","levemir","tresiba","xultophy", "novomix")

drugptn <- paste(drug_name, collapse = "|")

# filter to leave those with insulin prescription

temp <- gp_scr %>% filter(grepl(rd2ptn, read_2) | grepl(drugptn, drug_name, ignore.case = T))

# transform to date format
temp <- transform(temp, insulin_date=as.Date(as.character(issue_date), "%d/%m/%Y"))

# keep the earliest insulin use records for each id

temp <- temp %>% group_by(eid) %>% filter(insulin_date==min(insulin_date)) %>% select(eid, insulin_date) %>% distinct(eid, .keep_all = T)

# merge to UKB_a
UKB_a <- merge(UKB_a, temp, by.x = "id", by.y = "eid", all.x = T)

# non metformin anti-diabetic drugs
# the same method as above
read2 <- c("^f3", "^f5", "^f6", "^f8", "^ft")
# drug name: import from mannually generated list
drug_name <- read.csv(file.path(work_data, "non-metformin.csv"), header = F)
drug_name <- drug_name$V1

rd2ptn <- paste(read2, collapse = "|")
drugptn <- paste(drug_name, collapse = "|")

temp <- gp_scr %>% filter(grepl(rd2ptn, read_2) | grepl(drugptn, drug_name, ignore.case = T))
temp <- transform(temp, nonmetformin_date=as.Date(as.character(issue_date), "%d/%m/%Y"))
temp <- temp %>% group_by(eid) %>% filter(nonmetformin_date==min(nonmetformin_date)) %>% select(eid, nonmetformin_date) %>% distinct(eid, .keep_all = T)

UKB_a <- merge(UKB_a, temp, by.x = "id", by.y = "eid", all.x = T)

###


# use first occurrence data only
# all inpatient and death diabetes are included in first occurrence
# a very few baseline diabetes are not included 
# check by hand, those missing in first occurrence but available in baseline interview is because the date for the interview data is -1, representing "uncertain or unknown"
# for them, we cannot calculate duration, so will also abandon them

# use first occurrence data 

# five date.field for diabetes
# E10: bd$f.130706.0.0
# E11: bd$f.130708.0.0
# E12: bd$f.130710.0.0 # exclude
# E13: bd$f.130712.0.0
# E14: bd$f.130714.0.0

temp <- data.frame(E10=bd$f.130706.0.0, E11=bd$f.130708.0.0, 
                   E13=bd$f.130712.0.0, E14=bd$f.130714.0.0)

UKB_a$diabetes.fo <- 0
UKB_a$diabetes.fo[rowSums(!is.na(temp)) > 0] <- 1
# given one person can have more than one dates for diabetes, we should choose the first one
# before this, bd$f.130706.0.0 and bd$f.130708.0.0 have 1 and 3 "1902-02-02" respectively, 
# which means code has event date matching date of birth. recode them as NA
temp$E10[temp$E10=="1902-02-02"] <- NA
temp$E11[temp$E11=="1902-02-02"] <- NA
# after recoding, sum(rowSums(!is.na(temp)) > 0) has unchanged number
# it means no diabetes record will lose date
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(UKB_a), but keeps the same row names as names
# join to diabetes.date using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 

UKB_a$diabetes.fo.date <- temp$date
UKB_a$diabetes.fo.date <- as.Date(UKB_a$diabetes.fo.date)

# so far we have diabetes record and linked date as the first diagnosis
UKB_a$diabetes.fo.pre <- ifelse(is.na(UKB_a$diabetes.fo.date), 0,
                                ifelse(UKB_a$diabetes.fo.date<UKB_a$recruit.date, 1, 0))
UKB_a$diabetes.fo.pre.date <- UKB_a$diabetes.fo.date
UKB_a$diabetes.fo.pre.date[UKB_a$diabetes.fo.pre==0] <- NA

UKB_a$diabetes.fo.post <- ifelse(is.na(UKB_a$diabetes.fo.date), 0,
                                 ifelse(UKB_a$diabetes.fo.date>=UKB_a$recruit.date, 1, 0))
UKB_a$diabetes.fo.post.date <- UKB_a$diabetes.fo.date
UKB_a$diabetes.fo.post.date[UKB_a$diabetes.fo.post==0] <- NA

# # duration of baseline diabetes until recruitment
# UKB_a$diabetes.fo.pre.duration <- UKB_a$recruit.date-UKB_a$diabetes.fo.pre.date
# UKB_a$diabetes.fo.pre.duration <- as.numeric(UKB_a$diabetes.fo.pre.duration)

######## add medication-identified diabetes
# RW 2021-01-07

# pre

# keep the earlier from the above two prescription date
UKB_a <- transform(UKB_a, medidate_pre = pmin(insulin_date, nonmetformin_date, na.rm = T))
# T1, T2 could use this later as well.
# add medication-identified diabetes
UKB_a$diabetes.fo.pre[UKB_a$diabetes.fo.pre==0 & 
                        !is.na(UKB_a$medidate_pre) & UKB_a$medidate_pre<UKB_a$recruit.date] <- 1                  

# use the one above as the date of diagnosis  
UKB_a$diabetes.fo.pre.date[is.na(UKB_a$diabetes.fo.pre.date) & UKB_a$diabetes.fo.pre==1] <- 
  UKB_a$medidate_pre[is.na(UKB_a$diabetes.fo.pre.date) & UKB_a$diabetes.fo.pre==1]

# post
UKB_a$diabetes.fo.post[UKB_a$diabetes.fo.post==0 &
                         ((!is.na(UKB_a$insulin_date) & UKB_a$insulin_date>=UKB_a$recruit.date) | 
                            (!is.na(UKB_a$nonmetformin_date) & UKB_a$nonmetformin_date>=UKB_a$recruit.date))] <- 1  

UKB_a$diabetes.fo.post[UKB_a$diabetes.fo.pre==1] <- 0

# keep insulin issue date after recruitment 
UKB_a$medidate1 <- UKB_a$insulin_date
UKB_a$medidate1[!is.na(UKB_a$insulin_date) & UKB_a$insulin_date < UKB_a$recruit.date] <- NA
# keep nonmetformin diabetic drug issue date after recruitment
UKB_a$medidate2 <- UKB_a$nonmetformin_date
UKB_a$medidate2[!is.na(UKB_a$nonmetformin_date) & UKB_a$nonmetformin_date < UKB_a$recruit.date] <- NA
# keep the earlier from the above two
UKB_a <- transform(UKB_a, medidate = pmin(medidate1, medidate2, na.rm = T))
# T1, T2 could use this later as well.

# use the one above as the date of diagnosis  
UKB_a$diabetes.fo.post.date[is.na(UKB_a$diabetes.fo.post.date) & UKB_a$diabetes.fo.post==1] <- 
  UKB_a$medidate[is.na(UKB_a$diabetes.fo.post.date) & UKB_a$diabetes.fo.post==1]

UKB_a$diabetes.fo.post.date[UKB_a$diabetes.fo.post==0] <- NA

# -------------------------------------------------------------------------------------- # 
# split T1 T2 diabetes
# leave E13 and E14 unspecified diabetes at the end  

# T1 first
# E10: bd$f.130706.0.0
UKB_a$diabetes.T1.fo.date <- bd$f.130706.0.0
UKB_a$diabetes.T1.fo.date[UKB_a$diabetes.T1.fo.date=="1902-02-02"] <- NA

UKB_a$diabetes.T1.fo.date <- as.Date(UKB_a$diabetes.T1.fo.date)
UKB_a$diabetes.T1.fo <- ifelse(is.na(UKB_a$diabetes.T1.fo.date), 0 , 1)

# then T2

# E11: bd$f.130708.0.0
UKB_a$diabetes.T2.fo.date <- bd$f.130708.0.0
UKB_a$diabetes.T2.fo.date[UKB_a$diabetes.T2.fo.date=="1902-02-02"] <- NA

UKB_a$diabetes.T2.fo.date <- as.Date(UKB_a$diabetes.T2.fo.date)
UKB_a$diabetes.T2.fo <- ifelse(is.na(UKB_a$diabetes.T2.fo.date), 0 , 1)

# other/unspecified
# E13: bd$f.130712.0.0
# E14: bd$f.130714.0.0

temp <- data.frame(E13=bd$f.130712.0.0, E14=bd$f.130714.0.0)

UKB_a$diabetes.ukn.fo <- 0
UKB_a$diabetes.ukn.fo[rowSums(!is.na(temp)) > 0] <- 1
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(UKB_a), but keeps the same row names as names
# join to diabetes.date using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 

UKB_a$diabetes.ukn.fo.date <- temp$date
UKB_a$diabetes.ukn.fo.date <- as.Date(UKB_a$diabetes.ukn.fo.date)

# for those with both T1/T2 and unknown code
# replace unknown with T1/T2, but use the first diagnosis date
UKB_a$diabetes.T1.fo.date[UKB_a$diabetes.T1.fo==1 & UKB_a$diabetes.ukn.fo==1 & 
                            UKB_a$diabetes.ukn.fo.date < UKB_a$diabetes.T1.fo.date] <-
  UKB_a$diabetes.ukn.fo.date[UKB_a$diabetes.T1.fo==1 & UKB_a$diabetes.ukn.fo==1 & 
                               UKB_a$diabetes.ukn.fo.date < UKB_a$diabetes.T1.fo.date]

UKB_a$diabetes.T2.fo.date[UKB_a$diabetes.T2.fo==1 & UKB_a$diabetes.ukn.fo==1 & 
                            UKB_a$diabetes.ukn.fo.date < UKB_a$diabetes.T2.fo.date] <-
  UKB_a$diabetes.ukn.fo.date[UKB_a$diabetes.T2.fo==1 & UKB_a$diabetes.ukn.fo==1 & 
                               UKB_a$diabetes.ukn.fo.date < UKB_a$diabetes.T2.fo.date]

UKB_a$diabetes.ukn.fo[UKB_a$diabetes.T1.fo==1 | UKB_a$diabetes.T2.fo==1] <- 0
UKB_a$diabetes.ukn.fo.date[UKB_a$diabetes.ukn.fo==0] <- NA


# finally divide T1, T2 to pre and post recruitment
# pre
# T1
UKB_a$diabetes.T1.fo.pre <- ifelse(is.na(UKB_a$diabetes.T1.fo.date), 0,
                                   ifelse(UKB_a$diabetes.T1.fo.date<UKB_a$recruit.date, 1, 0))
UKB_a$diabetes.T1.fo.pre.date <- UKB_a$diabetes.T1.fo.date
UKB_a$diabetes.T1.fo.pre.date[UKB_a$diabetes.T1.fo.pre==0] <- NA
# # duration of baseline diabetes until recruitment
# UKB_a$diabetes.T1.fo.pre.duration <- UKB_a$recruit.date-UKB_a$diabetes.T1.fo.pre.date
# UKB_a$diabetes.T1.fo.pre.duration <- as.numeric(UKB_a$diabetes.T1.fo.pre.duration)

# T2
UKB_a$diabetes.T2.fo.pre <- ifelse(is.na(UKB_a$diabetes.T2.fo.date), 0,
                                   ifelse(UKB_a$diabetes.T2.fo.date<UKB_a$recruit.date, 1, 0))
UKB_a$diabetes.T2.fo.pre.date <- UKB_a$diabetes.T2.fo.date
UKB_a$diabetes.T2.fo.pre.date[UKB_a$diabetes.T2.fo.pre==0] <- NA
# # duration of baseline diabetes until recruitment
# UKB_a$diabetes.T2.fo.pre.duration <- UKB_a$recruit.date-UKB_a$diabetes.T2.fo.pre.date
# UKB_a$diabetes.T2.fo.pre.duration <- as.numeric(UKB_a$diabetes.T2.fo.pre.duration)

# unkown
UKB_a$diabetes.ukn.fo.pre <- ifelse(is.na(UKB_a$diabetes.ukn.fo.date), 0,
                                    ifelse(UKB_a$diabetes.ukn.fo.date<UKB_a$recruit.date, 1, 0))
UKB_a$diabetes.ukn.fo.pre.date <- UKB_a$diabetes.ukn.fo.date
UKB_a$diabetes.ukn.fo.pre.date[UKB_a$diabetes.ukn.fo.pre==0] <- NA

# apply the algorithm for pre-cruit dm 2021-01-07
# add unknown type diabetes have insulin within 12 month post-diagnosis to T1
UKB_a$diabetes.T1.fo.pre[UKB_a$diabetes.ukn.fo.pre==1 & !is.na(UKB_a$insulin_date) & 
                           UKB_a$insulin_date - UKB_a$diabetes.ukn.fo.pre.date <= 365 & 
                           UKB_a$insulin_date - UKB_a$diabetes.ukn.fo.pre.date >= 0] <- 1
UKB_a$diabetes.T1.fo.pre.date[is.na(UKB_a$diabetes.T1.fo.pre.date) & 
                                UKB_a$diabetes.T1.fo.pre==1] <- 
  UKB_a$diabetes.ukn.fo.pre.date[is.na(UKB_a$diabetes.T1.fo.pre.date) & 
                                   UKB_a$diabetes.T1.fo.pre==1]

# others T2
UKB_a$diabetes.T2.fo.pre[UKB_a$diabetes.ukn.fo.pre==1 & UKB_a$diabetes.T1.fo.pre==0] <- 1
UKB_a$diabetes.T2.fo.pre.date[is.na(UKB_a$diabetes.T2.fo.pre.date) & 
                                UKB_a$diabetes.T2.fo.pre==1] <- 
  UKB_a$diabetes.ukn.fo.pre.date[is.na(UKB_a$diabetes.T2.fo.pre.date) & 
                                   UKB_a$diabetes.T2.fo.pre==1]

# add people without any diabetes record (survey  + GP + HES + death)
# but have insulin or nonmetformin prescription after recruitment
# we have already add medication-identified dm for diabetes.fo.pre, so
# T1
# we don't know the diagnosis date for those with medication information only,
# so use insulin + prescription date<=20 yr

# the as.Date function assume the date and month of today for those with only year 
# so we minus the difference between the current date XXXX-XX-XX and the XXXX-01-01
# then add 182 to the mid of the year
dif <- as.numeric(Sys.Date() - as.Date(paste0(substr(as.character(Sys.Date()), 1, 4), "-01-01")))
UKB_a$birth.date.assumed <- as.Date(as.character(UKB_a$birth.year), format = "%Y") - dif +182

UKB_a$diabetes.T1.fo.pre[UKB_a$diabetes.T1.fo.pre==0 & !is.na(UKB_a$insulin_date) & 
                           UKB_a$insulin_date<UKB_a$recruit.date &  
                           UKB_a$insulin_date-UKB_a$birth.date.assumed<=365.25*20] <- 1 

UKB_a$diabetes.T1.fo.pre.date[is.na(UKB_a$diabetes.T1.fo.pre.date) & UKB_a$diabetes.T1.fo.pre==1] <- 
  UKB_a$medidate_pre[is.na(UKB_a$diabetes.T1.fo.pre.date) & UKB_a$diabetes.T1.fo.pre==1]
# above, not necessarily insulin date, but whatever date earilier

# T2
UKB_a$diabetes.T2.fo.pre[UKB_a$diabetes.T2.fo.pre==0 & UKB_a$diabetes.T1.fo.pre==0 & 
                           UKB_a$diabetes.fo.pre==1] <- 1

UKB_a$diabetes.T2.fo.pre.date[is.na(UKB_a$diabetes.T2.fo.pre.date) & UKB_a$diabetes.T2.fo.pre==1] <- 
  UKB_a$medidate_pre[is.na(UKB_a$diabetes.T2.fo.pre.date) & UKB_a$diabetes.T2.fo.pre==1]


# post
UKB_a$diabetes.T1.fo.post <- ifelse(UKB_a$diabetes.fo.pre==1 | is.na(UKB_a$diabetes.T1.fo.date), 0,  
                                    ifelse(UKB_a$diabetes.T1.fo.date>=UKB_a$recruit.date, 1, 0))
UKB_a$diabetes.T1.fo.post.date <- UKB_a$diabetes.T1.fo.date
UKB_a$diabetes.T1.fo.post.date[UKB_a$diabetes.T1.fo.post==0] <- NA

UKB_a$diabetes.T2.fo.post <- ifelse(UKB_a$diabetes.fo.pre==1 | is.na(UKB_a$diabetes.T2.fo.date), 0,
                                    ifelse(UKB_a$diabetes.T2.fo.date>=UKB_a$recruit.date, 1, 0))
UKB_a$diabetes.T2.fo.post.date <- UKB_a$diabetes.T2.fo.date
UKB_a$diabetes.T2.fo.post.date[UKB_a$diabetes.T2.fo.post==0] <- NA

UKB_a$diabetes.ukn.fo.post <- ifelse(UKB_a$diabetes.fo.pre==1 | is.na(UKB_a$diabetes.ukn.fo.date), 0,
                                     ifelse(UKB_a$diabetes.ukn.fo.date>=UKB_a$recruit.date, 1, 0))
UKB_a$diabetes.ukn.fo.post.date <- UKB_a$diabetes.ukn.fo.date
UKB_a$diabetes.ukn.fo.post.date[UKB_a$diabetes.ukn.fo.post==0] <- NA

# apply the algorithm 2020-12-01
# add unknown type diabetes have insulin within 12 month post-diagnosis to T1
UKB_a$diabetes.T1.fo.post[UKB_a$diabetes.ukn.fo.post==1 & !is.na(UKB_a$insulin_date) & 
                            UKB_a$insulin_date - UKB_a$diabetes.ukn.fo.post.date <= 365 & 
                            UKB_a$insulin_date - UKB_a$diabetes.ukn.fo.post.date >= 0] <- 1
UKB_a$diabetes.T1.fo.post.date[is.na(UKB_a$diabetes.T1.fo.post.date) & 
                                 UKB_a$diabetes.T1.fo.post==1] <- 
  UKB_a$diabetes.ukn.fo.post.date[is.na(UKB_a$diabetes.T1.fo.post.date) & 
                                    UKB_a$diabetes.T1.fo.post==1]

# others T2
UKB_a$diabetes.T2.fo.post[UKB_a$diabetes.ukn.fo.post==1 & UKB_a$diabetes.T1.fo.post==0] <- 1
UKB_a$diabetes.T2.fo.post.date[is.na(UKB_a$diabetes.T2.fo.post.date) & 
                                 UKB_a$diabetes.T2.fo.post==1] <- 
  UKB_a$diabetes.ukn.fo.post.date[is.na(UKB_a$diabetes.T2.fo.post.date) & 
                                    UKB_a$diabetes.T2.fo.post==1]

# add people without any diabetes record (survey  + GP + HES + death)
# but have insulin or nonmetformin prescription after recruitment
# T1
UKB_a$diabetes.T1.fo.post[UKB_a$diabetes.fo==0 & !is.na(UKB_a$insulin_date) &   
                            UKB_a$insulin_date - UKB_a$recruit.date <= 365 & 
                            UKB_a$insulin_date - UKB_a$recruit.date >= 0] <- 1

UKB_a$diabetes.T1.fo.post.date[is.na(UKB_a$diabetes.T1.fo.post.date) & UKB_a$diabetes.T1.fo.post==1] <- 
  UKB_a$medidate[is.na(UKB_a$diabetes.T1.fo.post.date) & UKB_a$diabetes.T1.fo.post==1]

# others T2
UKB_a$diabetes.T2.fo.post[UKB_a$diabetes.fo==0 & 
                            ((!is.na(UKB_a$insulin_date) & UKB_a$insulin_date>UKB_a$recruit.date) | 
                               (!is.na(UKB_a$nonmetformin_date) & UKB_a$nonmetformin_date>UKB_a$recruit.date)) &
                            UKB_a$diabetes.T1.fo.post==0] <- 1

UKB_a$diabetes.T2.fo.post.date[is.na(UKB_a$diabetes.T2.fo.post.date) & UKB_a$diabetes.T2.fo.post==1] <- 
  UKB_a$medidate[is.na(UKB_a$diabetes.T2.fo.post.date) & UKB_a$diabetes.T2.fo.post==1]

UKB_a$medidate <- NULL
UKB_a$medidate1 <- NULL
UKB_a$medidate2 <- NULL
UKB_a$medidate_pre <- NULL

# update T1 and T2 post, because medication information may add new pre dm as post dm exist
UKB_a$diabetes.T1.fo.post[UKB_a$diabetes.T1.fo.pre==1] <- 0
UKB_a$diabetes.T1.fo.post.date[UKB_a$diabetes.T1.fo.post==0] <- NA

UKB_a$diabetes.T2.fo.post[UKB_a$diabetes.T2.fo.pre==1] <- 0
UKB_a$diabetes.T2.fo.post.date[UKB_a$diabetes.T2.fo.post==0] <- NA

# update the combination of post and pre
UKB_a$diabetes.T1.fo <- UKB_a$diabetes.T1.fo.pre
UKB_a$diabetes.T1.fo[UKB_a$diabetes.T1.fo.post==1] <- 1

UKB_a$diabetes.T1.fo.date <- UKB_a$diabetes.T1.fo.pre.date
UKB_a$diabetes.T1.fo.date[is.na(UKB_a$diabetes.T1.fo.date) & UKB_a$diabetes.T1.fo==1] <- 
  UKB_a$diabetes.T1.fo.post.date[is.na(UKB_a$diabetes.T1.fo.date) & UKB_a$diabetes.T1.fo==1]

UKB_a$diabetes.T2.fo <- UKB_a$diabetes.T2.fo.pre
UKB_a$diabetes.T2.fo[UKB_a$diabetes.T2.fo.post==1] <- 1

UKB_a$diabetes.T2.fo.date <- UKB_a$diabetes.T2.fo.pre.date
UKB_a$diabetes.T2.fo.date[is.na(UKB_a$diabetes.T2.fo.date) & UKB_a$diabetes.T2.fo==1] <- 
  UKB_a$diabetes.T2.fo.post.date[is.na(UKB_a$diabetes.T2.fo.date) & UKB_a$diabetes.T2.fo==1]



##### check death date earlier than events----

# found one case with stroke.all.date=2014-07-29; death.vascular.date=2014-07-03
# change stroke.all.date=2014-07-03
# rowname=27773
UKB_a[UKB_a$death.vascular.date<UKB_a$stroke.all.date & !is.na(UKB_a$death.vascular.date) & 
        !is.na(UKB_a$stroke.all.date), "stroke.all.date"] <- 
  UKB_a[UKB_a$death.vascular.date<UKB_a$stroke.all.date & !is.na(UKB_a$death.vascular.date) & 
          !is.na(UKB_a$stroke.all.date), "death.vascular.date"]

# diabetes.fo.post.date=2017-01-08
# death.date=2016-12.31
# rowname=371444
UKB_a[UKB_a$death.date<UKB_a$diabetes.fo.post.date & !is.na(UKB_a$death.date) & 
        !is.na(UKB_a$diabetes.fo.post.date), c("diabetes.fo.post.date")] <- 
  UKB_a[UKB_a$death.date<UKB_a$diabetes.fo.post.date & !is.na(UKB_a$death.date) & 
          !is.na(UKB_a$diabetes.fo.post.date), c("death.date")]

# death.date = 2010-08-28
# cancer.all.date = 2010-08-31
# case name = 81073
# recode cancer.all.date = death.date
UKB_a[UKB_a$death.date<UKB_a$cancer.incident.only.date & !is.na(UKB_a$death.date) & 
        !is.na(UKB_a$cancer.incident.only.date), "cancer.incident.only.date"] <- 
  UKB_a[UKB_a$death.date<UKB_a$cancer.incident.only.date & !is.na(UKB_a$death.date) & 
          !is.na(UKB_a$cancer.incident.only.date), "death.date"]



# -------------------------------------------------------------------------------------- # 
# Other baseline disease algorithms #####
# -------------------------------------------------------------------------------------- # 


### hypertension ------------------------------------------------------------

# verbal interview data
# code 1065: hypertension
# code 1072: essential hypertension
text1 <- "hpt_sr"
UKB_a[[text1]] <- 0
for (j in 0:33) { # gather all records from multiple arrays
  text2 <- paste0("f.20002.0.", j)
  UKB_a[[text1]][bd[[text2]]== 1065 |
                   bd[[text2]]== 1072] <- 1
}

# ICD10: 
# first occurrence data
# I10 essential (primary) hypertension: f.131286.0.0
# I11 hypertensive heart disease: f.131288.0.0
# I12 hypertensive renal disease: f.131290.0.0
# I13 hypertensive heart and renal disease: f.131292.0.0
# I15 secondary hypertension: f.131294.0.0
UKB_a$hpt_fo <- 0
for (i in c("f.131286.0.0","f.131288.0.0","f.131290.0.0","f.131292.0.0", "f.131294.0.0")) {
  UKB_a$hpt_fo[bd[[i]]< UKB_a$recruit.date] <- 1
}

UKB_a$hypertension <- ifelse(UKB_a$hpt_fo==1 | UKB_a$hpt_sr==1, 1, 0)

UKB_a$hpt_fo <- NULL
UKB_a$hpt_sr <- NULL

# bp medication from verbal interview
# reference: Alice Carter. et al. Education inequalities in statin treatment  
bpmed <- read.csv(file.path(work_data, "HBPmed.csv"), header = F, encoding = "UTF-8")

# the first element is coded incorrectly
# revise it manually 
bpmed$V1[1] <- 1140860332

bpmed <- as.numeric(bpmed$V1)

UKB_a$BPMed_vi <- 0
for (j in 0:47) { # gather all records from multiple arrays
  text2 <- paste0("f.20003.0.", j)
  UKB_a$BPMed_vi[bd[[text2]] %in% bpmed] <- 1
}

# hypertension treatment
text1 <- "HBPtx"
text2 <- "hypertension"
text3 <- "BPMed_vi" 
UKB_a[[text1]] <- NA
UKB_a[[text1]][UKB_a[[text2]]==1 & UKB_a[[text3]]==1] <- 1
UKB_a[[text1]][UKB_a[[text2]]==0 | UKB_a[[text3]]==0] <- 0  




### peripheral arterial disease ---------------------------------------------
# PAD: verbal interview data
# code 1067: PVD
# code 1087: leg claudication/intermittent claudication
# code 1088: arterial embolism
# code 1492: aortic aneurysm
# code 1591: aortic aneurysm rupture
# code 1592: aortic dissection
# further add codes from operations

text1 <- "PVD.sr"
UKB_a[[text1]] <- 0
for (j in 0:33) { # gather all records from multiple arrays
  text2 <- paste0("f.20002.0.", j)
  UKB_a[[text1]][bd[[text2]]== 1067 |
                   bd[[text2]]== 1087 |
                   bd[[text2]]== 1088 |
                   bd[[text2]]== 1492 |
                   bd[[text2]]== 1591 |
                   bd[[text2]]== 1592] <- 1
}
for (k in 0:31) { # gather all records from multiple arrays
  text3 <- paste0("f.20004.0.", k)
  UKB_a[[text1]][bd[[text3]]== 1071 |
                   bd[[text3]]== 1102 |
                   bd[[text3]]== 1103 |
                   bd[[text3]]== 1555 |
                   bd[[text3]]== 1104 |
                   bd[[text3]]== 1105 |
                   # bd[[text3]]== 1106 |
                   bd[[text3]]== 1107 |
                   bd[[text3]]== 1108 |
                   bd[[text3]]== 1109 |
                   bd[[text3]]== 1110 |
                   bd[[text3]]== 1440 |
                   bd[[text3]]== 1441 |
                   bd[[text3]]== 1442 |
                   bd[[text3]]== 1443] <- 1
}

# first occurrence
# I71, I72, I73, I74, I77
UKB_a$PVD_fo <- 0
for (i in c("f.131382.0.0","f.131384.0.0","f.131386.0.0", "f.131388.0.0","f.131390.0.0")) {
  UKB_a$PVD_fo_date <- bd[[i]]
  UKB_a$PVD_fo_date[UKB_a$PVD_fo_date >= UKB_a$recruit.date] <- NA
  UKB_a$PVD_fo[!is.na(UKB_a$PVD_fo_date)] <- 1
}

# PAD OPCS code
# OPCS 4
UKB_a$PVD.opcs <- 0
for (i in 0:116) { # f.41272 has 117 arrays
  text2 <- paste0("f.41272.0.",i)
  UKB_a$PVD.opcs[grepl("^L16", bd[[text2]]) |
                   grepl("^L18", bd[[text2]]) |
                   grepl("^L19", bd[[text2]]) |
                   grepl("^L2", bd[[text2]]) |
                   grepl("^L30", bd[[text2]]) |
                   grepl("^L31", bd[[text2]]) |
                   grepl("^L37", bd[[text2]]) |
                   grepl("^L38", bd[[text2]]) |
                   grepl("^L39", bd[[text2]]) |
                   grepl("^L4", bd[[text2]]) |
                   grepl("^L5", bd[[text2]]) |
                   grepl("^L60", bd[[text2]]) |
                   grepl("^L62", bd[[text2]]) |
                   grepl("^L63", bd[[text2]]) |
                   grepl("^L65", bd[[text2]]) |
                   grepl("^L66", bd[[text2]]) |
                   grepl("^L67", bd[[text2]]) |
                   grepl("^L68", bd[[text2]]) |
                   grepl("^L70", bd[[text2]]) |
                   grepl("^L71", bd[[text2]]) |
                   grepl("^L74", bd[[text2]]) |
                   grepl("^L75", bd[[text2]]) |
                   grepl("^L76", bd[[text2]]) |
                   grepl("^L89", bd[[text2]])] <- 1
}
# OPCS 3
for (i in 0:15) {
  text2 <- paste0("f.41273.0.",i)
  UKB_a$PVD.opcs[grepl("88", bd[[text2]]) & !grepl("^888", bd[[text2]])] <- 1
}

temp <- as.data.frame(UKB_a$PVD.opcs)
# above is an easy way to create a date frame with the same structure with UKB_a
for (i in 0:116) {
  text1 <- paste0("PVD.opcs.date.", i)
  text2 <- paste0("f.41272.0.",i)
  text3 <- paste0("f.41282.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][(grepl("^L16", bd[[text2]]) |
                   grepl("^L18", bd[[text2]]) |
                   grepl("^L19", bd[[text2]]) |
                   grepl("^L2", bd[[text2]]) |
                   grepl("^L30", bd[[text2]]) |
                   grepl("^L31", bd[[text2]]) |
                   grepl("^L37", bd[[text2]]) |
                   grepl("^L38", bd[[text2]]) |
                   grepl("^L39", bd[[text2]]) |
                   grepl("^L4", bd[[text2]]) |
                   grepl("^L5", bd[[text2]]) |
                   grepl("^L60", bd[[text2]]) |
                   grepl("^L62", bd[[text2]]) |
                   grepl("^L63", bd[[text2]]) |
                   grepl("^L65", bd[[text2]]) |
                   grepl("^L66", bd[[text2]]) |
                   grepl("^L67", bd[[text2]]) |
                   grepl("^L68", bd[[text2]]) |
                   grepl("^L70", bd[[text2]]) |
                   grepl("^L71", bd[[text2]]) |
                   grepl("^L74", bd[[text2]]) |
                   grepl("^L75", bd[[text2]]) |
                   grepl("^L76", bd[[text2]]) |
                   grepl("^L89", bd[[text2]])) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][(grepl("^L16", bd[[text2]]) |
                   grepl("^L18", bd[[text2]]) |
                   grepl("^L19", bd[[text2]]) |
                   grepl("^L2", bd[[text2]]) |
                   grepl("^L30", bd[[text2]]) |
                   grepl("^L31", bd[[text2]]) |
                   grepl("^L37", bd[[text2]]) |
                   grepl("^L38", bd[[text2]]) |
                   grepl("^L39", bd[[text2]]) |
                   grepl("^L4", bd[[text2]]) |
                   grepl("^L5", bd[[text2]]) |
                   grepl("^L60", bd[[text2]]) |
                   grepl("^L62", bd[[text2]]) |
                   grepl("^L63", bd[[text2]]) |
                   grepl("^L65", bd[[text2]]) |
                   grepl("^L66", bd[[text2]]) |
                   grepl("^L67", bd[[text2]]) |
                   grepl("^L68", bd[[text2]]) |
                   grepl("^L70", bd[[text2]]) |
                   grepl("^L71", bd[[text2]]) |
                   grepl("^L74", bd[[text2]]) |
                   grepl("^L75", bd[[text2]]) |
                   grepl("^L76", bd[[text2]]) |
                   grepl("^L89", bd[[text2]])) & 
                  !is.na(bd[[text2]])]
}
# ICD-9
for (i in 0:15) {
  text1 <- paste0("PVD.opcs.date.", i+117)
  text2 <- paste0("f.41273.0.",i)
  text3 <- paste0("f.41283.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("88", bd[[text2]]) & !grepl("^888", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("88", bd[[text2]]) & !grepl("^888", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}

# if there are multiple date use the earliest one
# get rid of date after recruitment date and first column
for (i in 0:132){
  text1 <- paste0("PVD.opcs.date.", i)
  temp[[text1]][temp[[text1]] >= UKB_a$recruit.date] <- NA
}
temp$`UKB_a$PVD.opcs` <- NULL
# code PVD.opcs=0 if all dates are NAs, which means no PVD before recruitment
UKB_a$PVD.opcs[rowSums(is.na(temp)) == ncol(temp)] <- 0

# UKB_a$PVD_2 <- ifelse(UKB_a$PVD_fo==1 | UKB_a$PVD.opcs==1, 1, 0)

UKB_a$PVD <- ifelse(UKB_a$PVD.sr==1 | UKB_a$PVD_fo==1 | UKB_a$PVD.opcs==1, 1, 0)

UKB_a$PVD_fo <- NULL
UKB_a$PVD.sr <- NULL
UKB_a$PVD_fo_date <- NULL
UKB_a$PVD.opcs <- NULL


### CVD history (primary/secondary) -----------------------------------------

UKB_a$CVhist <- UKB_a$MI.baseline + UKB_a$stroke.baseline + UKB_a$PVD + UKB_a$othCHD

UKB_a$CVD <- ifelse(UKB_a$CVhist==0, "None",
                    ifelse(UKB_a$CVhist==1 & UKB_a$MI.baseline==1, "MI only",  
                           ifelse(UKB_a$CVhist==1 & UKB_a$stroke.baseline==1, "Stroke only", 
                                  ifelse(UKB_a$CVhist==1 & UKB_a$PVD==1, "PVD only", 
                                         ifelse(UKB_a$CVhist==1 & UKB_a$othCHD==1, "other CHD only", "Two or more")))))
UKB_a$CVD <- relevel(as.factor(UKB_a$CVD), ref = "None")



# -------------------------------------------------------------------------------------- # 
# Algorithms of other variables #####
# -------------------------------------------------------------------------------------- # 

### smoke status ####
UKB_a$smoke <- bd$f.20116.0.0
UKB_a$smoke[UKB_a$smoke=="Prefer not to answer"] <- NA
UKB_a$smoke <- factor(UKB_a$smoke)

### blood pressure ####
# 93,94 (manual, alternative), 4080,4079 (automatic, preferred), they rarely overlap
# two readings 
# systolic BP
text1 <- "manuSBP"
# two close readings, take average 
text2 <- "f.93.0.0"
text3 <- "f.93.0.1"
UKB_a[[text1]] <- ifelse(!is.na(bd[[text2]]) & !is.na(bd[[text3]]), 
                         (bd[[text2]]+bd[[text3]])/2, 
                         ifelse(is.na(bd[[text2]]), bd[[text3]], bd[[text2]]))
text1 <- "autoSBP"
# two close readings, take average 
text2 <- "f.4080.0.0"
text3 <- "f.4080.0.1"
UKB_a[[text1]] <- ifelse(!is.na(bd[[text2]]) & !is.na(bd[[text3]]), 
                         (bd[[text2]]+bd[[text3]])/2, 
                         ifelse(is.na(bd[[text2]]), bd[[text3]], bd[[text2]]))
# automatic reading is the priority
text1 <- "SBP"
text2 <- "autoSBP"
text3 <- "manuSBP"
UKB_a[[text1]] <- ifelse(!is.na(UKB_a[[text2]]), UKB_a[[text2]], UKB_a[[text3]])
# diastolic BP
text1 <- "manuDBP"
# two close readings, take average 
text2 <- "f.94.0.0"
text3 <- "f.94.0.1"
UKB_a[[text1]] <- ifelse(!is.na(bd[[text2]]) & !is.na(bd[[text3]]), 
                         (bd[[text2]]+bd[[text3]])/2, 
                         ifelse(is.na(bd[[text2]]), bd[[text3]], bd[[text2]]))
text1 <- "autoDBP"
# two close readings, take average 
text2 <- "f.4079.0.0"
text3 <- "f.4079.0.1"
UKB_a[[text1]] <- ifelse(!is.na(bd[[text2]]) & !is.na(bd[[text3]]), 
                         (bd[[text2]]+bd[[text3]])/2, 
                         ifelse(is.na(bd[[text2]]), bd[[text3]], bd[[text2]]))
text1 <- "DBP"
text2 <- "autoDBP"
text3 <- "manuDBP"
UKB_a[[text1]] <- ifelse(!is.na(UKB_a[[text2]]), UKB_a[[text2]], UKB_a[[text3]])
# delete the intermedium variables
UKB_a$autoSBP <- NULL
UKB_a$autoDBP <- NULL
UKB_a$manuDBP <- NULL
UKB_a$manuSBP <- NULL

### ethnicity ####
white <- c("White", "British", "Irish", "Any other white background")
black <- c("Black or Black British", "Caribbean", "African", "Any other Black background")
s.asian <- c("Asian or Asian British", "Indian", "Pakistani", "Bangladeshi")
other <- c("Mixed", "Chinese", "Other ethnic group", "White and Black Caribbean",
           "White and Black African", "White and Asian", "Any other mixed background",
           "Any other Asian background")
norecord <- c("Prefer not to answer", "Do not know")

text1 <- "ethnicity"
text2 <- "f.21000.0.0"
UKB_a[[text1]] <- NA
UKB_a[[text1]][bd[[text2]]%in% white] <- "White"
UKB_a[[text1]][bd[[text2]]%in% black] <- "Black"
UKB_a[[text1]][bd[[text2]]%in% s.asian] <- "South Asian"
UKB_a[[text1]][bd[[text2]]%in% other] <- "Others"
UKB_a[[text1]][bd[[text2]]%in% norecord] <- NA

UKB_a$ethn2 <- ifelse(UKB_a$ethnicity=="Black", "Black", 
                      ifelse(UKB_a$ethnicity=="White", "White", "Others"))
UKB_a$ethn2 <- relevel(as.factor(UKB_a$ethn2), ref = "White")

### BMI ####
text1 <- "BMI"
text2 <- "f.21001.0.0"
UKB_a[[text1]] <- bd[[text2]]
# BMI category
text1 <- "BMI_cat"
text2 <- "BMI"
UKB_a[[text1]] <- NA
UKB_a[[text1]][UKB_a[[text2]]<18.5] <- "<18.5"
UKB_a[[text1]][UKB_a[[text2]]>=18.5 & UKB_a[[text2]]<25] <- "18.5-25"
UKB_a[[text1]][UKB_a[[text2]]>=25 & UKB_a[[text2]]<30] <- "25-30"
UKB_a[[text1]][UKB_a[[text2]]>=30 & UKB_a[[text2]]<35] <- "30-35"
UKB_a[[text1]][UKB_a[[text2]]>=35 & UKB_a[[text2]]<40] <- "35-40"
UKB_a[[text1]][UKB_a[[text2]]>=40] <- "40+"



### IMD: England, Wales, and Scotland separately ####
# generate residential country using availability of IMD
UKB_a$England <- ifelse(is.na(bd$f.26410.0.0), 0, 1)
UKB_a$Wales <- ifelse(is.na(bd$f.26426.0.0), 0, 1)
UKB_a$Scotland <- ifelse(is.na(bd$f.26427.0.0), 0, 1)

# according to the imd_baseline.pdf from UKB, assign the IMD sources
# IMD2004: 2004, 2005, 2006
# IMD2007: 2007, 2008, 2009
# IMD2010: 2010
# WIMD2005: 2006, 2007
# WIMD2008: 2008, 2009, 2010
# SIMD2006: 2006, 2007, 2008
# SIMD2009: 2009, 2010

UKB_a$IMD.source <- NA
UKB_a$IMD.source[UKB_a$England==1 & UKB_a$recruit.date < "2007-01-01"] <- "IMD2004"
UKB_a$IMD.source[UKB_a$England==1 & UKB_a$recruit.date >= "2007-01-01" 
                 & UKB_a$recruit.date< "2010-01-01"] <- "IMD2007"
UKB_a$IMD.source[UKB_a$England==1 & UKB_a$recruit.date >= "2010-01-01"] <- "IMD2010"

UKB_a$IMD.source[UKB_a$Scotland==1 & UKB_a$recruit.date < "2009-01-01"] <- "SIMD2006"
UKB_a$IMD.source[UKB_a$Scotland==1 & UKB_a$recruit.date >= "2009-01-01"] <- "SIMD2009"

UKB_a$IMD.source[UKB_a$Wales==1 & UKB_a$recruit.date < "2008-01-01"] <- "WIMD2005"
UKB_a$IMD.source[UKB_a$Wales==1 & UKB_a$recruit.date >= "2008-01-01"] <- "WIMD2008"

# Quintile cutoffs of IMD of different sources
# IMD2004: 8.35; 13.72; 21.15; 34.20
# IMD2007: 8.32; 13.74; 21.22; 34.42
# IMD2010: 8.49; 13.79; 21.35; 34.17
# SIMD2006: 7.75; 13.56; 21.05; 33.70
# SIMD2009: 7.76; 13.76; 21.02; 33.72
# WIMD2005: 9.96; 14.94; 21.16; 32.70
# WIMD2008: 9.8; 14.8; 21.2; 32.5

# generate IMD quintile
UKB_a$IMD.Q5 <- NA
UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2004" & bd$f.26410.0.0 <= 8.35] <- 1 # wealthiest
UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2004" & bd$f.26410.0.0 > 8.35 & bd$f.26410.0.0 <= 13.72] <- 2
UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2004" & bd$f.26410.0.0 > 13.72 & bd$f.26410.0.0 <= 21.15] <- 3
UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2004" & bd$f.26410.0.0 > 21.15 & bd$f.26410.0.0 <= 34.20] <- 4
UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2004" & bd$f.26410.0.0 > 34.20] <- 5

UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2007" & bd$f.26410.0.0 <= 8.32] <- 1 # wealthiest
UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2007" & bd$f.26410.0.0 > 8.32 & bd$f.26410.0.0 <= 13.74] <- 2
UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2007" & bd$f.26410.0.0 > 13.74 & bd$f.26410.0.0 <= 21.22] <- 3
UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2007" & bd$f.26410.0.0 > 21.22 & bd$f.26410.0.0 <= 34.42] <- 4
UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2007" & bd$f.26410.0.0 > 34.42] <- 5

UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2010" & bd$f.26410.0.0 <= 8.49] <- 1 # wealthiest
UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2010" & bd$f.26410.0.0 > 8.49 & bd$f.26410.0.0 <= 13.79] <- 2
UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2010" & bd$f.26410.0.0 > 13.79 & bd$f.26410.0.0 <= 21.35] <- 3
UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2010" & bd$f.26410.0.0 > 21.35 & bd$f.26410.0.0 <= 34.17] <- 4
UKB_a$IMD.Q5[UKB_a$IMD.source=="IMD2010" & bd$f.26410.0.0 > 34.17] <- 5

UKB_a$IMD.Q5[UKB_a$IMD.source=="SIMD2006" & bd$f.26427.0.0 <= 7.75] <- 1 # wealthiest
UKB_a$IMD.Q5[UKB_a$IMD.source=="SIMD2006" & bd$f.26427.0.0 > 7.75 & bd$f.26427.0.0 <= 13.56] <- 2
UKB_a$IMD.Q5[UKB_a$IMD.source=="SIMD2006" & bd$f.26427.0.0 > 13.56 & bd$f.26427.0.0 <= 21.05] <- 3
UKB_a$IMD.Q5[UKB_a$IMD.source=="SIMD2006" & bd$f.26427.0.0 > 21.05 & bd$f.26427.0.0 <= 33.70] <- 4
UKB_a$IMD.Q5[UKB_a$IMD.source=="SIMD2006" & bd$f.26427.0.0 > 33.70] <- 5

UKB_a$IMD.Q5[UKB_a$IMD.source=="SIMD2009" & bd$f.26427.0.0 <= 7.76] <- 1 # wealthiest
UKB_a$IMD.Q5[UKB_a$IMD.source=="SIMD2009" & bd$f.26427.0.0 > 7.76 & bd$f.26427.0.0 <= 13.76] <- 2
UKB_a$IMD.Q5[UKB_a$IMD.source=="SIMD2009" & bd$f.26427.0.0 > 13.76 & bd$f.26427.0.0 <= 21.02] <- 3
UKB_a$IMD.Q5[UKB_a$IMD.source=="SIMD2009" & bd$f.26427.0.0 > 21.02 & bd$f.26427.0.0 <= 33.72] <- 4
UKB_a$IMD.Q5[UKB_a$IMD.source=="SIMD2009" & bd$f.26427.0.0 > 33.72] <- 5

UKB_a$IMD.Q5[UKB_a$IMD.source=="WIMD2005" & bd$f.26426.0.0 <= 9.96] <- 1 # wealthiest
UKB_a$IMD.Q5[UKB_a$IMD.source=="WIMD2005" & bd$f.26426.0.0 > 9.96 & bd$f.26426.0.0 <= 14.94] <- 2
UKB_a$IMD.Q5[UKB_a$IMD.source=="WIMD2005" & bd$f.26426.0.0 > 14.94 & bd$f.26426.0.0 <= 21.16] <- 3
UKB_a$IMD.Q5[UKB_a$IMD.source=="WIMD2005" & bd$f.26426.0.0 > 21.16 & bd$f.26426.0.0 <= 32.70] <- 4
UKB_a$IMD.Q5[UKB_a$IMD.source=="WIMD2005" & bd$f.26426.0.0 > 32.70] <- 5

UKB_a$IMD.Q5[UKB_a$IMD.source=="WIMD2008" & bd$f.26426.0.0 <= 9.8] <- 1 # wealthiest
UKB_a$IMD.Q5[UKB_a$IMD.source=="WIMD2008" & bd$f.26426.0.0 > 9.8 & bd$f.26426.0.0 <= 14.8] <- 2
UKB_a$IMD.Q5[UKB_a$IMD.source=="WIMD2008" & bd$f.26426.0.0 > 14.8 & bd$f.26426.0.0 <= 21.2] <- 3
UKB_a$IMD.Q5[UKB_a$IMD.source=="WIMD2008" & bd$f.26426.0.0 > 21.2 & bd$f.26426.0.0 <= 32.5] <- 4
UKB_a$IMD.Q5[UKB_a$IMD.source=="WIMD2008" & bd$f.26426.0.0 > 32.5] <- 5


### add gp record tag into UKB_a####
gp_clin <- read.delim(file.path(raw_data, "gp_clinical.txt"))

# keep unique id
gp <- data.frame(id=unique(gp_clin$eid))
gp$gp_record <- 1

UKB_a <- left_join(UKB_a, gp, by = "id")
UKB_a$gp_record[is.na(UKB_a$gp_record)] <- 0

# add loss of follow-up date
UKB_a$lossfudate <- as.Date(bd$f.191.0.0)




### CKD ----
# use ICD10 N18: chronic renal failure
UKB_a$CKD <- ifelse(bd$f.132032.0.0<UKB_a$recruit.date, 1, 0)
UKB_a$CKD[is.na(UKB_a$CKD)] <- 0

# Hba1c
UKB_a$hba1c <- bd$f.30750.0.0

# according to UKB showcase serum_hb1ac.pdf,
# the machine's analytical range is 15-184 mmol/mol
# 5 case >184 removed
# plus the 5 cases defined as no baseline diabetes

UKB_a$hba1c[UKB_a$hba1c>184] <- NA

# copy to UKB_a2
UKB_a2 <- merge(UKB_a2, UKB_a[, c("id", "hba1c")])
saveRDS(UKB_a2, file = file.path(work_data, "UKB_a2.rds"), compress = F)



### COPD ####

# baseline COPD based on f.20002; self-reported non-cancer illness
# UKB self-report codes within field 20002: 1112, 1113 & 1472

UKB_a$COPD.baseline <- 0

text1 <- "COPD.baseline"
UKB_a[[text1]] <- 0
for (j in 0:33) { # gather all records from multiple arrays
  text2 <- paste0("f.20002.0.", j)
  UKB_a[[text1]][bd[[text2]]== 1112 |
                   bd[[text2]]== 1113 |
                   bd[[text2]]== 1472] <- 1
}

subsetCOPD.baseline <- UKB_a[UKB_a$COPD.baseline == 1, ]
nrow(subsetCOPD.baseline) # 8,314


# ICD10-diagnosed COPD
UKB_a$COPD.inpatient <- 0
for (i in 0:212) { # f.41270 has 213 arrays
  text2 <- paste0("f.41270.0.",i)
  UKB_a$COPD.inpatient[grepl("^J43|^J44 ", bd[[text2]])] <- 1
}


# Diagnosis date
# ICD-10
# create a new date frame to collect all dates linked to COPD diagnoses
# one individual may have more than one diagnoses and linked dates
# create a date frame with the same structure with UKB_a:
temp <- as.data.frame(UKB_a$COPD.inpatient)
for (i in 0:212) { # f.41280 has 213 arrays, corresponding to 41270
  text1 <- paste0("COPD.inpatient", i)
  text2 <- paste0("f.41270.0.",i)
  text3 <- paste0("f.41280.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^J43|^J44", bd[[text2]]) &
                  !is.na(bd[[text2]])] <-
    bd[[text3]][grepl("^J43|^J44", bd[[text2]]) &
                  !is.na(bd[[text2]])]
}

# if there are multiple date for COPD ICD, use the earliest one after recruitment date
# get rid of date before recruitment date and first column
for (i in 0:212){
  text1 <- paste0("COPD.inpatient", i)
  temp[[text1]][temp[[text1]] < UKB_a$recruit.date] <- NA
}
temp$`UKB_a$COPD.inpatient` <- NULL
# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# select the earliest date
# the following command only apply to rows with at least one non-missing values
earliest <- apply(temp[rowSums(is.na(temp)) <
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(UKB_a), but keeps the same row
# names as names
# join to temp using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T)
# give the value to UKB_a
UKB_a$COPD.inpatient.date <- temp$date
UKB_a$COPD.inpatient.date <- as.Date(UKB_a$COPD.inpatient.date)

UKB_a$COPD.inpatient.post <- UKB_a$COPD.inpatient
UKB_a$COPD.inpatient.post[is.na(UKB_a$COPD.inpatient.date)] <- 0

# make inpatient record before recruitment is.na(UKB_a$COPD.inpatient.date) means the date precede recruit date
UKB_a$COPD.baseline[UKB_a$COPD.inpatient==1 & is.na(UKB_a$COPD.inpatient.date)] <- 1

subsetCOPD.baseline <- UKB_a[UKB_a$COPD.baseline == 1, ]
nrow(subsetCOPD.baseline) # 8,441





### sleep apnoea ####

# self-reported sleep apnoea:

UKB_a$sleep_apnoea.baseline <- 0

text1 <- "sleep_apnoea.baseline"
UKB_a[[text1]] <- 0
for (j in 0:33) { # gather all records from multiple arrays
  text2 <- paste0("f.20002.0.", j)
  UKB_a[[text1]][bd[[text2]]== 1123] <- 1
}

subset_sleep_apnoea <- UKB_a[UKB_a$sleep_apnoea.baseline == 1, ]
nrow(subset_sleep_apnoea) # 1,615 have self reported sleep apnoea in the entire UKB cohort


# ICD10-diagnosed sleep apnoea using ICD10 code G47.3
UKB_a$sleep_apnoea.inpatient <- 0
for (i in 0:212) { # f.41270 has 213 arrays
  text2 <- paste0("f.41270.0.",i)
  UKB_a$sleep_apnoea.inpatient[grepl("G473 ", bd[[text2]])] <- 1
}

subsetsleep_apnoea.inpatient <- UKB_a[UKB_a$sleep_apnoea.inpatient == 1, ]
nrow(subsetsleep_apnoea.inpatient) # 8,009 have diagnosed sleep apnoea in entire cohort


# ICD10 G47 first occurrence:
# data field 131060 = date G47 first reported (sleep disorders)

UKB_a$sleep_apnoea_fo <- 0
for (i in "f.131060.0.0") {
  UKB_a$hpt_fo[bd[[i]]< UKB_a$recruit.date] <- 1
}

UKB_a$sleep_apnoea_fo[is.na(UKB_a$sleep_apnoea_fo)] <- 0



# make inpatient record before recruitment is.na(UKB_a$sleep_apnoea.inpatient.date) means the date precede recruit date
UKB_a$sleep_apnoea.baseline[UKB_a$sleep_apnoea.inpatient==1 & is.na(UKB_a$sleep_apnoea_fo)] <- 1

subsetsleep_apnoea.baseline <- UKB_a[UKB_a$sleep_apnoea.baseline == 1, ]
nrow(subsetsleep_apnoea.baseline) # 1,615 have sleep apnoea BEFORE recruitment 






### Townsend deprivation score ####

# merge to get townsend scores (because bd = 1 row more than UKB_a) 
UKB_a <- merge(UKB_a, bd[,c("f.eid", "f.189.0.0")], 
               by.x = "eid", by.y = "f.eid", all.x = T) 

# rename
# UKB_a <- UKB_a %>% rename(townsend = f.189.0.0)

UKB_a$townsend[UKB_a$f.189.0.0 <= -2.8952442] <- 1 # least deprived
UKB_a$townsend[UKB_a$f.189.0.0 > -2.8952442 & UKB_a$f.189.0.0 <= -1.5053373] <- 2
UKB_a$townsend[UKB_a$f.189.0.0 > -1.5053373 & UKB_a$f.189.0.0 <= 0.3733443] <- 3
UKB_a$townsend[UKB_a$f.189.0.0 > 0.3733443 & UKB_a$f.189.0.0 <= 2.9642767] <- 4
UKB_a$townsend[UKB_a$f.189.0.0 > 2.9642767] <- 5 # most deprived

summary(UKB_a$townsend) # There are 623 NA's



### Depression ####

# Diagnosis of ICD-10 depressive disorder: (F32 & F33) in field 41270

# ICD10-diagnosed depression
UKB_a <- merge(UKB_a, bd[,c("f.eid", "f.41270.0.0")], 
                          by.x = "eid", by.y = "f.eid", all.x = T)


bd$depression <- 0
for (i in 0:212) { # f.41270 has 213 arrays
  text2 <- paste0("f.41270.0.",i)
  bd$depression[grepl("^F32|^F33 ", bd[[text2]])] <- 1
}

table(bd$depression)  # only 18,769 people in UKB have a depression ICD code


# Use Fields 41202 and 41204 instead:
UKB_a <- merge(UKB_a, bd[,c("f.eid", "f.41202.0.0", "f.41204.0.0")], 
                          by.x = "eid", by.y = "f.eid", all.x = T)

UKB_a$depression <- 0
UKB_a$depression[grepl("^F32|^F33", UKB_a$f.41202.0.0) |
                              grepl("^F32|^F33", UKB_a$f.41204.0.0) |
                              grepl("^F32|^F33", UKB_a$f.41270.0.0)] <- 1

table(UKB_a$depression) # 5,712 have depression using these codes



### Anxiety ####
# F40 and F41   ~  

# ICD-anxiety:
bd$icd_anxiety <- 0
for (i in 0:212) { # f.41270 has 213 arrays
  text2 <- paste0("f.41270.0.",i)
  bd$icd_anxiety[grepl("^F40|^F41", bd[[text2]])] <- 1
}
table(bd$icd_anxiety) # 11,536 have ICD-10 diagnosed anxiety


# How about 41202 AND 41204 ?
UKB_a$anxiety <- 0
UKB_a$anxiety[grepl("^F40|^F41", UKB_a$f.41202.0.0) |
                           grepl("^F40|^F41", UKB_a$f.41204.0.0) |
                           grepl("^F40|^F41", UKB_a$f.41270.0.0)] <- 1
table(UKB_a$anxiety) # 2,809





### Mental health illness ####

# (define a broad mental health category) 
# Data field 20002 (Non-cancer illness code, self reported) -> neurology/eye/psychiatry -> psychological/psychiatric problem
# self-reported illness at baseline or the subsequent 2 repeat assessments (3 instances)

# data field 20002 has 33 arrays

text1 <- "mental"
bd[[text1]] <- 0
for (j in 0:33) {
  text2 <- paste0("f.20002.0.", j) 
  bd[[text1]][bd[[text2]]== 1286 | 
                bd[[text2]]== 1531 | 
                bd[[text2]]== 1287 | 
                bd[[text2]]== 1288 | 
                bd[[text2]]== 1289 | 
                bd[[text2]]== 1290 |
                bd[[text2]]== 1291 |
                bd[[text2]]== 1469 |
                bd[[text2]]== 1470 |
                bd[[text2]]== 1614 |
                bd[[text2]]== 1615 |
                bd[[text2]]== 1616 |
                bd[[text2]]== 1408 |
                bd[[text2]]== 1409 |
                bd[[text2]]== 1410] <- 1  
}

table(bd$mental) # 37,109

# merge mental variable:
UKB_a <- merge(UKB_a, bd[,c("f.eid", "mental")], 
                          by.x = "eid", by.y = "f.eid", all.x = T)

table(UKB_a$mental)

#UKB_a <- subset(UKB_a, select = -c(f.20002.0.0, mental))

text1 <- "depanx"
bd[[text1]] <- 0
for (j in 0:33) {
  text2 <- paste0("f.20002.0.", j) 
  bd[[text1]][bd[[text2]]== 1286 | bd[[text2]]== 1287] <- 1 
}
table(bd$depanx) # depression AND/ OR anxiety = 33,036

text1 <- "dep"
bd[[text1]] <- 0
for (j in 0:33) {
  text2 <- paste0("f.20002.0.", j) 
  bd[[text1]][bd[[text2]]== 1286] <- 1 
}   
table(bd$dep) # depression ONLY = 28,206

text1 <- "anx"
bd[[text1]] <- 0
for (j in 0:33) {
  text2 <- paste0("f.20002.0.", j) 
  bd[[text1]][bd[[text2]]== 1287] <- 1 
} 
table(bd$anx) # anxiety ONLY = 6,723


# merge field 20002 into UKB_a:
UKB_a <- merge(UKB_a, bd[,c("f.eid", "f.20002.0.0")], 
                          by.x = "eid", by.y = "f.eid", all.x = T)

text1 <- "mental"
UKB_a[[text1]] <- 0
for (j in 0:33) {
  text2 <- paste0("f.20002.0.", j) 
  UKB_a[[text1]][UKB_a[[text2]]== 1286 | 
                              UKB_a[[text2]]== 1531 | 
                              UKB_a[[text2]]== 1287 | 
                              UKB_a[[text2]]== 1288 | 
                              UKB_a[[text2]]== 1289 | 
                              UKB_a[[text2]]== 1290 |
                              UKB_a[[text2]]== 1291 |
                              UKB_a[[text2]]== 1469 |
                              UKB_a[[text2]]== 1470 |
                              UKB_a[[text2]]== 1614 |
                              UKB_a[[text2]]== 1615 |
                              UKB_a[[text2]]== 1616 |
                              UKB_a[[text2]]== 1408 |
                              UKB_a[[text2]]== 1409 |
                              UKB_a[[text2]]== 1410] <- 1
}



### Assessment centre ####

UKB_a$centre <- bd$f.54.0.0
UKB_a$centre[is.na(UKB_a$centre)] <- 0


# Create censor date and fu_completion ##### ----------------------------------------

# censor date: ##

# merge the source variable (dsource) from hesin into UKB_a, by id and instance
# maintain the structure of UKB_a

hesin <- read.delim(file.path(raw_data, "hesin.txt"))

# but first rename eid so both datasets are called id
UKB_a <- UKB_a %>% rename(eid = id)

# merge dsource from hesin
UKB_a <- merge(UKB_a, hesin[,c("eid", "dsource")], by=c("eid"), all.x = T)

# generate a variable for censoring date: censordate

# first, use date of losing fu, death date or the earlier one if both exist
UKB_a <- transform(UKB_a, censordate = pmin(lossfudate, death.date, na.rm = T))

# check
head(UKB_a[!is.na(UKB_a$lossfudate) & !is.na(UKB_a$death.date), c("lossfudate", "death.date", "censordate")])

head(UKB_a[!is.na(UKB_a$lossfudate) | !is.na(UKB_a$death.date), c("lossfudate", "death.date", "censordate")])

# impute missing censordates with dsource: 
# censoring dates:
# HES: 2020-02-29
# SMR: 2016-09-30
# PEDW: 2016-03-31

UKB_a$censordate[is.na(UKB_a$censordate) & UKB_a$dsource=="HES"] <- 
  as.Date("2020-02-29")

UKB_a$censordate[is.na(UKB_a$censordate) & UKB_a$dsource=="SMR"] <- 
  as.Date("2016-09-30")

UKB_a$censordate[is.na(UKB_a$censordate) & UKB_a$dsource=="PEDW"] <- 
  as.Date("2016-03-31")

# check NA of dsource:
sum(is.na(UKB_a$dsource)) # 63,039 NA's! - Use England, Scotland, Wales variables for these instead

UKB_a$censordate[is.na(UKB_a$censordate) & is.na(UKB_a$dsource) & (UKB_a$England==1)] <- 
  as.Date("2020-02-29")

UKB_a$censordate[is.na(UKB_a$censordate) & is.na(UKB_a$dsource) & (UKB_a$Scotland==1)] <- 
  as.Date("2016-09-30")

UKB_a$censordate[is.na(UKB_a$censordate) & is.na(UKB_a$dsource) & (UKB_a$Wales==1)] <- 
  as.Date("2016-03-31")

sum(is.na(UKB_a$censordate)) # 1,665 NA's left for censordate. 

# assumption: those with NA censordate are from England
UKB_a$censordate[is.na(UKB_a$censordate)] <- 
  as.Date("2020-02-29")

# now we have a censoring date
# generate duration of follow up:
UKB_a$fu_duration <- as.numeric(UKB_a$censordate - UKB_a$recruit.date)

summary(UKB_a$fu_duration/365.25) # mean follow-up is 10.1 years

UKB_a <- UKB_a %>% 
  mutate(fu_year = 
           case_when(is.na(fu_duration) ~ "nofollowup",
                     fu_duration %/% (365.25) == -3 ~ "year_3_ago",
                     fu_duration %/% (365.25) == -2 ~ "year_2_ago",
                     fu_duration %/% (365.25) == -1 ~ "year_1_ago",
                     fu_duration %/% (365.25) == 0 ~ "year_1",
                     fu_duration %/% (365.25) == 1 ~ "year_2",
                     fu_duration %/% (365.25) == 2 ~ "year_3",
                     fu_duration %/% (365.25) == 3 ~ "year_4",
                     fu_duration %/% (365.25) == 4 ~ "year_5",
                     fu_duration %/% (365.25) == 5 ~ "year_6",
                     fu_duration %/% (365.25) == 6 ~ "year_7",
                     fu_duration %/% (365.25) == 7 ~ "year_8",
                     fu_duration %/% (365.25) == 8 ~ "year_9",
                     fu_duration %/% (365.25) == 9 ~ "year_10",
                     fu_duration %/% (365.25) == 10 ~ "year_11",
                     fu_duration %/% (365.25) == 11 ~ "year_12",
                     fu_duration %/% (365.25) == 12 ~ "year_13",
                     fu_duration %/% (365.25) == 13 ~ "year_14",
                     fu_duration %/% (365.25) == 14 ~ "year_15"
           ))

head(UKB_a)

# check the counts:
table(UKB_a$fu_year, useNA = "ifany")

# change the data structure:
# unique ids - i.e., 1 row per person
# make consentdate, censordate and fu_duration unique to each id
# keep the longest duration only

UKB_a <- UKB_a %>% group_by(eid) %>% filter(fu_duration==max(fu_duration)) %>% distinct(eid, .keep_all = T)
UKB_a_followup <- subset(UKB_a, 
                         select=c("eid", "recruit.date", "censordate", "fu_duration", "fu_year"))
head(UKB_a_followup)
## now we have a dataset with unique id and follow-up duration

# summarise longest year of followup for each id:
table(UKB_a$fu_year, useNA = "ifany") 

# because the longest fu year is 15.05
# generate 16 year columns:

for (i in 1:16) {
  
  # first create one column with all 0
  UKB_a[[paste0("year_", i)]] <- 0
  
  # change the year with full follow up to be 1
  UKB_a[[paste0("year_", i)]][ UKB_a$fu_duration / 365.25 >= i] <- 1
  
  # change the year with partial follow up to a proportion of the year
  UKB_a[[paste0("year_", i)]][ UKB_a$fu_duration / 365.25 < i &  UKB_a$fu_duration / 365.25 > i-1] <- 
    ( UKB_a$fu_duration[ UKB_a$fu_duration / 365.25 < i &  UKB_a$fu_duration / 365.25 > i-1] %% 365.25)/365.25 # %%: modulus ex. 366%%365.25=0.75
}


## Long version of UKB_a dataset (+ fu_completion) ####
# transform wide to long
UKB_a_2ndlong <- UKB_a %>% gather(paste0("year_", 1:16), key=fu_year, value = fu_completion) %>% arrange(eid)

head(UKB_a_2ndlong, 17)

# now each individual has 16 rows from year 1 to 16
head(UKB_a_2ndlong[,c("eid", "fu_year", "fu_completion")],25)



# -------------------------------------------------------------------------------------- # 
# Asthma diagnosis ####
# -------------------------------------------------------------------------------------- # 

# 42014 - binary variable - code asthma as 1 and no asthma as 0
# Date of asthma should be earlier than recruitment date

UKB_a$asthma <- 0

UKB_a$asthma <- ifelse(!is.na(bd$f.42014) & 
                         bd$f.42014<UKB_a$recruit.date, 1, 0)

table(UKB_a$asthma) # 59,998 have asthma at baseline, 442,507 = no asthma.

# identifying asthma medication codes from nurse interview:
biobank_med_codes <- readRDS("medcode.rds")

# in biobank showcase: Treatment/ medication code = 20003
asthma_all_med <- c(1140855302, 1140855304, 1140855308, 1140855372, 1140855374, 1140855376, 1140855380, 1140855384, 1140855390, 1140855466, 1140855500, 1140855508, 1140855520, 1140855524, 1140855528, 1140855534, 1140855536, 1140855538, 1140861996, 1140861998, 1140862008, 1140862016, 1140862060, 1140862066, 1140862070, 1140862086, 1140862110, 1140862118, 1140862120, 1140862124, 1140862134, 1140862140, 1140862380, 1140862382, 1140862406, 1140862476, 1140862526, 1140862532, 1140862560, 1140862572, 1140862574, 1140862584, 1140862600, 1140862610, 1140881856, 1141157418, 1141157486, 1141195280, 1140883548, 1141182628, 1140855320, 1140855322, 1140855328, 1140855330, 1140855332, 1140855358, 1140855360, 1140855366, 1140855400, 1140855424, 1140855426, 1140855442, 1140855496, 1140855504, 1140855506, 1140855530, 1140855540, 1140855542, 1140862092, 1140862144, 1140862148, 1140862162, 1140862168, 1140862222, 1140862224, 1140862236, 1140862260, 1140862266, 1140862274, 1140862280, 1140862290, 1140862292, 1140862294, 1140862306, 1140862310, 1140862320, 1140862336, 1140862346, 1140862348, 1140862362, 1140862364, 1140862374, 1140862412, 1140862414, 1140862418, 1140862424, 1140862432, 1140862438, 1140862474, 1140868364, 1140868370)
UKB_a$asthma_all <- 0
for (j in 0:47) { # gather all records from multiple arrays
  text2 <- paste0("f.20003.0.", j)
  UKB_a$asthma_all[bd[[text2]] %in% asthma_all_med] <- 1
}
# population with asthma code (42014) AND any medication code:
UKB_a$diagnosis_all_med <- 0
UKB_a$diagnosis_all_med[UKB_a$asthma==1 & UKB_a$asthma_all==1] <-1
table(UKB_a$diagnosis_all_med) #' THIS IS THE FINAL ASTHMA POPULATION. DOCTOR DIAGNOSED PLUS MEDICATION = 25,194

# (check no. of rows etc are exactly the same in both datasets:)
identical(bd$f.eid, UKB_a$id)
# = TRUE

# population with asthma code (42014) AND mild medication code
asthma_mild_med <- c(1140855302, 1140855304, 1140855308, 1140855372, 1140855374, 1140855376, 1140855380, 1140855384, 1140855390, 1140855466, 1140855500, 1140855508, 1140855520, 1140855524, 1140855528, 1140855534, 1140855536, 1140855538, 1140861996, 1140861998, 1140862008, 1140862016, 1140862060, 1140862066, 1140862070, 1140862086, 1140862110, 1140862118, 1140862120, 1140862124, 1140862134, 1140862140, 1140862380, 1140862382, 1140862406, 1140862476, 1140862526, 1140862532, 1140862560, 1140862572, 1140862574, 1140862584, 1140862600, 1140862610, 1140881856, 1141157418, 1141157486, 1141195280)
UKB_a$asthma_mild <- 0
for (j in 0:47) { # gather all records from multiple arrays
  text2 <- paste0("f.20003.0.", j)
  UKB_a$asthma_mild[bd[[text2]] %in% asthma_mild_med] <- 1
}
UKB_a$diagnosis_mild_med <- 0
UKB_a$diagnosis_mild_med[UKB_a$asthma==1 & UKB_a$asthma_mild==1] <-1
table(UKB_a$diagnosis_mild_med) 
# ^ asthma_mild = anyone with a mild medication. 
# diagnosis_mild_med = anyone with diagnosed asthma AND a mild medication


# population with asthma code (42014) AND moderate-severe medication code
asthma_moderate_severe <- c(1140883548, 1141182628, 1140855320, 1140855322, 1140855328, 1140855330, 1140855332, 1140855358, 1140855360, 1140855366, 1140855400, 1140855424, 1140855426, 1140855442, 1140855496, 1140855504, 1140855506, 1140855530, 1140855540, 1140855542, 1140862092, 1140862144, 1140862148, 1140862162, 1140862168, 1140862222, 1140862224, 1140862236, 1140862260, 1140862266, 1140862274, 1140862280, 1140862290, 1140862292, 1140862294, 1140862306, 1140862310, 1140862320, 1140862336, 1140862346, 1140862348, 1140862362, 1140862364, 1140862374, 1140862412, 1140862414, 1140862418, 1140862424, 1140862432, 1140862438, 1140862474, 1140868364, 1140868370)
UKB_a$asthma_mod_severe <- 0
for (j in 0:47) { # gather all records from multiple arrays
  text2 <- paste0("f.20003.0.", j)
  UKB_a$asthma_mod_severe[bd[[text2]] %in% asthma_moderate_severe] <- 1
}
UKB_a$diagnosis_severe_med <- 0
UKB_a$diagnosis_severe_med[UKB_a$asthma==1 & UKB_a$asthma_mod_severe==1] <-1
table(UKB_a$diagnosis_severe_med)
# ^ asthma_mod_severe = anyone with a moderate-severe medication. 
# diagnosis_severe_med = anyone with diagnosed asthma AND a mod-severe medication

## compare 
table(UKB_a$asthma_all) # 28,707 have an asthma medication code, 473,798 don't
table(UKB_a$diagnosis_all_med) # 25,194 have an asthma med code AND an asthma diagnosis code (exclude 3,513 people for not having doctor-diagnosed asthma)
## compare:    {medication categories in those with doctor-diagnosed asthma}
table(UKB_a$diagnosis_mild_med) #' 23,168 of asthma population are taking mild medication
table(UKB_a$diagnosis_all_med, UKB_a$diagnosis_mild_med, useNA = "ifany") #' However, 20,078 have mild medication ONLY (i.e., not alongside severe medication)*
table(UKB_a$diagnosis_severe_med) #' 5,116 of asthma population are taking severe medication
table(UKB_a$diagnosis_mild_med, UKB_a$diagnosis_severe_med, useNA = "ifany") #' 3,090 have BOTH a mild and mod-severe medication, while 2,026 have mod-severe medication ONLY

# create variable to identify mild population (doctor-diagnosed) with mild medication ONLY (i.e., no mild & severe medication overlap)
UKB_a$asthma_mild_ONLY <- 0
UKB_a$asthma_mild_ONLY[UKB_a$diagnosis_severe_med==0 & UKB_a$diagnosis_mild_med==1] <- 1
table(UKB_a$asthma_mild_ONLY) # no. = 20,078

### THEREFORE, the final asthma severity groups are:
#  table(UKB_a$diagnosis_severe_med)    = moderate-severe population (n = 5,116)
#  UKB_a$asthma_mild_ONLY        = mild population (n = 20,078)

# subset the population who have UKB asthma code but no medication code
asthma_no_med <- UKB_a[with(UKB_a, asthma == 1 & diagnosis_all_med == 0), ] 
nrow(asthma_no_med) # no. = 34,804 


# save ----
# save UKB_a dataframe
saveRDS(UKB_a, file.path(work_data, "UKB_a.rds"))
