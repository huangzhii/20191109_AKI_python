
###### mimic

setwd("C:/Users/congf/Desktop/AKI/AKI11.11/MIMIC11.11/non_aki")
mimic_nonaki_cohort = read.csv("non_aki_statistic/mimic_mimiciii_aki_kidgo_negative_stastic5.CSV", header = F, stringsAsFactors = F)


setwd("C:/Users/congf/Desktop/AKI/AKI11.11/MIMIC11.11/aki")
mimic_aki_cohort = read.csv("aki_statistic/mimic_mimiciii_aki_kidgo_positive_stastic5.CSV", header = F, stringsAsFactors = F)

label = levels(unlist(read.table("aki_statistic/mimic_label.tsv.txt", sep = "\t")))
colnames(mimic_aki_cohort) = label
colnames(mimic_nonaki_cohort) = label
mimic_aki_cohort$AKI = 1
mimic_nonaki_cohort$AKI = 0
mimic_aki_cohort$GENDER[is.na(mimic_aki_cohort$GENDER)] = 1
mimic_nonaki_cohort$GENDER[is.na(mimic_nonaki_cohort$GENDER)] = 1
mimic_aki_cohort$ETHNICITY[is.na(mimic_aki_cohort$ETHNICITY)] = 1
mimic_nonaki_cohort$ETHNICITY[is.na(mimic_nonaki_cohort$ETHNICITY)] = 1


setwd("C:/Users/congf/Desktop/AKI/AKI11.14")
alldata.imputed.shrinked = read.csv("MIMIC_final_ICUSTAY_ID_model_6_6.csv", header = F, stringsAsFactors = F) 
label2 = levels(unlist(read.table("final_LABLE.TSV.txt", sep = "\t")))
colnames(alldata.imputed.shrinked) = label2

finalid = unique(alldata.imputed.shrinked$ICUSTAY_ID)




mimic_cohort = rbind(mimic_aki_cohort, mimic_nonaki_cohort)
finaldata= mimic_cohort[mimic_cohort$ICUSTAY_ID %in% finalid,]
finaldata[finaldata$AGE == 300,] = 90
head(finaldata)

####???????????????????????????

install.packages("gmodels")
library(gmodels)

CrossTable(finaldata$HOS_DEATH,finaldata$AKI)
CrossTable(finaldata$AKI_RRT,finaldata$AKI)
#####?????????????????????????????????
install.packages("psych")
library(psych)

feature=c("LOS","ISOFA")
describeBy(finaldata[feature], list(am=finaldata$AKI))
describe(finaldata)
####???????????????????????????
colnames(finaldata)

finaldata2=finaldata[,c(6,10,24,27)]

corr.test(finaldata2,use="complete")



####???????????????????????????
cor.test(finaldata[,29],finaldata[,30])


mytable = with(finaldata,table(ETHNICITY))
mytable
prop.table(mytable)

mytable = with(finaldata,table(GENDER))
mytable
prop.table(mytable)

mytable = with(finaldata,table(SEP))
mytable
prop.table(mytable)

mytable = with(finaldata,table(CAR))
mytable
prop.table(mytable)

mytable = with(finaldata,table(RES))
mytable
prop.table(mytable)

mytable = with(finaldata,table(NEU))
mytable
prop.table(mytable)

mytable = with(finaldata,table(OD))
mytable
prop.table(mytable)

mytable = with(finaldata,table(CKD))
mytable
prop.table(mytable)

mytable = with(finaldata,table(DIA))
mytable
prop.table(mytable)

mytable = with(finaldata,table(CHF))
mytable
prop.table(mytable)

mytable = with(finaldata,table(CLD))
mytable
prop.table(mytable)

mytable = with(finaldata,table(CPD))
mytable
prop.table(mytable)

mytable = with(finaldata,table(HYP))
mytable
prop.table(mytable)

mytable = with(finaldata,table(AKI))
mytable
prop.table(mytable)


mytable = with(finaldata,table(MECH))
mytable
prop.table(mytable)

mytable = with(finaldata,table(VASO))
mytable
prop.table(mytable)

mytable = with(finaldata,table(HOS_DEATH))
mytable
prop.table(mytable)

mytable = with(finaldata,table(AKI_RRT))
mytable
prop.table(mytable)

mytable = with(finaldata,table(AKI_UO))
mytable
prop.table(mytable)

mytable = with(finaldata,table(CST))
mytable
prop.table(mytable)


install.packages("pastecs")
library(pastecs)
mimic_desc<- c("AGE","LOS","ISOFA")
stat.desc(finaldata[mimic_desc])

install.packages("Hmisc")
library(Hmisc)
describe(finaldata[mimic_desc])


head(finaldata)


###### eicu

setwd("C:/Users/congf/Desktop/AKI/AKI11.14/EICU-ICM-11.14/non_aki")
mimic_nonaki_cohort = read.csv("non_aki_statistic2/eicu_eicu_crd_aki_negative_statistic6.CSV", header = F, stringsAsFactors = F)


setwd("C:/Users/congf/Desktop/AKI/AKI11.14/EICU-ICM-11.14/aki")
mimic_aki_cohort = read.csv("aki_statistic/eicu_eicu_crd_aki_positive_statistic6.CSV", header = F, stringsAsFactors = F)

label = levels(unlist(read.table("aki_statistic/eicu_label.tsv.txt", sep = "\t")))
colnames(mimic_aki_cohort) = label
colnames(mimic_nonaki_cohort) = label
mimic_aki_cohort$AKI = 1
mimic_nonaki_cohort$AKI = 0
mimic_aki_cohort$GENDER[is.na(mimic_aki_cohort$GENDER)] = 1
mimic_nonaki_cohort$GENDER[is.na(mimic_nonaki_cohort$GENDER)] = 1
mimic_aki_cohort$ETHNICITY[is.na(mimic_aki_cohort$ETHNICITY)] = 1
mimic_nonaki_cohort$ETHNICITY[is.na(mimic_nonaki_cohort$ETHNICITY)] = 1



setwd("C:/Users/congf/Desktop/AKI/AKI11.14")
alldata.imputed.shrinked = read.csv("EICU_final_ICUSTAY_ID_model_6_6.csv", header = F, stringsAsFactors = F) 
label2 = levels(unlist(read.table("final_LABLE.TSV.txt", sep = "\t")))
colnames(alldata.imputed.shrinked) = label2

finalid = unique(alldata.imputed.shrinked$ICUSTAY_ID)



mimic_cohort = rbind(mimic_aki_cohort, mimic_nonaki_cohort)
finaldata= mimic_cohort[mimic_cohort$ICUSTAY_ID %in% finalid,]

colnames(finaldata)

####???????????????????????????

install.packages("gmodels")
library(gmodels)

CrossTable(finaldata$HOSP_DEATH,finaldata$AKI)
CrossTable(finaldata$AKI_RRT,finaldata$AKI)
#####?????????????????????????????????
install.packages("psych")
library(psych)

feature=c("LOS","ISOFA")
describeBy(finaldata[feature], list(am=finaldata$AKI))

####???????????????????????????
finaldata2=finaldata[,c(8,13,25,28)]

corr.test(finaldata2,use="complete")

cor.test(finaldata[,30],finaldata[,31])



mytable = with(finaldata,table(ETHNICITY))
mytable
prop.table(mytable)

mytable = with(finaldata,table(GENDER))
mytable
prop.table(mytable)

mytable = with(finaldata,table(SEP))
mytable
prop.table(mytable)

mytable = with(finaldata,table(CAR))
mytable
prop.table(mytable)

mytable = with(finaldata,table(RES))
mytable
prop.table(mytable)

mytable = with(finaldata,table(NEU))
mytable
prop.table(mytable)

mytable = with(finaldata,table(OD))
mytable
prop.table(mytable)

mytable = with(finaldata,table(CKD))
mytable
prop.table(mytable)

mytable = with(finaldata,table(DIA))
mytable
prop.table(mytable)

mytable = with(finaldata,table(CHF))
mytable
prop.table(mytable)

mytable = with(finaldata,table(CLD))
mytable
prop.table(mytable)

mytable = with(finaldata,table(CPD))
mytable
prop.table(mytable)

mytable = with(finaldata,table(HYP))
mytable
prop.table(mytable)

mytable = with(finaldata,table(AKI))
mytable
prop.table(mytable)


mytable = with(finaldata,table(MECH))
mytable
prop.table(mytable)

mytable = with(finaldata,table(VASO))
mytable
prop.table(mytable)

mytable = with(finaldata,table(HOSP_DEATH))
mytable
prop.table(mytable)


mytable = with(finaldata,table(AKI_RRT))
mytable
prop.table(mytable)

mytable = with(finaldata,table(AKI_UO))
mytable
prop.table(mytable)

mytable = with(finaldata,table(CST))
mytable
prop.table(mytable)


install.packages("pastecs")
library(pastecs)
mimic_desc<- c("AGE","LOS","ISOFA")
stat.desc(finaldata[mimic_desc])

install.packages("Hmisc")
library(Hmisc)
describe(finaldata[mimic_desc])
