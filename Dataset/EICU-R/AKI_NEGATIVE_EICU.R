setwd("C:/Users/congf/Desktop/AKI/aki_kidgo_EICU_2019.3.13/non_aki")
install.packages(c("anytime", "tictoc","lubridate","dplyr"))
library("anytime")
library("tictoc")
library("lubridate")
help(anytime)options(stringsAsFactors = F)

####################################################################
####################################################################
############################ non_aki ###############################
####################################################################
####################################################################

###########################################################
###            Expand data
###########################################################
alldata = read.csv("non_aki_cohort/eicu_eicu_crd_aki_negative_demo.csv", header = F, stringsAsFactors = F)
label = unlist(read.table("non_aki_cohort/eicu_label.tsv.txt", sep = "\t"))
colnames(alldata) = label
head(alldata)
alldata$AGE = as.numeric(gsub("[^0-9]", "", alldata$AGE))
alldata = alldata[alldata$ADMI_TIME < alldata$OUT_TIME, ] # this is to remove rows of which admi time greater than out time.

# remove patients that age < 15
#alldata = alldata[alldata$AGE >= 15,]
#alldata = alldata[complete.cases(alldata),]
#alldata$HEIGHT[alldata$HEIGHT<3] = alldata$HEIGHT[alldata$HEIGHT<3] * 100 # some heights looks like calculated in meter.
#alldata$BMI = alldata$WEIGHT / (alldata$HEIGHT/100)^2
#alldata = alldata[complete.cases(alldata),]
#min(alldata$OUT_TIME) # ADMI_TIME is useless
# ADMI_TIME and OUT_TIME is in hour.
#sum(alldata$OUT_TIME/24 != alldata$LOS) # make sure length of stay (LOS) is exactly the ICU stay time


# start expanding data
expand.data = NULL
start_time <- proc.time()

allhours = 0
hourslist = NULL
for (i in 1:dim(alldata)[1]){
  if (i %% 2500 == 0){
    message(i, "  -  ", allhours, "    Elapsed: ", round((proc.time() - start_time)[3], 2), " secs")
}
row = alldata[i,]
time.in = row$ADMI_TIME
time.out = row$OUT_TIME
hours = (time.out - time.in) + 1
hourslist = c(hourslist, hours)
allhours = allhours + hours
}
alldata$hours = hourslist
# expand
alldata.expanded <- alldata[rep(row.names(alldata), alldata$hours), 1:dim(alldata)[2]]
head(alldata.expanded)

start_time <- proc.time()
currtime_list = NULL
allhours = 0

for (i in 1:dim(alldata)[1]){
  if (i %% 2500 == 0){
    message(i, "  -  ", allhours, "    Elapsed: ", round((proc.time() - start_time)[3], 2), " secs")
  }
  row = alldata[i,]
  time.in = row$ADMI_TIME
  time.out = row$OUT_TIME
  hours = time.out - time.in
  allhours = allhours + hours
  
  numerical.times = seq(time.in, time.in + hours, by = 1)
  currtime_list[[i]] = numerical.times
}
alldata.expanded$CURR_TIME = unlist(currtime_list)
head(alldata.expanded)



save(alldata.expanded, file = "expanded.all.eicu.nonaki.data.Rdata")
load("expanded.all.eicu.nonaki.data.Rdata")
###########################################################
###             merge vit
###########################################################
library(dplyr)
library(anytime)

hy_vit = read.csv("non_aki_vit/eicu_eicu_crd_aki_negative_vit2.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("non_aki_vit/EICU_vit_label.tsv.txt", sep = "\t"))))
colnames(hy_vit) = label
hy_vit$mergekey = paste0(hy_vit$ICUSTAY_ID, sprintf("_%06d", hy_vit$VIT_TIME))
head(hy_vit)
#take average of every hours (since hours are converted from minutes, it exist duplications)
hy_vit = aggregate(hy_vit[, 1:9], list(hy_vit$mergekey), mean, na.rm=TRUE)
colnames(hy_vit)[colnames(hy_vit) == "Group.1"] = "mergekey"
head(hy_vit)

alldata.expanded$mergekey = paste0(alldata.expanded$ICUSTAY_ID, sprintf("_%06d", alldata.expanded$CURR_TIME))
alldata.merged = full_join(alldata.expanded, hy_vit, by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
head(alldata.merged)
###########################################################
###             merge lab
###########################################################
hy_lab3 = read.csv("non_aki_lab/eicu_eicu_crd_aki_negative_lab1.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("non_aki_lab/EICU_LAB_LABEL.TSV.txt", header = F, sep = "\t"))))
colnames(hy_lab3) = label
##hy_lab3$LAB_TIME = round(hy_lab3$LAB_TIME/60)
##hy_lab3 = hy_lab3[hy_lab3$LAB_TIME >= 0,]
hy_lab3$mergekey = paste0(hy_lab3$ICUSTAY_ID, sprintf("_%06d", hy_lab3$LAB_TIME))
head(hy_lab3)
#take average of every hours (since hours are converted from minutes, it exist duplications)
hy_lab3 = aggregate(hy_lab3[, 1:35], list(hy_lab3$mergekey), mean, na.rm=TRUE)
colnames(hy_lab3)[colnames(hy_lab3) == "Group.1"] = "mergekey"
head(hy_lab3)

alldata.merged = full_join(alldata.merged, hy_lab3, by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
head(alldata.merged)
###########################################################
###             merge INPUT_6hr
###########################################################

hy_input_6hr = read.csv("non_aki_input/6hr/eicu_eicu_crd_aki_negative_input_6hr.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("non_aki_input/6hr/EICU_inputPUT6_label.tsv.txt", sep = "\t"))))
colnames(hy_input_6hr) = label

hy_input_6hr$mergekey = paste0(hy_input_6hr$ICUSTAY_ID, sprintf("_%06d", hy_input_6hr$INPUT_6HR_TIME))
head(hy_input_6hr)
#take average of every hours (since hours are converted from minutes, it exist duplications)
hy_input_6hr = aggregate(hy_input_6hr[, 1:3], list(hy_input_6hr$mergekey), mean, na.rm=TRUE)
colnames(hy_input_6hr)[colnames(hy_input_6hr) == "Group.1"] = "mergekey"
head(hy_input_6hr)

alldata.merged = full_join(alldata.merged, hy_input_6hr, by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
head(alldata.merged)
###########################################################
###             merge INPUT_12hr
###########################################################

hy_input_12hr = read.csv("non_aki_input/12hr/eicu_eicu_crd_aki_negative_input_12hr.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("non_aki_input/12hr/EICU_inputPUT12_label.tsv.txt", sep = "\t"))))
colnames(hy_input_12hr) = label

hy_input_12hr$mergekey = paste0(hy_input_12hr$ICUSTAY_ID, sprintf("_%06d", hy_input_12hr$INPUT_12HR_TIME))
head(hy_input_12hr)
#take average of every hours (since hours are converted from minutes, it exist duplications)
hy_input_12hr = aggregate(hy_input_12hr[, 1:3], list(hy_input_12hr$mergekey), mean, na.rm=TRUE)
colnames(hy_input_12hr)[colnames(hy_input_12hr) == "Group.1"] = "mergekey"
head(hy_input_12hr)

alldata.merged = full_join(alldata.merged, hy_input_12hr, by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
head(alldata.merged)
###########################################################
###             merge INPUT_24hr
###########################################################

hy_input_24hr = read.csv("non_aki_input/24hr/eicu_eicu_crd_aki_negative_input_24hr.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("non_aki_input/24hr/EICU_inputPUT24_label.tsv.txt", sep = "\t"))))
colnames(hy_input_24hr) = label

hy_input_24hr$mergekey = paste0(hy_input_24hr$ICUSTAY_ID, sprintf("_%06d", hy_input_24hr$INPUT_24HR_TIME))
head(hy_input_24hr)
#take average of every hours (since hours are converted from minutes, it exist duplications)
hy_input_24hr = aggregate(hy_input_24hr[, 1:3], list(hy_input_24hr$mergekey), mean, na.rm=TRUE)
colnames(hy_input_24hr)[colnames(hy_input_24hr) == "Group.1"] = "mergekey"
head(hy_input_24hr)

alldata.merged = full_join(alldata.merged, hy_input_24hr, by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
head(alldata.merged)
###########################################################
###             merge OUTPUT_6hr
###########################################################

hy_output_6hr = read.csv("non_aki_output/6hr/eicu_eicu_crd_aki_negative_output_6hr.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("non_aki_output/6hr/EICU_OUTPUT6_label.tsv.txt", sep = "\t"))))
colnames(hy_output_6hr) = label

hy_output_6hr$mergekey = paste0(hy_output_6hr$ICUSTAY_ID, sprintf("_%06d", hy_output_6hr$OUTPUT_6HR_TIME))
head(hy_output_6hr)
#take average of every hours (since hours are converted from minutes, it exist duplications)
hy_output_6hr = aggregate(hy_output_6hr[, 1:3], list(hy_output_6hr$mergekey), mean, na.rm=TRUE)
colnames(hy_output_6hr)[colnames(hy_output_6hr) == "Group.1"] = "mergekey"
head(hy_output_6hr)

alldata.merged = full_join(alldata.merged, hy_output_6hr, by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
head(alldata.merged)
###########################################################
###             merge OUTPUT_12hr
###########################################################

hy_output_12hr = read.csv("non_aki_output/12hr/eicu_eicu_crd_aki_negative_output_12hr.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("non_aki_output/12hr/EICU_OUTPUT12_label.tsv.txt", sep = "\t"))))
colnames(hy_output_12hr) = label

hy_output_12hr$mergekey = paste0(hy_output_12hr$ICUSTAY_ID, sprintf("_%06d", hy_output_12hr$OUTPUT_12HR_TIME))
head(hy_output_12hr)
#take average of every hours (since hours are converted from minutes, it exist duplications)
hy_output_12hr = aggregate(hy_output_12hr[, 1:3], list(hy_output_12hr$mergekey), mean, na.rm=TRUE)
colnames(hy_output_12hr)[colnames(hy_output_12hr) == "Group.1"] = "mergekey"
head(hy_output_12hr)

alldata.merged = full_join(alldata.merged, hy_output_12hr, by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
head(alldata.merged)

###########################################################
###             merge OUTPUT_24hr
###########################################################

hy_output_24hr = read.csv("non_aki_output/24hr/eicu_eicu_crd_aki_negative_output_24hr.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("non_aki_output/24hr/EICU_OUTPUT24_label.tsv.txt", sep = "\t"))))
colnames(hy_output_24hr) = label

hy_output_24hr$mergekey = paste0(hy_output_24hr$ICUSTAY_ID, sprintf("_%06d", hy_output_24hr$OUTPUT_24HR_TIME))
head(hy_output_24hr)
#take average of every hours (since hours are converted from minutes, it exist duplications)
hy_output_24hr = aggregate(hy_output_24hr[, 1:3], list(hy_output_24hr$mergekey), mean, na.rm=TRUE)
colnames(hy_output_24hr)[colnames(hy_output_24hr) == "Group.1"] = "mergekey"
head(hy_output_24hr)

alldata.merged = full_join(alldata.merged, hy_output_24hr, by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
head(alldata.merged)

###########################################################
###             merge ABG
###########################################################
hy_abg = read.csv("non_aki_bg/eicu_eicu_crd_aki_negative_bg2.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("non_aki_bg/eicu_ABG_LABLE.TSV.txt", sep = "\t"))))
colnames(hy_abg) = label

hy_abg$mergekey = paste0(hy_abg$ICUSTAY_ID, sprintf("_%06d", hy_abg$ABG_TIME))
head(hy_abg)
#take average of every hours (since hours are converted from minutes, it exist duplications)
hy_abg = aggregate(hy_abg[, 1:6], list(hy_abg$mergekey), mean, na.rm=TRUE)
colnames(hy_abg)[colnames(hy_abg) == "Group.1"] = "mergekey"
head(hy_abg)

alldata.merged = full_join(alldata.merged, hy_abg, by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
head(alldata.merged)
###########################################################
###             merge mechvent_starttime
###########################################################
install.packages(c("dplyr"))
library(dplyr)
hy_mechvent = read.csv("non_aki_mechvent/starttime/eicu_eicu_crd_aki_negative_mechvent_starttime.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("non_aki_mechvent/starttime/EICU_AKI_mechvent_start.txt", sep = "\t"))))
colnames(hy_mechvent) = label
hy_mechvent$MECHVENT_STARTTIME = ceiling(hy_mechvent$MECHVENT_STARTTIME)
hy_mechvent$mergekey = paste0(hy_mechvent$ICUSTAY_ID, sprintf("_%06d", hy_mechvent$MECHVENT_STARTTIME))
head(hy_mechvent)
#take average of every hours (since hours are converted from minutes, it exist duplications)
hy_mechvent = aggregate(hy_mechvent[, 1:2], list(hy_mechvent$mergekey), mean, na.rm=TRUE)
colnames(hy_mechvent)[colnames(hy_mechvent) == "Group.1"] = "mergekey"
head(hy_mechvent)

alldata.merged = full_join(alldata.merged, hy_mechvent, by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
head(alldata.merged)
###########################################################
###             merge mechvent_endtime
###########################################################

hy_mechventend = read.csv("non_aki_mechvent/endtime/eicu_eicu_crd_aki_negative_mechvent_endtime.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("non_aki_mechvent/endtime/EICU_AKI_mechvent_end.txt", sep = "\t"))))
colnames(hy_mechventend) = label
hy_mechventend$MECHVENT_ENDTIME = ceiling(hy_mechventend$MECHVENT_ENDTIME)
hy_mechventend$mergekey = paste0(hy_mechventend$ICUSTAY_ID, sprintf("_%06d", hy_mechventend$MECHVENT_ENDTIME))
head(hy_mechventend)
#take average of every hours (since hours are converted from minutes, it exist duplications)
hy_mechventend = aggregate(hy_mechventend[, 1:2], list(hy_mechventend$mergekey), mean, na.rm=TRUE)
colnames(hy_mechventend)[colnames(hy_mechventend) == "Group.1"] = "mergekey"
head(hy_mechventend)

alldata.merged = full_join(alldata.merged, hy_mechventend, by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
head(alldata.merged)
###########################################################
###             merge vaso_STARTTIME
###########################################################
hy_VASO = read.csv("non_aki_vaso/eicu_eicu_crd_aki_negative_vaso.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("non_aki_vaso/EICU_vaso_label.tsv.txt", sep = "\t"))))
colnames(hy_VASO) = label
hy_VASO$VASO_TIME = as.character(anytime(hy_VASO$VASO_TIME))
hy_VASO$mergekey = paste0(hy_VASO$ICUSTAY_ID, "_", hy_VASO$VASO_TIME)

alldata.merged = full_join(alldata.merged, hy_VASO,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
head(alldata.merged)
colnames(alldata.merged)
###########################################################
###             Remove redundant data
###########################################################


  alldata.merged2 = alldata.merged[, c("ICUSTAY_ID.x","ICU_CLASS","ETHNICITY",
                                       "AGE","GENDER","LOS","HEIGHT","WEIGHT",
                                       "BMI","ISOFA","SEP","CAR","RES","NEU.x","OD","CKD","DIA","CHF",
                                       "CLD","CPD","HYP","ADMI_TIME","OUT_TIME","CURR_TIME","DIAS_BP","HR","SYS_BP","MEAN_BP",
                                       "RR","TEM","SPO2","PH","CA_ION","HGB",
                                       "WBC","RBC","HCT","PLT","RDW","CRP","HCO3",
                                       "ALT","AST","ALB","TBB","TNT",
                                       "CK","CKMB","CR","UN","AMI","LIP",
                                       "BNP","CL_ION","GLU","K_ION","NA_ION","APTT",
                                       "INR","FIB","LAC","AG","P_ION","MG_ION",
                                       "INPUT_6HR","INPUT_12HR","INPUT_24HR","OUTPUT_6HR","OUTPUT_12HR","OUTPUT_24HR",
                                       "FIO2","PCO2","PO2","MECHVENT_STARTTIME","MECHVENT_ENDTIME","VASO_TIME","VASO_LABEL")]
colnames(alldata.merged2)[1] = c("ICUSTAY_ID")
colnames(alldata.merged2)
alldata.merged2 = alldata.merged2[!is.na(alldata.merged2$ICUSTAY_ID), ]
save(alldata.merged2, file = "expanded.all.data.nonaki.eicu.merged0.Rdata")
load("expanded.all.data.nonaki.eicu.merged0.Rdata")
###########################################################
###       Sort data based on ICUSTAY_ID and TIME
###########################################################
library(anytime)
index2sort = paste0(sprintf("%07d",alldata.merged2$ICUSTAY_ID), "_", sprintf("%06d",alldata.merged2$CURR_TIME + abs(min(alldata.merged2$CURR_TIME))+1))
order = sort.int(index2sort, index.return = T)$ix
alldata.merged2 = alldata.merged2[order,]
save(alldata.merged2, file = "expanded.all.data.nonaki.eicu.merged.Rdata")


###########################################################
###             Fill missing values
###########################################################
##setwd("/home/zhihuan/Documents/Cong_Feng/20180908_Hypoxemia/Hypoxemia - LSTM/PFéçåµ/is_vent")
library(zoo)
library(anytime)
load("expanded.all.data.nonaki.eicu.merged.Rdata")
colnames(alldata.merged2)
alldata.merged2[is.na(alldata.merged2$INPUT_6HR), "INPUT_6HR"] = 0
alldata.merged2[is.na(alldata.merged2$INPUT_12HR), "INPUT_12HR"] = 0
alldata.merged2[is.na(alldata.merged2$INPUT_24HR), "INPUT_24HR"] = 0
alldata.merged2[is.na(alldata.merged2$OUTPUT_6HR), "OUTPUT_6HR"] = 0
alldata.merged2[is.na(alldata.merged2$OUTPUT_12HR), "OUTPUT_12HR"] = 0
alldata.merged2[is.na(alldata.merged2$OUTPUT_24HR), "OUTPUT_24HR"] = 0
alldata.merged2$INPUT_MINUS_OUTPUT_6HR = alldata.merged2$INPUT_6HR - alldata.merged2$OUTPUT_6HR
alldata.merged2$INPUT_MINUS_OUTPUT_12HR = alldata.merged2$INPUT_12HR - alldata.merged2$OUTPUT_12HR
alldata.merged2$INPUT_MINUS_OUTPUT_24HR = alldata.merged2$INPUT_24HR - alldata.merged2$OUTPUT_24HR


alldata.merged2$HEIGHT[is.na(alldata.merged2$HEIGHT)] = mean(alldata.merged2$HEIGHT, na.rm=TRUE)
alldata.merged2$WEIGHT[is.na(alldata.merged2$WEIGHT)] = mean(alldata.merged2$WEIGHT, na.rm=TRUE)
alldata.merged2$BMI = alldata.merged2$WEIGHT / (alldata.merged2$HEIGHT/100)^2


# is.na(alldata.merged2$FIO2)
icustayIDlist = unique(alldata.merged2$ICUSTAY_ID)
### Expanding original alldata
newrows1 = data.frame(matrix(-Inf, nrow=(length(icustayIDlist)-1), dim(alldata.merged2)[2]))
newrows2 = data.frame(matrix(-Inf, nrow=(length(icustayIDlist)-1), dim(alldata.merged2)[2]))
colnames(newrows1) = colnames(alldata.merged2)
colnames(newrows2) = colnames(alldata.merged2)
newrows1$ICUSTAY_ID = unique(alldata.merged2$ICUSTAY_ID)[1:length(unique(alldata.merged2$ICUSTAY_ID))-1]
newrows2$ICUSTAY_ID = unique(alldata.merged2$ICUSTAY_ID)[2:length(unique(alldata.merged2$ICUSTAY_ID))]
newrows1$CURR_TIME = 999999
newrows2$CURR_TIME = 0
newrows1 = data.frame(newrows1)
newrows2 = data.frame(newrows2)



alldata.merged2$LIP = as.numeric(alldata.merged2$LIP)
alldata.merged2$BNP = as.numeric(alldata.merged2$BNP)
##alldata.merged2$DD = as.numeric(alldata.merged2$DD)
##alldata.merged2$FIB = as.numeric(alldata.merged2$FIB)
alldata.merged2$MECHVENT_STARTTIME = as.numeric(anytime(alldata.merged2$MECHVENT_STARTTIME))
alldata.merged2$MECHVENT_ENDTIME = as.numeric(anytime(alldata.merged2$MECHVENT_ENDTIME))
alldata.merged2$VASO_TIME = as.numeric(anytime(alldata.merged2$VASO_TIME))
alldata.merged2$ADMI_TIME = as.numeric(anytime(alldata.merged2$ADMI_TIME))
alldata.merged2$OUT_TIME = as.numeric(anytime(alldata.merged2$OUT_TIME))


tobind = rbind(newrows1, newrows2)
alldata.merged3 = rbind(alldata.merged2, tobind)
alldata.merged2$CURR_TIME_temp = alldata.merged2$CURR_TIME + abs(min(alldata.merged2$CURR_TIME))+1
index2sort.1 = paste0(sprintf("%07d",alldata.merged2$ICUSTAY_ID), sprintf("_%06d", alldata.merged2$CURR_TIME_temp))
index2sort.2 = paste0(sprintf("%07d",tobind$ICUSTAY_ID), sprintf("_%06d", tobind$CURR_TIME))
index2sort = c(index2sort.1, index2sort.2)
order = sort.int(index2sort, index.return = T)$ix
alldata.merged3 = alldata.merged3[order,]

save(alldata.merged3, file = "expanded.all.inputoutput.missingfill1.Rdata")
# load("expanded.all.inputoutput.missingfill1.Rdata")

colnames(alldata.merged3)

############## lab vit fio2????????????
columns2impute = colnames(alldata.merged3)[c(25:64,72:73)]
linear.impute = na.approx(zoo(alldata.merged3[, columns2impute]), na.rm = F)
linear.impute[is.nan(linear.impute)] = NA
alldata.merged3[, columns2impute] = as.data.frame(linear.impute)
alldata.merged3 = alldata.merged3[alldata.merged3$ICU_CLASS != -Inf,]
alldata.merged3[alldata.merged3 == -Inf] = NA




alldata.merged3 = rbind(alldata.merged3, tobind)
index2sort.1 = paste0(sprintf("%07d",alldata.merged2$ICUSTAY_ID), sprintf("_%06d", alldata.merged2$CURR_TIME_temp))
index2sort.2 = paste0(sprintf("%07d",tobind$ICUSTAY_ID), sprintf("_%06d", tobind$CURR_TIME))
index2sort = c(index2sort.1, index2sort.2)
order = sort.int(index2sort, index.return = T)$ix
alldata.merged3 = alldata.merged3[order,]
save(alldata.merged3, file = "expanded.all.labvitbglinear.missingfill.Rdata")



#######################ABG LAB VIT FIO2?????????
columns4impute = colnames(alldata.merged3)[c(25:64,71:73)]
front.impute = data.frame(na.locf(zoo(alldata.merged3[, columns4impute]), na.rm = F))
alldata.merged3[, columns4impute] = front.impute
alldata.merged3 = alldata.merged3[alldata.merged3$ICU_CLASS != -Inf,]
alldata.merged3[alldata.merged3 == -Inf] = NA

alldata.merged3 = rbind(alldata.merged3, tobind)
index2sort.1 = paste0(sprintf("%07d",alldata.merged2$ICUSTAY_ID), sprintf("_%06d", alldata.merged2$CURR_TIME_temp))
index2sort.2 = paste0(sprintf("%07d",tobind$ICUSTAY_ID), sprintf("_%06d", tobind$CURR_TIME))
index2sort = c(index2sort.1, index2sort.2)
order = sort.int(index2sort, index.return = T)$ix
alldata.merged3 = alldata.merged3[order,]
save(alldata.merged3, file = "expanded.all.ABGLABVITFIO2frontimpute.missingfill.Rdata")



#######################ABG LAB VIT FIO2 ?????????
columns3impute = colnames(alldata.merged3)[c(25:64,71:73)]
back.impute = data.frame(na.locf(zoo(alldata.merged3[, columns3impute]), na.rm = F, fromLast = T))
alldata.merged3[, columns3impute] = back.impute
alldata.merged3 = alldata.merged3[alldata.merged3$ICU_CLASS != -Inf,]
alldata.merged3[alldata.merged3 == -Inf] = NA

alldata.merged3 = rbind(alldata.merged3, tobind)
index2sort.1 = paste0(sprintf("%07d",alldata.merged2$ICUSTAY_ID), sprintf("_%06d", alldata.merged2$CURR_TIME_temp))
index2sort.2 = paste0(sprintf("%07d",tobind$ICUSTAY_ID), sprintf("_%06d", tobind$CURR_TIME))
index2sort = c(index2sort.1, index2sort.2)
order = sort.int(index2sort, index.return = T)$ix
alldata.merged3 = alldata.merged3[order,]



save(alldata.merged3, file = "expanded.all.ABGLABVITFIO2backimpute.missingfill.nonaki.eicu.Rdata")
load("expanded.all.ABGLABVITFIO2backimpute.missingfill.nonaki.eicu.Rdata")
install.packages(c("zoo"))
library(zoo)

# mechvent time
alldata.merged3$MECH = NA
alldata.merged3$MECH = alldata.merged3$MECHVENT_STARTTIME
alldata.merged3$MECH[!is.na(alldata.merged3$MECHVENT_ENDTIME)] = alldata.merged3$MECHVENT_ENDTIME[!is.na(alldata.merged3$MECHVENT_ENDTIME)]
front.impute = data.frame(na.locf(zoo(alldata.merged3$MECH), na.rm = F))
back.impute = data.frame(na.locf(zoo(alldata.merged3$MECH), na.rm = F, fromLast = T))
mask = front.impute * back.impute
mask[mask == Inf] = 0
mask[mask == -Inf] = 0
mask[mask != 0] = 1
alldata.merged3$MECH = mask


alldata.merged3 = alldata.merged3[alldata.merged3$ICU_CLASS != -Inf,]
alldata.merged3[alldata.merged3 == -Inf] = NA

# ADD P_F_ratio
alldata.merged3$PF = alldata.merged3$PO2/alldata.merged3$FIO2

save(alldata.merged3, file = "expanded.all.vasomechventimpute.nonaki.eicu.Rdata")
load("expanded.all.vasomechventimpute.nonaki.eicu.Rdata")
colnames(alldata.merged3)

### Add new column: hours
currtime = as.numeric(anytime(alldata.merged3$CURR_TIME))
icustayIDtable = table(alldata.merged3$ICUSTAY_ID)
icustayIDtable_cumsum = cumsum(table(alldata.merged3$ICUSTAY_ID))
hours = unlist(sapply(icustayIDtable, function(x) rep(1:x)))
alldata.merged3$HOURS = hours
alldata.merged3$AKI = 0
alldata.merged3$AKI_PF = 0
alldata.merged3$AKI_PF[alldata.merged3$P_F_ratio < 3] = 1
alldata.merged3$AKI_PF[alldata.merged3$P_F_ratio < 2] = 2
alldata.merged3$AKI_PF[alldata.merged3$P_F_ratio < 1] = 3
alldata.merged3$AKI_BMI = 0
alldata.merged3$AKI_BMI[alldata.merged3$BMI >= 40] = 4
alldata.merged3$AKI_BMI[alldata.merged3$BMI < 40] = 3
alldata.merged3$AKI_BMI[alldata.merged3$BMI < 25] = 2
alldata.merged3$AKI_BMI[alldata.merged3$BMI < 23] = 1
alldata.merged3$AKI_BMI[alldata.merged3$BMI < 18.5] = 0
save(alldata.merged3, file = "expanded.all.data.merged.imputed.calculated.eicu.nonAKI.Rdata")

###########################################################
###########################################################
###########################################################
###  Merge aki and non_aki
###########################################################
###########################################################
###########################################################
#temp = alldata.imputed.aki$MECHVENT_TIME
#alldata.imputed.aki$MECHVENT_TIME = as.numeric(alldata.imputed.aki$MECHVENT_TIME)
#alldata.imputed.aki$MECHVENT_TIME = unlist(alldata.imputed.aki$MECHVENT_TIME)

###typeof(rownames(alldata.merged3)[1])
##????????????0??????8??????
##rownames(alldata.merged3) = str_pad(rownames(alldata.merged3), 8, pad = "0")
##????????????????????????1
##rownames(alldata.merged3)= paste0("1",rownames(alldata.merged3))
##rownames(nalldata.merged3) = str_pad(rownames(nalldata.merged3), 8, pad = "0")
##rownames(nalldata.merged3)= paste0("2",rownames(nalldata.merged3))
library(stringr)

setwd("C:/Users/congf/Desktop/AKI/aki_kidgo_EICU_2019.3.13/non_aki")
load("expanded.all.data.merged.imputed.calculated.eicu.nonAKI.Rdata")
alldata.imputed.naki = alldata.merged3



setwd("C:/Users/congf/Desktop/AKI/aki_kidgo_EICU_2019.3.13/aki")
load("expanded.all.data.merged.imputed.calculated.eicu.AKI.Rdata")
alldata.imputed.aki = alldata.merged3


# If rbind fails, it means one of the attributes of the data frame was a list
alldata.imputed.aki$MECH = unlist(alldata.imputed.aki$MECH)
alldata.imputed.aki$MECH[is.na(alldata.imputed.aki$MECH)] = 0
alldata.imputed.aki$VASO_TIME = unlist(alldata.imputed.aki$VASO_TIME)
alldata.imputed.aki$VASO_LABEL[is.na(alldata.imputed.aki$VASO_LABEL)] = 0
alldata.imputed.aki$VASO = unlist(alldata.imputed.aki$VASO_LABEL)
alldata.imputed.aki$GENDER[is.na(alldata.imputed.aki$GENDER)] = 1
alldata.imputed.aki$ETHNICITY[is.na(alldata.imputed.aki$ETHNICITY)] = 1


alldata.imputed.naki$MECH = unlist(alldata.imputed.naki$MECH)
alldata.imputed.naki$MECH[is.na(alldata.imputed.naki$MECH)] = 0
alldata.imputed.naki$VASO_TIME = unlist(alldata.imputed.naki$VASO_TIME)
alldata.imputed.naki$VASO_LABEL[is.na(alldata.imputed.naki$VASO_LABEL)] = 0
alldata.imputed.naki$VASO = unlist(alldata.imputed.naki$VASO_LABEL)
alldata.imputed.naki$GENDER[is.na(alldata.imputed.naki$GENDER)] = 1
alldata.imputed.naki$ETHNICITY[is.na(alldata.imputed.naki$ETHNICITY)] = 1



alldata.imputed = rbind(alldata.imputed.aki, alldata.imputed.naki)
alldata.imputed$NEU=alldata.imputed$NEU.x
save(alldata.imputed, file = "expanded.all.data.merged.imputed.calculated.merged.eicu.Rdata")
load("expanded.all.data.merged.imputed.calculated.merged.eicu.Rdata")

###########################################################
###             Remove some colnames data which data is all missing
###########################################################
colnames(alldata.imputed)



alldata.final = alldata.imputed[, c("ICUSTAY_ID","ICU_CLASS","ETHNICITY",
                                    "AGE","GENDER","HEIGHT","WEIGHT",
                                    "AKI_BMI","ISOFA","SEP","CAR","RES","NEU","OD","CKD","DIA","CHF",
                                    "CLD","CPD","HYP","DIAS_BP","HR","SYS_BP","MEAN_BP",
                                    "RR","TEM","SPO2","PH","CA_ION","HGB",
                                    "WBC","RBC","HCT","PLT","HCO3",
                                    "ALT","AST","ALB","TBB",
                                    "CR","UN",
                                    "CL_ION","GLU","K_ION","NA_ION","APTT",
                                    "INR","LAC","AG","P_ION","MG_ION",
                                    "INPUT_6HR","INPUT_12HR","INPUT_24HR","OUTPUT_6HR","OUTPUT_12HR","OUTPUT_24HR",
                                    "AKI_PF",
                                    "INPUT_MINUS_OUTPUT_6HR","INPUT_MINUS_OUTPUT_12HR","INPUT_MINUS_OUTPUT_24HR",
                                    "MECH","VASO","AKI","HOURS")]
colnames(alldata.final)
save(alldata.final, file = "expanded.Remove.some.colnames.missingdata.aki.eicu.Rdata")
###########################################################
### Drop some columns, make complete table
###########################################################

load("expanded.Remove.some.colnames.missingdata.aki.eicu.Rdata")
colnames(alldata.final)

droplist = NULL
#droplist = c("FIO2","PO2","PCO2","SUBJECT_ID","HADM_ID","LOS")
# droplist = c("HEIGHT","WEIGHT","BMI")
remains = NULL
for (i in 1:dim(alldata.final)[2]){
  message(colnames(alldata.final)[i], "\t", sum(is.na(alldata.final[,i])) )
  if (sum(is.na(alldata.final[,i])) > 300000){
    droplist = c(droplist, colnames(alldata.final)[i])
  }
  else {
    remains = c(remains, colnames(alldata.final)[i])
  }
}

droplist = c("PH","LAC","P_ION","APTT","INR","MG_ION","ALB")

#remains = c(remains, "FIO2", "PO2", "PCO2", "P_F_ratio", "IS_VENT", "IS_VENT_P_F_ratio_target")
todrop = which(names(alldata.final) %in% droplist)
if (length(todrop) == 0){
  alldata.imputed.shrinked = alldata.final[complete.cases(alldata.final), remains]
} else {
  alldata.imputed.shrinked = alldata.final[complete.cases(alldata.final[,-todrop]), -todrop]
}
sum(is.na(alldata.imputed.shrinked)) == 0
length(unique(alldata.imputed.shrinked$ICUSTAY_ID))




save(alldata.imputed.shrinked, file="AKI.EICU.FINAL.Rdata")
load("AKI.EICU.FINAL.Rdata")
# write.csv(alldata.imputed.shrinked, file = "expanded.all.data.merged.imputed.calculated.shrinked.csv", row.names = F)
write.csv(alldata.imputed.shrinked, file = "AKI.EICU.FINAL.csv", row.names = F)


colnames(alldata.imputed.shrinked)
alldata.ed.eicu.aki = alldata.imputed.shrinked[,c("ICUSTAY_ID","ICU_CLASS","ETHNICITY",
                                                  "AGE","GENDER","HEIGHT","WEIGHT",
                                                  "BMI","DIAS_BP","HR","SYS_BP","MEAN_BP",
                                                  "RESPRATE","Temprate","SPO2",
                                                  "MECHVENT_TIME","VASO_TIME","AKI","HOURS")]

save(alldata.ed.eicu.aki, file="alldata.ed.eicu.aki.Rdata")

# write.csv(alldata.imputed.shrinked, file = "expanded.all.data.merged.imputed.calculated.shrinked.csv", row.names = F)
write.csv(alldata.ed.eicu.aki, file = "alldata.ed.eicu.aki.csv", row.names = F)



