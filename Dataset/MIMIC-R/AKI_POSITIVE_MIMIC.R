setwd("/media/zhihuan/Drive3/20191109_AKI_python/Dataset/MIMIC11.11/aki/")
install.packages(c("anytime", "tictoc","lubridate"))
library("anytime")
library("tictoc")
library("lubridate")

####################################################################
####################################################################
############################ MIMIC_AKI###############################
####################################################################
####################################################################

###########################################################
###            Expand data
###########################################################
alldata = read.csv("aki_cohort/mimic_mimiciii_aki_kidgo_positive_demo3.csv", header = F, stringsAsFactors = F)
label = levels(unlist(read.table("aki_cohort/mimic_label.tsv.txt", sep = "\t")))
colnames(alldata) = label

alldata$ADMI_TIME = anytime(alldata$ADMI_TIME)
alldata$OUT_TIME = anytime(alldata$OUT_TIME)
alldata$ADMI_TIME_numeric = as.numeric(alldata$ADMI_TIME)
alldata$OUT_TIME_numeric = as.numeric(alldata$OUT_TIME)
alldata = alldata[alldata$AGE < 150, ]
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
  time.in = row$ADMI_TIME_numeric
  time.out = row$OUT_TIME_numeric
  hours = (time.out - time.in)/3600 + 1
  hourslist = c(hourslist, hours)
  allhours = allhours + hours
}
alldata$hours = hourslist
# expand
alldata.expanded <- alldata[rep(row.names(alldata), alldata$hours), 1:25]

start_time <- proc.time()
currtime_list = NULL
allhours = 0
for (i in 1:dim(alldata)[1]){
  if (i %% 2500 == 0){
    message(i, "  -  ", allhours, "    Elapsed: ", round((proc.time() - start_time)[3], 2), " secs")
  }
  row = alldata[i,]
  time.in = row$ADMI_TIME_numeric
  time.out = row$OUT_TIME_numeric
  hours = (time.out - time.in)/3600
  allhours = allhours + hours
  
  numerical.times = seq(time.in, time.in + 3600*hours, by=3600)
  currtime_list[[i]] = numerical.times
}
currtime_list2 = unlist(currtime_list)#do.call(c, unlist(currtime_list, recursive=FALSE))
times = as.character(anytime(currtime_list2))

alldata.expanded$CURR_TIME = times

save(alldata.expanded, file = "expanded.all.data.Rdata")

load("expanded.all.data.Rdata")

###########################################################
###             merge vit
###########################################################
library(dplyr)
library(anytime)

hy_vit = read.csv("aki_vit/mimic_mimiciii_aki_kidgo_positive_vit.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("aki_vit/MIMIC_vit_label.tsv.txt", sep = "\t"))))
colnames(hy_vit) = label
hy_vit$VIT_TIME = as.character(anytime(hy_vit$VIT_TIME))

alldata.expanded$mergekey = paste0(alldata.expanded$ICUSTAY_ID, "_", alldata.expanded$CURR_TIME)
hy_vit$mergekey = paste0(hy_vit$ICUSTAY_ID, "_", hy_vit$VIT_TIME)
alldata.merged = full_join(alldata.expanded, hy_vit,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]

###########################################################
###             merge aki_lab
###########################################################
hy_lab3 = read.csv("aki_lab/mimic_mimiciii_aki_kidgo_positive_lab.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("aki_lab/mimic_LAB_LABEL.TSV.txt", sep = "\t"))))
colnames(hy_lab3) = label
hy_lab3$LAB_TIME = as.character(anytime(hy_lab3$LAB_TIME))
hy_lab3$mergekey = paste0(hy_lab3$ICUSTAY_ID, "_", hy_lab3$LAB_TIME)

alldata.merged = full_join(alldata.merged, hy_lab3,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
###########################################################
###             merge INPUT_6hr
###########################################################
hy_input_6hr = read.csv("aki_input/input_6hr/mimic_mimiciii_aki_kidgo_positive_input_6h.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("aki_input/input_6hr/input_6hr.tsv.txt", sep = "\t"))))
colnames(hy_input_6hr) = label
hy_input_6hr$INPUT_6HR_CHARTTIME = as.character(anytime(hy_input_6hr$INPUT_6HR_CHARTTIME))
hy_input_6hr$mergekey = paste0(hy_input_6hr$ICUSTAY_ID, "_", hy_input_6hr$INPUT_6HR_CHARTTIME)

alldata.merged = full_join(alldata.merged, hy_input_6hr,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
###########################################################
###             merge INPUT_12hr
###########################################################
hy_input_12hr = read.csv("aki_input/input_12hr/mimic_mimiciii_aki_kidgo_positive_input_12h.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("aki_input/input_12hr/input_12hr.tsv.txt", sep = "\t"))))
colnames(hy_input_12hr) = label
hy_input_12hr$INPUT_12HR_CHARTTIME = as.character(anytime(hy_input_12hr$INPUT_12HR_CHARTTIME))
hy_input_12hr$mergekey = paste0(hy_input_12hr$ICUSTAY_ID, "_", hy_input_12hr$INPUT_12HR_CHARTTIME)

alldata.merged = full_join(alldata.merged, hy_input_12hr,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
###########################################################
###             merge INPUT_24hr
###########################################################
hy_input_24hr = read.csv("aki_input/input_24hr/mimic_mimiciii_aki_kidgo_positive_input_24h.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("aki_input/input_24hr/input_24hr.tsv.txt", sep = "\t"))))
colnames(hy_input_24hr) = label
hy_input_24hr$INPUT_24HR_CHARTTIME = as.character(anytime(hy_input_24hr$INPUT_24HR_CHARTTIME))
hy_input_24hr$mergekey = paste0(hy_input_24hr$ICUSTAY_ID, "_", hy_input_24hr$INPUT_24HR_CHARTTIME)

alldata.merged = full_join(alldata.merged, hy_input_24hr,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]

###########################################################
###             merge outPUT_6hr
###########################################################
hy_output_6hr = read.csv("aki_output/output_6hr/mimic_mimiciii_aki_kidgo_positive_output_6h.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("aki_output/output_6hr/output_6hr.tsv.txt", sep = "\t"))))
colnames(hy_output_6hr) = label
hy_output_6hr$OUT_6HR_TIME = as.character(anytime(hy_output_6hr$OUTPUT_6HR_TIME))
hy_output_6hr$mergekey = paste0(hy_output_6hr$ICUSTAY_ID, "_", hy_output_6hr$OUTPUT_6HR_TIME)

alldata.merged = full_join(alldata.merged, hy_output_6hr,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
###########################################################
###             merge outPUT_12hr
###########################################################
hy_output_12hr = read.csv("aki_output/output_12hr/mimic_mimiciii_aki_kidgo_positive_output_12h.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("aki_output/output_12hr/output_12hr.tsv.txt", sep = "\t"))))
colnames(hy_output_12hr) = label
hy_output_12hr$OUT_12HR_TIME = as.character(anytime(hy_output_12hr$OUTPUT_12HR_TIME))
hy_output_12hr$mergekey = paste0(hy_output_12hr$ICUSTAY_ID, "_", hy_output_12hr$OUTPUT_12HR_TIME)

alldata.merged = full_join(alldata.merged, hy_output_12hr,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
###########################################################
###             merge outPUT_24hr
###########################################################
hy_output_24hr = read.csv("aki_output/output_24hr/mimic_mimiciii_aki_kidgo_positive_output_24h.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("aki_output/output_24hr/output_24hr.tsv.txt", sep = "\t"))))
colnames(hy_output_24hr) = label
hy_output_24hr$OUT_24HR_TIME = as.character(anytime(hy_output_24hr$OUTPUT_24HR_TIME))
hy_output_24hr$mergekey = paste0(hy_output_24hr$ICUSTAY_ID, "_", hy_output_24hr$OUTPUT_24HR_TIME)

alldata.merged = full_join(alldata.merged, hy_output_24hr,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
###########################################################
###             merge ABG
###########################################################
hy_abg = read.csv("aki_bg/mimic_mimiciii_aki_kidgo_positive_abg.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("aki_bg/mimic_ABG_LABLE.TSV.txt", sep = "\t"))))
colnames(hy_abg) = label
hy_abg$ABG_TIME = as.character(anytime(hy_abg$ABG_TIME))
hy_abg$mergekey = paste0(hy_abg$ICUSTAY_ID, "_", hy_abg$ABG_TIME)

alldata.merged = full_join(alldata.merged, hy_abg,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
colnames(alldata.merged)
###########################################################
###             merge mechvent_starttime
###########################################################
install.packages(c("dplyr"))
library(dplyr)
hy_mechvent = read.csv("aki_mechvent/mechvent_starttime/mimic_mimiciii_aki_kidgo_positive_mechvent_starttime.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("aki_mechvent/mechvent_starttime/MIMIC_AKI_mechvent.txt", sep = "\t"))))
colnames(hy_mechvent) = label
hy_mechvent$MECHVENT_STARTTIME = as.character(anytime(hy_mechvent$MECHVENT_STARTTIME))

hy_mechvent$mergekey = paste0(hy_mechvent$ICUSTAY_ID, "_", hy_mechvent$MECHVENT_STARTTIME)

alldata.merged = full_join(alldata.merged, hy_mechvent,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]

###########################################################
###             merge mechvent_endtime
###########################################################
hy_mechvent2 = read.csv("aki_mechvent/mechvent_endtime/mimic_mimiciii_aki_kidgo_positive_mechvent_endtime.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("aki_mechvent/mechvent_endtime/MIMIC_AKI_mechvent.txt", sep = "\t"))))
colnames(hy_mechvent2) = label
hy_mechvent2$MECHVENT_ENDTIME = as.character(anytime(hy_mechvent2$MECHVENT_ENDTIME))

hy_mechvent2$mergekey = paste0(hy_mechvent2$ICUSTAY_ID, "_", hy_mechvent2$MECHVENT_ENDTIME)

alldata.merged = full_join(alldata.merged, hy_mechvent2,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
###########################################################
###             merge vaso_STARTTIME
###########################################################
hy_VASO = read.csv("aki_vaso/vaso_start/mimic_mimiciii_aki_kidgo_positive_vasopressor_starttime.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("aki_vaso/vaso_start/MIMIC_vaso_label.tsv.txt", sep = "\t"))))
colnames(hy_VASO) = label
hy_VASO$VASO_STARTTIME = as.character(anytime(hy_VASO$VASO_STARTTIME))
hy_VASO$mergekey = paste0(hy_VASO$ICUSTAY_ID, "_", hy_VASO$VASO_STARTTIME)

alldata.merged = full_join(alldata.merged, hy_VASO,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]
###########################################################
###             merge vaso_ENDTIME
###########################################################
hy_VASO2 = read.csv("aki_vaso/vaso_end/mimic_mimiciii_aki_kidgo_positive_vasopressor_endtime.csv", header = F, stringsAsFactors = F)
label = as.character(unlist(t(read.table("aki_vaso/vaso_end/MIMIC_vaso_label.tsv.txt", sep = "\t"))))
colnames(hy_VASO2) = label
hy_VASO2$VASO_ENDTIME = as.character(anytime(hy_VASO2$VASO_ENDTIME))
hy_VASO2$mergekey = paste0(hy_VASO2$ICUSTAY_ID, "_", hy_VASO2$VASO_ENDTIME)

alldata.merged = full_join(alldata.merged, hy_VASO2,
                           by = 'mergekey')
alldata.merged <- alldata.merged[!duplicated(alldata.merged$mergekey),]


#alldata.merged$OUT_6HR=alldata.merged$OUTPUT_6HR
#alldata.merged$OUT_12HR=alldata.merged$OUTPUT_12HR
#alldata.merged$OUT_24HR=alldata.merged$OUTPUT_24HR
colnames(alldata.merged)

###########################################################
###             Remove redundant data
###########################################################

alldata.merged2 = alldata.merged[, c("ICUSTAY_ID.x","ICU_CLASS","ETHNICITY",
                                     "AGE","GENDER","LOS","HEIGHT","WEIGHT",
                                     "BMI","ISOFA","SEP","CAR","RES","NEU","OD","CKD","DIA","CHF",
                                     "CLD","CPD","HYP",
                                     "ADMI_TIME","OUT_TIME","CURR_TIME","DIAS_BP","HR","SYS_BP","MEAN_BP",
                                     "RR","TEM","SPO2","PH","CA_ION","HGB",
                                     "WBC","RBC","NEU_PER","HCT","PLT","RDW","CRP","HCO3",
                                     "ALT","AST","ALB","TBB","TNT",
                                     "CK","CKMB","CR","UN","AMI","LIP",
                                     "BNP","CL_ION","GLU","K_ION","NA_ION","APTT",
                                     "PT","INR","DD","FIB","LAC","AG","P_ION","MG_ION",
                                     "INPUT_6HR","INPUT_12HR","INPUT_24HR","OUTPUT_6HR","OUTPUT_12HR","OUTPUT_24HR",
                                     "FIO2","PCO2","PO2","MECHVENT_STARTTIME","MECHVENT_ENDTIME","VASO_STARTTIME","VASO_ENDTIME")]
colnames(alldata.merged2)[1] = c("ICUSTAY_ID")
colnames(alldata.merged2)
alldata.merged2 = alldata.merged2[!is.na(alldata.merged2$ICUSTAY_ID), ]
###########################################################
###       Sort data based on ICUSTAY_ID and TIME
###########################################################
library(anytime)
index2sort = paste0(alldata.merged2$ICUSTAY_ID, "_", alldata.merged2$CURR_TIME)
order = sort.int(index2sort, index.return = T)$ix
alldata.merged2 = alldata.merged2[order,]
save(alldata.merged2, file = "expanded.all.data.merged.Rdata")


###########################################################
###             Fill missing values
###########################################################
##setwd("/home/zhihuan/Documents/Cong_Feng/20180908_Hypoxemia/Hypoxemia - LSTM/PFéçåµ/is_vent")
library(zoo)
library(anytime)
load("expanded.all.data.merged.Rdata")
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
newrows1$CURR_TIME = "9999-12-31 23:59:59"
newrows2$CURR_TIME = "0001-01-01 00:00:00"
newrows1 = data.frame(newrows1)
newrows2 = data.frame(newrows2)


alldata.merged2$NEU_PER = as.numeric(alldata.merged2$NEU_PER)
alldata.merged2$LIP = as.numeric(alldata.merged2$LIP)
alldata.merged2$BNP = as.numeric(alldata.merged2$BNP)
alldata.merged2$DD = as.numeric(alldata.merged2$DD)
alldata.merged2$FIB = as.numeric(alldata.merged2$FIB)
alldata.merged2$MECHVENT_STARTTIME = as.numeric(anytime(alldata.merged2$MECHVENT_STARTTIME))
alldata.merged2$MECHVENT_ENDTIME = as.numeric(anytime(alldata.merged2$MECHVENT_ENDTIME))
alldata.merged2$VASO_STARTTIME = as.numeric(anytime(alldata.merged2$VASO_STARTTIME))
alldata.merged2$VASO_ENDTIME = as.numeric(anytime(alldata.merged2$VASO_ENDTIME))
alldata.merged2$ADMI_TIME = as.numeric(anytime(alldata.merged2$ADMI_TIME))
alldata.merged2$OUT_TIME = as.numeric(anytime(alldata.merged2$OUT_TIME))



tobind = rbind(newrows1, newrows2)
alldata.merged3 = rbind(alldata.merged2, tobind)
index2sort = paste0(alldata.merged3$ICUSTAY_ID, "_", alldata.merged3$CURR_TIME)
order = sort.int(index2sort, index.return = T)$ix
alldata.merged3 = alldata.merged3[order,]
save(alldata.merged3, file = "expanded.all.inputoutput.missingfill1.Rdata")
colnames(alldata.merged3)
############## lab vit fio2????????????
columns2impute = colnames(alldata.merged3)[c(25:67,75:76)]
linear.impute = na.approx(zoo(alldata.merged3[, columns2impute]), na.rm = F)
linear.impute[is.nan(linear.impute)] = NA
alldata.merged3[, columns2impute] = as.data.frame(linear.impute)
alldata.merged3 = alldata.merged3[alldata.merged3$ICU_CLASS != -Inf,]
alldata.merged3[alldata.merged3 == -Inf] = NA

alldata.merged3 = rbind(alldata.merged3, rbind(newrows1, newrows2))
index2sort = paste0(alldata.merged3$ICUSTAY_ID, "_", alldata.merged3$CURR_TIME)
order = sort.int(index2sort, index.return = T)$ix
alldata.merged3 = alldata.merged3[order,]
save(alldata.merged3, file = "expanded.all.labvitbglinear.missingfill.Rdata")

#######################ABG LAB VIT FIO2?????????
columns4impute = colnames(alldata.merged3)[c(25:67,74:76)]
front.impute = data.frame(na.locf(zoo(alldata.merged3[, columns4impute]), na.rm = F))
alldata.merged3[, columns4impute] = front.impute
alldata.merged3 = alldata.merged3[alldata.merged3$ICU_CLASS != -Inf,]
alldata.merged3[alldata.merged3 == -Inf] = NA

alldata.merged3 = rbind(alldata.merged3, rbind(newrows1, newrows2))
index2sort = paste0(alldata.merged3$ICUSTAY_ID, "_", alldata.merged3$CURR_TIME)
order = sort.int(index2sort, index.return = T)$ix
alldata.merged3 = alldata.merged3[order,]
save(alldata.merged3, file = "expanded.all.ABGLABVITFIO2frontimpute.missingfill.Rdata")



#######################ABG LAB VIT FIO2 ?????????
columns3impute = colnames(alldata.merged3)[c(25:67,74:76)]
back.impute = data.frame(na.locf(zoo(alldata.merged3[, columns3impute]), na.rm = F, fromLast = T))
alldata.merged3[, columns3impute] = back.impute
alldata.merged3 = alldata.merged3[alldata.merged3$ICU_CLASS != -Inf,]
alldata.merged3[alldata.merged3 == -Inf] = NA

alldata.merged3 = rbind(alldata.merged3, rbind(newrows1, newrows2))
index2sort = paste0(alldata.merged3$ICUSTAY_ID, "_", alldata.merged3$CURR_TIME)
order = sort.int(index2sort, index.return = T)$ix
alldata.merged3 = alldata.merged3[order,]
save(alldata.merged3, file = "expanded.all.ABGLABVITFIO2backimpute.missingfill.Rdata")
load("expanded.all.ABGLABVITFIO2backimpute.missingfill.Rdata")
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




# vaso time
alldata.merged3$VASO = NA
alldata.merged3$VASO = alldata.merged3$VASO_STARTTIME
alldata.merged3$VASO[!is.na(alldata.merged3$VASO_ENDTIME)] = alldata.merged3$VASO_ENDTIME[!is.na(alldata.merged3$VASO_ENDTIME)]
front.impute = data.frame(na.locf(zoo(alldata.merged3$VASO), na.rm = F))
back.impute = data.frame(na.locf(zoo(alldata.merged3$VASO), na.rm = F, fromLast = T))
mask = front.impute * back.impute
mask[mask == Inf] = 0
mask[mask == -Inf] = 0
mask[mask != 0] = 1
alldata.merged3$VASO = mask

alldata.merged3 = alldata.merged3[alldata.merged3$ICU_CLASS != -Inf,]
alldata.merged3[alldata.merged3 == -Inf] = NA




alldata.merged3$PF = alldata.merged3$PO2/alldata.merged3$FIO2

save(alldata.merged3, file = "expanded.all.vasomechventimpute.missingfill.Rdata")
load("expanded.all.vasomechventimpute.missingfill.Rdata")
colnames(alldata.merged3)

### Add new column: hours
currtime = as.numeric(anytime(alldata.merged3$CURR_TIME))/3600
icustayIDtable = table(alldata.merged3$ICUSTAY_ID)
icustayIDtable_cumsum = cumsum(table(alldata.merged3$ICUSTAY_ID))
hours = unlist(sapply(icustayIDtable, function(x) rep(1:x)))
alldata.merged3$HOURS = hours
alldata.merged3$AKI = 1
alldata.merged3$AKI_PF = 0
alldata.merged3$AKI_PF[alldata.merged3$PF < 3] = 1
alldata.merged3$AKI_PF[alldata.merged3$PF < 2] = 2
alldata.merged3$AKI_PF[alldata.merged3$PF < 1] = 3
alldata.merged3$AKI_BMI = 0
alldata.merged3$AKI_BMI[alldata.merged3$BMI >= 40] = 4
alldata.merged3$AKI_BMI[alldata.merged3$BMI < 40] = 3
alldata.merged3$AKI_BMI[alldata.merged3$BMI < 25] = 2
alldata.merged3$AKI_BMI[alldata.merged3$BMI < 23] = 1
alldata.merged3$AKI_BMI[alldata.merged3$BMI < 18.5] = 0
save(alldata.merged3, file = "expanded.all.data.merged.imputed.calculated.AKI.Rdata")


length(unique(alldata.merged3$ICUSTAY_ID))

