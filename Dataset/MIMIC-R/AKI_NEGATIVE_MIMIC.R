setwd("C:/Users/congf/Desktop/AKI/AKI_KIDGO_MIMIC_2019.3.13/non_aki")
install.packages(c("anytime", "tictoc","lubridate","dplyr"))
library("anytime")
library("tictoc")
library("lubridate")

####################################################################
####################################################################
############################ MIMIC_non_AKI###############################
####################################################################
####################################################################

###########################################################
###            Expand data
###########################################################
nalldata = read.csv("non_aki_cohort/mimic_mimiciii_aki_kidgo_negative_demo2.CSV", header = F, stringsAsFactors = F)
nlabel = levels(unlist(read.table("non_aki_cohort/mimic_label.tsv.txt", sep = "\t")))
colnames(nalldata) = nlabel

nalldata$ADMI_TIME = anytime(nalldata$ADMI_TIME)
nalldata$OUT_TIME = anytime(nalldata$OUT_TIME)
nalldata$ADMI_TIME_numeric = as.numeric(nalldata$ADMI_TIME)
nalldata$OUT_TIME_numeric = as.numeric(nalldata$OUT_TIME)
nalldata = nalldata[nalldata$AGE < 150, ]
# start expanding data
expand.data = NULL
start_time <- proc.time()

allhours = 0
hourslist = NULL
for (i in 1:dim(nalldata)[1]){
  if (i %% 2500 == 0){
    message(i, "  -  ", allhours, "    Elapsed: ", round((proc.time() - start_time)[3], 2), " secs")
  }
  row = nalldata[i,]
  time.in = row$ADMI_TIME_numeric
  time.out = row$OUT_TIME_numeric
  hours = (time.out - time.in)/3600 + 1
  hourslist = c(hourslist, hours)
  allhours = allhours + hours
}
nalldata$hours = hourslist
# expand
nalldata.expanded <- nalldata[rep(row.names(nalldata), nalldata$hours), 1:25]

start_time <- proc.time()
currtime_list = NULL
allhours = 0
for (i in 1:dim(nalldata)[1]){
  if (i %% 2500 == 0){
    message(i, "  -  ", allhours, "    Elapsed: ", round((proc.time() - start_time)[3], 2), " secs")
  }
  row = nalldata[i,]
  time.in = row$ADMI_TIME_numeric
  time.out = row$OUT_TIME_numeric
  hours = (time.out - time.in)/3600
  allhours = allhours + hours
  
  numerical.times = seq(time.in, time.in + 3600*hours, by=3600)
  currtime_list[[i]] = numerical.times
}
currtime_list2 = unlist(currtime_list)#do.call(c, unlist(currtime_list, recursive=FALSE))
times = as.character(anytime(currtime_list2))

nalldata.expanded$CURR_TIME = times

save(nalldata.expanded, file = "expanded.nall.data.Rdata")

load("expanded.nall.data.Rdata")

###########################################################
###             merge non aki vit
###########################################################
library(dplyr)
library(anytime)

nhy_vit = read.csv("non_aki_vit/mimic_mimiciii_aki_kidgo_negative_vit.csv", header = F, stringsAsFactors = F)
nlabel = as.character(unlist(t(read.table("non_aki_vit/MIMIC_vit_label.tsv.txt", sep = "\t"))))
colnames(nhy_vit) = nlabel
nhy_vit$VIT_TIME = as.character(anytime(nhy_vit$VIT_TIME))

nalldata.expanded$mergekey = paste0(nalldata.expanded$ICUSTAY_ID, "_", nalldata.expanded$CURR_TIME)
nhy_vit$mergekey = paste0(nhy_vit$ICUSTAY_ID, "_", nhy_vit$VIT_TIME)
nalldata.merged = full_join(nalldata.expanded, nhy_vit,
                            by = 'mergekey')
nalldata.merged <- nalldata.merged[!duplicated(nalldata.merged$mergekey),]

###########################################################
###             merge non_aki_lab
###########################################################
nhy_lab3 = read.csv("non_aki_lab/mimic_mimiciii_aki_kidgo_negative_lab.csv", header = F, stringsAsFactors = F)
nlabel = as.character(unlist(t(read.table("non_aki_lab/mimic_LAB_LABEL.TSV.txt", sep = "\t"))))
colnames(nhy_lab3) = nlabel
nhy_lab3$LAB_TIME = as.character(anytime(nhy_lab3$LAB_TIME))
nhy_lab3$mergekey = paste0(nhy_lab3$ICUSTAY_ID, "_", nhy_lab3$LAB_TIME)

nalldata.merged = full_join(nalldata.merged, nhy_lab3,
                            by = 'mergekey')
nalldata.merged <- nalldata.merged[!duplicated(nalldata.merged$mergekey),]
###########################################################
###             merge NON INPUT_6hr
###########################################################
nhy_input_6hr = read.csv("non_aki_input/input_6hr/mimic_mimiciii_aki_kidgo_negative_input_6h.csv", header = F, stringsAsFactors = F)
nlabel = as.character(unlist(t(read.table("non_aki_input/input_6hr/input_6hr.tsv.txt", sep = "\t"))))
colnames(nhy_input_6hr) = nlabel
nhy_input_6hr$INPUT_6HR_CHARTTIME = as.character(anytime(nhy_input_6hr$INPUT_6HR_CHARTTIME))
nhy_input_6hr$mergekey = paste0(nhy_input_6hr$ICUSTAY_ID, "_", nhy_input_6hr$INPUT_6HR_CHARTTIME)

nalldata.merged = full_join(nalldata.merged, nhy_input_6hr,
                            by = 'mergekey')
nalldata.merged <- nalldata.merged[!duplicated(nalldata.merged$mergekey),]
###########################################################
###             merge non INPUT_12hr
###########################################################
nhy_input_12hr = read.csv("non_aki_input/input_12hr/mimic_mimiciii_aki_kidgo_negative_input_12h.csv", header = F, stringsAsFactors = F)
nlabel = as.character(unlist(t(read.table("non_aki_input/input_12hr/input_12hr.tsv.txt", sep = "\t"))))
colnames(nhy_input_12hr) = nlabel
nhy_input_12hr$INPUT_12HR_CHARTTIME = as.character(anytime(nhy_input_12hr$INPUT_12HR_CHARTTIME))
nhy_input_12hr$mergekey = paste0(nhy_input_12hr$ICUSTAY_ID, "_", nhy_input_12hr$INPUT_12HR_CHARTTIME)

nalldata.merged = full_join(nalldata.merged, nhy_input_12hr,
                            by = 'mergekey')
nalldata.merged <- nalldata.merged[!duplicated(nalldata.merged$mergekey),]
###########################################################
###             merge non INPUT_24hr
###########################################################
nhy_input_24hr = read.csv("non_aki_input/input_24hr/mimic_mimiciii_aki_kidgo_negative_input_24h.csv", header = F, stringsAsFactors = F)
nlabel = as.character(unlist(t(read.table("non_aki_input/input_24hr/input_24hr.tsv.txt", sep = "\t"))))
colnames(nhy_input_24hr) = nlabel
nhy_input_24hr$INPUT_24HR_CHARTTIME = as.character(anytime(nhy_input_24hr$INPUT_24HR_CHARTTIME))
nhy_input_24hr$mergekey = paste0(nhy_input_24hr$ICUSTAY_ID, "_", nhy_input_24hr$INPUT_24HR_CHARTTIME)

nalldata.merged = full_join(nalldata.merged, nhy_input_24hr,
                            by = 'mergekey')
nalldata.merged <- nalldata.merged[!duplicated(nalldata.merged$mergekey),]

###########################################################
###             merge non outPUT_6hr
###########################################################
nhy_output_6hr = read.csv("non_aki_output/output_6hr/mimic_mimiciii_aki_kidgo_negative_output_6h.csv", header = F, stringsAsFactors = F)
nlabel = as.character(unlist(t(read.table("non_aki_output/output_6hr/output_6hr.tsv.txt", sep = "\t"))))
colnames(nhy_output_6hr) = nlabel
nhy_output_6hr$OUTPUT_6HR_TIME = as.character(anytime(nhy_output_6hr$OUTPUT_6HR_TIME))
nhy_output_6hr$mergekey = paste0(nhy_output_6hr$ICUSTAY_ID, "_", nhy_output_6hr$OUTPUT_6HR_TIME)

nalldata.merged = full_join(nalldata.merged, nhy_output_6hr,
                            by = 'mergekey')
nalldata.merged <- nalldata.merged[!duplicated(nalldata.merged$mergekey),]
###########################################################
###             merge non outPUT_12hr
###########################################################
nhy_output_12hr = read.csv("non_aki_output/output_12hr/mimic_mimiciii_aki_kidgo_negative_output_12h.csv", header = F, stringsAsFactors = F)
nlabel = as.character(unlist(t(read.table("non_aki_output/output_12hr/output_12hr.tsv.txt", sep = "\t"))))
colnames(nhy_output_12hr) = nlabel
nhy_output_12hr$OUTPUT_12HR_TIME = as.character(anytime(nhy_output_12hr$OUTPUT_12HR_TIME))
nhy_output_12hr$mergekey = paste0(nhy_output_12hr$ICUSTAY_ID, "_", nhy_output_12hr$OUTPUT_12HR_TIME)

nalldata.merged = full_join(nalldata.merged, nhy_output_12hr,
                            by = 'mergekey')
nalldata.merged <- nalldata.merged[!duplicated(nalldata.merged$mergekey),]
###########################################################
###             merge non outPUT_24hr
###########################################################
nhy_output_24hr = read.csv("non_aki_output/output_24hr/mimic_mimiciii_aki_kidgo_negative_output_24h.csv", header = F, stringsAsFactors = F)
nlabel = as.character(unlist(t(read.table("non_aki_output/output_24hr/output_24hr.tsv.txt", sep = "\t"))))
colnames(nhy_output_24hr) = nlabel
nhy_output_24hr$OUTPUT_24HR_TIME = as.character(anytime(nhy_output_24hr$OUTPUT_24HR_TIME))
nhy_output_24hr$mergekey = paste0(nhy_output_24hr$ICUSTAY_ID, "_", nhy_output_24hr$OUTPUT_24HR_TIME)

nalldata.merged = full_join(nalldata.merged, nhy_output_24hr,
                            by = 'mergekey')
nalldata.merged <- nalldata.merged[!duplicated(nalldata.merged$mergekey),]
###########################################################
###             merge non ABG
###########################################################
nhy_abg = read.csv("non_aki_bg/mimic_mimiciii_aki_kidgo_negative_abg.csv", header = F, stringsAsFactors = F)
nlabel = as.character(unlist(t(read.table("non_aki_bg/mimic_ABG_LABLE.TSV.txt", sep = "\t"))))
colnames(nhy_abg) = nlabel
nhy_abg$ABG_TIME = as.character(anytime(nhy_abg$ABG_TIME))
nhy_abg$mergekey = paste0(nhy_abg$ICUSTAY_ID, "_", nhy_abg$ABG_TIME)

nalldata.merged = full_join(nalldata.merged, nhy_abg,
                            by = 'mergekey')
nalldata.merged <- nalldata.merged[!duplicated(nalldata.merged$mergekey),]
###########################################################
###             merge non mechvent_starttime
###########################################################
nhy_mechvent = read.csv("non_aki_mechvent/mechvent_starttime/mimic_mimiciii_aki_kidgo_negative_mechvent_starttime.csv", header = F, stringsAsFactors = F)
nlabel = as.character(unlist(t(read.table("non_aki_mechvent/mechvent_starttime/MIMIC_AKI_mechvent.txt", sep = "\t"))))
colnames(nhy_mechvent) = nlabel
nhy_mechvent$MECHVENT_STARTTIME = as.character(anytime(nhy_mechvent$MECHVENT_STARTTIME))

nhy_mechvent$mergekey = paste0(nhy_mechvent$ICUSTAY_ID, "_", nhy_mechvent$MECHVENT_STARTTIME)

nalldata.merged = full_join(nalldata.merged, nhy_mechvent,
                            by = 'mergekey')
nalldata.merged <- nalldata.merged[!duplicated(nalldata.merged$mergekey),]

###########################################################
###             merge non mechvent_endtime
###########################################################
nhy_mechvent2 = read.csv("non_aki_mechvent/mechvent_endtime/mimic_mimiciii_aki_kidgo_negative_mechvent_endtime.csv", header = F, stringsAsFactors = F)
nlabel = as.character(unlist(t(read.table("non_aki_mechvent/mechvent_endtime/MIMIC_AKI_mechvent.txt", sep = "\t"))))
colnames(nhy_mechvent2) = nlabel
nhy_mechvent2$MECHVENT_ENDTIME = as.character(anytime(nhy_mechvent2$MECHVENT_ENDTIME))

nhy_mechvent2$mergekey = paste0(nhy_mechvent2$ICUSTAY_ID, "_", nhy_mechvent2$MECHVENT_ENDTIME)

nalldata.merged = full_join(nalldata.merged, nhy_mechvent2,
                            by = 'mergekey')
nalldata.merged <- nalldata.merged[!duplicated(nalldata.merged$mergekey),]
###########################################################
###             merge non vaso_STARTTIME
###########################################################
nhy_VASO = read.csv("non_aki_vaso/vaso_start/mimic_mimiciii_aki_kidgo_negative_vasopressor_starttime.csv", header = F, stringsAsFactors = F)
nlabel = as.character(unlist(t(read.table("non_aki_vaso/vaso_start/MIMIC_vaso_label.tsv.txt", sep = "\t"))))
colnames(nhy_VASO) = nlabel
nhy_VASO$VASO_STARTTIME = as.character(anytime(nhy_VASO$VASO_STARTTIME))
nhy_VASO$mergekey = paste0(nhy_VASO$ICUSTAY_ID, "_", nhy_VASO$VASO_STARTTIME)

nalldata.merged = full_join(nalldata.merged, nhy_VASO,
                            by = 'mergekey')
nalldata.merged <- nalldata.merged[!duplicated(nalldata.merged$mergekey),]
###########################################################
###             merge non vaso_ENDTIME
###########################################################
nhy_VASO2 = read.csv("non_aki_vaso/vaso_end/mimic_mimiciii_aki_kidgo_negative_vasopressor_endtime.csv", header = F, stringsAsFactors = F)
nlabel = as.character(unlist(t(read.table("non_aki_vaso/vaso_end/MIMIC_vaso_label.tsv.txt", sep = "\t"))))
colnames(nhy_VASO2) = nlabel
nhy_VASO2$VASO_ENDTIME = as.character(anytime(nhy_VASO2$VASO_ENDTIME))
nhy_VASO2$mergekey = paste0(nhy_VASO2$ICUSTAY_ID, "_", nhy_VASO2$VASO_ENDTIME)

nalldata.merged = full_join(nalldata.merged, nhy_VASO2,
                            by = 'mergekey')
nalldata.merged <- nalldata.merged[!duplicated(nalldata.merged$mergekey),]


colnames(nalldata.merged)
###########################################################
###             Remove redundant data
###########################################################

nalldata.merged2 = nalldata.merged[, c("ICUSTAY_ID.x","ICU_CLASS","ETHNICITY",
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
colnames(nalldata.merged2)[1] = c("ICUSTAY_ID")
colnames(nalldata.merged2)
nalldata.merged2 = nalldata.merged2[!is.na(nalldata.merged2$ICUSTAY_ID), ]
###########################################################
###       Sort data based on ICUSTAY_ID and TIME
###########################################################
library(anytime)
index2sort = paste0(nalldata.merged2$ICUSTAY_ID, "_", nalldata.merged2$CURR_TIME)
order = sort.int(index2sort, index.return = T)$ix
nalldata.merged2 = nalldata.merged2[order,]
save(nalldata.merged2, file = "expanded.nall.data.merged.Rdata")
load("expanded.nall.data.merged.Rdata")

###########################################################
###             Fill missing values
###########################################################
##setwd("/home/zhihuan/Documents/Cong_Feng/20180908_Hypoxemia/Hypoxemia - LSTM/PFéçåµ?/is_vent")
library(zoo)
library(anytime)
load("expanded.nall.data.merged.Rdata")
colnames(nalldata.merged2)
nalldata.merged2[is.na(nalldata.merged2$INPUT_6HR), "INPUT_6HR"] = 0
nalldata.merged2[is.na(nalldata.merged2$INPUT_12HR), "INPUT_12HR"] = 0
nalldata.merged2[is.na(nalldata.merged2$INPUT_24HR), "INPUT_24HR"] = 0
nalldata.merged2[is.na(nalldata.merged2$OUTPUT_6HR), "OUTPUT_6HR"] = 0
nalldata.merged2[is.na(nalldata.merged2$OUTPUT_12HR), "OUTPUT_12HR"] = 0
nalldata.merged2[is.na(nalldata.merged2$OUTPUT_24HR), "OUTPUT_24HR"] = 0
nalldata.merged2$INPUT_MINUS_OUTPUT_6HR = nalldata.merged2$INPUT_6HR - nalldata.merged2$OUTPUT_6HR
nalldata.merged2$INPUT_MINUS_OUTPUT_12HR = nalldata.merged2$INPUT_12HR - nalldata.merged2$OUTPUT_12HR
nalldata.merged2$INPUT_MINUS_OUTPUT_24HR = nalldata.merged2$INPUT_24HR - nalldata.merged2$OUTPUT_24HR


nalldata.merged2$HEIGHT[is.na(nalldata.merged2$HEIGHT)] = mean(nalldata.merged2$HEIGHT, na.rm=TRUE)
nalldata.merged2$WEIGHT[is.na(nalldata.merged2$WEIGHT)] = mean(nalldata.merged2$WEIGHT, na.rm=TRUE)
nalldata.merged2$BMI = nalldata.merged2$WEIGHT / (nalldata.merged2$HEIGHT/100)^2


# is.na(alldata.merged2$FIO2)
icustayIDlist = unique(nalldata.merged2$ICUSTAY_ID)
### Expanding original alldata
newrows1 = data.frame(matrix(-Inf, nrow=(length(icustayIDlist)-1), dim(nalldata.merged2)[2]))
newrows2 = data.frame(matrix(-Inf, nrow=(length(icustayIDlist)-1), dim(nalldata.merged2)[2]))
colnames(newrows1) = colnames(nalldata.merged2)
colnames(newrows2) = colnames(nalldata.merged2)
newrows1$ICUSTAY_ID = unique(nalldata.merged2$ICUSTAY_ID)[1:length(unique(nalldata.merged2$ICUSTAY_ID))-1]
newrows2$ICUSTAY_ID = unique(nalldata.merged2$ICUSTAY_ID)[2:length(unique(nalldata.merged2$ICUSTAY_ID))]
newrows1$CURR_TIME = "9999-12-31 23:59:59"
newrows2$CURR_TIME = "0001-01-01 00:00:00"
newrows1 = data.frame(newrows1)
newrows2 = data.frame(newrows2)



nalldata.merged2$NEU_PER = as.numeric(nalldata.merged2$NEU_PER)
nalldata.merged2$LIP = as.numeric(nalldata.merged2$LIP)
nalldata.merged2$BNP = as.numeric(nalldata.merged2$BNP)
nalldata.merged2$DD = as.numeric(nalldata.merged2$DD)
nalldata.merged2$FIB = as.numeric(nalldata.merged2$FIB)
nalldata.merged2$MECHVENT_STARTTIME = as.numeric(anytime(nalldata.merged2$MECHVENT_STARTTIME))
nalldata.merged2$MECHVENT_ENDTIME = as.numeric(anytime(nalldata.merged2$MECHVENT_ENDTIME))
nalldata.merged2$VASO_STARTTIME = as.numeric(anytime(nalldata.merged2$VASO_STARTTIME))
nalldata.merged2$VASO_ENDTIME = as.numeric(anytime(nalldata.merged2$VASO_ENDTIME))
nalldata.merged2$ADMI_TIME = as.numeric(anytime(nalldata.merged2$ADMI_TIME))
nalldata.merged2$OUT_TIME = as.numeric(anytime(nalldata.merged2$OUT_TIME))


# for (i in 1:dim(nalldata.merged2)[2]){
#   message(colnames(nalldata.merged2)[i])
#   print(typeof(nalldata.merged2[1,i]))
#   print(typeof(tobind[1,i]))
#   tmp = rbind(nalldata.merged2[1:10,1:i], tobind[1:10,1:i])
# }

tobind = rbind(newrows1, newrows2)
nalldata.merged3 = rbind(nalldata.merged2, tobind)
index2sort = paste0(nalldata.merged3$ICUSTAY_ID, "_", nalldata.merged3$CURR_TIME)
order = sort.int(index2sort, index.return = T)$ix
nalldata.merged3 = nalldata.merged3[order,]
save(nalldata.merged3, file = "expanded.nall.inputoutput.missingfill1.Rdata")
colnames(nalldata.merged3)
############## ????????????
ncolumns2impute = colnames(nalldata.merged3)[c(25:67,75:76)]
nlinear.impute = na.approx(zoo(nalldata.merged3[, ncolumns2impute]), na.rm = F)
nlinear.impute[is.nan(nlinear.impute)] = NA
nalldata.merged3[, ncolumns2impute] = as.data.frame(nlinear.impute)
nalldata.merged3 = nalldata.merged3[nalldata.merged3$ICU_CLASS != -Inf,]
nalldata.merged3[nalldata.merged3 == -Inf] = NA

nalldata.merged3 = rbind(nalldata.merged3, rbind(newrows1, newrows2))
index2sort = paste0(nalldata.merged3$ICUSTAY_ID, "_", nalldata.merged3$CURR_TIME)
order = sort.int(index2sort, index.return = T)$ix
nalldata.merged3 = nalldata.merged3[order,]
save(nalldata.merged3, file = "expanded.nall.labvitbglinear.missingfill.Rdata")

#######################ABG LAB VIT FIO2?????????
ncolumns4impute = colnames(nalldata.merged3)[c(25:67,74:76)]
nfront.impute = data.frame(na.locf(zoo(nalldata.merged3[, ncolumns4impute]), na.rm = F))
nalldata.merged3[, ncolumns4impute] = nfront.impute
nalldata.merged3 = nalldata.merged3[nalldata.merged3$ICU_CLASS != -Inf,]
nalldata.merged3[nalldata.merged3 == -Inf] = NA

nalldata.merged3 = rbind(nalldata.merged3, rbind(newrows1, newrows2))
index2sort = paste0(nalldata.merged3$ICUSTAY_ID, "_", nalldata.merged3$CURR_TIME)
order = sort.int(index2sort, index.return = T)$ix
nalldata.merged3 = nalldata.merged3[order,]
save(nalldata.merged3, file = "expanded.nall.ABGLABVITFIO2frontimpute.missingfill.Rdata")



#######################ABG LAB VIT FIO2 ?????????
ncolumns3impute = colnames(nalldata.merged3)[c(25:67,74:76)]
nback.impute = data.frame(na.locf(zoo(nalldata.merged3[, ncolumns3impute]), na.rm = F, fromLast = T))
nalldata.merged3[, ncolumns3impute] = nback.impute
nalldata.merged3 = nalldata.merged3[nalldata.merged3$ICU_CLASS != -Inf,]
nalldata.merged3[nalldata.merged3 == -Inf] = NA

nalldata.merged3 = rbind(nalldata.merged3, rbind(newrows1, newrows2))
index2sort = paste0(nalldata.merged3$ICUSTAY_ID, "_", nalldata.merged3$CURR_TIME)
order = sort.int(index2sort, index.return = T)$ix
nalldata.merged3 = nalldata.merged3[order,]
save(nalldata.merged3, file = "expanded.nall.ABGLABVITFIO2backimpute.missingfill.Rdata")
load("expanded.nall.ABGLABVITFIO2backimpute.missingfill.Rdata")
install.packages(c("zoo"))
library(zoo)

# mechvent time
nalldata.merged3$MECH = NA
nalldata.merged3$MECH = nalldata.merged3$MECHVENT_STARTTIME
nalldata.merged3$MECH[!is.na(nalldata.merged3$MECHVENT_ENDTIME)] = nalldata.merged3$MECHVENT_ENDTIME[!is.na(nalldata.merged3$MECHVENT_ENDTIME)]
front.impute = data.frame(na.locf(zoo(nalldata.merged3$MECH), na.rm = F))
back.impute = data.frame(na.locf(zoo(nalldata.merged3$MECH), na.rm = F, fromLast = T))
mask = front.impute * back.impute
mask[mask == Inf] = 0
mask[mask == -Inf] = 0
mask[mask != 0] = 1
nalldata.merged3$MECH = mask


# vaso time
nalldata.merged3$VASO = NA
nalldata.merged3$VASO = nalldata.merged3$VASO_STARTTIME
nalldata.merged3$VASO[!is.na(nalldata.merged3$VASO_ENDTIME)] = nalldata.merged3$VASO_ENDTIME[!is.na(nalldata.merged3$VASO_ENDTIME)]
front.impute = data.frame(na.locf(zoo(nalldata.merged3$VASO), na.rm = F))
back.impute = data.frame(na.locf(zoo(nalldata.merged3$VASO), na.rm = F, fromLast = T))
mask = front.impute * back.impute
mask[mask == Inf] = 0
mask[mask == -Inf] = 0
mask[mask != 0] = 1
nalldata.merged3$VASO = mask

nalldata.merged3 = nalldata.merged3[nalldata.merged3$ICU_CLASS != -Inf,]
nalldata.merged3[nalldata.merged3 == -Inf] = NA




save(nalldata.merged3, file = "expanded.nall.vasomechventimpute.missingfill.Rdata")
load("expanded.nall.vasomechventimpute.missingfill.Rdata")
colnames(nalldata.merged3)

### Add new column: hours
ncurrtime = as.numeric(anytime(nalldata.merged3$CURR_TIME))/3600
nicustayIDtable = table(nalldata.merged3$ICUSTAY_ID)
nicustayIDtable_cumsum = cumsum(table(nalldata.merged3$ICUSTAY_ID))
nhours = unlist(sapply(nicustayIDtable, function(x) rep(1:x)))
nalldata.merged3$HOURS = nhours
nalldata.merged3$PF = nalldata.merged3$PO2/nalldata.merged3$FIO2
nalldata.merged3$AKI = 0
nalldata.merged3$AKI_PF = 0
nalldata.merged3$AKI_PF[nalldata.merged3$P_F_ratio < 3] = 1
nalldata.merged3$AKI_PF[nalldata.merged3$P_F_ratio < 2] = 2
nalldata.merged3$AKI_PF[nalldata.merged3$P_F_ratio < 1] = 3
nalldata.merged3$AKI_BMI = 0
nalldata.merged3$AKI_BMI[nalldata.merged3$BMI >= 40] = 4
nalldata.merged3$AKI_BMI[nalldata.merged3$BMI < 40] = 3
nalldata.merged3$AKI_BMI[nalldata.merged3$BMI < 25] = 2
nalldata.merged3$AKI_BMI[nalldata.merged3$BMI < 23] = 1
nalldata.merged3$AKI_BMI[nalldata.merged3$BMI < 18.5] = 0
length(unique(nalldata.merged3$ICUSTAY_ID))

save(nalldata.merged3, file = "expanded.nall.data.merged.imputed.calculated.AKI.Rdata")

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
setwd("C:/Users/congf/Desktop/AKI/AKI_KIDGO_MIMIC_2019.3.13/non_aki")
load("expanded.nall.data.merged.imputed.calculated.AKI.Rdata")

setwd("C:/Users/congf/Desktop/AKI/AKI_KIDGO_MIMIC_2019.3.13/aki")
load("expanded.all.data.merged.imputed.calculated.AKI.Rdata")


alldata.imputed.aki = alldata.merged3
alldata.imputed.naki = nalldata.merged3
colnames(alldata.imputed.aki)
colnames(alldata.imputed.naki)
# If rbind fails, it means one of the attributes of the data frame was a list
alldata.imputed.aki$MECH = unlist(alldata.imputed.aki$MECH)
alldata.imputed.aki$MECH[is.na(alldata.imputed.aki$MECH)] = 0
alldata.imputed.aki$VASO = unlist(alldata.imputed.aki$VASO)
alldata.imputed.aki$VASO[is.na(alldata.imputed.aki$VASO)] = 0
alldata.imputed.aki$GENDER[is.na(alldata.imputed.aki$GENDER)] = 1
alldata.imputed.aki$ETHNICITY[is.na(alldata.imputed.aki$ETHNICITY)] = 1

alldata.imputed.naki$MECH = unlist(alldata.imputed.naki$MECH)
alldata.imputed.naki$MECH[is.na(alldata.imputed.naki$MECH)] = 0
alldata.imputed.naki$VASO = unlist(alldata.imputed.naki$VASO)
alldata.imputed.naki$VASO[is.na(alldata.imputed.naki$VASO)] = 0
alldata.imputed.naki$GENDER[is.na(alldata.imputed.naki$GENDER)] = 1
alldata.imputed.naki$ETHNICITY[is.na(alldata.imputed.naki$ETHNICITY)] = 1

alldata.imputed = rbind(alldata.imputed.aki, alldata.imputed.naki)
save(alldata.imputed, file = "expanded.all.data.merged.imputed.calculated.merged.Rdata")
load("expanded.all.data.merged.imputed.calculated.merged.Rdata")

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
save(alldata.final, file = "expanded.Remove.some.colnames.missingdata.aki.mimic.Rdata")
###########################################################
### Drop some columns, make complete table
###########################################################
load("expanded.Remove.some.colnames.missingdata.aki.mimic.Rdata")


colnames(alldata.final)

droplist = NULL
# droplist = c("FIO2","PO2","PCO2","SUBJECT_ID","HADM_ID","LOS")
# droplist = c("HEIGHT","WEIGHT","BMI")
remains = NULL
for (i in 1:dim(alldata.final)[2]){
  message(colnames(alldata.final)[i], "\t", sum(is.na(alldata.final[,i])) )
  if (sum(is.na(alldata.final[,i])) > 600000){
    droplist = c(droplist, colnames(alldata.final)[i])
  }
  else {
    remains = c(remains, colnames(alldata.final)[i])
  }
}

droplist = c("PH","LAC","P_ION","APTT","INR","MG_ION","ALB")

# remains = c(remains, "FIO2", "PO2", "PCO2", "P_F_ratio", "IS_VENT", "IS_VENT_P_F_ratio_target")
todrop = which(names(alldata.final) %in% droplist)
if (length(todrop) == 0){
  alldata.imputed.shrinked = alldata.final[complete.cases(alldata.final), remains]
} else {
  alldata.imputed.shrinked = alldata.final[complete.cases(alldata.final[,-todrop]), -todrop]
}

sum(is.na(alldata.imputed.shrinked)) == 0
length(unique(alldata.imputed.shrinked$ICUSTAY_ID))

colnames(alldata.imputed.shrinked)


save(alldata.imputed.shrinked, file="AKI.MIMIC.FINAL.Rdata")
load("AKI.MIMIC.FINAL.Rdata")
# load("expanded.all.data.merged.imputed.calculated.merged.shrinked.Rdata")
write.csv(alldata.imputed.shrinked, file = "AKI.MIMIC.FINAL.csv", row.names = F)

































alldata.ed.mimic.aki = alldata.imputed.shrinked[,c("ICUSTAY_ID","ICU_CLASS","ETHNICITY",
                                                   "AGE","GENDER","HEIGHT","WEIGHT",
                                                   "BMI","DIAS_BP","HR","SYS_BP","MEAN_BP",
                                                   "RESPRATE","Temprate","SPO2",
                                                   "MECHVENT_TIME","VASO_TIME","AKI","HOURS")]

save(alldata.ed.mimic.aki, file="alldata.ed.mimic.aki")

# write.csv(alldata.imputed.shrinked, file = "expanded.all.data.merged.imputed.calculated.shrinked.csv", row.names = F)
write.csv(alldata.ed.mimic.aki, file = "alldata.ed.mimic.aki.csv", row.names = F)


