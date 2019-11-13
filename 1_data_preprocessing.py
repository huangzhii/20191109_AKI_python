#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 20:43:30 2019

@author: zhihuan
"""


import re, copy
import time
import pandas as pd
import numpy as np
import sys, os, platform
import pickle
import itertools
import argparse
import urllib
import json
from tqdm import tqdm
tqdm.pandas()


print("Currently working on " + list(platform.uname())[1] + " Machine")
if 'Y710' in list(platform.uname())[1]:
    workdir = '/media/zhihuan/Drive3/20191109_AKI_python/'
else:
    workdir = ''
    
def expand_data(x, group):
    start_row = copy.deepcopy(x)
    start_row['CURR_TIME'] = pd.to_datetime(x['ADMI_TIME'])
    end_row = copy.deepcopy(x)
    end_row['CURR_TIME'] = pd.to_datetime(x['OUT_TIME'])
    
    before_expand = pd.concat((start_row, end_row))
    before_expand.index = before_expand['CURR_TIME']
    after_expand = before_expand.resample('1H').pad()
    after_expand['CURR_TIME'] = after_expand.index
    after_expand['HOURS'] = range(1, len(after_expand) + 1)
    after_expand['AKI'] = 0
    if group == 'aki': after_expand['AKI'][-1] = 1
    return after_expand


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', default='MIMIC', type=str)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    
    datadir = {}
    datadir['MIMIC'] = workdir + 'Dataset/' + 'MIMIC11.11/'
    datadir['EICU'] = workdir + 'Dataset/' + 'EICU11.11/'
    
    datagroup = ['bg', 'cohort', 'input', 'lab', 'mechvent', 'output', 'statistic', 'vaso', 'vit']
    data_expand = {}
    
    
    for group in ['aki', 'non_aki']:
        if group == 'aki': positivity = 'positive'
        if group == 'non_aki': positivity = 'negative'
    # =============================================================================
    #      Cohort
    # =============================================================================
        colnames = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_cohort/mimic_label.tsv.txt', sep = '\t').columns.values.astype(str)
        data = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_cohort/mimic_mimiciii_aki_kidgo_'+positivity+'_demo3.csv', header=None)
        data.columns = colnames
        data = data.loc[data.AGE < 150,:]
        if len(data.ICUSTAY_ID.unique()) < len(data): print('Error: ICUSTAY_ID is not unique')
        # expand data
        data_expand[group] = data.groupby(by = 'ICUSTAY_ID').progress_apply(lambda x: expand_data(x, group))
        data_expand[group].reset_index(drop = True, inplace = True)
        print('length of '+ group + ' data (expanded):', len(data_expand[group]))
    # =============================================================================
    #      Vit (vital)
    # =============================================================================
        colnames = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_vit/MIMIC_vit_label.tsv.txt', sep = '\t').columns.values.astype(str)
        data = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_vit/mimic_mimiciii_aki_kidgo_'+positivity+'_vit.csv', header=None)
        data.columns = colnames
        data['CURR_TIME'] = pd.to_datetime(data['VIT_TIME'])
        data.drop_duplicates(subset = ['ICUSTAY_ID', 'CURR_TIME'], keep = 'first', inplace = True)
        data_expand[group] = data_expand[group].merge(data, on = ['ICUSTAY_ID', 'CURR_TIME'], how = 'left')
        
    # =============================================================================
    #      LAB
    # =============================================================================
        colnames = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_lab/mimic_LAB_LABEL.TSV.txt', sep = '\t').columns.values.astype(str)
        data = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_lab/mimic_mimiciii_aki_kidgo_'+positivity+'_lab.csv', header=None)
        data.columns = colnames
        data['CURR_TIME'] = pd.to_datetime(data['LAB_TIME'])
        data.drop_duplicates(subset = ['ICUSTAY_ID', 'CURR_TIME'], keep = 'first', inplace = True)
        data_expand[group] = data_expand[group].merge(data, on = ['ICUSTAY_ID', 'CURR_TIME'], how = 'left')
        
    # =============================================================================
    #      Input 6, 12, 24 hr
    # =============================================================================
        hours = [6, 12, 24]
        for h in hours:
            colnames = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_input/input_' + str(h) +'hr/input_' + str(h) +'hr.tsv.txt', sep = '\t').columns.values.astype(str)
            data = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_input/input_' + str(h) +'hr/mimic_mimiciii_aki_kidgo_'+positivity+'_input_' + str(h) +'h.csv', header=None)
            data.columns = colnames
            data['CURR_TIME'] = pd.to_datetime(data['INPUT_' + str(h) +'HR_CHARTTIME'])
            data.drop_duplicates(subset = ['ICUSTAY_ID', 'CURR_TIME'], keep = 'first', inplace = True)
            data_expand[group] = pd.merge(data_expand[group], data, on = ['ICUSTAY_ID', 'CURR_TIME'], how = 'left', validate = 'm:1')
        
    # =============================================================================
    #      Output 6, 12, 24 hr
    # =============================================================================
        hours = [6, 12, 24]
        for h in hours:
            colnames = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_output/output_' + str(h) +'hr/output_' + str(h) +'hr.tsv.txt', sep = '\t').columns.values.astype(str)
            data = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_output/output_' + str(h) +'hr/mimic_mimiciii_aki_kidgo_'+positivity+'_output_' + str(h) +'h.csv', header=None)
            data.columns = colnames
            data['CURR_TIME'] = pd.to_datetime(data['OUTPUT_' + str(h) +'HR_TIME'])
            data.drop_duplicates(subset = ['ICUSTAY_ID', 'CURR_TIME'], keep = 'first', inplace = True)
            data_expand[group] = pd.merge(data_expand[group], data, on = ['ICUSTAY_ID', 'CURR_TIME'], how = 'left', validate = 'm:1')
        
    # =============================================================================
    #      ABG
    # =============================================================================
        colnames = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_bg/mimic_ABG_LABLE.TSV.txt', sep = '\t').columns.values.astype(str)
        data = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_bg/mimic_mimiciii_aki_kidgo_'+positivity+'_abg.csv', header=None)
        data.columns = colnames
        data['CURR_TIME'] = pd.to_datetime(data['ABG_TIME'])
        data.drop_duplicates(subset = ['ICUSTAY_ID', 'CURR_TIME'], keep = 'first', inplace = True)
        data_expand[group] = data_expand[group].merge(data, on = ['ICUSTAY_ID', 'CURR_TIME'], how = 'left')
    
    # =============================================================================
    #      mechvent_starttime mechvent_endtime
    # =============================================================================
        for t in ['starttime', 'endtime']:
            colnames = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_mechvent/mechvent_'+ t +'/MIMIC_AKI_mechvent.txt', sep = '\t').columns.values.astype(str)
            data = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_mechvent/mechvent_'+ t +'/mimic_mimiciii_aki_kidgo_'+positivity+'_mechvent_'+ t +'.csv', header=None)
            data.columns = colnames
            data['CURR_TIME'] = pd.to_datetime(data['MECHVENT_' + t.upper()])
            data.drop_duplicates(subset = ['ICUSTAY_ID', 'CURR_TIME'], keep = 'first', inplace = True)
            data_expand[group] = data_expand[group].merge(data, on = ['ICUSTAY_ID', 'CURR_TIME'], how = 'left')
            
    # =============================================================================
    #      vaso_starttime vaso_endtime
    # =============================================================================
        for t in ['start', 'end']:
            colnames = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_vaso/vaso_'+ t +'/MIMIC_vaso_label.tsv.txt', sep = '\t').columns.values.astype(str)
            data = pd.read_csv(datadir['MIMIC'] + group +'/' + group + '_vaso/vaso_'+ t +'/mimic_mimiciii_aki_kidgo_'+positivity+'_vasopressor_'+t+'time.csv', header=None)
            data.columns = colnames
            data['CURR_TIME'] = pd.to_datetime(data['VASO_' + t.upper() + 'TIME'])
            data.drop_duplicates(subset = ['ICUSTAY_ID', 'CURR_TIME'], keep = 'first', inplace = True)
            data_expand[group] = data_expand[group].merge(data, on = ['ICUSTAY_ID', 'CURR_TIME'], how = 'left')
            
            
            
    data_expand_all = pd.concat((data_expand['aki'], data_expand['non_aki']))

#    colnames_we_need = ["ICUSTAY_ID","ICU_CLASS","ETHNICITY","AGE","GENDER","LOS",\
#                        "HEIGHT","WEIGHT","BMI","ISOFA","SEP","CAR","RES","OD",\
#                        "CKD","DIA","CHF","CLD","CPD","HYP","ADMI_TIME","OUT_TIME",\
#                        "CURR_TIME","DIAS_BP","HR","SYS_BP","MEAN_BP","RR","TEM","SPO2",\
#                        "PH","CA_ION","HGB","WBC","RBC","HCT","PLT","RDW",\
#                        "CRP","HCO3","ALT","AST","ALB","TBB","TNT","CK","CKMB","CR",\
#                        "UN","AMI","LIP","BNP","CL_ION","GLU","K_ION","NA_ION",\
#                        "APTT","PT","INR","DD","FIB","LAC","AG","P_ION","MG_ION",\
#                        "INPUT_6HR","INPUT_12HR","INPUT_24HR","OUTPUT_6HR","OUTPUT_12HR",\
#                        "OUTPUT_24HR","FIO2","PCO2","PO2","MECHVENT_STARTTIME",\
#                        "MECHVENT_ENDTIME","VASO_STARTTIME","VASO_ENDTIME"] # ['NEU', 'NEU_PER'] not in index
    
    demotraphic_information = ["ICUSTAY_ID","AGE","GENDER","ETHNICITY","HEIGHT","WEIGHT","BMI"]
    vital_signs = ["TEM","HR","RR","DIAS_BP","SYS_BP","MEAN_BP","SPO2"]
    laboratory_values = ["NA_ION","CA_ION","K_ION","CL_ION","AG","GLU","HCO3","WBC","RBC",\
                         "HGB","HCT","PLT","PF","ALT","AST","TBB","UN","CR"] #CR is SCr
    fluid_balance = ["INPUT_6HR","INPUT_12HR","INPUT_24HR","OUTPUT_6HR","OUTPUT_12HR","OUTPUT_24HR"]
    addtional_respiratory_and_hemodynamic_support = ["VASO", "MECH"]
    primary_diagnosis = ["SEP", "CAR", "NEU", "RES", "OD"]
    comorbidities = ["HYP","DIA","CHF", "CPD", "CKD", "CLD"]
    others = ["HOURS", "ISOFA", "AKI"]
    
    colnames_we_need = demotraphic_information + vital_signs + laboratory_values + fluid_balance + \
                        addtional_respiratory_and_hemodynamic_support + primary_diagnosis + comorbidities + others
    
    data_expand = data_expand[colnames_we_need]
