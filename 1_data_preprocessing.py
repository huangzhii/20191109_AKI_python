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
elif 'Zhi-G7-7790' in list(platform.uname())[1]:
    workdir = '/media/zhihuan/DATA/20191109_AKI_python/'
elif 'DESKTOP-05QACO1' in list(platform.uname())[1]:
    workdir = 'D:/20191109_AKI_python/'
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
#    if group == 'aki': after_expand['AKI'][-1] = 1
    if group == 'aki': after_expand['AKI'] = 1 # should all be 1
    return after_expand

def merge_mech_vaso(x, data, feature_name, cohort):
    keys = list(data.keys())
    icuid = x['ICUSTAY_ID'].values[0]
    
    # Start time
    startkey, endkey = [k for k in keys if 'start' in k][0], [k for k in keys if 'end' in k][0]
    data_starts = data[startkey].loc[data[startkey].ICUSTAY_ID == icuid,:]
    data_ends = data[endkey].loc[data[endkey].ICUSTAY_ID == icuid,:]
    if len(data_starts) < 1 and len(data_ends) < 1: return
    if len(data_starts) != len(data_ends) < 1:
        print('Error: start time and end time inconsistent')
        
    data_starts.reset_index(drop = True, inplace = True)
    data_ends.reset_index(drop = True, inplace = True)
    x[feature_name] = 0
    for i in range(len(data_starts)):
        if cohort == 'MIMIC':
            starttime = pd.to_datetime(data_starts.iloc[i, 1])
            endtime = pd.to_datetime(data_ends.iloc[i, 1])
        if cohort == 'EICU':
            starttime = pd.to_datetime(data_starts.iloc[i, 1], unit = 'h')
            endtime = pd.to_datetime(data_ends.iloc[i, 1], unit = 'h')
        x.loc[x.CURR_TIME.between(starttime, endtime), feature_name] = 1
    return x


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', default='MIMIC', type=str)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    
    datadir = {}
    datadir['MIMIC'] = workdir + 'Dataset/AKI without overt AKI patients/' + 'MIMIC-ICM-11.14/'
    datadir['EICU'] = workdir + 'Dataset/AKI without overt AKI patients/' + 'EICU-ICM-11.14/'
    
    for cohort in datadir.keys():
        print('\nCurrently processing cohort', cohort)
        datagroup = ['bg', 'cohort', 'input', 'lab', 'mechvent', 'output', 'statistic', 'vaso', 'vit']
        data_expand = {}
        
        
        for group in ['aki', 'non_aki']:
            if group == 'aki': positivity = 'positive'
            if group == 'non_aki': positivity = 'negative'
        # =============================================================================
        #      Cohort
        # =============================================================================
            if cohort == 'MIMIC':
                colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_cohort/mimic_label.tsv.txt', sep = '\t').columns.values.astype(str)
                data = pd.read_csv(datadir[cohort] + group +'/' + group + '_cohort/mimic_mimiciii_aki_kidgo_'+positivity+'_demo3.csv', header=None)
            elif cohort == 'EICU':
                colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_cohort/eicu_label.tsv.txt', sep = '\t').columns.values.astype(str)
                data = pd.read_csv(datadir[cohort] + group +'/' + group + '_cohort/eicu_eicu_crd_aki_'+positivity+'_demo2.csv', header=None)
            data.columns = colnames
            data = data.loc[data.AGE < 150,:]
            
            if cohort == 'EICU':
                data['ADMI_TIME'] = pd.to_datetime(data['ADMI_TIME'], unit = 'h')
                data['OUT_TIME'] = pd.to_datetime(data['OUT_TIME'], unit = 'h')
                
            if len(data.ICUSTAY_ID.unique()) < len(data): print('Error: ICUSTAY_ID is not unique')
            # expand data
            data_expand[group] = data.groupby(by = 'ICUSTAY_ID').progress_apply(lambda x: expand_data(x, group))
            data_expand[group].reset_index(drop = True, inplace = True)
            print('length of '+ group + ' data (expanded):', len(data_expand[group]))
        # =============================================================================
        #      Vit (vital)
        # =============================================================================
            if cohort == 'MIMIC':
                colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_vit/MIMIC_vit_label.tsv.txt', sep = '\t').columns.values.astype(str)
                data = pd.read_csv(datadir[cohort] + group +'/' + group + '_vit/mimic_mimiciii_aki_kidgo_'+positivity+'_vit.csv', header=None)
                data.columns = colnames
                data['CURR_TIME'] = pd.to_datetime(data['VIT_TIME'])
            if cohort == 'EICU':
                colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_vit/EICU_vit_label.tsv.txt', sep = '\t').columns.values.astype(str)
                data = pd.read_csv(datadir[cohort] + group +'/' + group + '_vit/eicu_eicu_crd_aki_'+positivity+'_vit2.csv', header=None)
                data.columns = colnames
                data['CURR_TIME'] = pd.to_datetime(data['VIT_TIME'], unit = 'h')
            data.drop_duplicates(subset = ['ICUSTAY_ID', 'CURR_TIME'], keep = 'first', inplace = True)
            data_expand[group] = data_expand[group].merge(data, on = ['ICUSTAY_ID', 'CURR_TIME'], how = 'left')
            
        # =============================================================================
        #      LAB
        # =============================================================================
            if cohort == 'MIMIC':
                colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_lab/mimic_LAB_LABEL.TSV.txt', sep = '\t').columns.values.astype(str)
                data = pd.read_csv(datadir[cohort] + group +'/' + group + '_lab/mimic_mimiciii_aki_kidgo_'+positivity+'_lab.csv', header=None)
                data.columns = colnames
                data['CURR_TIME'] = pd.to_datetime(data['LAB_TIME'])
            if cohort == 'EICU':
                colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_lab/EICU_LAB_LABEL.TSV.txt', sep = '\t').columns.values.astype(str)
                data = pd.read_csv(datadir[cohort] + group +'/' + group + '_lab/eicu_eicu_crd_aki_'+positivity+'_lab1.csv', header=None)
                data.columns = colnames
                data['CURR_TIME'] = pd.to_datetime(data['LAB_TIME'], unit = 'h')
            data.drop_duplicates(subset = ['ICUSTAY_ID', 'CURR_TIME'], keep = 'first', inplace = True)
            data_expand[group] = data_expand[group].merge(data, on = ['ICUSTAY_ID', 'CURR_TIME'], how = 'left')
            
        # =============================================================================
        #      Input 6, 12, 24 hr
        # =============================================================================
            hours = [6, 12, 24]
            for h in hours:
                if cohort == 'MIMIC':
                    colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_input/input_' + str(h) +'hr/input_' + str(h) +'hr.tsv.txt', sep = '\t').columns.values.astype(str)
                    data = pd.read_csv(datadir[cohort] + group +'/' + group + '_input/input_' + str(h) +'hr/mimic_mimiciii_aki_kidgo_'+positivity+'_input_' + str(h) +'h.csv', header=None)
                    data.columns = colnames
                    data['CURR_TIME'] = pd.to_datetime(data['INPUT_' + str(h) +'HR_CHARTTIME'])
                if cohort == 'EICU':
                    colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_input/' + str(h) +'hr/EICU_inputPUT' + str(h) +'_label.tsv.txt', sep = '\t').columns.values.astype(str)
                    data = pd.read_csv(datadir[cohort] + group +'/' + group + '_input/' + str(h) +'hr/eicu_eicu_crd_aki_'+positivity+'_input_' + str(h) +'hr.csv', header=None)
                    data.columns = colnames
                    data['CURR_TIME'] = pd.to_datetime(data['INPUT_' + str(h) +'HR_TIME'], unit = 'h')
                data.drop_duplicates(subset = ['ICUSTAY_ID', 'CURR_TIME'], keep = 'first', inplace = True)
                data_expand[group] = pd.merge(data_expand[group], data, on = ['ICUSTAY_ID', 'CURR_TIME'], how = 'left', validate = 'm:1')
            
        # =============================================================================
        #      Output 6, 12, 24 hr
        # =============================================================================
            hours = [6, 12, 24]
            for h in hours:
                if cohort == 'MIMIC':
                    colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_output/output_' + str(h) +'hr/output_' + str(h) +'hr.tsv.txt', sep = '\t').columns.values.astype(str)
                    data = pd.read_csv(datadir[cohort] + group +'/' + group + '_output/output_' + str(h) +'hr/mimic_mimiciii_aki_kidgo_'+positivity+'_output_' + str(h) +'h.csv', header=None)
                    data.columns = colnames
                    data['CURR_TIME'] = pd.to_datetime(data['OUTPUT_' + str(h) +'HR_TIME'])
                if cohort == 'EICU':
                    colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_output/' + str(h) + 'hr/EICU_OUTPUT' + str(h) +'_label.tsv.txt', sep = '\t').columns.values.astype(str)
                    data = pd.read_csv(datadir[cohort] + group +'/' + group + '_output/' + str(h) + 'hr/eicu_eicu_crd_aki_'+positivity+'_output_' + str(h) +'hr.csv', header=None)
                    data.columns = colnames
                    data['CURR_TIME'] = pd.to_datetime(data['OUTPUT_' + str(h) +'HR_TIME'], unit = 'h')
                data.drop_duplicates(subset = ['ICUSTAY_ID', 'CURR_TIME'], keep = 'first', inplace = True)
                data_expand[group] = pd.merge(data_expand[group], data, on = ['ICUSTAY_ID', 'CURR_TIME'], how = 'left', validate = 'm:1')
            
        # =============================================================================
        #      ABG
        # =============================================================================
            if cohort == 'MIMIC':
                colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_bg/mimic_ABG_LABLE.TSV.txt', sep = '\t').columns.values.astype(str)
                data = pd.read_csv(datadir[cohort] + group +'/' + group + '_bg/mimic_mimiciii_aki_kidgo_'+positivity+'_abg.csv', header=None)
                data.columns = colnames
                data['CURR_TIME'] = pd.to_datetime(data['ABG_TIME'])
            if cohort == 'EICU':
                colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_bg/eicu_ABG_LABLE.TSV.txt', sep = '\t').columns.values.astype(str)
                data = pd.read_csv(datadir[cohort] + group +'/' + group + '_bg/eicu_eicu_crd_aki_'+positivity+'_bg2.csv', header=None)
                data.columns = colnames
                data['CURR_TIME'] = pd.to_datetime(data['ABG_TIME'], unit = 'h')
            data.drop_duplicates(subset = ['ICUSTAY_ID', 'CURR_TIME'], keep = 'first', inplace = True)
            data_expand[group] = data_expand[group].merge(data, on = ['ICUSTAY_ID', 'CURR_TIME'], how = 'left')
        
        # =============================================================================
        #      mechvent_starttime mechvent_endtime
        # =============================================================================
            data = {}    
            for t in ['starttime', 'endtime']:
                if cohort == 'MIMIC':
                    colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_mechvent/mechvent_'+ t +'/MIMIC_AKI_mechvent.txt', sep = '\t').columns.values.astype(str)
                    data[t] = pd.read_csv(datadir[cohort] + group +'/' + group + '_mechvent/mechvent_'+ t +'/mimic_mimiciii_aki_kidgo_'+positivity+'_mechvent_'+ t +'.csv', header=None)
                if cohort == 'EICU':
                    colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_mechvent/' + t +'/EICU_AKI_mechvent_'+t[:-4]+'.txt', sep = '\t').columns.values.astype(str)
                    data[t] = pd.read_csv(datadir[cohort] + group +'/' + group + '_mechvent/'+ t +'/eicu_eicu_crd_aki_'+positivity+'_mechvent_'+ t +'.csv', header=None)
                data[t].columns = colnames
                
                
            data_expand[group] = data_expand[group].groupby(by = 'ICUSTAY_ID').progress_apply(lambda x: merge_mech_vaso(x, data, 'MECH', cohort))
            data_expand[group].reset_index(drop = True, inplace = True)
        # =============================================================================
        #      vaso_starttime vaso_endtime
        # =============================================================================
            data = {}    
            if cohort == 'MIMIC':
                for t in ['start', 'end']:
                    colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_vaso/vaso_'+ t +'/MIMIC_vaso_label.tsv.txt', sep = '\t').columns.values.astype(str)
                    data[t] = pd.read_csv(datadir[cohort] + group +'/' + group + '_vaso/vaso_'+ t +'/mimic_mimiciii_aki_kidgo_'+positivity+'_vasopressor_'+t+'time.csv', header=None)
                    data[t].columns = colnames
            
            if cohort == 'EICU':
                t = 'start'
                colnames = pd.read_csv(datadir[cohort] + group +'/' + group + '_vaso/EICU_vaso_label.tsv.txt', sep = '\t').columns.values.astype(str)
                data['start'] = pd.read_csv(datadir[cohort] + group +'/' + group + '_vaso/eicu_eicu_crd_aki_'+positivity+'_vaso.csv', header=None)
                data['start'].columns = colnames
                data['end'] = pd.read_csv(datadir[cohort] + group +'/' + group + '_vaso/eicu_eicu_crd_aki_'+positivity+'_vaso.csv', header=None)
                data['end'].columns = colnames
                
            data_expand[group] = data_expand[group].groupby(by = 'ICUSTAY_ID').progress_apply(lambda x: merge_mech_vaso(x, data, 'VASO', cohort))
            data_expand[group].reset_index(drop = True, inplace = True)
                
                
        data_expand_all = pd.concat((data_expand['aki'], data_expand['non_aki']))
        
        unique_ID = {}
        unique_ID['aki'] = data_expand['aki'].ICUSTAY_ID.unique()
        unique_ID['non_aki'] = data_expand['non_aki'].ICUSTAY_ID.unique()
        
        print('Unique ICUSTAY_ID: AKI: %d, None-AKI: %d' % (len(unique_ID['aki']), len(unique_ID['non_aki'])))
        
        data_expand_all['PF'] = data_expand_all['PO2'] / data_expand_all['FIO2']
        data_expand_all['BMI'] = data_expand_all['WEIGHT'] / ((data_expand_all['HEIGHT']/100) ** 2)
        
        # =============================================================================
        #     Get columns we needed
        # =============================================================================
        
        demotraphic_information = ["ICUSTAY_ID","AGE","GENDER","ETHNICITY","HEIGHT","WEIGHT","BMI"]
        vital_signs = ["TEM","HR","RR","DIAS_BP","SYS_BP","MEAN_BP","SPO2"]
        laboratory_values = ["NA_ION","CA_ION","K_ION","CL_ION","AG","GLU","HCO3","WBC","RBC",\
                             "HGB","HCT","PLT","PF","ALT","AST","TBB","UN","CR"] #CR is SCr
        fluid_balance = ["INPUT_6HR","INPUT_12HR","INPUT_24HR","OUTPUT_6HR","OUTPUT_12HR","OUTPUT_24HR"]
        addtional_respiratory_and_hemodynamic_support = ["VASO", "MECH"]
        primary_diagnosis = ["SEP", "CAR", "NEU", "RES", "OD"] # we don't need NEU
        primary_diagnosis = ["SEP", "CAR", "RES", "OD"] # we don't need NEU
        comorbidities = ["HYP","DIA","CHF", "CPD", "CKD", "CLD"]
        others = ["HOURS", "ISOFA", "AKI"]
        
        colnames_we_need = demotraphic_information + vital_signs + laboratory_values + fluid_balance + \
                            addtional_respiratory_and_hemodynamic_support + primary_diagnosis + comorbidities + others
        
        data_expand_all = data_expand_all[colnames_we_need]
        
        data_expand_all.reset_index(drop = True, inplace = True)
        data_expand_all.to_csv(workdir + 'Processed_Data/data_expand_' + cohort + '.csv')
