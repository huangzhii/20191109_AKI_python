#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 01:34:39 2019

@author: zhihuan
"""


import re, copy
import time
import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
import sys, os, platform
import pickle
import itertools
import argparse
import urllib
import json
from tqdm import tqdm
tqdm.pandas()


#print("Currently working on " + list(platform.uname())[1] + " Machine")
#if 'Y710' in list(platform.uname())[1]:
#    workdir = '/media/zhihuan/Drive3/20191109_AKI_python/'
#elif 'Zhi-G7-7790' in list(platform.uname())[1]:
#    workdir = '/home/zhihuan/Documents/20191109_AKI_python/'
#else:
#    workdir = ''
    

def processing(data, series = 6, gap = 6):
    Y = data[['AKI', 'HOURS']]
    for col in ['INPUT_MINUS_OUTPUT_6HR','INPUT_MINUS_OUTPUT_12HR','INPUT_MINUS_OUTPUT_24HR', 'AKI']:
        if col in data.columns:
            data.drop(col, 1, inplace = True)
    X = data
#    if 'SUBJECT_ID' in X.columns or 'HADM_ID' in X.columns or 'PT' in X.columns:
#        X = X.drop(['SUBJECT_ID','HADM_ID','PT'], 1) # PT is not showing up in EICU cohort
        
    if not all(np.isreal(item) == True for item in X['GENDER'].values):
        X['GENDER'][X['GENDER'] == 'F'] = 0
        X['GENDER'][X['GENDER'] == 'M'] = 1
    
    X_is_vent = X.loc[Y['AKI'] == 1]
    X_non_vent = X.loc[Y['AKI'] == 0]
#    X_is_vent.reset_index(drop = True, inplace = True)
#    X_non_vent.reset_index(drop = True, inplace = True)
    X_non_vent.index = range(0,X_non_vent.shape[0])
    diff = [y-x for x, y in zip(Y.loc[:,'HOURS'], Y.loc[1:,'HOURS'])]
    diff.append(-1)
    diff = np.asarray(diff)
    
    y_is_vent = Y.loc[[a and b for a, b in zip((Y['AKI'] == 1).values, diff < 0)],:]
    y_non_vent = Y.loc[Y['AKI'] == 0, ]
#    y_is_vent.reset_index(drop = True, inplace = True)
#    y_non_vent.reset_index(drop = True, inplace = True)
    y_non_vent.index = range(0,X_non_vent.shape[0])
    
    y_is_vent_index = y_is_vent.index.values
    y_non_vent_index = y_non_vent.index.values
    
    
    X_is_vent_reshaped, columns, ICUSTAY_ID_is_vent = crop_reshape(X_is_vent, y_is_vent_index, series, gap)
    y_is_vent = [1]*X_is_vent_reshaped.shape[0]
    X_non_vent_reshaped, y_non_vent, ICUSTAY_ID_non_vent = None, None, None
    X, y = X_is_vent_reshaped, y_is_vent
    ICUSTAY_ID = ICUSTAY_ID_is_vent
    
#    if len(y_non_vent_index) > 10: # not EICU data
    X_non_vent_reshaped, columns, ICUSTAY_ID_non_vent = crop_reshape(X_non_vent, y_non_vent_index, series, gap)
    y_non_vent = [0]*X_non_vent_reshaped.shape[0]
    X = np.concatenate((X, X_non_vent_reshaped), axis=0)
    y = np.concatenate((y, y_non_vent))
    ICUSTAY_ID = np.concatenate((ICUSTAY_ID,ICUSTAY_ID_non_vent)).astype(int)
    
    columns = columns[[c not in ['AKI'] for c in columns]]
    return X, y, columns, ICUSTAY_ID

def crop_reshape(X, y, series, gap):
    data_index = y - gap - 1
    data_index = data_index.reshape(1, len(data_index))
    data_index_concat = y.reshape(1, len(y)) # also add the label row
    for i in range(series):
        data_index_concat = np.concatenate((data_index-i, data_index_concat), axis = 0)
    X_ICU = X['ICUSTAY_ID'].values
    X_ICU_series = X_ICU[data_index_concat]
    index2keep1 = data_index_concat[0,:] >= 0 # remove columns that have negative index
    index2keep2 = np.all(X_ICU_series == X_ICU_series[0,:], axis = 0) # is all the column with the same ICUSTAY_ID?
    index2keep = np.logical_and(index2keep1, index2keep2)
    data_index_concat = data_index_concat[:, index2keep]
    data_index_concat = data_index_concat[0:data_index_concat.shape[0]-1,:] #remove the label row
    y = y[index2keep]
    
    # before reshape, let's remove some irrelevant columns
    cols2keep = ['ICUSTAY_ID', 'SUBJECT_ID', 'HADM_ID', 'LOS']
    cols2keep = [X.columns.get_loc(c) for c in X.columns if c not in cols2keep]
    cols2keep_names = X.columns[cols2keep].values
    X_np = X.values.astype(float)
    
# =============================================================================
#     Reshape
# =============================================================================
    X_reshaped = X_np[data_index_concat, :]
    ICUSTAY_ID = X_reshaped[0,:,0]
    X_reshaped = X_reshaped[:,:, cols2keep]
    X_reshaped = np.swapaxes(X_reshaped,0,1)
    return X_reshaped, cols2keep_names, ICUSTAY_ID


def imputation_on_the_fly(data_expand, series = 6, gap = 6):
    # Load expanded data
    if 'PF' in data_expand.columns:
        print('Remove PF_ratio (so many missing)')
        data_expand.drop('PF', axis = 1, inplace = True)
    if 'ETHNICITY' in data_expand.columns:
        print('Remove ETHNICITY')
        data_expand.drop('ETHNICITY', axis = 1, inplace = True)
        
        
    # =============================================================================
    #       Data interpolation
    # =============================================================================
    print('Now perform data interpolation...')
    def fillNAinterpolation(x):
        x2 = x.fillna(method='ffill', axis = 0)
        x2 = x2.fillna(method='bfill', axis = 0)
        return x2
    data_expand_fillNA = data_expand.groupby(by = 'ICUSTAY_ID').progress_apply(lambda x: fillNAinterpolation(x))
    
    print('Before drop NA - data shape:', data_expand_fillNA.shape)
    data_expand_fillNA.dropna(inplace = True)
    data_expand_fillNA.reset_index(inplace = True, drop = True)
    print('After drop NA - data shape:', data_expand_fillNA.shape)
        
    # =============================================================================
    #       Data reshape
    # =============================================================================
    X, y, columns, ICUSTAY_ID = processing(data = copy.deepcopy(data_expand_fillNA),\
                                              series = series, gap = gap)
    print('After data reshape, AKI:', sum(y == 1), '( unique:',len(np.unique(ICUSTAY_ID[y == 1])),')', \
          ' |  None AKI:', sum(y == 0), '( unique:',len(np.unique(ICUSTAY_ID[y == 0])),')')
        
    return X, y, columns, ICUSTAY_ID
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    