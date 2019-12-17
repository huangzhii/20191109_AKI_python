#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 13:19:23 2019

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
import numpy as np
import pandas as pd
from collections import Counter
from sklearn.metrics import auc, roc_curve, f1_score, recall_score, precision_score
from tqdm import tqdm
import gc, logging, copy, pickle, math, random, argparse, time
import matplotlib.pyplot as plt
from sklearn import svm
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier, GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.metrics import matthews_corrcoef
np.warnings.filterwarnings('ignore')
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
sys.path.append(workdir)
from utils_imputation import imputation_on_the_fly
    
def check_missing(x):
    series = pd.Series([int(x[c].isnull().all()) for c in x.columns])
    series.index = x.columns
    return series
    
if __name__ == '__main__':
    with open(workdir + 'Processed_Data/data_expand_MIMIC.pkl', 'rb') as f:
        data_expand_MIMIC = pickle.load(f)
    with open(workdir + 'Processed_Data/data_expand_EICU.pkl', 'rb') as f:
        data_expand_EICU = pickle.load(f)
        
    columns_missing = data_expand_EICU.groupby(by = 'ICUSTAY_ID').progress_apply(lambda x: check_missing(x))
    columns_missing.sum(axis = 0)
    
    columns_missing = data_expand_MIMIC.groupby(by = 'ICUSTAY_ID').progress_apply(lambda x: check_missing(x))
    columns_missing.sum(axis = 0)
    
    
    col_2b_removed = ['ETHNICITY','PF','ALT','AST','TBB']
    print('Remove features that has so many missing data:', col_2b_removed)
    data_expand_MIMIC.drop(col_2b_removed, axis = 1, inplace = True)
    data_expand_EICU.drop(col_2b_removed, axis = 1, inplace = True)
    _ = imputation_on_the_fly(data_expand_MIMIC, series = 6, gap = 6)
    _ = imputation_on_the_fly(data_expand_EICU, series = 6, gap = 6)
