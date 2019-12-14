#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 15:46:59 2019

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
from tqdm import tqdm
import gc, logging, copy, pickle, math, random, argparse, time
import matplotlib.pyplot as plt
np.warnings.filterwarnings('ignore')
tqdm.pandas()


print("Currently working on " + list(platform.uname())[1] + " Machine")
if 'Y710' in list(platform.uname())[1]:
    workdir = '/media/zhihuan/Drive3/20191109_AKI_python/'
elif 'Zhi-G7-7790' in list(platform.uname())[1]:
    workdir = '/home/zhihuan/Documents/20191109_AKI_python/'
else:
    workdir = ''
    
    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--result_dir', default=workdir+'Results/All_features/', type=str)
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    result_dirs = [args.result_dir + s + '/' for s in os.listdir(args.result_dir) if 'Serie' in s]
    
    series, gaps = set(), set()
    for result_dir in result_dirs:
        series.add(int(result_dir.split('Serie_')[1].split('_Gap')[0]))
        gaps.add(int(result_dir.split('Gap_')[1].split('/')[0]))
    
    serie_gap_list = []
    for s in sorted(series):
        for g in sorted(gaps):
            serie_gap_list.append('Serie_%d_Gap_%d' % (s, g))
        
    result = {}
    metrics = ['auc_train', 'auc_test', 'f1_train', 'f1_test', 'precision_test', \
       'recall_test', 'sensitivity', 'specificity', 'ppv', 'npv', 'hitrate', 'MCC']
    methods = [mtd for mtd in os.listdir(args.result_dir + serie_gap_list[0])]
    methods = ['logit_l1', 'logit_l2']
    for c in metrics:
        result[c] = pd.DataFrame(columns = serie_gap_list, index = methods)
    
    for c in metrics:
        for mtd in methods:
            for sg in serie_gap_list:
                result[c].loc[mtd,sg] = pd.read_csv(args.result_dir + sg + '/' + mtd + '/performance_MIMIC.csv', index_col = 0).loc[c,].values[0]
    
    
    for c in metrics:
        result[c].to_csv(args.result_dir + 'performance_' + c + '.csv')
    
    