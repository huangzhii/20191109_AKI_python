#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 21:36:43 2019

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
import aki_utils
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
from sklearn.preprocessing import OneHotEncoder

sys.path.append(workdir)
from utils_imputation import imputation_on_the_fly
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
    parser.add_argument('--series', default=6, type=int)
    parser.add_argument('--gap', default=6, type=int)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    series = args.series
    gap = args.gap
    data_expand = pd.read_csv(workdir+'Processed_Data/data_expand_MIMIC.csv', index_col = 0)
    X, y, columns, ICUSTAY_ID = imputation_on_the_fly(data_expand, series = 6, gap = 6)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    