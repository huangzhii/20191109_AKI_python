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
from sklearn.tree import export_graphviz
import matplotlib.cm as mplcm
import matplotlib.colors as colors
np.warnings.filterwarnings('ignore')
tqdm.pandas()


print("Currently working on " + list(platform.uname())[1] + " Machine")
if 'Y710' in list(platform.uname())[1]:
    workdir = '/media/zhihuan/Drive3/20191109_AKI_python/'
elif 'Zhi-G7-7790' in list(platform.uname())[1]:
    workdir = '/media/zhihuan/DATA/20191109_AKI_python/'
elif 'DESKTOP-05QACO1' in list(platform.uname())[1]:
    workdir = 'D:/20191109_AKI_python/'
elif 'dl' in list(platform.uname())[1]: # IU deep learning server
    workdir = '/gpfs/home/z/h/zhihuan/Carbonate/Desktop/20191109_AKI_python/'
elif 'uits' in list(platform.uname())[1]: # IU deep learning server
    workdir = '/gpfs/home/z/h/zhihuan/Carbonate/Desktop/20191109_AKI_python/'
else:
    workdir = ''
    
    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--result_dir', default=workdir+'Results/All_features/', type=str)
#    parser.add_argument('--result_dir', default=workdir+'Results/No_input_output/', type=str)
#    parser.add_argument('--result_dir', default=workdir+'Results/No_input_output_SCr/', type=str)
#    parser.add_argument('--result_dir', default=workdir+'Results/No_SCr/', type=str)
#    parser.add_argument('--result_dir', default=workdir+'Results_limit=2/All_features/', type=str)
#    parser.add_argument('--result_dir', default=workdir+'Results_limit=4/All_features/', type=str)
#    parser.add_argument('--result_dir', default=workdir+'Results_limit=6/All_features/', type=str)
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
    methods = [mtd for mtd in os.listdir(args.result_dir + serie_gap_list[0])]
    methods = ['logit_l1', 'logit_l2', 'DT', 'RF', 'AdaBoost', 'GBM', 'NN', 'NN_l2']
#    methods = ['logit_l1', 'logit_l2']
    metrics = ['auc_train', 'auc_test', 'f1_train', 'f1_test', 'precision_test', \
       'recall_test', 'sensitivity', 'specificity', 'ppv', 'npv', 'hitrate', 'MCC']
    for c in metrics:
        result[c] = pd.DataFrame(columns = serie_gap_list, index = methods)
    for c in metrics:
        for mtd in methods:
            for sg in serie_gap_list:
                result[c].loc[mtd,sg] = pd.read_csv(args.result_dir + sg + '/' + mtd + '/performance_MIMIC.csv', index_col = 0).loc[c,].values[0]
        result[c].to_csv(args.result_dir + 'performance_' + c + '.csv')
        
    metrics = ['auc_test', 'f1_test', 'precision_test', \
       'recall_test', 'sensitivity', 'specificity', 'ppv', 'npv', 'hitrate', 'MCC']
    for c in metrics:
        result[c] = pd.DataFrame(columns = serie_gap_list, index = methods)
    for c in metrics:
        for mtd in methods:
            for sg in serie_gap_list:
                result[c].loc[mtd,sg] = pd.read_csv(args.result_dir + sg + '/' + mtd + '/performance_EICU.csv', index_col = 0).loc[c,].values[0]
        result[c].to_csv(args.result_dir + 'eICU_performance_' + c + '.csv')
    
## =============================================================================
##     Feature ranking
## =============================================================================
#    mdl_dir = workdir + 'Results/All_features/Serie_6_Gap_6/logit_l1/'
#    with open(mdl_dir + 'regr_model.pickle', 'rb') as f:
#        model = pickle.load(f)
#    with open(mdl_dir + 'column_names.pickle', 'rb') as f:
#        colnames = pickle.load(f)
#    
#    rank = pd.DataFrame(index = colnames, columns = ['coefficient'])
#    rank['coefficient'] = model.coef_.reshape(-1)
#    rank.sort_values(by = 'coefficient', inplace = True, ascending = False)
#    rank.to_csv(mdl_dir + 'feature_ranking.csv')
        
        
    mdl_dir = workdir + 'Results/All_features/Serie_6_Gap_6/GBM/'
    with open(mdl_dir + 'regr_model.pickle', 'rb') as f:
        model = pickle.load(f)
    with open(mdl_dir + 'column_names.pickle', 'rb') as f:
        colnames = pickle.load(f)
        colnames[colnames == 'CR'] = 'SCr'
    rank = pd.DataFrame(index = colnames, columns = ['coefficient'])
    rank['coefficient'] = np.sum(model.feature_importances_.reshape(len(colnames),-1), 1)
    rank.sort_values(by = 'coefficient', inplace = True, ascending = False)
    rank.to_csv(mdl_dir + 'feature_importances_.csv')



# =============================================================================
# cm = plt.get_cmap('OrRd_r')
# possible colormap are: Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, 
# BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r, Greens, Greens_r, Greys,
# Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1,
# Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r,
# PuRd, PuRd_r, Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu,
# RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1, Set1_r, Set2, Set2_r, Set3, Set3_r,
# Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r,
# YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r, brg,
# brg_r, bwr, bwr_r, cividis, cividis_r, cool, cool_r, coolwarm, coolwarm_r, copper,
# copper_r, cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r, gist_gray,
# gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r,
# gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r,
# gray, gray_r, hot, hot_r, hsv, hsv_r, inferno, inferno_r, jet, jet_r, magma, magma_r,
# nipy_spectral, nipy_spectral_r, ocean, ocean_r, pink, pink_r, plasma, plasma_r, prism,
# prism_r, rainbow, rainbow_r, seismic, seismic_r, spring, spring_r, summer, summer_r,
# tab10, tab10_r, tab20, tab20_r, tab20b, tab20b_r, tab20c, tab20c_r, terrain, terrain_r,
# twilight, twilight_r, twilight_shifted, twilight_shifted_r, viridis, viridis_r, winter, winter_r
# =============================================================================
    top = 17
    cm = plt.get_cmap('OrRd_r')
    cNorm  = colors.Normalize(vmin=0, vmax=top-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    rgba_color = [scalarMap.to_rgba(i) for i in range(top)]
    plt.figure(figsize=(12,8))
    plt.barh(range(len(rank))[:top], rank.coefficient.values[:top], height=0.8, align='center', \
             color = rgba_color, edgecolor = 'grey')
    for i, v in enumerate(rank.coefficient.values[:top]):
        plt.text(v + 0.002, i + .25, str(round(v,4)), color='black', fontname="Arial", fontsize=16)
    plt.yticks(np.arange(top), rank.index.values[:top], fontname="Arial", fontsize=16)
    plt.xticks(fontname="Arial", fontsize=16)
    plt.xlabel('Gini importance', fontname="Arial", fontsize=20)
    plt.ylabel('Important features (Gini importance >= 0.01)', fontname="Arial", fontsize=20)
    plt.ylim((-1, top))
    plt.xlim((0, max(rank.coefficient.values) + 0.02))
    plt.gca().invert_yaxis()
    plt.suptitle('Feature importance derived from Gradient boosting classifier', fontname="Arial", fontsize=20)
#    plt.grid(on)
    plt.savefig(fname = mdl_dir + 'feature_importances_top_' + str(top) +'.png', dpi = 600, bbox_inches='tight')



# =============================================================================
# Check optimal hyper-parameters
# =============================================================================
    result_dirs = [args.result_dir + s + '/' for s in os.listdir(args.result_dir) if 'Serie' in s]
    methods = ['logit_l2', 'DT', 'RF', 'AdaBoost', 'GBM', 'NN', 'NN_l2']
    for r_dir in result_dirs:
        print('--------------------------------------------------------------')
        print(r_dir.split('/')[-2])
        for mtd in methods:
            print(mtd)
            with open(r_dir + mtd + '/regr_model.pickle', 'rb') as f:
                regr_model = pickle.load(f)
                print(regr_model)
