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
    

def perf_measure(y_actual, y_hat):
    TP = 0
    FP = 0
    TN = 0
    FN = 0
    for i in range(len(y_hat)): 
        if y_actual[i]==y_hat[i]==1:
           TP += 1
        if y_hat[i]==1 and y_actual[i]!=y_hat[i]:
           FP += 1
        if y_actual[i]==y_hat[i]==0:
           TN += 1
        if y_hat[i]==0 and y_actual[i]!=y_hat[i]:
           FN += 1
    return(TP, FP, TN, FN)
    
def get_evaluation_res(tp, fp, tn, fn):
    sensitivity = tp/(tp+fn) # recall
    specificity = tn/(tn+fp)
    if (tp+fp) == 0:
        ppv = np.nan
    else:
        ppv = tp/(tp+fp) # precision or positive predictive value (PPV)
    npv = tn/(tn+fn) # negative predictive value (NPV)
    hitrate = (tp+tn)/(tp+tn+fp+fn) # accuracy (ACC)
    return sensitivity, specificity, ppv, npv, hitrate
    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--series', default=1, type=int)
    parser.add_argument('--gap', default=6, type=int)
    parser.add_argument('--no_input_output', default=False, action='store_true')
    parser.add_argument('--no_SCr', default=False, action='store_true')
    parser.add_argument('--result_dir', default=workdir+'Results/', type=str)
    return parser.parse_args()


if __name__ == '__main__':
    random.seed(0)
    np.random.seed(0)
    args = parse_args()
    series = args.series
    gap = args.gap
    data_expand_MIMIC = pd.read_csv(workdir+'Processed_Data/data_expand_MIMIC.csv', index_col = 0)
    data_expand_EICU = pd.read_csv(workdir+'Processed_Data/data_expand_EICU.csv', index_col = 0)
    if args.no_input_output:
        print('Remove input output feature ...')
        data_expand_MIMIC = data_expand_MIMIC[[col for col in data_expand_MIMIC.columns if 'INPUT' not in col and 'OUTPUT' not in col]]
        data_expand_EICU = data_expand_EICU[[col for col in data_expand_MIMIC.columns if 'INPUT' not in col and 'OUTPUT' not in col]]
    if args.no_SCr:
        data_expand_MIMIC = data_expand_MIMIC[[col for col in data_expand_MIMIC.columns if 'CR' not in col]]
        data_expand_EICU = data_expand_EICU[[col for col in data_expand_MIMIC.columns if 'CR' not in col]]
        print('Remove SCr feature ...')
        
    col_2b_removed = ['ETHNICITY','PF','ALT','AST','TBB']
    print('Remove features that has so many missing data:', col_2b_removed)
    data_expand_MIMIC.drop(col_2b_removed, axis = 1, inplace = True)
    data_expand_EICU.drop(col_2b_removed, axis = 1, inplace = True)

# =============================================================================
#     Data normalization
# =============================================================================
    print('Performs normalization from: MIMIC; to: MIMIC and EICU')
    scaler = StandardScaler().fit(data_expand_MIMIC)
    data_expand_MIMIC_scaled = pd.DataFrame(scaler.transform(data_expand_MIMIC))
    data_expand_EICU_scaled = pd.DataFrame(scaler.transform(data_expand_EICU))
    data_expand_MIMIC_scaled.columns, data_expand_MIMIC_scaled.index = data_expand_MIMIC.columns, data_expand_MIMIC.index
    data_expand_EICU_scaled.columns, data_expand_EICU_scaled.index = data_expand_EICU.columns, data_expand_EICU.index
    data_expand_MIMIC_scaled['ICUSTAY_ID'], data_expand_EICU_scaled['ICUSTAY_ID'] = data_expand_MIMIC['ICUSTAY_ID'], data_expand_EICU['ICUSTAY_ID']
    data_expand_MIMIC_scaled['AKI'], data_expand_EICU_scaled['AKI'] = data_expand_MIMIC['AKI'], data_expand_EICU['AKI']
# =============================================================================
#     Imputation on-the-fly
# =============================================================================
    X_MIMIC, y_MIMIC, columns_MIMIC, ICUSTAY_ID_MIMIC = imputation_on_the_fly(data_expand_MIMIC_scaled, series = series, gap = gap)
    X_EICU, y_EICU, columns_EICU, ICUSTAY_ID_EICU = imputation_on_the_fly(data_expand_EICU_scaled, series = series, gap = gap)
    
    
# =============================================================================
#     Prepare dataset
# =============================================================================
    print('forming dataset ... [Training]: MIMIC; [External validation]: EICU.')
    dataset = {}
    dataset['MIMIC'] = {}
    dataset['MIMIC']['column names'] = columns_MIMIC
    dataset['MIMIC']['ICUSTAY_ID'] = ICUSTAY_ID_MIMIC
    
    dataset['EICU'] = {}
    dataset['EICU']['column names'] = columns_EICU
    dataset['EICU']['ICUSTAY_ID'] = ICUSTAY_ID_EICU
    dataset['EICU']['X'] = X_EICU
    dataset['EICU']['y'] = y_EICU
    
    testing_uniqueID = random.sample(set(ICUSTAY_ID_MIMIC), round(len(np.unique(ICUSTAY_ID_MIMIC))/4))
    trainval_uniqueID = [i for i in np.unique(ICUSTAY_ID_MIMIC) if i not in testing_uniqueID]
    dataset['MIMIC']['ICUSTAY_ID_unique_test'] = testing_uniqueID
    dataset['MIMIC']['ICUSTAY_ID_unique_trainval'] = trainval_uniqueID
    
    dataset['MIMIC']['test'] = {}
    dataset['MIMIC']['test']['X'] = X_MIMIC[[v in testing_uniqueID for v in ICUSTAY_ID_MIMIC],:,:]
    dataset['MIMIC']['test']['y'] = y_MIMIC[[v in testing_uniqueID for v in ICUSTAY_ID_MIMIC]]
    
    nfolds = 5
    kf = KFold(n_splits=nfolds, random_state=0, shuffle=True)
    for i, (train_index, val_index) in zip(range(1,nfolds+1), kf.split(trainval_uniqueID)):
        dataset['MIMIC']['fold_' + str(i)] = {}
        dataset['MIMIC']['fold_' + str(i)]['train'] = {}
        dataset['MIMIC']['fold_' + str(i)]['train']['ICUSTAY_ID_unique'] = np.array(trainval_uniqueID)[train_index]
        dataset['MIMIC']['fold_' + str(i)]['train']['X'] = X_MIMIC[[v in np.array(trainval_uniqueID)[train_index] for v in ICUSTAY_ID_MIMIC],:,:]
        dataset['MIMIC']['fold_' + str(i)]['train']['y'] = y_MIMIC[[v in np.array(trainval_uniqueID)[train_index] for v in ICUSTAY_ID_MIMIC]]
        dataset['MIMIC']['fold_' + str(i)]['val'] = {}
        dataset['MIMIC']['fold_' + str(i)]['val']['ICUSTAY_ID_unique'] = np.array(trainval_uniqueID)[val_index]
        dataset['MIMIC']['fold_' + str(i)]['val']['X'] = X_MIMIC[[v in np.array(trainval_uniqueID)[val_index] for v in ICUSTAY_ID_MIMIC],:,:]
        dataset['MIMIC']['fold_' + str(i)]['val']['y'] = y_MIMIC[[v in np.array(trainval_uniqueID)[val_index] for v in ICUSTAY_ID_MIMIC]]
    
# =============================================================================
#     Machine Learning
# =============================================================================
        
    class_weight = {0:0.5, 1:0.5}
    class_weight = dict(Counter(1-y_MIMIC))
    
    for i in range(1, nfolds+1):
        print("%d fold CV -- %d/%d" % (nfolds, i, nfolds))
        mtds = ["logit_l1", "logit_l2", "NN", "NN_l2", "DT", "RF", "AdaBoost", "GBM"]
        for mtd in mtds:
            print(mtd)
            if mtd == "logit_l1": # around 8 mins for all folds
                regr_list = []
                hyperparam_list = [0.5, 1, 1.5, 2, 2.5, 3]
                for c in hyperparam_list:
                    regr_list.append(LogisticRegression(penalty='l1', C=c, solver='saga',class_weight = class_weight,max_iter=2000))
            if mtd == "logit_l2": # around 4 mins for all folds
                regr_list = []
                hyperparam_list = [0.5, 1, 1.5, 2, 2.5, 3]
                for c in hyperparam_list:
                    regr_list.append(LogisticRegression(penalty='l2', C=c, solver='sag',class_weight = class_weight,max_iter=2000))
            if mtd == "NN": # around 2 mins for all folds
                regr_list = []
                hyperparam_list = [16, 32, 64, 128]
                for c in hyperparam_list:
                    regr_list.append(MLPClassifier(solver='adam', alpha=0, max_iter=2000, # alpha is L2 reg
                                     hidden_layer_sizes=(c), random_state=1))
            if mtd == "NN_l2": # around 2 mins for all folds
                regr_list = []
                hyperparam_list = [16, 32, 64, 128]
                for c in hyperparam_list:
                    regr_list.append(MLPClassifier(solver='adam', alpha=1e-5, max_iter=2000, # alpha is L2 reg
                                     hidden_layer_sizes=(c), random_state=1))
            if mtd == "DT":
                regr_list = []
                hyperparam_list = [1,2,3,4,5,6,7,8,9,10]
                for c in hyperparam_list:
                    regr_list.append(DecisionTreeClassifier(max_depth = c, class_weight = class_weight))
            if mtd == "RF":
                regr_list = []
                hyperparam_list = [16, 32, 64, 128, 256]
                for c in hyperparam_list:
                    regr_list.append(RandomForestClassifier(n_estimators = c, class_weight = class_weight))
            if mtd == "AdaBoost":
                regr_list = []
                hyperparam_list = [16, 32, 64, 128, 256]
                for c in hyperparam_list:
                    regr_list.append(AdaBoostClassifier(n_estimators = c))
            if mtd == "GBM":
                regr_list = []
                hyperparam_list = [16, 32, 64, 128, 256]
                for c in hyperparam_list:
                    regr_list.append(GradientBoostingClassifier(n_estimators = c))
                    
    
            TIMESTRING  = time.strftime("%Y%m%d-%H.%M.%S", time.localtime())
            if args.no_input_output and not args.no_SCr :
                results_dir_dataset = args.result_dir + 'No_input_output/Serie_' + str(args.series) + '_Gap_' + str(args.gap) + '/' + mtd + '/'
            elif args.no_SCr and not args.no_input_output :
                results_dir_dataset = args.result_dir + 'No_SCr/Serie_' + str(args.series) + '_Gap_' + str(args.gap) + '/' + mtd + '/'
            elif args.no_input_output and args.no_SCr:
                results_dir_dataset = args.result_dir + 'No_input_output_SCr/Serie_' + str(args.series) + '_Gap_' + str(args.gap) + '/' + mtd + '/'
            else:
                results_dir_dataset = args.result_dir + 'All_features/Serie_' + str(args.series) + '_Gap_' + str(args.gap) + '/' + mtd + '/'
            if not os.path.exists(results_dir_dataset):
                os.makedirs(results_dir_dataset)
    
            # create logger
            logger = logging.getLogger(TIMESTRING)
            logger.setLevel(logging.DEBUG)
            # create file handler which logs even debug messages
            fh = logging.FileHandler(results_dir_dataset+'mainlog.log', mode='w')
            fh.setLevel(logging.DEBUG)
            logger.addHandler(fh)
            
            curr_data = dataset['MIMIC']['fold_' + str(i)]
            
            X_train = curr_data['train']['X']
            y_train = curr_data['train']['y']
            X_val = curr_data['val']['X']
            y_val = curr_data['val']['y']
            
            X_train_reshaped = X_train.reshape(X_train.shape[0], -1)
            X_val_reshaped = X_val.reshape(X_val.shape[0], -1)
        # =============================================================================
        #             Find Hyper-parameters
        # =============================================================================
            f1_val_list = []
            for idx, regr in enumerate(regr_list):
                regr.fit(X_train_reshaped, y_train)
                y_pred = regr.predict(X_val_reshaped)
                f1_val = f1_score(y_pred, y_val, average = 'macro')
                f1_val_list.append(f1_val)
                print("Current model: %s, current hyper-parameter: %s, current F1 score: %.8f" % (mtd, str(hyperparam_list[idx]), f1_val) )
                logger.log(logging.INFO, "Current model: %s, current hyper-parameter: %s, current F1 score: %.8f" % (mtd, str(hyperparam_list[idx]), f1_val) )
            
            # Choose the optimal regr
            print("Best hyper-parameter: %s" % str(hyperparam_list[np.argmax(f1_val_list)]))
            logger.log(logging.INFO, "Best hyper-parameter: %s" % str(hyperparam_list[np.argmax(f1_val_list)]))
            print("Model:")
            logger.log(logging.INFO, "Model:")
            regr = regr_list[np.argmax(f1_val_list)]
            print(regr)
            logger.log(logging.INFO, regr)
            
        # =============================================================================
        #             Training after searching hyper-parameters
        # =============================================================================
            
            X_train = np.concatenate((curr_data['train']['X'], curr_data['val']['X']), 0)
            y_train = np.concatenate((curr_data['train']['y'], curr_data['val']['y']), 0)
            X_test = dataset['MIMIC']['test']['X']
            y_test = dataset['MIMIC']['test']['y']
            X_train_reshaped = X_train.reshape(X_train.shape[0], -1)
            X_test_reshaped = X_test.reshape(X_test.shape[0], -1)
            
            regr.fit(X_train_reshaped, y_train)
            
            # AUC F1 train
            y_pred = regr.predict(X_train_reshaped)
            y_pred_proba = regr.predict_proba(X_train_reshaped)[:,1]
            fpr, tpr, thresholds = roc_curve(y_train, y_pred_proba)
            auc_train = auc(fpr, tpr)
            f1_train = f1_score(y_pred, y_train, average = 'macro')
            
            # AUC F1 test
            y_pred = regr.predict(X_test_reshaped)
            y_pred_proba = regr.predict_proba(X_test_reshaped)[:,1]
            fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
            auc_test = auc(fpr, tpr)
            f1_test = f1_score(y_pred, y_test, average = 'macro')
            
            # P, R test
            precision_test = precision_score(y_test, y_pred, average = 'macro')
            recall_test = recall_score(y_test, y_pred, average = 'macro')
            
            # other statistics
            TP, FP, TN, FN = perf_measure(y_test, y_pred)
            mcc = matthews_corrcoef(y_test, y_pred)
            sensitivity, specificity, ppv, npv, hitrate = get_evaluation_res(TP, FP, TN, FN)
            
            print("[MIMIC] train AUC: %.8f, test AUC: %.8f, train F1: %.8f, test F1: %.8f, test Precision: %.8f, test Recall: %.8f" \
                  % (auc_train, auc_test, f1_train, f1_test, precision_test, recall_test))
            logger.log(logging.INFO, "[MIMIC] train AUC: %.8f, test AUC: %.8f, train F1: %.8f, test F1: %.8f, test Precision: %.8f, test Recall: %.8f" \
                  % (auc_train, auc_test, f1_train, f1_test, precision_test, recall_test))
        
            with open(results_dir_dataset + 'regr_model.pickle', 'wb') as handle:
                pickle.dump(regr, handle, protocol=pickle.HIGHEST_PROTOCOL)
            with open(results_dir_dataset + 'outputs_test_proba_MIMIC.pickle', 'wb') as handle:
                pickle.dump(y_pred_proba, handle, protocol=pickle.HIGHEST_PROTOCOL)
            with open(results_dir_dataset + 'outputs_test_bin_MIMIC.pickle', 'wb') as handle:
                pickle.dump(y_pred, handle, protocol=pickle.HIGHEST_PROTOCOL)
            with open(results_dir_dataset + 'column_names.pickle', 'wb') as handle:
                pickle.dump(dataset['MIMIC']['column names'], handle, protocol=pickle.HIGHEST_PROTOCOL)
                
            res_table = pd.DataFrame([auc_train, auc_test, f1_train, f1_test, precision_test, recall_test, \
                                      sensitivity, specificity, ppv, npv, hitrate, mcc])
            res_table.index = ['auc_train', 'auc_test', 'f1_train', 'f1_test', 'precision_test', 'recall_test', \
                               'sensitivity', 'specificity', 'ppv', 'npv', 'hitrate', 'MCC']
            res_table.to_csv(results_dir_dataset + 'performance_MIMIC.csv')
            
            plt.figure(figsize=(8,4))
            plt.plot(fpr, tpr, color='darkorange',
                     lw=2, label='%s test set (AUC = %0.4f%%)' % (mtd, 100*auc(fpr, tpr)))
            plt.axes().set_aspect('equal')
            plt.plot([0, 1], [0, 1], color='navy', lw=1, linestyle='--')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('False Positive Rate (FPR)')
            plt.ylabel('True Positive Rate (TPR)')
            plt.title('Receiver Operating Characteristic Curve')
            plt.legend(loc="lower right")
            plt.savefig(results_dir_dataset + "AUC_test_MIMIC.png",dpi=300)
        
            # =============================================================================
            #             External validation on EICU
            # =============================================================================
            
            X_test = dataset['EICU']['X']
            y_test = dataset['EICU']['y']
            
            X_test_reshaped = X_test.reshape(X_test.shape[0], -1)
            
            # AUC F1 test
            y_pred = regr.predict(X_test_reshaped)
            y_pred_proba = regr.predict_proba(X_test_reshaped)[:,1]
            fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
            auc_test = auc(fpr, tpr)
            f1_test = f1_score(y_pred, y_test, average = 'macro')
            
            # P, R test
            precision_test = precision_score(y_test, y_pred, average = 'macro')
            recall_test = recall_score(y_test, y_pred, average = 'macro')
            
            # other statistics
            TP, FP, TN, FN = perf_measure(y_test, y_pred)
            mcc = matthews_corrcoef(y_test, y_pred)
            sensitivity, specificity, ppv, npv, hitrate = get_evaluation_res(TP, FP, TN, FN)
            
            print("[EICU] test AUC: %.8f, test F1: %.8f, test Precision: %.8f, test Recall: %.8f" \
                  % (auc_test, f1_test, precision_test, recall_test))
            logger.log(logging.INFO, "[EICU] test AUC: %.8f, test F1: %.8f, test Precision: %.8f, test Recall: %.8f" \
                  % (auc_test, f1_test, precision_test, recall_test))
        
            with open(results_dir_dataset + 'regr_model.pickle', 'wb') as handle:
                pickle.dump(regr, handle, protocol=pickle.HIGHEST_PROTOCOL)
            with open(results_dir_dataset + 'outputs_test_proba_EICU.pickle', 'wb') as handle:
                pickle.dump(y_pred_proba, handle, protocol=pickle.HIGHEST_PROTOCOL)
            with open(results_dir_dataset + 'outputs_test_bin_EICU.pickle', 'wb') as handle:
                pickle.dump(y_pred, handle, protocol=pickle.HIGHEST_PROTOCOL)
                
            res_table = pd.DataFrame([auc_test, f1_test, precision_test, recall_test, \
                                      sensitivity, specificity, ppv, npv, hitrate, mcc])
            res_table.index = ['auc_test', 'f1_test', 'precision_test', 'recall_test', \
                               'sensitivity', 'specificity', 'ppv', 'npv', 'hitrate', 'MCC']
            res_table.to_csv(results_dir_dataset + 'performance_EICU.csv')
            
            plt.figure(figsize=(8,4))
            plt.plot(fpr, tpr, color='darkorange',
                     lw=2, label='%s test set (AUC = %0.4f%%)' % (mtd, 100*auc(fpr, tpr)))
            plt.axes().set_aspect('equal')
            plt.plot([0, 1], [0, 1], color='navy', lw=1, linestyle='--')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('False Positive Rate (FPR)')
            plt.ylabel('True Positive Rate (TPR)')
            plt.title('Receiver Operating Characteristic Curve')
            plt.legend(loc="lower right")
            plt.savefig(results_dir_dataset + "AUC_test_EICU.png",dpi=300)
        
    
    
    
    
    
    
    
    
    
    
    
    
