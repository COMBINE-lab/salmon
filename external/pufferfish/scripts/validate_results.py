#!/usr/bin/env python

import argparse
import scipy.stats
import scipy.spatial
import sklearn.metrics
import math
import pandas as pd
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(__file__, description = "")
    parser.add_argument('-truer', '-tr', help = 'True counts file', default = None)
    parser.add_argument('-est', '-s', help = 'out.tab file', default = None)
    parser.add_argument('-truthType', '-tt', help = 'what is type of your truth, transcriptome, microbiome, or taxanomy', default=None)

args = parser.parse_args()

print(args.est)

print("\n")

def mard(df, c1, c2):
    df['ard'] = abs(df[c1]-df[c2])/(df[c1]+df[c2])
    df.loc[~np.isfinite(df['ard']), 'ard'] = 0
    return df['ard'].mean()#sum()/df['ard'].count()

def calc_ref_metrics(truer, quants,  est_col, tru_col='cnt'):
    truth = pd.read_table(truer)
    truth.columns = ['id', 'cnt']
    estimated = pd.read_table(quants)
    estimated['id'] = estimated['Name'].str.split('|').str[0]
    est = estimated.groupby('id')[est_col].sum().reset_index()
    metrics = { "R^2" : sklearn.metrics.r2_score, \
                "Explained Var." : sklearn.metrics.explained_variance_score, \
                "Mean Abs Error" : sklearn.metrics.mean_absolute_error, \
                "Mean Sq. Error" : sklearn.metrics.mean_squared_error, \
                "Mean Sq. Log. Error" : sklearn.metrics.mean_squared_log_error, \
                "Med Abs Error" : sklearn.metrics.median_absolute_error, \
                "Bray Curtis" : scipy.spatial.distance.braycurtis, \
                "Kendall Tau" : scipy.stats.kendalltau, \
                "Cosine Similarity" : scipy.spatial.distance.cosine, \
#                "Minkowski Distance" : scipy.spatial.distance.minkowski, \
                "Canberra Distance": scipy.spatial.distance.canberra}

    metric_res = []
    merged=pd.merge(truth, est, on='id', how='outer').fillna(0)
    MARD = mard(merged, tru_col, est_col)
    pcc = merged[[tru_col, est_col]].corr(method='pearson')[tru_col][est_col]
    sp = merged[[tru_col, est_col]].corr(method='spearman')[tru_col][est_col]
    metric_dict = {'mard': MARD, 'pcc': pcc, 'sp': sp}
    for k, v in metrics.items():
        val = v(merged[tru_col], merged[est_col])
        if k == "Kendall Tau":
            metric_dict[k] = val.correlation
        else:
            metric_dict[k] = val
    df = pd.DataFrame(metric_dict, index=[0])
    merged.columns=['refid','truth', 'pred', 'ard']
    merged['abs_err'] = abs(merged['truth']-merged['pred'])
    return merged, df.transpose()

def calc_taxlevel_metrics(truer, est):
    true = pd.read_table(truer)
    estimated = pd.read_table(est)
    metrics = { "R^2" : sklearn.metrics.r2_score, \
                "Explained Var." : sklearn.metrics.explained_variance_score, \
                "Mean Abs Error" : sklearn.metrics.mean_absolute_error, \
                "Mean Sq. Error" : sklearn.metrics.mean_squared_error, \
                "Mean Sq. Log. Error" : sklearn.metrics.mean_squared_log_error, \
                "Med Abs Error" : sklearn.metrics.median_absolute_error, \
                "Bray Curtis" : scipy.spatial.distance.braycurtis, \
                "Kendall Tau" : scipy.stats.kendalltau, \
                "Cosine Similarity" : scipy.spatial.distance.cosine, \
#                "Minkowski Distance" : scipy.spatial.distance.minkowski, \
                "Canberra Distance": scipy.spatial.distance.canberra}


    metric_res = {}
    for r in ["Phylum", "Genus", "Species", "Scietific_Name", "TaxID"]:
        metric_res[r] = {}
        tr_taxid_counts = true.groupby(r)['NumCounts'].sum().reset_index()
        es_taxid_counts = estimated.groupby(r)['NumCounts'].sum().reset_index()
        merged=pd.merge(tr_taxid_counts, es_taxid_counts, on=r, how='outer').fillna(0)
        MARD = mard(merged, 'NumCounts_x', 'NumCounts_y')
        pcc = merged[['NumCounts_x', 'NumCounts_y']].corr(method='pearson')['NumCounts_x']['NumCounts_y']
        sp = merged[['NumCounts_x', 'NumCounts_y']].corr(method='spearman')['NumCounts_x']['NumCounts_y']
        metric_res[r]['mard'] = MARD
        metric_res[r]['pcc'] = pcc
        metric_res[r]['sp'] = sp
        #metric_res[r] += [MARD, pcc, sp]
        for k, v in metrics.items():
            val = v(merged['NumCounts_x'], merged['NumCounts_y'])
            if k == "Kendall Tau":
                metric_res[r][k] = val.correlation
            else:
                metric_res[r][k] = val
        df = pd.DataFrame(metric_res)
        #df.columns=['mard', 'pcc', 'sp', 'r2', 'ex_var', 'mae', 'mse', 'msle', 'medae']
        merged.columns=['taxid','truth', 'pred', 'ard']
        merged['abs_err'] = abs(merged['truth']-merged['pred'])
    return merged, df[['Phylum', 'Genus', 'Species', 'Scietific_Name', 'TaxID']]

if args.truthType == 'transcriptome':
    d, m = calc_ref_metrics(args.truer, args.est, 'TPM')
    print(m)
elif args.truthType == 'microbiome':
    d, m = calc_ref_metrics(args.truer, args.est, 'NumReads')
    print(m)
elif args.truthType == 'taxonomy':
    d, m = calc_taxlevel_metrics(args.truer, args.est)
    print(m)
print("\n")
