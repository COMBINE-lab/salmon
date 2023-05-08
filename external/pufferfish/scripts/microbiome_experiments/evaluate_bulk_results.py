import sys
import numpy as np
import pandas as pd
from sklearn import metrics

def mard(df, c1, c2):
    df['ard'] = abs(df[c1]-df[c2])/(df[c1]+df[c2])
    df.loc[~np.isfinite(df['ard']), 'ard'] = 0#np.nan
    return df['ard'].sum()/df['ard'].count()

def evaluate_results(sample, readName="simulatedReads"):
    print("========================")
    print("\n\nResults for sample: {}\n".format(sample))
    print(truthFile)
    truth = pd.read_csv('{}/{}/{}.profile'.format(truthFile, sample, readName), sep=" ", header=None, names=['count', 'ref'])
    totalReads = truth['count'].sum()
    print("\t==> # references: {}".format(len(truth)))
    print("\t==> # reads: {}".format(totalReads))
    print("\n\n")
    file = '{}/quants/salmon/{}'.format(root, sample)
    results = []
    fileList = [
                '.pufferfish.primary',
                '.bowtie2.primary',
                '.deBGA.primary',
                '.pufferfish.maxNumHits20',
                '.bowtie2.maxNumHits20',
                '.deBGA.maxNumHits20',
                '.pufferfish.maxNumHits200',
                '.bowtie2.maxNumHits200',
                '.deBGA.maxNumHits200'        
    ]
    for fs in zip([file]*len(fileList), fileList):
        filename = fs[0]+fs[1]
        df = pd.read_csv(filename+'/quant.sf', sep="\t")
        df_truth = pd.merge(df, truth, left_on="Name", right_on="ref", how="outer").fillna(0)
        sp = df_truth[['NumReads', 'count']].corr(method='spearman')['NumReads']['count']
        pr = df_truth[['NumReads', 'count']].corr(method='pearson')['NumReads']['count']
        m = mard(df_truth,'NumReads', 'count')
        splited = fs[1].split('.')
        results += [['{} [refCnt({}) -- readCnt({})]'.format(sample,len(truth), totalReads),
                     splited[1], splited[2], np.round(df_truth['NumReads'].sum()*100/totalReads, 1), 
                     np.round(sp, 3), 
#                      np.round(pr, 3), 
#                      np.round(metrics.explained_variance_score(df_truth['NumReads'], df_truth['count']), 3),
#                      np.round(metrics.r2_score(df_truth['NumReads'], df_truth['count']), 3),
                     np.round(m, 3), 
#                      np.round(metrics.median_absolute_error(df_truth['NumReads'], df_truth['count']), 3),
                     np.round(metrics.mean_absolute_error(df_truth['NumReads'], df_truth['count']), 3),
                     np.round(metrics.mean_squared_error(df_truth['NumReads'], df_truth['count']), 1),
                     np.round(metrics.mean_squared_log_error(df_truth['NumReads'], df_truth['count']), 3)]]
    res_df = pd.DataFrame.from_records(results, columns=['sample', 'Tool', 'option',
                                                         'mapping rate(%)',
                                                         'spearman', 
#                                                          'pearson',
#                                                          'explainedVariance', 
#                                                          'R2',
                                                         'mard',
#                                                          'medAE', 
                                                         'MAE',
                                                         'MSE',
                                                         'MSLE'])
    print(res_df)
    return res_df


if __name__ == "__main__":
    if len(sys.argv) < 5 or len(sys.argv) > 5:
        print("""\
        
        This script is part of the metagenomic bulk experiment for Puffaligner and follows the required directory structure for it.
        It compares Salmon quantifications vs truth for a series of samples.

            Required packages:
                numpy
                pandas
                scikit-learn
                openpyxl

        Usage:  python <script> <root_address> <truthFile_address_relative2root> <sample1,sample2,...> <truthFileName_without_the_profile_extension>
        """)
        sys.exit(0)


    root = sys.argv[1]
    truthFile_relative2root = sys.argv[2]
    truthFile = '{}/{}'.format(root, truthFile_relative2root)
    samples = sys.argv[3].split(',')
    with pd.option_context('expand_frame_repr', False):
        l = []
        for s in samples:
            s = s.strip()
            if s == "":
                continue
            print(s)
            res = evaluate_results(s)
            l += [res]
        res_df = pd.concat(l)
        res_dfg = res_df.groupby(['sample','option', 'Tool']).first()
        res_dfg.to_excel('{}/results/bulk_accuracy.xlsx'.format(root))
