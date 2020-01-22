# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 11:14:45 2020

@author: spborder

Viewing and cleaning up GTEx phenotype data

"""

#import numpy as np
#import glob
import pandas as pd
import os

#cwd = os.getcwd()
cwd = 'C:\\Users\\spborder\\Desktop\\'
print(cwd)

Pheno_file = cwd + 'GTEx_Project\\PhenotypeFiles\\'

filenames = os.listdir(Pheno_file)

feat_list = []
col_starts = []
file_list = []
for file in filenames:
    if file[-1] == 't':
        print("reading file:{}".format(file))
        
        ##Extracting the column headers for each file
        with open(Pheno_file+file,"r") as f:
            col_start = 0
            for line in f:
                #print(line)
                if "dbGaP_Subject_ID" in line:
                    feat_list.append(line.split("\t"))
                    file_list.append(file)
                    col_starts.append(col_start)
                    continue
                elif "dbGaP_Sample_ID" in line:
                    feat_list.append(line.split("\t"))
                    file_list.append(file)
                    col_starts.append(col_start)
                    continue
                
                col_start+=1
                
filedata = []            
for j in range(0,len(file_list)):
        
    df = pd.read_csv(Pheno_file+file_list[j],sep="\t",skiprows = col_starts[j]+1, header = None)
   
    df.columns = feat_list[j]
    ###################Comment this line out on multiple run throughs, permissions error########################################
    df.to_csv(Pheno_file+file_list[j].replace('.txt','.csv'))
    ############################################################################################################################
    filedata.append(df)
    


################### Filtering data ##########################################
#############################################################################
    
#Step 1, remove samples with consent = 0
consent = filedata[0]

no_consent = consent[consent['CONSENT']==0]
consent_filt = consent[consent['CONSENT']==1]

#Step 2, remove completely empty columns in Subject_Phenotypes
phenos = filedata[2]

phenos_filt = phenos.dropna(how = 'all',axis = 'columns')

#Step 3, count number of nan's per column to see if there are sparsely recorded features
sum_na = phenos_filt.isna().sum()

## dropping columns with more than half nan's
for index, value in sum_na.iteritems():
    if value>=(0.5)*phenos_filt.shape[0]:
        print('dropping :{}, #NA :{}'.format(index, value))
        phenos_filt = phenos_filt.drop(columns = index)

#Step 4, histograms of every phenotype feature
for feat in feat_list[2]:
    try:
        try:
            
            phenos_filt.hist(column = feat)
            print("feature is :{}".format(feat))
            
        except ValueError:
            print("feature is :{}".format(feat))
            phenos_filt[feat].value_counts().plot(kind='bar',title=feat)
    except KeyError:
        print('Nixed Feature')
        continue

#Step 5, Remove features where >(some percentage) have the same value
for feat in feat_list[2]:
    feat_count = phenos_filt[feat].value_counts()
    
    #Removing single value features
    if feat_count.shape[0]<2:
        phenos_filt = phenos_filt.drop(columns=feat)
        continue
    
def diversity_percentage(df, columns):
    """
    This function returns the number of different elements in each column as a percentage of the total elements in the group.
    A low value indicates there are many repeated elements.
    Example 1: a value of 0 indicates all values are the same.
    Example 2: a value of 100 indicates all values are different.
    """
    diversity = dict()

    for col in columns:
        diversity[col] = len(df[col].unique())

    diversity_series = pd.Series(diversity)
    return (100*diversity_series/len(df)).sort_values()

dp = diversity_percentage(phenos_filt,list(phenos_filt.columns.values))
filt_feats = dp[dp<1] + dp[dp==100]

phenos_filt = phenos_filt.drop(columns = list(filt_feats.index))

phenos_filt.to_csv(Pheno_file+"phenotype_data_filtered.csv")











