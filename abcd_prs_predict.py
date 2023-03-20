import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import pyreadr
import sklearn
from sklearn.naive_bayes import MultinomialNB
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
import scipy.stats as stats
import statsmodels.api as sm
from pathlib import Path
import random
import glob
from functools import reduce
import os.path
from os import path
from sklearn import preprocessing


# ============= WIP script for genetic prediction models in ABCD =================

# first input class data from mplus output > merge with original NDAR IDs
# read in, select, clean and normalise all variables [0,1]
# build design matrix
# impute missing data
# split test and train > run ML models
# evaluate feature importance

# class data
bpm_4k = pd.read_table('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/gmm/gmm_abcd/'
                      'mplus_data/4k_probs_bpm.txt', delim_whitespace=True, header=None)
bpm_4k.columns = ['y0.5','y1','y1.5','y2','y2.5','y3','v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','class','id']
bpm_4k = bpm_4k[['id','class']] # subset
bpm_4k.replace('*', np.nan) # fill missing

# class data has unique id's only > merge with anth to align NDAR id's
os.chdir('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/gmm/gmm_abcd')

# if running script on eddie
# from abcd_PREP.py import get_dataframes
# df1, df2 = get_dataframes()

# if local
# anth cols: [id time age sex og_id]
anth = pd.read_csv("/Volumes/igmm/GenScotDepression/users/poppy/abcd/abcd_anth.csv")
# append class enumeration
bpm_4k = pd.merge(bpm_4k, anth, on='id')

print('class data appended')

# ================================VARIABLES=====================================
# make df with all vars, clean all vars (X), append class col (Y)

data = bpm_4k
data = data.rename(columns={'og_id': 'src_subject_id', 'sex': 'SEX', 'time': 'eventname'})
print(data) # this should include mapped numeric and og id's and class enumeration

# =================== add bpm score ===================

# change to var storage area
os.chdir('/Volumes/igmm/GenScotDepression/data/abcd/release4.0/iii.data/')
# all files containing VOIs
rds = ['abcd_yssbpm01.rds']

datasets = [data] # list of datasets with the base data in
count = 1 # 1 df already in list

# get all dataframes from rds list and append using os.walk
for root, dirs, files in os.walk(".", topdown=False):
   for file in files:
      if file in rds:
          df = pyreadr.read_r(os.path.join(root, file))
          df = df[None]  # changes from dict
          datasets.append(df)
          count = count + 1

# to prevent duplicate cols
for df in datasets:
    df.drop(['interview_date','interview_age','sex'], axis=1, inplace=True, errors='ignore')
    df['eventname'] = df['eventname'].replace(['baseline_year_1_arm_1', '6_month_follow_up_arm_1', '1_year_follow_up_y_arm_1',
                                               '18_month_follow_up_arm_1', '2_year_follow_up_y_arm_1', '30_month_follow_up_arm_1',
                                               '3_year_follow_up_y_arm_1'], [0, 0.5, 1, 1.5, 2, 2.5, 3],)
    df['eventname']=df['eventname'].astype(float) # all floats for merging
    df['src_subject_id']=df['src_subject_id'].astype(object) # all objects for merging

print('done - datasets is a list of dataframes')

# merging all datasets on id and time point
all = reduce(lambda left,right: pd.merge(left,right,on=['src_subject_id','eventname'], how='outer'), datasets)

# filling NaNs in data cols
# Create a dictionary for each
id_class_dict = all.dropna(subset=['class']).set_index('src_subject_id')['class'].to_dict()
id_sex_dict = all.dropna(subset=['SEX']).set_index('src_subject_id')['SEX'].to_dict()
id_event_age_dict = all.dropna(subset=['age']).set_index(['src_subject_id', 'eventname'])['age'].to_dict() # age changes with time

# Use the dictionary to fill NaN values in specified column
all['class'] = all['class'].fillna(all['src_subject_id'].map(id_class_dict))
all['SEX'] = all['SEX'].fillna(all['src_subject_id'].map(id_sex_dict))
all['age'] = all.apply(lambda row: id_event_age_dict.get((row['src_subject_id'], row['eventname']), row['age']), axis=1)

print('all is a single df')

# remove the numeric id col used in mplus
all = all.drop(columns=['id'])
column_to_move = all.pop("src_subject_id") # moving id to first col
data = all.copy() # to prevent fragmented dataframe
data.insert(0, "src_subject_id", column_to_move) # insert column with insert(location, column_name, column_value)

df_bkup = data

data = data.rename(columns={'src_subject_id':'IID', 'bpm_y_scr_internal_r':'bpm'})

print('data is the current working df')

# ================================GENETIC VARS=====================================

# change to var storage area
os.chdir('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/prs_bpm_OUT/t2')
# all files containing PRSs [change from alspac to abcd]
csvs = ['abcd_anx_prs_0320.best','abcd_neu_prs_0320.best',
        'abcd_mdd_prs_0320.best','abcd_scz_prs_0320.best',
        'abcd_asd_prs_0320.best','abcd_bip_prs_0320.best',
        'abcd_adhd_prs_0320.best']


# iterate over each file
for csv_file in csvs:
    df = pd.read_csv(csv_file, sep='\s+')
    df = df[['IID', 'PRS']]
    # extract the prefix of the column name from the file name
    col_prefix = csv_file.split('_')[1]
    # rename the second column to the appropriate prefix
    df = df.rename(columns={'PRS': col_prefix + '_prs'})
    # merge the DataFrame with the 'als' DataFrame on the 'IID' column
    data = pd.merge(data, df, on='IID', how='left')

print('genetic data appended')

# extract vars we want
b_vars = ['SEX']
c_vars = ['bpm', 'age']
g_vars = ['scz_prs', 'neu_prs', 'mdd_prs','bip_prs', 'asd_prs', 'anx_prs', 'adhd_prs']
dem_vars = ['IID','class','eventname']
all_vars = dem_vars + b_vars + c_vars + g_vars

# ================================= MAKING DESIGN MATRIX ===================================

# useful code to make table
data['class'].fillna('', inplace=True)  # fill empty values in 'class' column with empty string
grouped = data.groupby(['eventname', 'class'])['class'].count().unstack(fill_value=0)
result = grouped.astype(int) # shows 0 NaNs at baseline measure
print(result)

# can idenitfy which class is which trajec
bpm_table = pd.pivot_table(data, values='bpm', index='eventname', columns='class', aggfunc='mean')
print(bpm_table)


# select only vars we want
data = data[all_vars]

# binary vars restrict to [0,1,NaN]
data[b_vars] = data[b_vars].mask(~data[b_vars].isin([0,1]))


# standardise PRS
for g_var in g_vars:
    data[g_var]  = preprocessing.scale(data[g_var])


# ============================= MN log reg =====================================

# 'scz_prs', 'neu_prs', 'mdd_prs','bip_prs', 'asd_prs', 'anx_prs', 'adhd_prs'

# Filter out rows with missing values in 'prs_mdd' and 'class' columns
df2 = data.dropna(subset=['anx_prs', 'class'])
df3 = df2.drop_duplicates(subset=["IID"])  # as data is in long format

# class 3 is stable low for bpm > change to 0 for reference
# c1 high, c3 decreasing, c4 increasing

df3.loc[:, 'class'] = df3['class'].replace(3.0, 0.0)  # replace with 0 for ref level
x = df3['anx_prs']
y = df3['class']
X = sm.add_constant(x)
model = sm.MNLogit(y, X)
result = model.fit()

# Format the coefficients and standard errors in exponential notation
params_exp = np.exp(result.params)
conf_int_exp = np.exp(result.conf_int())
print(params_exp) # result of this is c1,c2,c4,blank class
print(conf_int_exp)