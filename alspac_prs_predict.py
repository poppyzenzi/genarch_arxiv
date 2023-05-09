import numpy as np
import os
import pandas as pd
import pyreadr
import statsmodels.api as sm
from sklearn import preprocessing

# ================================ APPEND CLASS ENUMERATION =================================

# class data from mplus
alspac_4k = pd.read_table('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/gmm/gmm_alspac/mplus_data/4k_gmm_smfq.txt', delim_whitespace=True, header=None)  # this is just test will need to change
alspac_4k.columns = ['y0','y1','y2','y3','v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','class','id']
alspac_4k = alspac_4k[['id','class']] # subsetting
alspac_4k.replace('*', np.nan) # replace missing

# smfq data - scores calculated in R (alspac_smfq_scores.Rmd)
df = pyreadr.read_r('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/gmm/gmm_alspac/alspac_dep_long.rds')
df = df[None]
df['id'].nunique() # check 15,645
df = df[['id','time','dep']]

data = pd.merge(df, alspac_4k, on=["id"]) # merge class with smfq data
data = data[data["time"] < 5] # exclude after 4 time points as this is what model was built on

print('class data merged')

# =================================== EXTRACTING VARS ===================================

# takes a while to read in
alspac = pd.read_stata("/Volumes/cmvm/scs/groups/ALSPAC/data/B3421_Whalley_04Nov2021.dta")

# select only cols we want
als_cols = ['cidB3421','qlet','kz021','c804']
als = alspac[als_cols]

# make alspac have same ids as originally coded in R
als['cidB3421'] = als['cidB3421'].astype(int)
als['IID'] = als['cidB3421'].astype(str).str.cat(als['qlet'], sep='')
als['id'] = als['cidB3421'].astype(str) + als['qlet']
als = als.drop(['cidB3421', 'qlet'], axis=1)
als['id'] = pd.factorize(als['id'])[0] + 1  # makes ids unique and numeric, should be 15,645
als = als.rename(columns={'kz021':'sex', 'c804':'ethnicity'}) # rename some cols
column_to_move = als.pop("id") # moving id to first col
als.insert(0, "id", column_to_move) # insert column with insert(location, column_name, column_value)

print('alspac data read in as als')

# ============================ GENETIC SCORES ===================================

os.chdir('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/prs_alspac_OUT')
csvs = ['alspac_scz_prs_0317.best','alspac_neu_prs_0316.best','alspac_mdd_prs_0320.best','alspac_bip_prs_0317.best',
       'alspac_asd_prs_0317.best','alspac_anx_prs_0316.best', 'alspac_adhd_prs_0317.best', 'alspac_meta_anx_prs_0405.best',
        'alspac_cf_prs_0418.best', 'mood_prs0501.best']

# iterate over each file
for csv_file in csvs:
    df = pd.read_csv(csv_file, sep='\s+')
    df = df[['IID', 'PRS']]
    col_prefix = csv_file.split('_')[1]  # extract the prefix of the column name from the file name
    df = df.rename(columns={'PRS': col_prefix + '_prs'})  # rename the second column to the appropriate prefix
    als = pd.merge(als, df, on='IID', how='left') # merge on IID

print(als.head())

als = als.rename(columns={'prs0501.best_prs':'mood_prs'}) # rename parcel cols

# ==============================================================================================

# merge als with alspac_4k by id [but need to make sure these IDs are the same
alspac_4k['id'] = alspac_4k['id'].astype(int) # make ids same object type (int)
df = pd.merge(alspac_4k, als, on=["id"])

# select and clean variables
dem_vars = ['id', 'IID', 'ethnicity', 'class'] # demographic
g_vars = ['scz_prs', 'neu_prs', 'mdd_prs','bip_prs', 'asd_prs', 'anx_prs', 'adhd_prs', 'meta_prs','cf_prs', 'mood_prs'] # genetic
b_vars = ['sex'] # binary
all_vars = dem_vars + g_vars + b_vars
X_vars = g_vars + b_vars

# ================================= MAKING DESIGN MATRIX ===================================

df = df[all_vars] # select only vars we want
df['sex'] = df['sex'].replace(['Female','Male'], [1,0]) # first recode sex
df[b_vars] = df[b_vars].mask(~df[b_vars].isin([0,1])) # binary vars restrict to [0,1,NaN]

for g_var in g_vars:
    df[g_var]  = preprocessing.scale(df[g_var])  # standardise PRS

df['mood_prs'].nunique()

# ============================= MN log reg =====================================

df2 = df.dropna(subset=['neu_prs', 'class']) # Filter out rows with missing vals in prs and class
df3 = df2.drop_duplicates(subset=["id"]) # as data is in long format
df3.loc[:, 'class'] = df3['class'].replace(2.0, 0.0) # replace with 0 for ref level

# fit model
x = df3['neu_prs']
y = df3['class']
X = sm.add_constant(x)
model = sm.MNLogit(y, X)
result = model.fit()

# Format the coefficients and standard errors in exponential notation
params_exp = np.exp(result.params)
conf_int_exp = np.exp(result.conf_int())
print(params_exp)
print(conf_int_exp)