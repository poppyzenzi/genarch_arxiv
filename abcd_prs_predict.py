import numpy as np
import pandas as pd
import statsmodels.api as sm
import os.path
from sklearn import preprocessing

# ================================ APPEND CLASS ENUMERATION =================================

# read in class data from mplus output
bpm_4k = pd.read_table('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/gmm/gmm_abcd/'
                      'mplus_data/4k_gmm_bpm.txt', delim_whitespace=True, header=None)
bpm_4k.columns = ['y0.5','y1','y1.5','y2','y2.5','y3','v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','class','id']
bpm_4k = bpm_4k[['id','class']] # subset
bpm_4k.replace('*', np.nan) # fill missing

# class data has unique id's only > merge with anth to align NDAR id's
anth = pd.read_csv("/Volumes/igmm/GenScotDepression/users/poppy/abcd/abcd_anth.csv") # 11,733 nunique

data = pd.merge(bpm_4k, anth, on='id') # append class enumeration
data['id'].nunique() # check unique id's
print('class data merged')

# ================================ VARIABLES =====================================

data = data.rename(columns={'og_id': 'src_subject_id', 'sex': 'SEX', 'time': 'eventname'})
print(data) # this should include mapped numeric and og id's and class enumeration

df = data.drop(columns=['id']) # remove the numeric id col used in mplus
column_to_move = df.pop("src_subject_id") # moving id to first col
data = df.copy() # to prevent fragmented dataframe
data.insert(0, "src_subject_id", column_to_move) # insert column with insert(location, column_name, column_value)
data = data.rename(columns={'src_subject_id':'IID', 'bpm_y_scr_internal_r':'bpm'})

# ================================ GENETIC VARS =====================================

os.chdir('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/prs_bpm_OUT/t2')
csvs = ['abcd_anx_prs_0320.best','abcd_neu_prs_0320.best',
        'abcd_mdd_prs_0320.best','abcd_scz_prs_0320.best',
        'abcd_asd_prs_0320.best','abcd_bip_prs_0320.best',
        'abcd_adhd_prs_0320.best', 'abcd_meta_anx_prs_0405.best',
        'mood_prs0501.best', 'bpm_cf_prs_0418.best']

# iterate over each file
for csv_file in csvs:
    df = pd.read_csv(csv_file, sep='\s+')
    df = df[['IID', 'PRS']]
    col_prefix = csv_file.split('_')[1]  # extract the prefix of the column name from the file name
    df = df.rename(columns={'PRS': col_prefix + '_prs'})  # rename the second column to the appropriate prefix
    data = pd.merge(data, df, on='IID', how='left')  # merge the DataFrame with the 'data' DataFrame on the 'IID' column

data = data.sort_values('IID') # sort the data by the 'IID' column
print('genetic data appended')

data = data.rename(columns={'prs0501.best_prs':'mood_prs'}) # rename parcel cols

# extract vars
b_vars = []
c_vars = ['int_r', 'age']
g_vars = ['scz_prs', 'neu_prs', 'mdd_prs','bip_prs', 'asd_prs', 'anx_prs', 'adhd_prs', 'meta_prs', 'mood_prs', 'cf_prs']
dem_vars = ['IID','class','eventname']
all_vars = dem_vars + b_vars + c_vars + g_vars

# ================================= MAKING DESIGN MATRIX ===================================

# useful code to make table
grouped = data.groupby(['eventname', 'class'])['class'].count().unstack(fill_value=0)
result = grouped.astype(int) # shows 0 NaNs at baseline measure
print(result)

# check which class is which trajectory
bpm_table = pd.pivot_table(data, values='int_r', index='eventname', columns='class', aggfunc='mean')
print(bpm_table)

data = data[all_vars] # select only vars we want
data[b_vars] = data[b_vars].mask(~data[b_vars].isin([0,1])) # restrict binary vars to [0,1,NaN]
for g_var in g_vars:
    data[g_var] = preprocessing.scale(data[g_var])  # standardise PRS

# =============================== CLASSIFICATION =====================================

df2 = data.dropna(subset=['cf_prs', 'class']) # Filter out rows with missing values in vat and class cols
df3 = df2.drop_duplicates(subset=["IID"])  # as data is in long format
df3.loc[:, 'class'] = df3['class'].replace(1.0, 0.0)  # change reference class (non-depressed) to 0

# fit mn logistic regression model
x = df3['cf_prs']
y = df3['class']
X = sm.add_constant(x)
model = sm.MNLogit(y, X)
result = model.fit()

# Format the coefficients and standard errors in exponential notation
params_exp = np.exp(result.params)
conf_int_exp = np.exp(result.conf_int())
print(params_exp) # result of this is c1,c2,c4,blank class
print(conf_int_exp)


############ linear regression of MDD PRS and BPM at 6 time points ############
'''df = data
df = df.dropna(subset=['mdd_prs', 'bpm'])
event_names = df['eventname'].unique()
regression_results = {}
# Iterate over the unique eventnames and perform regression for each group
for event_name in event_names:
    df_group = df[df['eventname'] == event_name]
    X = df_group['mdd_prs']
    y = df_group['bpm']
    X = sm.add_constant(X)
    model = sm.OLS(y, X)
    results = model.fit()

    regression_results[event_name] = results

    print(f'Regression results for eventname: {event_name}')
    print(results.summary())
    print('---------------------------------------------------------------------')'''