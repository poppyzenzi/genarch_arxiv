import numpy as np
import pandas as pd
import statsmodels.api as sm
import os.path
from sklearn import preprocessing

# ================================ APPEND CLASS ENUMERATION =================================

# read in class data from mplus output
bpm_4k = pd.read_csv('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/gmm/gmm_abcd/mplus_data/4k_gmm_bpm.txt',
                     delim_whitespace=True, header=None, names=['y0.5','y1','y1.5','y2','y2.5','y3','v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','class','id'])
bpm_4k = bpm_4k[['id','class']] # subset
bpm_4k = bpm_4k.replace('*', np.nan) # fill missing


# class data has unique id's only > merge with anth to align NDAR id's
anth = pd.read_csv("/Volumes/igmm/GenScotDepression/users/poppy/abcd/abcd_anth.csv") # 11,733 nunique

merged = pd.merge(bpm_4k, anth, on='id') # append class enumeration
merged['id'].nunique() # check unique id's
print('class data merged')

# ================================ VARIABLES =====================================

abcd = merged.rename(columns={'og_id': 'src_subject_id', 'sex': 'SEX', 'time': 'eventname'})
df = abcd.drop(columns=['id']) # remove the numeric id col used in mplus
column_to_move = df.pop("src_subject_id") # moving id to first col
df.insert(0, "src_subject_id", column_to_move) # insert column with insert(location, column_name, column_value)
#data = abcd.copy() # to prevent fragmented dataframe
data = df.rename(columns={'src_subject_id':'IID', 'bpm_y_scr_internal_r':'bpm'})

bpm_means = data.groupby('eventname')['int_r'].std()
print(bpm_means)

# ================================ GENETIC VARS =====================================

os.chdir('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/prs_bpm_OUT/t2')
csvs = ['abcd_anx_prs_0320.best','abcd_neu_prs_0320.best',
        'abcd_mdd_prs_0320.best','abcd_scz_prs_0320.best',
        'abcd_asd_prs_0320.best','abcd_bip_prs_0320.best',
        'abcd_adhd_prs_0320.best', 'abcd_meta_anx_prs_0405.best',
        'mood_prs0501.best', 'bpm_cf_prs_0418.best',
        'abcd_adhd_23_prs_0512.best']

# iterate over each file
for csv_file in csvs:
    df = pd.read_csv(csv_file, sep='\s+')
    df = df[['IID', 'PRS']]
    col_prefix = csv_file.split('_')[1]  # extract the prefix of the column name from the file name
    df = df.rename(columns={'PRS': col_prefix + '_prs'})  # rename the second column to the appropriate prefix
    data = pd.merge(data, df, on='IID', how='left')  # merge the DataFrame with the 'data' DataFrame on the 'IID' column

data = data.sort_values('IID') # sort the data by the 'IID' column
print('genetic data appended')

data = data.rename(columns={'prs0501.best_prs':'mood_prs', 'prs0511.best':'high_prs', 'adhd_prs_y':'adhd23_prs'}) # rename parcel cols

# extract vars
b_vars = []
c_vars = ['int_r', 'age']
g_vars = ['scz_prs', 'neu_prs', 'mdd_prs','bip_prs', 'asd_prs',
          'meta_prs', 'mood_prs', 'cf_prs', 'adhd23_prs']
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

df = data.drop(columns=['eventname','int_r']) # remove the numeric id col used in mplus
df = df.drop_duplicates(subset=["IID"])  # as data is in long format

#data['high_prs'].nunique() # check how many have risk scores, can perform checks with .best file length

# =============================== CLASSIFICATION =====================================
df.loc[:, 'class'] = df['class'].replace(1.0, 0.0)  # change reference class (non-depressed) to 0

results_list = []

for var in g_vars:
    df = df.dropna(subset=[var, 'class'])     # drop rows with NA's in the variable and class columns
    x = df[var]
    y = df['class']
    X = sm.add_constant(x)
    model = sm.MNLogit(y, X)
    result = model.fit()

    odds_ratios = np.exp(result.params.iloc[1])
    odds_ratios.columns = ['Acute', 'Increasing', 'Decreasing']  # set index names
    conf_int = np.exp(result.conf_int().iloc[[1,3,5],:])
    p_values = result.pvalues.iloc[1]
    p_values.columns = ['Acute', 'Increasing', 'Decreasing']
    results_df = pd.DataFrame({'Variable': [var],
                               #'Class': result.model.endog_names,
                               'Acute OR (95% CI)': f"{odds_ratios.iloc[0]:.2f} ({conf_int.iloc[0, 0]:.2f}, {conf_int.iloc[0, 1]:.2f})",
                               'Increasing OR (95% CI)': f"{odds_ratios.iloc[1]:.2f} ({conf_int.iloc[1, 0]:.2f}, {conf_int.iloc[1, 1]:.2f})",
                               'Decreasing OR (95% CI)': f"{odds_ratios.iloc[2]:.2f} ({conf_int.iloc[2, 0]:.2f}, {conf_int.iloc[2, 1]:.2f})"
                               #'p-value (Acute)': p_values.iloc[0],
                               #'p-value (Increasing)': p_values.iloc[1],
                               #'p-value (Decreasing)': p_values.iloc[2]
                               })

    results_list.append(results_df)

results_table = pd.concat(results_list, ignore_index=True)
pd.set_option('display.max_columns', None)
print(results_table)
