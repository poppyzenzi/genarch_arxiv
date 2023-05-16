import numpy as np
import os
import pandas as pd
import pyreadr
import statsmodels.api as sm
from sklearn import preprocessing
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

# select only cols we want, re-run from here in console if needed
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

als.to_csv("/Volumes/igmm/GenScotDepression/users/poppy/alspac/alspac_demog.csv", header=True, sep='\t')

# ============================ GENETIC SCORES ===================================

als = pd.read_csv("/Volumes/igmm/GenScotDepression/users/poppy/alspac/alspac_demog.csv",  sep='\t')

os.chdir('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/prs_alspac_OUT')
csvs = ['alspac_scz_prs_0317.best','alspac_neu_prs_0316.best','alspac_mdd_prs_0320.best','alspac_bip_prs_0317.best',
       'alspac_asd_prs_0317.best', 'alspac_meta_anx_prs_0405.best',
        'alspac_cf_prs_0418.best', 'mood_prs0501.best', 'alspac_adhd_23_prs_0512.best', 'alspac_highfac_prs0516.best',
        'alspac_neurodev_corr_prs0601.best']

# iterate over each file
for csv_file in csvs:
    df = pd.read_csv(csv_file, sep='\s+')
    df = df[['IID', 'PRS']]
    col_prefix = csv_file.split('_')[1]  # extract the prefix of the column name from the file name
    df = df.rename(columns={'PRS': col_prefix + '_prs'})  # rename the second column to the appropriate prefix
    als = pd.merge(als, df, on='IID', how='left') # merge on IID

print(als.head())

als = als.rename(columns={'prs0501.best_prs':'mood_prs',
                          'meta_prs':'anx_prs', 'cf_prs':'common_fac',
                          'highfac_prs':'high_fac'}) # rename parcel cols

# ==============================================================================================

# merge als with alspac_4k by id [but need to make sure these IDs are the same
alspac_4k['id'] = alspac_4k['id'].astype(int) # make ids same object type (int)
df = pd.merge(alspac_4k, als, on=["id"])


#smfq_means = df.groupby('time')['age'].mean()

# select and clean variables
dem_vars = ['id', 'IID', 'ethnicity', 'class'] # demographic
g_vars = ['scz_prs', 'neu_prs', 'mdd_prs','bip_prs', 'neurodev_prs', 'asd_prs', 'anx_prs',
          'common_fac', 'mood_prs', 'adhd_prs', 'high_fac'] # genetic
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

# =============================== CLASSIFICATION =====================================
#df.loc[:, 'class'] = df['class'].replace(2.0, 0.0)  # change reference class (non-depressed) 2.0 to 0
df['class'] = df['class'].replace({2: 0, 1: 3, 3: 4, 4: 2})

results_list = []
plot_result_list = []

for var in g_vars:
    df = df.dropna(subset=[var, 'class'])     # drop rows with NA's in the variable and class columns
    x = df[var]
    y = df['class']
    X = sm.add_constant(x)
    model = sm.MNLogit(y, X)
    result = model.fit()

    odds_ratios = np.exp(result.params.iloc[1])
    odds_ratios.columns = ['Decreasing', 'Increasing', 'Acute']  # set index names
    conf_int = np.exp(result.conf_int().iloc[[1,3,5],:])
    p_values = result.pvalues.iloc[1]
    p_values.columns = ['Decreasing', 'Increasing', 'Acute']
    results_df = pd.DataFrame({'Variable': [var],
                               #'Class': result.model.endog_names,
                               'Decreasing OR (95% CI)': f"{odds_ratios.iloc[0]:.2f} ({conf_int.iloc[0, 0]:.2f}, {conf_int.iloc[0, 1]:.2f})",
                               'Increasing OR (95% CI)': f"{odds_ratios.iloc[1]:.2f} ({conf_int.iloc[1, 0]:.2f}, {conf_int.iloc[1, 1]:.2f})",
                               'Acute OR (95% CI)': f"{odds_ratios.iloc[2]:.2f} ({conf_int.iloc[2, 0]:.2f}, {conf_int.iloc[2, 1]:.2f})"
                               #'p-value (Increasing)': p_values.iloc[0],
                               #'p-value (Acute)': p_values.iloc[1],
                               #'p-value (Decreasing)': p_values.iloc[2]
                               })
    plot_df = pd.DataFrame({'Variable': [var] * 3,
                            'Odds Ratio': odds_ratios.values.tolist(),
                            'Lower CI': conf_int.iloc[:, 0].values.tolist(),
                            'Upper CI': conf_int.iloc[:, 1].values.tolist(),
                            'Category': ['Decreasing', 'Increasing', 'Acute']
                            })

    plot_result_list.append(plot_df)
    results_list.append(results_df)

results_table = pd.concat(results_list, ignore_index=True)
pd.set_option('display.max_columns', None)
print(results_table)

odds = pd.concat(plot_result_list, ignore_index=True)
#odds = odds.loc[odds['Variable'].isin(['mdd_prs', 'cf_prs', 'meta_prs', 'adhd23_prs', 'asd_prs'])]
odds = odds.set_index('Variable')
groups = odds.groupby('Variable')

fig, ax = plt.subplots(figsize=(8, 10))
width = 0.06
y_pos = 0.5
variable_order = ['common_fac', 'high_fac', 'mood_prs', 'mdd_prs', 'neu_prs', 'anx_prs', 'scz_prs','bip_prs', 'neurodev_prs', 'adhd_prs', 'asd_prs']
green_cm = plt.cm.get_cmap('jet', len(odds.index.unique()))
green_cm = green_cm.reversed()
purple_cm = plt.cm.get_cmap('viridis', len(odds.index.unique()))
markers = ['o', 's', 'D', '>', 'P', '^', 'p', '*', 'X', 'h', '+']

def sorting_key(name_group):
    name, group = name_group
    return variable_order.index(name)

sorted_groups = sorted(groups, key=sorting_key) # Sort the groups based on the custom sorting key function

for i, (name, group) in enumerate(sorted_groups):
    if i < 3:  # Cluster the bars 1, 2, 3 together
        x_pos = np.arange(len(group.index)) + i * width
    else:  # Add white space between the third and fourth error bars
        x_pos = np.arange(len(group.index)) + (i + 3) * width
    ax.errorbar(group['Odds Ratio'], x_pos+width/2, xerr=[group['Odds Ratio']-group['Lower CI'],     # plot the odds ratios with error bars
                                                  group['Upper CI']-group['Odds Ratio']], fmt=markers[i],
                                                capsize=0.3,capthick=0.8, linewidth=0.8,
                color = green_cm(i / len(sorted_groups)) if i < len(sorted_groups) - 7 else purple_cm((i - (len(sorted_groups) - 7)) / 7)
                )


plt.subplots_adjust(left=0.15, bottom=0.08, right=0.95, top=0.95) # adjust the margins
categories = list(odds['Category'].unique()) # get the unique categories
ax.set_yticklabels(odds['Category'].unique(), fontsize=8)
ax.set_yticks(np.arange(len(categories))+y_pos)
ax.set_xlabel('95% CI', fontsize=8)
ax.set_xlim(0.8, max(odds['Upper CI'])+0.1)
legend_labels = [name for name, _ in sorted_groups]
ax.legend(legend_labels)
ax.axvline(x=1, linestyle='--', color='black', linewidth=0.8)
plt.rcParams['legend.fontsize'] = 8
plt.show()