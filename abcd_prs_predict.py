import numpy as np
import pandas as pd
import statsmodels.api as sm
import os.path
from sklearn import preprocessing
import matplotlib.pyplot as plt

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
csvs = ['abcd_neu_prs_0320.best',
        'abcd_mdd_prs_0320.best','abcd_scz_prs_0320.best',
        'abcd_asd_prs_0320.best','abcd_bip_prs_0320.best',
        'abcd_meta_anx_prs_0405.best',
        'mood_prs0501.best', 'bpm_cf_prs_0418.best',
        'abcd_adhd_23_prs_0512.best', 'abcd_highfac_prs_0516.best',
        'ABCD3.0_imputed_blackonly_mdd_prs_230529.best', 'abcd_neurodev_corr_prs0529.best']

# iterate over each file
for csv_file in csvs:
    df = pd.read_csv(csv_file, sep='\s+')
    df = df[['IID', 'PRS']]
    col_prefix = csv_file.split('_')[1]  # extract the prefix of the column name from the file name
    df = df.rename(columns={'PRS': col_prefix + '_prs'})  # rename the second column to the appropriate prefix
    data = pd.merge(data, df, on='IID', how='left')  # merge the DataFrame with the 'data' DataFrame on the 'IID' column

data = data.sort_values('IID') # sort the data by the 'IID' column
print('genetic data appended')

data = data.rename(columns={'prs0501.best_prs':'mood_prs', 'prs0511.best':'high_prs',
                            'meta_prs':'anx_prs', 'cf_prs':'common_fac', 'highfac_prs':'high_fac', 'imputed_prs':'mdd_black'}) # rename parcel cols

# extract vars
b_vars = []
c_vars = ['int_r', 'age']
g_vars = ['common_fac', 'high_fac', 'mood_prs', 'mdd_prs', 'neu_prs', 'anx_prs', 'scz_prs',
                  'bip_prs', 'neurodev_prs', 'adhd_prs', 'asd_prs']
dem_vars = ['IID','class','eventname']
all_vars = dem_vars + b_vars + c_vars + g_vars
# mdd black not working for plotting?
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
#df.loc[:, 'class'] = df['class'].replace(1.0, 0.0)  # change reference class (non-depressed) to 0
# 1 = low, 2 = acute, 3 = increasing, 4 = decreasing
df['class'] = df['class'].replace({1: 0, 2: 4, 4: 2})
# 0 = low, 2 = decreasing, 3 = increasing, 4 = acute

results_list = []
plot_result_list = []

for var in g_vars:
    df.dropna(subset=[var, 'class'], inplace=True)  # Drop rows with NA's in var and class columns
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
                               #'p-value (Acute)': p_values.iloc[0],
                               #'p-value (Increasing)': p_values.iloc[1],
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
results_table.to_csv('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/logreg_results/results_table.csv', index=False)

odds = pd.concat(plot_result_list, ignore_index=True)
#odds = odds.loc[odds['Variable'].isin(['mdd_prs', 'cf_prs', 'meta_prs', 'adhd23_prs', 'asd_prs'])]
odds = odds.set_index('Variable')
groups = odds.groupby('Variable')

fig, ax = plt.subplots(figsize=(6, 8))
width = 0.06
y_pos = 0.5
variable_order = ['common_fac', 'high_fac', 'mood_prs', 'mdd_prs', 'neu_prs', 'anx_prs', 'scz_prs',
                  'bip_prs', 'neurodev_prs', 'adhd_prs', 'asd_prs']
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
ax.legend(legend_labels, loc='lower right')
#ax.legend(legend_labels, bbox_to_anchor=(1.05, 1), loc='upper right')
ax.axvline(x=1, linestyle='--', color='black', linewidth=0.8)
plt.rcParams['legend.fontsize'] = 5
plt.show()