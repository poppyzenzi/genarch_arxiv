import numpy as np
import pandas as pd
import statsmodels.api as sm
import os.path
from sklearn import preprocessing
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D


# ================================ APPEND CLASS ENUMERATION =================================
os.chdir('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/gmm/gmm_abcd')

# read in class data from mplus output
bpm_4k = pd.read_csv('mplus_data/4k_gmm_bpm_v5_gen.txt', delim_whitespace=True, header=None,
                     names=['y0.5','y1','y1.5','y2','y2.5','y3','y3.5','y4',
                            'v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','class','SubjectNumeric'])
bpm_4k = bpm_4k[['SubjectNumeric','class']] # subset
bpm_4k = bpm_4k.replace('*', np.nan) # fill missing

# class data has unique id's only > merge with anth to align NDAR id's
demo = pd.read_csv("genetic_subjects_only/abcd_dep_long_gen_only.csv")
merged = pd.merge(bpm_4k, demo, on='SubjectNumeric') # append class enumeration
merged['SubjectNumeric'].nunique() # check unique id's are 8016
print('class data merged')

# ================================ VARIABLES =====================================

data = merged.rename(columns={'src_subject_id':'IID', 'bpm_sum':'bpm'})

bpm_means = data.groupby('eventname')['bpm'].mean()
bpm_std = data.groupby('eventname')['bpm'].std()

print(bpm_means, bpm_std)

# ================================ GENETIC VARS =====================================

os.chdir('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/prs_bpm_OUT/all_thresholds')
csvs = ['abcd_adhd_23_prs_0813.t2.best',
        'abcd_asd_prs_0813.t2.best','abcd_bip_prs_0814.t1.best',
        'abcd_mddipsych_prs_0831.t2.5.best','abcd_meta_anx_prs_0813.t1.best',
        'abcd_neu_prs_0813.t1.best', 'abcd_scz_prs_0814.t3.best',
        'bpm_cf_prs_0817.t3.best', 'bpm_high_prs_0817.t2.best',
        'ABCD3.0_AFR_mdd_prs_230608.best',
        'ABCD3.0_AMR_mdd_prs_230608.best', 'ABCD3.0_EAS_mdd_prs_230609.best']

# iterate over each file
for csv_file in csvs:
    df = pd.read_csv(csv_file, sep='\s+')
    df = df[['IID', 'PRS']]
    col_prefix = csv_file.split('_')[1]  # extract the prefix of the column name from the file name
    df = df.rename(columns={'PRS': col_prefix + '_prs'})  # rename the second column to the appropriate prefix
    data = pd.merge(data, df, on='IID', how='left')  # merge the DataFrame with the 'data' DataFrame on the 'IID' column

print('genetic data appended')

# useful code to get genetic info counts

afr_count = data.drop_duplicates(subset='IID').dropna(subset=['AFR_prs'])

data = data.rename(columns={'prs0501.best_prs':'mood_prs', 'prs0511.best':'high_prs',
                            'meta_prs':'anx_prs', 'cf_prs':'common', 'high_prs':'hierarchical', 'mddipsych_prs':'mdd_prs'}) # rename parcel cols

# extract vars
b_vars = ['sex']
c_vars = ['bpm', 'interview_age']
g_vars = ['common', 'hierarchical', 'mdd_prs', 'neu_prs', 'anx_prs', 'scz_prs',
                  'bip_prs', 'adhd_prs', 'asd_prs']
dem_vars = ['IID','class','eventname']
all_vars = dem_vars + b_vars + c_vars + g_vars
# ================================= MAKING DESIGN MATRIX ===================================

# useful code to make table
grouped = data.groupby(['eventname', 'class'])['class'].count().unstack(fill_value=0)
result = grouped.astype(int) # shows 0 NaNs at baseline measure
print(result)

# check which class is which trajectory
bpm_table = pd.pivot_table(data, values='bpm', index='eventname', columns='class', aggfunc='mean')
print(bpm_table)

data = data[all_vars] # select only vars we want
# data[b_vars] = data[b_vars].mask(~data[b_vars].isin([0,1])) # restrict binary vars to [0,1,NaN]
for g_var in g_vars:
    data[g_var] = preprocessing.scale(data[g_var])  # standardise PRS

df = data.drop(columns=['eventname','bpm']) # remove bpm and time data as not needed and has dups
df = df.drop_duplicates(subset=["IID"])  # as data is in long format. Should be 8016 rows.

# ======== PRINCIPAL COMPONENTS ============
principal_components = pd.read_table('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/multiancestry/abcd.randomforest.ancestries.tsv')
pcs = principal_components[['IID', 'PC1_AVG', 'PC2_AVG', 'PC3_AVG', 'PC4_AVG', 'PC5_AVG','PC6_AVG']]
data_with_pcs = pd.merge(df, pcs, on='IID')

# =============================== CLASSIFICATION =====================================
# 1 = persistent, 2 = increasing, 3 = low, 4 = decreasing
data_with_pcs['class'] = data_with_pcs['class'].replace({3: 0, 4: 2, 2: 3, 1: 4})
# 0 = low, 2 = decreasing, 3 = increasing, 4 = Persistent

results_list = []
plot_result_list = []

covariate_columns = ['PC1_AVG', 'PC2_AVG', 'PC3_AVG', 'PC4_AVG', 'PC5_AVG', 'PC6_AVG']

for var in g_vars:
    df2 = data_with_pcs.dropna(subset=[var])  # Drop rows with NA's in var and class columns
    x = df2[[var] + ['sex'] + covariate_columns]
    y = df2['class']
    X = sm.add_constant(x)
    model = sm.MNLogit(y, X)
    result = model.fit()

    odds_ratios = np.exp(result.params.iloc[1])
    odds_ratios.columns = ['Decreasing', 'Increasing', 'Persistent']  # set index names
    conf_int = np.exp(result.conf_int().iloc[[1,10,19],:])# will need changing depending on #of covariates
    p_values = result.pvalues.iloc[1]
    p_values.columns = ['Decreasing', 'Increasing', 'Persistent']
    results_df = pd.DataFrame({'Variable': [var],
                               'Decreasing OR (95% CI)': f"{odds_ratios.iloc[0]:.2f} ({conf_int.iloc[0, 0]:.2f}, {conf_int.iloc[0, 1]:.2f})",
                               'Increasing OR (95% CI)': f"{odds_ratios.iloc[1]:.2f} ({conf_int.iloc[1, 0]:.2f}, {conf_int.iloc[1, 1]:.2f})",
                               'Persistent OR (95% CI)': f"{odds_ratios.iloc[2]:.2f} ({conf_int.iloc[2, 0]:.2f}, {conf_int.iloc[2, 1]:.2f})"
                               #'p-value (Persistent)': p_values.iloc[0],
                               #'p-value (Increasing)': p_values.iloc[1],
                               #'p-value (Decreasing)': p_values.iloc[2]
                               })
    plot_df = pd.DataFrame({'Variable': [var] * 3,
                            'Odds Ratio': odds_ratios.values.tolist(),
                            'Lower CI': conf_int.iloc[:, 0].values.tolist(),
                            'Upper CI': conf_int.iloc[:, 1].values.tolist(),
                            'Category': ['Decreasing', 'Increasing', 'Persistent']
                            })

    plot_result_list.append(plot_df)
    results_list.append(results_df)

results_table = pd.concat(results_list, ignore_index=True)
pd.set_option('display.max_columns', None)
print(results_table)

results_table.to_csv('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/logreg_results/abcd_results_table.csv', index=False)

odds = pd.concat(plot_result_list, ignore_index=True)
odds = odds.set_index('Variable')
groups = odds.groupby('Variable')

odds.to_csv('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/logreg_results/abcd_plot_results.csv')

fig, ax = plt.subplots(figsize=(6, 8))

width = 0.06
y_pos = 0.5
variable_order = ['common', 'hierarchical', 'mdd_prs', 'neu_prs', 'anx_prs', 'scz_prs',
                  'bip_prs', 'adhd_prs', 'asd_prs']

markers = ['o', 's', 'X', 'D', '>', 'P', '^', 'p', '*']

def sorting_key(name_group):
    name, group = name_group
    return variable_order.index(name)

sorted_groups = sorted(groups, key=sorting_key)  # Reverse the order of the groups

variables = {
    'common': {
        'legend_text': 'Common (PT 0.01)',
        'marker': 'o',
        'color': 'moccasin'
    },
    'hierarchical': {
        'legend_text': 'Hierarchical (PT 0.01)',
        'marker': 's',
        'color': 'orange'
    },
    'mdd_prs': {
        'legend_text': 'MDD (PT 0.05)',
        'marker': 'X',
        'color': 'darkblue'
    },
    'neu_prs': {
        'legend_text': 'NEU (PT 0.2)',
        'marker': 'D',
        'color': 'blue'
    },
    'anx_prs': {
        'legend_text': 'ANX (PT 0.05)',
        'marker': '>',
        'color': 'dodgerblue'
    },
    'scz_prs': {
        'legend_text': 'SCZ (PT 0.5)',
        'marker': 'P',
        'color': 'darkslategrey'
    },
    'bip_prs': {
        'legend_text': 'BIP (PT 0.0001)',
        'marker': '^',
        'color': 'slategrey'
    },
    'adhd_prs': {
        'legend_text': 'ADHD (PT 0.2)',
        'marker': 'p',
        'color': 'deeppink'
    },
    'asd_prs': {
        'legend_text': 'ASD (PT 0.001)',
        'marker': '*',
        'color': 'lightpink'
    }
}

legend_handles = []  # Create an empty list to store legend patches

for i, (name, group) in enumerate(sorted_groups):
    if i < 2: # Cluster the bars 1, 2, 3 together
        x_pos = np.arange(len(group.index)) + i * width
    else:  # Add white space between the third and fourth error bars
        x_pos = np.arange(len(group.index)) + (i + 2) * width

    ax.errorbar(
        group['Odds Ratio'], x_pos + width / 2,
        xerr=[group['Odds Ratio'] - group['Lower CI'], group['Upper CI'] - group['Odds Ratio']],
        fmt=variables[name]['marker'],
        capsize=0.3, capthick=0.8, linewidth=0.8, color=variables[name]['color']
    )

    legend_handles.append(
        Line2D([0], [0], marker=variables[name]['marker'], color='w',
               markerfacecolor=variables[name]['color'], markersize=8, label=variables[name]['legend_text'])
    )

plt.subplots_adjust(left=0.15, bottom=0.08, right=0.8, top=0.95)
categories = list(odds['Category'].unique())
ax.set_yticklabels(odds['Category'].unique(), fontsize=8)
ax.set_yticks(np.arange(len(categories)) + y_pos)
ax.set_xlabel('Odds Ratio (95% CI)', fontsize=8)
ax.set_xscale('log')

ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0))
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{x:.1f}"))
ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
ax.xaxis.set_minor_formatter(ticker.ScalarFormatter())
ax.tick_params(axis='x', which='both', labelsize=7)

ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=7)

ax.axvline(x=1, linestyle='--', color='black', linewidth=0.8)
plt.show()