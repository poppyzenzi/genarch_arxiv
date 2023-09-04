import numpy as np
import os
import pandas as pd
import statsmodels.api as sm
from sklearn import preprocessing
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D


# ================================ APPEND CLASS ENUMERATION =================================
os.chdir('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/gmm/gmm_alspac')
# class data from mplus
alspac_4k = pd.read_table('mplus_data/4k_gmm_smfq_gen_only.txt', delim_whitespace=True, header=None)  # this is just test will need to change
alspac_4k.columns = ['y0','y1','y2','y3','v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','class','SubjectNumeric']
alspac_4k = alspac_4k[['SubjectNumeric','class']] # subsetting
alspac_4k.replace('*', np.nan) # replace missing

# smfq data - scores calculated in R (alspac_smfq_scores.Rmd)
df = pd.read_csv("genetic_subjects_only/alspac_dep_long_gen_only.csv")
df['IID'].nunique() # check 7862
df = df[['IID','SubjectNumeric','time','dep', 'sex', 'ethnicity', 'age']]
class_data = pd.merge(df, alspac_4k, on=["SubjectNumeric"]) # merge class with smfq data
print('class data merged with IIDs')
# should = 6096 for genetic data (x4 time points = 24384)

# ============================ GENETIC SCORES ===================================

os.chdir('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/prs_alspac_OUT/all_thresholds')

csvs = ['alspac_adhd_23_prs_0813.t3.best','alspac_asd_prs_0813.t3.best', 'alspac_bip_prs_0813.t2.best',
        'alspac_mdd_ipsych_prs_0816.t2.best', 'alspac_meta_anx_prs_0813.t3.best', 'alspac_neu_prs_0813.t3.best',
        'alspac_scz_prs_0813.t3.best', 'alspac_cf_prs_0817.t2.best', 'alspac_high_prs_0817.t2.best']

# iterate over each file
for csv_file in csvs:
    df = pd.read_csv(csv_file, sep='\s+')
    df = df[['IID', 'PRS']]
    col_prefix = csv_file.split('_')[1]  # extract the prefix of the column name from the file name
    df = df.rename(columns={'PRS': col_prefix + '_prs'})  # rename the second column to the appropriate prefix
    class_data = pd.merge(class_data, df, on='IID', how='left')

# rename PRS cols
allData = class_data.rename(columns={'meta_prs':'anx_prs', 'cf_prs':'common','high_prs':'hierarchical'})

# ==============================================================================================

# get class counts
prop = allData.drop_duplicates(['IID']).groupby('class').size()
print(prop)
#smfq_means = df.groupby('time')['age'].mean()

# check which class is which trajectory
smfq_table = pd.pivot_table(allData, values='dep', index='time', columns='class', aggfunc='mean')
print(smfq_table)

# ================================= MAKING DESIGN MATRIX ===================================

df = allData

#df[b_vars] = df[b_vars].mask(~df[b_vars].isin([0,1])) # binary vars restrict to [0,1,NaN]

g_vars = ['scz_prs', 'neu_prs', 'mdd_prs','bip_prs', 'asd_prs', 'anx_prs',
          'common', 'adhd_prs', 'hierarchical']

for g_var in g_vars:
    df[g_var]  = preprocessing.scale(df[g_var])  # standardise PRS

df = df.drop(columns=['time','dep']) # remove bpm and time data as not needed and has dups
df = df.drop_duplicates(subset=["IID"])  # as data is in long format. Should be 8016 rows.

# ===================== PRINCIPAL COMPONENTS ============================

principal_components = pd.read_csv('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/multiancestry/alspac_pca/alspac_pcs')
data_with_pcs = pd.merge(df, principal_components, on='IID')
covariate_columns = ['pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'pc6']

# =============================== CLASSIFICATION =====================================
# check using smfq_table
# 1 dec, 2 low, 3 persistent, 4 inc
data_with_pcs['class'] = data_with_pcs['class'].replace({2: 0, 1: 2, 4: 3, 3: 4})
# 0 = low, 2 = decreasing, 3 = increasing, 4 = persistent

results_list = []
plot_result_list = []

for var in g_vars:
    df = data_with_pcs.dropna(subset=[var, 'class', 'sex'])     # drop rows with NA's in the variable and class columns
    x = df[[var] + ['sex'] + covariate_columns]
    y = df['class']
    X = sm.add_constant(x)
    model = sm.MNLogit(y, X)
    result = model.fit()

    odds_ratios = np.exp(result.params.iloc[1])
    odds_ratios.columns = ['Decreasing', 'Increasing', 'Persistent']  # set index names
    conf_int = np.exp(result.conf_int().iloc[[1,10,19],:]) # will need changing depending on #of covariates
    p_values = result.pvalues.iloc[1]
    p_values.columns = ['Decreasing', 'Increasing', 'Persistent']
    results_df = pd.DataFrame({'Variable': [var],
                               #'Class': result.model.endog_names,
                               'Decreasing OR (95% CI)': f"{odds_ratios.iloc[0]:.2f} ({conf_int.iloc[0, 0]:.2f}, {conf_int.iloc[0, 1]:.2f})",
                               'Increasing OR (95% CI)': f"{odds_ratios.iloc[1]:.2f} ({conf_int.iloc[1, 0]:.2f}, {conf_int.iloc[1, 1]:.2f})",
                               'Persistent OR (95% CI)': f"{odds_ratios.iloc[2]:.2f} ({conf_int.iloc[2, 0]:.2f}, {conf_int.iloc[2, 1]:.2f})"
                               #'p-value (Increasing)': p_values.iloc[0],
                               #'p-value (Acute)': p_values.iloc[1],
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

results_table.to_csv('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/logreg_results/alspac_results_table.csv', index=False)

odds = pd.concat(plot_result_list, ignore_index=True)
#odds = odds.loc[odds['Variable'].isin(['mdd_prs', 'cf_prs', 'meta_prs', 'adhd23_prs', 'asd_prs'])]
odds = odds.set_index('Variable')
groups = odds.groupby('Variable')

odds.to_csv('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/logreg_results/alspac_plot_results.csv')

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