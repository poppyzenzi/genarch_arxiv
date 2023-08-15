import numpy as np
import pandas as pd
import statsmodels.api as sm
import os.path
from sklearn import preprocessing
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ================================ APPEND CLASS ENUMERATION =================================

os.chdir('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/gmm/gmm_abcd')
# read in class data from mplus output
bpm_4k = pd.read_csv('mplus_data/4k_gmm_bpm_v5_gen.txt',
                     delim_whitespace=True, header=None,
                     names=['y0.5','y1','y1.5','y2','y2.5','y3','y3.5','y4','v1','v2','v3','v4','v5','v6','v7','v8',
                            'v9','v10','class','SubjectNumeric'])
bpm_4k = bpm_4k[['SubjectNumeric','class']] # subset
bpm_4k = bpm_4k.replace('*', np.nan) # fill missing

# merge with NDAR IDs
anth = pd.read_csv("genetic_subjects_only/abcd_dep_long_gen_only.csv")
merged = pd.merge(bpm_4k, anth, on='SubjectNumeric')
merged['SubjectNumeric'].nunique() # check unique id's
print('class data merged')

# ================================ VARIABLES =====================================

data = merged.rename(columns={'src_subject_id':'IID', 'bpm_sum':'bpm'})

bpm_means = data.groupby('eventname')['bpm'].mean()
print(bpm_means)

# ================================ GENETIC VARS =====================================
os.chdir('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/prs_bpm_OUT/t2')
csvs = ['abcd_mdd_prs_0320.best','ABCD3.0_AFR_mdd_prs_230608.best',
        'ABCD3.0_AMR_mdd_prs_230608.best', 'ABCD3.0_EAS_mdd_prs_230609.best']

# iterate over each file
for csv_file in csvs:
    df = pd.read_csv(csv_file, sep='\s+')
    df = df[['IID', 'PRS']]
    col_prefix = csv_file.split('_')[1]  # extract the prefix of the column name from the file name
    df = df.rename(columns={'PRS': col_prefix})  # rename the second column to the appropriate prefix
    data = pd.merge(data, df, on='IID', how='left')  # merge the DataFrame with the 'data' DataFrame on the 'IID' column

print('genetic data appended')

data = data.rename(columns={'mdd':'EUR'}) # rename cols

# ================================= MAKING DESIGN MATRIX ===================================
g_vars = ['EUR', 'AFR', 'AMR', 'EAS']
for var in g_vars:
    data[var] = preprocessing.scale(data[var])  # standardise PRS

df = data.drop(columns=['eventname']) # remove the numeric id col used in mplus
df = df.drop_duplicates(subset=["IID"])  # as data is in long format

#data['high_prs'].nunique() # check how many have risk scores, can perform checks with .best file length


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
    conf_int = np.exp(result.conf_int().iloc[[1,10,19],:])
    p_values = result.pvalues.iloc[1]
    p_values.columns = ['Decreasing', 'Increasing', 'Persistent']
    results_df = pd.DataFrame({'Variable': [var],
                               #'Class': result.model.endog_names,
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
results_table.to_csv('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/logreg_results/transancestry_results.csv', index=False)


odds = pd.concat(plot_result_list, ignore_index=True)
odds = odds.set_index('Variable')
groups = odds.groupby('Variable')

odds.to_csv('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/logreg_results/transancestry_plot_results.csv')


fig, ax = plt.subplots(figsize=(6, 8))
width = 0.1
y_pos = 0.3
variable_order = ['EUR', 'AFR', 'AMR', 'EAS']

color_map = {
    'EUR': 'black',
    'AFR': 'black',
    'AMR': 'black',
    'EAS': 'black'}
markers = ['d', 'X', 'h', 'v']

def sorting_key(name_group):
    name, group = name_group
    return variable_order.index(name)

sorted_groups = sorted(groups, key=sorting_key) # Sort the groups based on the custom sorting key function

for i, (name, group) in enumerate(sorted_groups):
    x_pos = np.arange(len(group.index)) + i * width

    ax.errorbar(group['Odds Ratio'], x_pos+width/2, xerr=[group['Odds Ratio']-group['Lower CI'],     # plot the odds ratios with error bars
                                                  group['Upper CI']-group['Odds Ratio']], fmt=markers[i],
                                                capsize=0.3,capthick=0.8, linewidth=0.8,color = color_map[name])

plt.subplots_adjust(left=0.15, bottom=0.08, right=0.8, top=0.95) # adjust the margins
categories = list(odds['Category'].unique()) # get the unique categories
ax.set_yticklabels(odds['Category'].unique(), fontsize=8)
ax.set_yticks(np.arange(len(categories))+y_pos)
ax.set_xlabel('Odds Ratio (95% CI)', fontsize=8)
#ax.set_xlim(min(odds['Lower CI'])-0.1, max(odds['Upper CI'])+0.1)
ax.set_xscale('log')  # Set the x-axis to a logarithmic scale

# Set the tick locator and formatter for the x-axis
ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0))
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{x:.1f}"))
#ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10, 4) * 0.1))
#ax.xaxis.set_minor_formatter(ticker.ScalarFormatter())
ax.tick_params(axis='x', which='both', labelsize=7)

legend_labels = [name for name, _ in sorted_groups]
#ax.legend(legend_labels, loc='lower right')
ax.legend(legend_labels, bbox_to_anchor=(1, 1), loc='upper right')
ax.axvline(x=1, linestyle='--', color='black', linewidth=0.8)
plt.rcParams['legend.fontsize'] = 6
plt.show()