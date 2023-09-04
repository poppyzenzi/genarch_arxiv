import numpy as np
import pandas as pd
import statsmodels.api as sm
import os.path
from sklearn import preprocessing
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

os.chdir('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/logreg_results')

abcd = pd.read_csv('abcd_plot_results.csv')
alspac = pd.read_csv('alspac_plot_results.csv')

cohorts = [abcd, alspac]

fig, axs = plt.subplots(1, 2, figsize=(12, 8))
width = 0.06
y_pos = 0.5
variable_order = ['common', 'hierarchical', 'mdd_prs',
                  'neu_prs', 'anx_prs', 'scz_prs',
                  'bip_prs', 'adhd_prs', 'asd_prs']

markers = ['o', 's', 'X', 'D', '>', 'P', '^', 'p', '*']

variables = {
    'common': {
        'legend_text': 'Common',
        'marker': 'o',
        'color': 'moccasin'
    },
    'hierarchical': {
        'legend_text': 'Hierarchical',
        'marker': 's',
        'color': 'orange'
    },
    'mdd_prs': {
        'legend_text': 'MDD',
        'marker': 'X',
        'color': 'darkblue'
    },
    'neu_prs': {
        'legend_text': 'NEU',
        'marker': 'D',
        'color': 'blue'
    },
    'anx_prs': {
        'legend_text': 'ANX',
        'marker': '>',
        'color': 'dodgerblue'
    },
    'scz_prs': {
        'legend_text': 'SCZ',
        'marker': 'P',
        'color': 'darkslategrey'
    },
    'bip_prs': {
        'legend_text': 'BIP',
        'marker': '^',
        'color': 'slategrey'
    },
    'adhd_prs': {
        'legend_text': 'ADHD',
        'marker': 'p',
        'color': 'deeppink'
    },
    'asd_prs': {
        'legend_text': 'ASD',
        'marker': '*',
        'color': 'lightpink'
    }
}

def sorting_key(name_group):
    name, group = name_group
    return variable_order.index(name)

cohort_names = ['ABCD', 'ALSPAC']

for cohort_idx, odds in enumerate(cohorts):
    groups = odds.groupby('Variable')
    sorted_groups = sorted(groups, key=sorting_key)
    ax = axs[cohort_idx]

    for i, (name, group) in enumerate(sorted_groups):
        if i < 2:
            x_pos = np.arange(len(group.index)) + i * width
        else:
            x_pos = np.arange(len(group.index)) + (i + 2) * width

        ax.errorbar(
            group['Odds Ratio'], x_pos + width / 2,
            xerr=[group['Odds Ratio'] - group['Lower CI'], group['Upper CI'] - group['Odds Ratio']],
            fmt=markers[i],
            capsize=0.3, capthick=0.8, linewidth=0.8, color=variables[name]['color']
        )

    categories = odds['Category'].unique()  # Extract categories from the current cohort
    ax.set_yticklabels(categories, fontsize=8)
    ax.set_yticks(np.arange(len(categories)) + y_pos)
    ax.set_xlabel('Odds Ratio (95% CI)', fontsize=8)
    ax.set_xscale('log')
    ax.set_xlim([0.85, 1.7])

    ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{x:.1f}"))
    ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax.xaxis.set_minor_formatter(ticker.ScalarFormatter())
    ax.tick_params(axis='x', which='both', labelsize=7)

    ax.set_title(cohort_names[cohort_idx], fontsize=10)

    if cohort_idx == 1:
        legend_handles = []
        for i, name in enumerate(variable_order):
            legend_handles.append(
                Line2D([0], [0], marker=markers[i], color='w',
                       markerfacecolor=variables[name]['color'], markersize=8,
                       label=variables[name]['legend_text'])
            )
        ax.legend(handles=legend_handles, loc='upper right', fontsize=7)

    ax.axvline(x=1, linestyle='--', color='black', linewidth=0.8)

plt.subplots_adjust(left=0.08, bottom=0.08, right=0.95, top=0.95, wspace=0.25)
plt.show()
