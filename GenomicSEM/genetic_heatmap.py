import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# add significance level?

# all GWAS
df = pd.read_table("/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/GenomicSEM/output_7_factors/gen_cor_matrix.txt")
df.index = ['ADHD', 'BIP', 'SCZ', 'ASD','MDD', 'NEU','ANX']

# create a mask for the upper diagonal half
mask = np.triu(np.ones_like(df, dtype=bool), k=1)

# set the diagonal values to 1
np.fill_diagonal(df.values, 1)

# plot the heatmap with the mask
sns.heatmap(df, vmin=-1, vmax=1, annot=True, cmap=sns.color_palette("viridis", as_cmap=True), mask=mask)

plt.show()



