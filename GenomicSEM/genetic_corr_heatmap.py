import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# all GWAS
df = pd.read_table("/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/GenomicSEM/gen_cor_matrix.txt")
df.index = ['ADHD', 'BIP', 'SCZ', 'ASD','MDD', 'NEU', 'ANX']
sns.heatmap(df,vmin=-1, vmax=1, annot=True, cmap=sns.color_palette("viridis", as_cmap=True) )
plt.show()


# excluding neuroticism
df2 = df.drop(columns=['NEU'], index=['NEU'])
mask = np.triu(np.ones_like(df2, dtype=bool))

sns.heatmap(df2,vmin=-1, vmax=1, annot=True, cmap=sns.color_palette("viridis", as_cmap=True), mask=mask)
plt.show()
