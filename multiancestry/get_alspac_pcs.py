import pandas as pd

#Read in .dta file (connected to VPN)
alspac = pd.read_stata("/Volumes/cmvm/scs/groups/ALSPAC/data/B3421_Whalley_04Nov2021.dta")

bkup = alspac
bkup = bkup[['cidB3421', 'qlet', 'pc1', 'pc2', 'pc3', 'pc4', 'pc5', 'pc6', 'pc7', 'pc8', 'pc9', 'pc10']]

bkup['cidB3421'] = bkup['cidB3421'].astype(int)
bkup['IID'] = bkup['cidB3421'].astype(str) + bkup['qlet']

bkup.drop(columns=['cidB3421','qlet'], inplace=True)
bkup

file_path = '/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/multiancestry/alspac_pca/alspac_pcs'

bkup.to_csv(file_path, index=False)