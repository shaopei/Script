import pandas as pd 
import numpy as np
hic_data = pd.read_csv("/local/workdir/Gilad/merged_nodups_HS_NHS_HIC069_HIC070_HIC071_HIC072_HIC073_HIC074_sorted_filtered.txt", sep = ' ', header = None, usecols = [1,2,3,4,5,6,7,8])
hic_data.columns = ['str1', 'chr1', 'pos1', 'frag1', 'str2', 'chr2', 'pos2', 'frag2']

print 'finished loading'

chr1 = []
for i in list(hic_data['chr1']):
  chr1.append(str(i))
hic_data['chr1'] = chr1

chr2 = []
for i in list(hic_data['chr2']):
  chr2.append(str(i))
hic_data['chr2'] = chr2

chromosomes = ['X', 'Y']
for i in range(1,23):
	chromosomes.append(str(i))

for chr in chromosomes:
	out = hic_data.loc[np.logical_and(hic_data['chr1'] == chr, hic_data['chr2'] == chr)]
	out.to_csv("/local/workdir/Gilad/K562_data/K562_Hi-C_interactions_chr_" + chr + ".txt", sep = ' ')
