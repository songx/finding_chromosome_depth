import pandas as pd

import numpy as np


#read data into dataframe

df1 = pd.read_csv('sample1_depths.txt.gz',compression='gzip', header=0, sep='\t', quotechar='"')

df2 = pd.read_csv('sample2_depths.txt.gz',compression='gzip', header=0, sep='\t', quotechar='"')

df3 = pd.read_csv('sample3_depths.txt.gz',compression='gzip', header=0, sep='\t', quotechar='"')


#convert position to sets and find unqiue and common positions

df1_pos = set(df1['Chromosome']+":"+df1['Position'].astype(str))

df2_pos = set(df2['Chromosome']+":"+df2['Position'].astype(str))

df3_pos = set(df3['Chromosome']+":"+df3['Position'].astype(str))

print('1.',len(df2_pos - df1_pos - df3_pos),'positions are unique to sample2_depths.txt')

print('2.',len(df1_pos.intersection(df2_pos).intersection(df3_pos)),'positions are common between all three depth files')


#concatenate by chromosome

dfs = pd.concat([df1,df2,df3])

chr1_df = dfs[dfs['Chromosome']=='chr1'].sort_values('Position')

chr2_df = dfs[dfs['Chromosome']=='chr2'].sort_values('Position')

chr3_df = dfs[dfs['Chromosome']=='chr3'].sort_values('Position')


#find average and max depth at each position of each chromosome

print('3.')

#chr1 = chr1_df.groupby(['Chromosome','Position']).apply(lambda x: x['Depth'].sum()/3).to_frame('Depth').reset_index()

#print('Chr1 average read depth:',chr1['Depth'].sum()/1000000)

#chr1_max = chr1[chr1['Depth']==chr1['Depth'].max()]

chr1 = chr1_df.groupby(['Chromosome','Position']).agg({'Depth':np.sum}).reset_index()

chr1['Depth'] /= 3

print('Chr1 average read depth:',chr1['Depth'].sum()/1000000)

chr1_max = chr1[chr1['Depth']==chr1['Depth'].max()]


#chr2 = chr2_df.groupby(['Chromosome','Position']).apply(lambda x: x['Depth'].sum()/3).to_frame('Depth').reset_index()

#print('Chr2 average read depth:',chr2['Depth'].sum()/1000000)

#chr2_max = chr2[chr2['Depth']==chr2['Depth'].max()]

chr2 = chr2_df.groupby(['Chromosome','Position']).agg({'Depth':np.sum}).reset_index()

chr2['Depth'] /= 3

print('Chr2 average read depth:',chr2['Depth'].sum()/1000000)

chr2_max = chr2[chr2['Depth']==chr2['Depth'].max()]


#chr3 = chr3_df.groupby(['Chromosome','Position']).apply(lambda x: x['Depth'].sum()/3).to_frame('Depth').reset_index()

#print('Chr3 average read depth:',chr3['Depth'].sum()/1000000)

#chr3_max = chr3[chr3['Depth']==chr3['Depth'].max()]

chr3 = chr3_df.groupby(['Chromosome','Position']).agg({'Depth':np.sum}).reset_index()

chr3['Depth'] /= 3

print('Chr3 average read depth:',chr3['Depth'].sum()/1000000)

chr3_max = chr3[chr3['Depth']==chr3['Depth'].max()]


#find chromosome and position with largest average depth

max_dfs = pd.concat([chr1_max,chr2_max,chr3_max])

print('4.',max_dfs[max_dfs['Depth']==max_dfs['Depth'].max()]['Chromosome'].iloc[0],max_dfs[max_dfs['Depth']==max_dfs['Depth'].max()]['Position'].iloc[0],'has the largest average depth')
