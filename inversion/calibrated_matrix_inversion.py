#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse 
from scipy.sparse import spdiags, hstack, vstack
from scipy.sparse.linalg import lsqr
import glob
from scipy.sparse.linalg import inv

import obspy
import datetime
import networkx as nx
import itertools


###################### USER INPUTS ##########################


catalog = pd.read_csv()
    

input_files = glob.glob()

CC_cut = 
scaling = 

min_num_sta = 

scaling_magnitude_tag = 

output_file = 

###################################################################3


file_list = []

for item in input_files: 
    
    frame = pd.read_csv(item)
    
    frame = frame[frame['dtm'].isin(catalog['datetime_ID']) == True]
    frame = frame[frame['dtn'].isin(catalog['datetime_ID']) == True]
    
    sta_label = item.split('_')[2]
    l = [sta_label]*frame.shape[0]
    frame['sta'] = l
    
    frame = frame[frame['cc'] >= CC_cut]
    
    file_list.append(frame)
    print(frame.shape[0])
    

df = pd.concat(file_list, axis = 0, ignore_index = True)

aggregation_functions = {'cc':'mean', 'AR':'mean', 'sta':list}
df_new = df.groupby(['dtm', 'dtn'], as_index=False).aggregate(aggregation_functions).reindex(columns=df.columns)



ix = []
for i in range(df_new.shape[0]):
    if len(df_new['sta'].iloc[i]) < min_num_sta:
        ix.append(i)
        
df_new = df_new.drop(index = ix)



L = df_new[['dtm', 'dtn']].values.tolist()
L = list(tuple(sub) for sub in L)

G=nx.from_edgelist(L)

l=list(nx.connected_components(G))
# after that we create the map dict , for get the unique id for each nodes
mapdict={z:x for x, y in enumerate(l) for z in y }
# then append the id back to original data for groupby 
newlist=[ x+(mapdict[x[0]],)for  x in L]

#using groupby make the same id into one sublist
newlist=sorted(newlist,key=lambda x : x[2])
yourlist=[list(y) for x , y in itertools.groupby(newlist,key=lambda x : x[2])]
yourlist[2]

l = 0
yourlist_index = 0

for i in range(len(yourlist)):

    if len(yourlist[i]) > l:
        l = len(yourlist[i])
        yourlist_index = i

x = []

for i in range(len(yourlist[yourlist_index])):
    x.append(yourlist[yourlist_index][i][0])
    x.append(yourlist[yourlist_index][i][1])

x = list(set(x))

df_new = df_new[df_new['dtm'].isin(x)]
df_new = df_new[df_new['dtn'].isin(x)]



AR_array = np.log10(abs(df_new['AR']))*scaling

print('remapping events to new indices')

# combine both lists of timestamps from dataframe, convert to numpy array and remove 
# the repeats
all_timestamps = pd.concat([df_new['dtm'], df_new['dtn']]).to_numpy()
unique_events = np.unique(all_timestamps)

number_of_events = len(unique_events)
map_label = np.arange(0, number_of_events, 1)

re_mapped_i = np.zeros(df_new.shape[0])
re_mapped_j = np.zeros(df_new.shape[0])

for idx in range(df_new.shape[0]):
    event = df_new['dtm'].iloc[idx]
    index = np.where(unique_events == event)
    re_mapped_i[idx] = index[0][0]

print("\t done with event i's") 

for idx in range(df_new.shape[0]):
    event = df_new['dtn'].iloc[idx]
    index = np.where(unique_events == event)
    re_mapped_j[idx] = index[0][0]

print("\t done with event j's")

df_new['re_mapped_i'] = re_mapped_i
df_new['re_mapped_j'] = re_mapped_j


catalog_mag = []
for item in unique_events:
    event = catalog[catalog['datetime_ID'] == item]

    if event.shape[0] < 1:
        catalog_mag.append(np.nan)

    else:

        event_mag = float(event[comparison_magnitude_tag].iloc[0])
        catalog_mag.append(event_mag)
        
data = {'event ID':unique_events, 're-mapped index':map_label, 'catalog mag': catalog_mag}
results = pd.DataFrame(data = data)


scaling_magnitudes = []
for item in unique_events:
    event = catalog[catalog['datetime_ID'] == item]

    if event.shape[0] < 1:
        scaling_magnitudes.append(np.nan)

    else:

        event_mag = float(event[scaling_magnitude_tag].iloc[0])
        scaling_magnitudes.append(event_mag)

results['scaling mag'] = scaling_magnitudes


matrix_shape = [df_new.shape[0], number_of_events]
matrix = np.zeros(matrix_shape)

# populate the matrix with -1 for the i events and 1 for the j events 
for m in range(df_new.shape[0]):
    n = int(re_mapped_i[m])
    matrix[m][n] = -1 

    p = int(re_mapped_j[m])
    matrix[m][p] = 1

scaling_results = results.dropna()

constraint_rows = np.empty(results.shape[0])

for i in range(scaling_results.shape[0]):

    row = np.zeros(results.shape[0])    
    row[scaling_results['re-mapped index'].iloc[i]] = 1 

    if i == 0:
        constraint_rows = row

    else:
        constraint_rows = np.vstack([constraint_rows,row])  

print('converting matrices to sparse format')
matrix_sparse = sparse.csr_matrix(matrix)
constraint_rows_sparse = sparse.csr_matrix(constraint_rows)

A_sparse = vstack([matrix_sparse, constraint_rows_sparse])

del matrix
del constraint_rows

b = np.append(AR_array, np.asarray(scaling_results['scaling mag']))
b_sparse = sparse.csr_matrix(b)
b_sparse = sparse.csr_matrix.transpose(b_sparse)
del b
print('done converting')

cc = df_new['cc']
cc = cc.tolist()

for i in range(scaling_results.shape[0]):
    cc.append(10)

N = A_sparse.shape[0]
data = np.array(cc)
diags = np.array([0])
W = spdiags(data, diags, N, N)

print('adding constant column')
constant_col = sparse.csr_matrix(np.ones((N,1)))

A_new = A_sparse

print('weighting A and b matrices')
Aw_new = W @ A_new
bw = W@b_sparse


bw = bw.toarray()

print('performing least squares inversion')
soln = lsqr(Aw_new, bw)

results['new mag'] = soln[0]

results.to_csv(output_file)



