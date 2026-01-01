import pandas as pd
import numpy as np
from scipy import sparse 
from scipy.sparse import spdiags, vstack
from scipy.sparse.linalg import lsqr
import glob
import datetime

################################ USER INPUTS ##############################

catalog = pd.read_csv()
prior_events = pd.read_csv()

input_files = glob.glob()

CC_cut = 
scaling = 

min_num_sta = 

############################################################################


timestamp = []
for i in range(catalog.shape[0]):
    timestamp.append(datetime.datetime.strptime(catalog['datetime_ID'].iloc[i], '%Y%m%dT%H%M%SZ'))

catalog['timestamp'] = timestamp

file_list = []

for item in input_files: 

    frame = pd.read_csv(item)


    sta_label = item[12:-4]
    l = [sta_label]*frame.shape[0]
    frame['sta'] = l

    frame = frame[frame['cc'] > CC_cut]

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


final_results = [np.nan] * catalog.shape[0]

for q in range(catalog.shape[0]):
    print(q)
    
    catalog_date = catalog['datetime_ID'].iloc[q]

    b1 = df_new[(df_new['dtm'].isin(prior_events['event ID'])) & (df_new['dtn'] == catalog_date)]
    b2 = df_new[(df_new['dtn'].isin(prior_events['event ID'])) & (df_new['dtm'] == catalog_date)]

    bb = pd.concat([b1,b2])
    print('remapping events to new indices')


    #unique_events = np.append(unique_events, catalog_date)
    all_timestamps = pd.concat([bb['dtm'], bb['dtn']]).to_numpy()
    unique_events = np.unique(all_timestamps)

    
    number_of_events = len(unique_events)
    map_label = np.arange(0, number_of_events, 1)

    re_mapped_i = np.zeros(bb.shape[0])
    re_mapped_j = np.zeros(bb.shape[0])

    for idx in range(bb.shape[0]):
        event = bb['dtm'].iloc[idx]
        index = np.where(unique_events == event)
        re_mapped_i[idx] = index[0][0]

    print("\t done with event i's") 

    for idx in range(bb.shape[0]):
        event = bb['dtn'].iloc[idx]
        index = np.where(unique_events == event)
        re_mapped_j[idx] = index[0][0]

    print("\t done with event j's")

    bb['re_mapped_i'] = re_mapped_i
    bb['re_mapped_j'] = re_mapped_j
    
    
    d = bb

    if bb.shape[0] < 1:
        continue
    
    
    AR_array = np.log10(abs(d['AR']))*scaling

    catalog_mag = []
    for item in unique_events:
        event = catalog[catalog['datetime_ID'] == item]

        if event.shape[0] < 1:
            catalog_mag.append(np.nan)

        else:

            event_mag = float(event['Local Magnitude'].iloc[0])
            catalog_mag.append(event_mag)

    data = {'event ID':unique_events, 're-mapped index':map_label, 'catalog mag': catalog_mag}
    results = pd.DataFrame(data = data)

    scaling_magnitudes = []
    for item in unique_events:
        event = prior_events[prior_events['event ID'] == item]

        if event.shape[0] < 1:
            scaling_magnitudes.append(np.nan)

        else:

            event_mag = float(event['new mag'].iloc[0])
            scaling_magnitudes.append(event_mag)

    results['scaling mag'] = scaling_magnitudes

    matrix_shape = [d.shape[0], number_of_events]
    matrix = np.zeros(matrix_shape)

    #populate the matrix with -1 for the i events and 1 for the j events 
    for m in range(d.shape[0]):
        n = int(d['re_mapped_i'].iloc[m])
        matrix[m][n] = -1 

        p = int(d['re_mapped_j'].iloc[m])
        matrix[m][p] = 1

    scaling_results = results[results['scaling mag'].isna() == False]


    constraint_rows = np.zeros((scaling_results.shape[0], results.shape[0]))

    for i in range(scaling_results.shape[0]):
        #print(i)
        j = scaling_results['re-mapped index'].iloc[i]

        constraint_rows[i,j] = 1

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

    cc = d['cc']
    cc = cc.tolist()

    for i in range(scaling_results.shape[0]):
        cc.append(10)

    N = A_sparse.shape[0]
    data = np.array(cc)
    diags = np.array([0])
    W = spdiags(data, diags, N, N)

    print('adding constant column')
    #constant_col = sparse.csr_matrix(np.ones((N,1)))

    #A_new = hstack((constant_col,A_sparse))
    A_new = A_sparse
    print('weighting A and b matrices')
    Aw_new = W @ A_new
    bw = W@b_sparse

    print('converting weighted sparse matrices to regular arrays')
    # Aw_new = Aw_new.toarray()
    bw = bw.toarray()

    print('performing least squares inversion')
    soln = lsqr(Aw_new, bw)

    results['new mag'] = soln[0]
    
    if results[results['event ID'] == catalog_date].shape[0] < 1:
        continue

    r = results[results['event ID'] == catalog_date]['new mag'].iloc[0]
    final_results[q] = r

    
catalog['final_results'] = final_results

catalog.to_csv('new_results_v2.csv')
