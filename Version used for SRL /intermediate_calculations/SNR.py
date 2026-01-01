#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 13:19:50 2023

@author: gablesydney
"""

import obspy 
import numpy as np
import pandas as pd
import glob
import sys


# In[2]:


windowLength = 8

#inputs = pd.read_csv('../example_data/SNR_inp.csv')
inputs = pd.read_csv(sys.argv[1], sep=',', dtype = str)

#idx = int(1)
idx = int(sys.argv[2])

catalog = pd.read_csv(inputs['catalog'].iloc[idx])

sta = inputs['station'].iloc[idx]
et = int(inputs['extra_time'].iloc[idx])
data_dir = inputs['data_dir'].iloc[idx]
SNR_thresh = float(inputs['SNR_thresh'].iloc[idx])
dist_thresh = float(inputs['dist_thresh'].iloc[idx])
output_dir = inputs['output_dir'].iloc[idx]
pairs_num = 1000000

catalog = catalog[np.isnan(catalog[sta+'_ptime']) == False]
catalog = catalog.reset_index(drop = True)

file_list = glob.glob(data_dir+'/*.SAC')
datafile_datetime = []
for i in range(len(file_list)):
    f = file_list[i]
    datafile_datetime.append(f[-20:-4])
    
unique_list = list(set(datafile_datetime))
print(unique_list)
########################################################################


# In[3]:


index_list = []
for i in range(catalog.shape[0]):
    ds = catalog['datetime_ID'].iloc[i]
    boolean = ds in unique_list
    
    if boolean == False:
        index_list.append(i)
    
catalog.drop(index = index_list, inplace = True)

SNR1 = np.zeros(catalog.shape[0])
SNR2 = np.zeros(catalog.shape[0])
SNRZ = np.zeros(catalog.shape[0])
pWinStart = np.zeros(catalog.shape[0])
pWinEnd = np.zeros(catalog.shape[0])
sWinStart = np.zeros(catalog.shape[0])
sWinEnd = np.zeros(catalog.shape[0])


for i in range(catalog.shape[0]):
    print('SNR', i)
    
    # read in data file 
    st = obspy.read(data_dir + '/*' + catalog['datetime_ID'].iloc[i] + '.SAC')
    
    # check that the event has all three channels, we could process without all three
    # but it is likely that if a channel is missing there is a problem 
    if len(st) != 3:
        SNR1[i] = np.nan
        SNR2[i] = np.nan
        SNRZ[i] = np.nan
        pWinStart[i] = np.nan
        pWinEnd[i] = np.nan
        sWinStart[i] = np.nan
        sWinEnd[i] = np.nan
        
    else:
        # get p and s arrival times from the catalog 
        p = catalog[sta+'_ptime'].iloc[i] + et
        s = catalog[sta+'_stime'].iloc[i] + et
        
        for trace in st:
            
            # set variables for the data and times of each channel based on the the 
            # component name
            if trace.stats.channel[-1] == 'Z':
                dataZ = trace.data
                timesZ = trace.times()
                
            elif trace.stats.channel[-1] == '1' or trace.stats.channel[-1] == 'E':
                data1 = trace.data
                times1 = trace.times()
                
            elif trace.stats.channel[-1] == '2' or trace.stats.channel[-1] == 'N':
                data2 = trace.data
                times2 = trace.times()
                

        
        numSamples = windowLength * trace.stats.sampling_rate
        
        pWindowStart = np.absolute(times1 - (p - 1)).argmin()
        pWindowEnd = np.absolute(times1 - (p + 7)).argmin()
        
        sWindowStart = np.absolute(times1 - (s - 1)).argmin()
        sWindowEnd = np.absolute(times1 - (s + 7)).argmin()

        nWindowStart = np.absolute(times1 - (p - 13)).argmin()
        nWindowEnd = np.absolute(times1 - (p - 5)).argmin()
        
        if pWindowEnd - pWindowStart < numSamples or sWindowEnd - sWindowStart < numSamples or nWindowEnd - nWindowStart < numSamples:
            
            SNR1[i] = np.nan
            SNR2[i] = np.nan
            SNRZ[i] = np.nan
            pWinStart[i] = np.nan
            pWinEnd[i] = np.nan
            sWinStart[i] = np.nan
            sWinEnd[i] = np.nan
        
        else:
            rms1 = np.sqrt(np.sum(data1[nWindowStart:nWindowEnd]**2)/len(data1[nWindowStart:nWindowEnd]))
            sig1 = max(data1[sWindowStart:sWindowEnd])
        
            rms2 = np.sqrt(np.sum(data2[nWindowStart:nWindowEnd]**2)/len(data2[nWindowStart:nWindowEnd]))
            sig2 = max(data2[sWindowStart:sWindowEnd])
    
            rmsZ = np.sqrt(np.sum(dataZ[nWindowStart:nWindowEnd]**2)/len(dataZ[nWindowStart:nWindowEnd]))
            sigZ = max(dataZ[pWindowStart:pWindowEnd])
            
            snr1 = sig1/rms1
            snr2 = sig2/rms2
            snrZ = sigZ/rmsZ
        
            SNR1[i] = snr1
            SNR2[i] = snr2
            SNRZ[i] = snrZ
            pWinStart[i] = pWindowStart
            pWinEnd[i] = pWindowEnd
            sWinStart[i] = sWindowStart
            sWinEnd[i] = sWindowEnd
        
        
 

data = {'datetime_ID' : catalog['datetime_ID'], 'SNR_1' : SNR1, 'SNR_2' : SNR2, 'SNR_Z' : SNRZ, 'p_start' : pWinStart, 'p_end' : pWinEnd, 's_start' : sWinStart, 's_end' : sWinEnd}
SNR = pd.DataFrame(data)

SNR = SNR.dropna()
SNR.to_csv(output_dir + '/snrEventResults_'+sta+'.csv', index = False)


# In[4]:



cha_ext = ['1','2','Z']

pairs_num = 1000000
pairs = []
counter = 0
output_counter = 0
for m in range(SNR.shape[0]):
    print('Creating Pairs', m)
    for n in range(counter+1):
        if m != n:
            
            dtm = SNR['datetime_ID'].iloc[m]
            dtn = SNR['datetime_ID'].iloc[n]

            m_cat = catalog[catalog['datetime_ID'] == dtm]
            n_cat = catalog[catalog['datetime_ID'] == dtn]

            dist = obspy.geodetics.base.locations2degrees(m_cat['Lat'].iloc[0], m_cat['Lon'].iloc[0], n_cat['Lat'].iloc[0], n_cat['Lon'].iloc[0])
            
            if dist <= dist_thresh:
            
                for i in [0,1,2]:
                    if i ==0:
                        cha = cha_ext[0]
                        snr = SNR['SNR_1']
                        if snr.iloc[m]>SNR_thresh and snr.iloc[n]>SNR_thresh:
                            pair = np.array([dtm,dtn,cha])
                            pairs.append(pair)
                            continue

                    elif i == 1:
                        cha = cha_ext[1]
                        snr = SNR['SNR_2']
                        if snr.iloc[m] > SNR_thresh and snr.iloc[n] > SNR_thresh:
                            pair = np.array([dtm,dtn,cha])
                            pairs.append(pair)
                            continue

                    elif i == 2:
                        cha = cha_ext[2]
                        snr = SNR['SNR_Z']
                        if snr.iloc[m] > SNR_thresh and snr.iloc[n] > SNR_thresh:
                            pair = np.array([dtm,dtn,cha])
                            pairs.append(pair)
                            continue 
                            
                if len(pairs) >= pairs_num:
                    pairs_df = pd.DataFrame(pairs, columns = ['dtm', 'dtn', 'cha'])
                    pairs_df.to_csv(output_dir+'/pairs_output_'+sta+'_'+str(output_counter))
                    
                    pairs = []
                    output_counter += 1 



    counter+=1
    
    
pairs_df = pd.DataFrame(pairs, columns = ['dtm', 'dtn', 'cha'])
pairs_df.to_csv(output_dir+'/pairs_output_'+sta+'_'+str(output_counter))

#pairs = np.asarray(pairs)
#np.save(output_dir + '/snrOutfile_'+sta+'.npy', pairs)


# In[ ]:




