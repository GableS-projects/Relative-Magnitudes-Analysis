#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 22:25:00 2023

@author: gablesydney
"""
import pandas as pd
import obspy
from obspy.signal.cross_correlation import correlate, xcorr_max
import numpy as np
from numpy import linalg as LA
import sys

#####################  INPUTS  #################################
inputs = pd.read_csv(sys.argv[1])

idx = int(sys.argv[2])

sta = inputs['station'].iloc[idx]
file_tag = inputs['file_tag'].iloc[idx]

data_dir = inputs['data_dir'].iloc[idx]
output_dir = inputs['output_dir'].iloc[idx]
pairs_file = inputs['pairs_file'].iloc[idx]
pairs = pd.read_csv(pairs_file)

catalog_file = inputs['catalog'].iloc[idx]
catalog = pd.read_csv(catalog_file)
                      
pairs = pairs[pairs['dtm'] != pairs['dtn']]                    
################################################################

final_results = [None]*pairs.shape[0]

for i in range(pairs.shape[0]):
    #if i%10000 == 0:
    #    print(i)
    
    pair = pairs.iloc[i]
    dtm = pair['dtm']
    dtn = pair['dtn']
    cha = pair['cha'][-1]
    
    if dtm == dtn:
        continue
    
    if cha == '1' or cha == '2':
        s1 = catalog[catalog['datetime_ID'] == dtm].iloc[0][sta+'_stime']
        s2 = catalog[catalog['datetime_ID'] == dtn].iloc[0][sta+'_stime']

    else:
        s1 = catalog[catalog['datetime_ID'] == dtm].iloc[0][sta+'_ptime']
        s2 = catalog[catalog['datetime_ID'] == dtn].iloc[0][sta+'_ptime']
    
    if cha == '1':
        cha = 'E'
    if cha == '2':
        cha = 'N'

    waveform1 = obspy.read(data_dir+file_tag+cha+'__'+dtm+'.SAC')
    waveform2 = obspy.read(data_dir+file_tag+cha+'__'+dtn+'.SAC')
  
    
    waveform1.trim(waveform1[0].stats.starttime +20+s1 -1, waveform1[0].stats.starttime +20+s1 +7)
    waveform2.trim(waveform2[0].stats.starttime +20+s2 -1, waveform2[0].stats.starttime +20+s2 +7)
    
    trace1 = waveform1[0]
    trace2 = waveform2[0] 

    if len(trace1.data) > len(trace2.data):
        trace1.data = trace1.data[:-1]
    elif len(trace1.data) < len(trace2.data):
        trace2.data = trace2.data[:-1]
        
    cc_func=correlate(trace1, trace2, 100)
    
    shift, cc = xcorr_max(cc_func)
    
    
    if shift == 0:
        shifted_trace1 = trace1.data
        shifted_trace2 = trace2.data

    elif shift != 0:

        if shift > 0:
            shifted_trace1 = trace1.data[shift:]
            shifted_trace2 = trace2.data[:shift*-1]

        elif shift < 0:
            shifted_trace1 = trace1.data[:shift]
            shifted_trace2 = trace2.data[shift*-1:]

    COV = np.cov(shifted_trace1, shifted_trace2)
    vals, vects = LA.eig(COV)
    maxcol = list(vals).index(max(vals))
    eigenvect = vects[:,maxcol]
    AR = eigenvect[1]/eigenvect[0]
    
    
    result = [dtm, dtn, cha, cc, AR]
    final_results[i] = result
    
result_df = pd.DataFrame(final_results, columns=['dtm', 'dtn', 'cha', 'cc', 'AR'])
result_df.to_csv(output_dir+'/results_'+sta+'_'+str(idx))  
    


