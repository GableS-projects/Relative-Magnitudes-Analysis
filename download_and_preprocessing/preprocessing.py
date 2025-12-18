#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 11:40:51 2023

@author: sydneygable
"""

# read in necessary modules 
import obspy
import glob
import os 
import sys 
import pandas as pd

# read in input file and subsequent input information

inputs = pd.read_csv(sys.argv[1], sep = ',', dtype = str)
idx = int(sys.argv[2])
print(idx)

seed_dir = str(inputs['seed_dir'].iloc[idx])
sac_dir = str(inputs['sac_dir'].iloc[idx])+'/'
freqmin = float(inputs['freqmin'].iloc[idx])
freqmax = float(inputs['freqmax'].iloc[idx])

extra_time = True
et = float(inputs['extra_time'].iloc[idx])



# check to see if the output directory already exists, if not create it

if not os.path.exists(sac_dir):
    os.mkdir(sac_dir)
    
# read in all of the miniseed files from the seed_dir, then detrend taper and filter the files, and save as sac file with standard filename

st_list = glob.glob(seed_dir+'/*.mseed')
print(seed_dir+'*.mseed')
print(len(st_list))

for i in range(len(st_list)):
    print(i)
    st = obspy.read(st_list[i])
    
    st.detrend('demean')
    st.taper(max_percentage = 0.05)
    st.filter('bandpass', freqmin = freqmin, freqmax = freqmax, corners = 4)

    stream = st[0]
    
    hdr = stream.stats
    
    start = hdr.starttime+et
    
    start = start.strftime('%Y%m%dT%H%M%SZ')
    new_filename = stream.id+'__'+start
    
    stream.write(sac_dir+new_filename+'.SAC', format = 'SAC')
