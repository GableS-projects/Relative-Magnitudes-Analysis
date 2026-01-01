#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 18:49:24 2025

@author: sydneygable
"""

import pandas as pd
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees

######################## USER INPUTS #####################################

catalog = pd.read_csv()
stations = pd.read_csv()

model = TauPyModel(model="ak135")
output_catalog = pd.read_csv()

##########################################################################
for j in range(stations.shape[0]):
    print(j)
    
    p_list = []
    s_list = []
    
    sta = stations['Station'].iloc[j]
    
    sta_lat = stations['Latitude'].iloc[j]
    sta_lon = stations['Longitude'].iloc[j]
    
    for i in range(catalog.shape[0]):
    
       #print(i, 'of', catalog.shape[0])
        
        evt_lat = catalog['lat'].iloc[i]
        evt_lon = catalog['lon'].iloc[i]
        
        deg_dist = locations2degrees(sta_lat, sta_lon, evt_lat, evt_lon) 
        
        p_arrivals = model.get_travel_times(source_depth_in_km=catalog['d(km)'].iloc[i], distance_in_degree=deg_dist, phase_list=['p','P'])
        s_arrivals = model.get_travel_times(source_depth_in_km=catalog['d(km)'].iloc[i], distance_in_degree=deg_dist, phase_list=['s','S'])
        
        p_list.append(p_arrivals[0].time)
        s_list.append(s_arrivals[0].time)
        
    catalog[sta+'_ptime'] = p_list
    catalog[sta+'_stime'] = s_list
    
    
catalog.to_csv(output_catalog)





