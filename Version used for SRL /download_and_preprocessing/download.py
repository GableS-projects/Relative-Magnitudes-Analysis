#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 12:09:23 2023@author: sydneygable
"""
# import the necessary modules and functions

import pandas as pd
import obspy
from obspy.clients.fdsn.mass_downloader import RectangularDomain, Restrictions, MassDownloader
import os
import sys


import logging
logger = logging.getLogger("obspy.clients.fdsn.mass_downloader")
logger.setLevel(logging.ERROR)

########################################################
# function that reads in the following information and downloads the data for a particular event 
# st: start time (origin time) of event  
# net: network
# sta: station
# loc: location code
# dc: data center 
# tbo & tao: time (in seconds) to download before provided origin time, tao is same but after origin time 
#              total timespan of downloaded waveform is tbo + ta0

def download(st, net, sta, loc, cha, dc, tbo, tao):
    st = st-tbo
    et = st+tbo+tao

    restrictions = Restrictions(
        starttime=st,
        endtime=et,
        network=net, station=sta, location=loc, channel=cha,
        reject_channels_with_gaps=True,
        minimum_length=0.0,
        minimum_interstation_distance_in_m=100.0)

    mdl = MassDownloader(providers=[dc])
    domain = RectangularDomain(minlatitude=-90, maxlatitude=90,
                           minlongitude=-180, maxlongitude=180)
    mdl.download(domain, restrictions, mseed_storage=mseed_directory, stationxml_storage=station_response_directory)

    return

################################################################
# read in csv file which should contain all of the following input values
inputs = pd.read_csv(sys.argv[1], sep = ',', dtype = str)

# specify which row to read in system inputs 
idx = int(sys.argv[2])

net = str(inputs['Network'].iloc[idx])
sta = str(inputs['Station'].iloc[idx])
loc = str(inputs['Location'].iloc[idx])
cha = str(inputs['Channels'].iloc[idx])
dc = str(inputs['Data_center'].iloc[idx])

#print(net, type(net))
#print(sta, type(sta))
#print(loc, type(loc))
#print(cha, type(cha))
#print(dc, type(dc))


if loc != '00' and loc != '01':
    loc = ''


time_before_origin = float(inputs['Time_before_origin'].iloc[idx])
time_after_origin = float(inputs['Time_after_origin'].iloc[idx])

catalog = pd.read_csv(str(inputs['catalog_file'].iloc[idx]))

mseed_directory = str(inputs['mseed_filepath'].iloc[idx])
station_response_directory = str(inputs['xml_filepath'].iloc[idx])


###########################################################
# create mseed and station response file directories if they do not already exist

if not os.path.exists(mseed_directory):
    os.mkdir(mseed_directory)

if not os.path.exists(station_response_directory):
    os.mkdir(station_response_directory)


datetime = catalog['datetime_ID']

# read start times from catalog datetime_IDs, convert to datetime objects and append in start_time array
# then loop over this array to download each event 

start_time = []
for i in range(catalog.shape[0]):
    st = obspy.core.utcdatetime.UTCDateTime(datetime.iloc[i])
    start_time.append(st)



for st in start_time:
    download(st, net, sta, loc, cha, dc, time_before_origin, time_after_origin)


