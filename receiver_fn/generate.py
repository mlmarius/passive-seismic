import os.path

import matplotlib.pyplot as plt
import numpy as np
from obspy import read_inventory, read_events, UTCDateTime as UTC
from obspy.clients.fdsn import Client
from rf import read_rf, RFStream
from rf import get_profile_boxes, iter_event_data, IterMultipleComponents
from rf.imaging import plot_profile_map
from rf.profile import profile
from tqdm import tqdm

import pyasdf
from pyasdf import ASDFDataSet
from obspy.core import Stream

import argparse

data = os.path.join('data', '')
invfile = data + '7D-inventory.xml'
catfile = data + '7D-catalog.xml'
datafile = data + '7D-rf_profile_data.h5'
rffile = data + '7D-rf_profile_rfs.h5'
profilefile = data + '7D-rf_profile_profile.h5'

parser = argparse.ArgumentParser(
    description="Calculation of receiver functions for all stations in a temporary survey.. ",
)

parser.add_argument(
    'asdffile',
    default='/g/data/ha3/Passive/_ANU/7D(2012-2013)/ASDF/7D(2012-2013).h5',
    type=str,
    help='filepath of the asdf input file .')

parser.add_argument(
    '--inv',
    dest='inventory',
    default=invfile,
    type=str,
    help='filepath of the station xml file.')

parser.add_argument(
    '--cat',
    dest='catalog',
    default=catfile,
    type=str,
    help='filepath of the event catalog xml file.')

parser.add_argument(
    '--bandpasslow',
    dest='bpl',
    default=0.5,
    type=float,
    help='lower threshold for the band pass filter.')

parser.add_argument(
    '--bandpasshigh',
    dest='bph',
    default=2,
    type=float,
    help='higher threshold for the band pass filter..')

parser.add_argument(
    '--outdata',
    dest='data',
    default=datafile,
    type=str,
    help='filepath of the desired output waveform file in asdf format.')

parser.add_argument(
    '--outrf',
    dest='rff',
    default=rffile,
    type=str,
    help='filepath of the desired output receiver function file in asdf format.')

args = parser.parse_args()

def custom_get_waveforms(network, station, location, channel, starttime,
                  endtime, quality=None, minimumlength=None,
                  longestonly=None, filename=None, attach_response=False,
                  **kwargs):
    with pyasdf.ASDFDataSet(args.asdffile, mode='r') as asdfDataSet:
        st = Stream()
        #ignoring channel for now as all the 7D network waveforms have only BH? channels
        filteredList = [i for i in asdfDataSet.waveforms[network+'.'+station].list() if
                        'raw_recording' in i and
                        UTC(i.split("__")[1]) < starttime and
                        UTC(i.split("__")[2]) > endtime]
        for t in filteredList:
            st += asdfDataSet.waveforms[network+'.'+station][t]
        return st

inventory = read_inventory(args.inventory)
catalog = read_events(args.catalog)

stream = RFStream()
with tqdm() as pbar:
    for s in iter_event_data(catalog, inventory, custom_get_waveforms, pbar=pbar):
        stream.extend(s)

stream.write(args.data, 'H5')

rfstream = RFStream()
for stream3c in tqdm(IterMultipleComponents(stream, 'onset', 3)):
    stream3c.filter('bandpass', freqmin=args.bpl, freqmax=args.bph)
    stream3c.trim2(-25, 75, 'onset')
    if len(stream3c) != 3:
        continue
    stream3c.rf()
    stream3c.moveout()
    rfstream.extend(stream3c)
rfstream.write(args.rff, 'H5')

