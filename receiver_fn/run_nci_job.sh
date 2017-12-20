#!/bin/bash

#PBS -P vy72
#PBS -q hugemem
#PBS -l walltime=15:00:00,mem=257GB
#PBS -l ncpus=7
#PBS -l jobfs=2GB
#PBS -l wd
#PBS -l software=python

# sample invocation of the generate.py script
# python generate.py --asdf '/g/data/ha3/Passive/_ANU/7A(2011-2012)/ASDF/AA.h5' --inv '/g/data/ha3/Passive/_ANU/7A(2011-2012)/network_metadata/stnXML/AA.xml' --cat '/g/data/ha3/Passive/_ANU/7A(2011-2012)/event_metadata/earthquake/quakeML/fdsnws-event_2016-11-15T07_34_21.xml' --outdata 'data/7A-rf_profile_data.h5' --outrf 'data/7A-rf_profile_rfs.h5'

python generate.py > output.txt

