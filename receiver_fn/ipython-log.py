import os.path
from obspy import read_inventory, read_events, UTCDateTime as UTC
from obspy.clients.fdsn import Client
from rf import read_rf, RFStream
from rf import get_profile_boxes, iter_event_data, IterMultipleComponents
from rf.imaging import plot_profile_map
from rf.profile import profile
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
invfile = 'rfprofile_stations.xml'
catfile = 'rfprofile_events.xml'
datafile = 'rfprofile_data.h5'
rffile = 'rfprofile_rfs.h5'
profilefile = 'rfprofile.h5'
inventory = read_inventory(invfile)
stream = read_rf(rffile, 'H5')
ppoints = stream.ppoints(70)
boxes = get_profile_boxes((-21.3, -70.7), 90, np.linspace(0, 180, 73), width=530)
plt.figure(figsize=(10, 10))
plot_profile_map(boxes, inventory=inventory, ppoints=ppoints)
pstream = profile(tqdm(stream), boxes) # this is having some issues
help(RFStream.profile)
pstream = RFStream.profile(stream, boxes) # instead this works
pstream.write(profilefile, 'H5')
print(pstream)
pstream = read_rf(profilefile)
pstream.trim2(-5, 20, 'onset')
pstream.select(channel='??Q').normalize().plot_profile(scale=1.5, top='hist') # this doesn't seem to plot
fig = pstream.select(channel='??Q').normalize().plot_profile(scale=1.5, top='hist') # instead assign to variable and then call .show()
fig.show()
