import pyasdf
import time
import os
from xcorqc.xcorqc import xcorr2
from obspy import Stream, Trace
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from collections import defaultdict
from netCDF4 import Dataset
import numpy as np

code_start_time = time.time()


# =========================== User Input Required =========================== #

#Path to the ASDF file for xcor
data_path = '/g/data1/ha3/Passive/_AusArray/OA/processing/test_xcor_ashby/one_stn_xcor_test.h5'
asdf_out = '/g/data1/ha3/Passive/_AusArray/OA/processing/test_xcor_ashby/one_stn_xcor_out.h5'

temp_net = "OA"

# Interval over which windows are stacked
intervalSeconds     = 3600*24 # 1 day
# Size of data-window to be stacked
windowSeconds       = 3600

# == Output Path ==
outputDir           = '/g/data1/ha3/Passive/_AusArray/OA/processing/test_xcor_ashby/'

#decimation factor for temp stations data and refernece station data
tempDecFactor       = 5
refDecFactor        = 1


# =========================================================================== #


if os.path.exists(asdf_out):
    os.remove(asdf_out)

ds = pyasdf.ASDFDataSet(data_path)
print(ds.waveforms.list())

# dictionary with unique ids as keys and stream objects (for permanent station data) as values
id_st_dict = {}

# first extract obspy stream data from the ASDF file that is from a permanent station (i.e. not the temporary network)
for net_sta in ds.waveforms.list():
    iter_net = net_sta.split('.')[0]

    if not temp_net == iter_net:
        # get the station accessor
        iter_sta_accessor = ds.waveforms[net_sta]

        # get a list of the waveforms
        iter_waveforms = iter_sta_accessor.list()

        for wave in iter_waveforms:
            if not wave == "StationXML":
                iter_st = iter_sta_accessor[wave]
                # print(iter_st)

                # add it to the dict
                tag = wave.split('__')[3]

                id_st_dict[tag] = iter_st


# print(id_st_dict)

print('')


# function to make a number into a 3 digit string with leading zeros
def make_threedig(a):
    if len(a) == 1:
        return '00' + a
    elif len(a) == 2:
        return '0' + a
    return a




def xcor_process(st, inv):

    xcor_st = Stream()


    for tr in st:
        temp_st = Stream(traces=[tr])
        print('')
        print(tr.stats.asdf.labels)

        # get the uid label
        uid_label = tr.stats.asdf.labels[1]

        # get the region name
        region_id = tr.stats.asdf.labels[0]

        print(uid_label)


        perm_st = id_st_dict[uid_label]

        temp_tr = temp_st[0]
        ref_tr = perm_st[0]

        # decimate the traces
        temp_tr.decimate(tempDecFactor)
        ref_tr.decimate(refDecFactor)

        stationPair = ref_tr.stats.station + '.' + temp_tr.stats.station

        print(temp_st)
        print(perm_st)
        print("DO XCOR......")

        xcl, winsPerInterval = xcorr2(ref_tr, temp_tr, interval_seconds=intervalSeconds, window_seconds=windowSeconds)

        if (xcl is None):
            print("\t\tWarning: no cross-correlation results returned for station-pair %s, " %
                  (stationPair) + " due to gaps in data.")
            continue




        print(xcl)
        print(xcl.shape)
        print(winsPerInterval)
        print(winsPerInterval.shape)


        for _i, day_xcor in enumerate(xcl):
            print(day_xcor)

            # fill in headers for new xcor trace
            stats = {'network': tr.stats.network, 'station': tr.stats.station, 'location': "",
                     'channel': make_threedig(str(_i)), 'npts': len(xcl[0]), 'sampling_rate': tr.stats.sampling_rate,
                     'mseed': {'dataquality': 'D'}, "asdf": {}}


            stats["starttime"] = tr.stats.starttime

            temp_tr = Trace(data=day_xcor, header=stats)

            added_labels = [region_id, stationPair, str(winsPerInterval[_i])]

            print(added_labels)

            temp_tr.stats.asdf.labels = added_labels

            xcor_st += temp_tr

    return xcor_st



ds.process(process_function=xcor_process,
           output_filename=asdf_out,
           tag_map={"raw_recording": "xcor"})



tmp_results_dict = defaultdict(list)  # Results dictionary indexed by station-pair string
tmp_windows_count_dict = defaultdict(list)  # Window-count dictionary indexed by station-pair string

# load in new ASDF that has been written
new_ds = pyasdf.ASDFDataSet(asdf_out)

# extract all of the xcor
net_sta_list = new_ds.waveforms.list()

for net_sta in net_sta_list:
    sta_accessor = new_ds.waveforms[net_sta]

    tags_list = sta_accessor.list()



    for wave_id in tags_list:
        if wave_id == "StationXML":
            # ignore the station xml file
            continue

        #open up the trace
        tmp_st = sta_accessor[wave_id]
        tmp_tr = tmp_st[0]

        #get labels
        tmp_labels = tmp_tr.stats.asdf.labels

        print(tmp_labels)

        # make new key
        window_key = tmp_labels[1]+'_'+tmp_labels[0]

        tmp_results_dict[window_key].append(tmp_tr.data)
        tmp_windows_count_dict[window_key].append(int(tmp_labels[2]))







del ds

print(tmp_results_dict)
print(tmp_windows_count_dict)



print('=== Collating Results ===')
x = None
# Concatenate results
for k in tmp_results_dict.keys():
    combinedXcorrResults = None
    combinedWindowCountResults = None

    print(k)


    print(len(tmp_results_dict[k]))

    for i in np.arange(len(tmp_results_dict[k])):
        print(i)
        if (i == 0):
            combinedXcorrResults = np.array([tmp_results_dict[k][0]])
            combinedWindowCountResults = np.array([tmp_windows_count_dict[k][0]])

            print(combinedXcorrResults.shape)

            # Generate time samples (only needs to be done once)
            if (x is None):
                x = np.linspace(-windowSeconds, windowSeconds,
                                len(tmp_results_dict[k][0]))
                # end if
        else:
            combinedXcorrResults = np.concatenate((combinedXcorrResults,
                                                   np.array([tmp_results_dict[k][i]])), axis=0)

            combinedWindowCountResults = np.concatenate((combinedWindowCountResults,
                                                         np.array([
                                                             tmp_windows_count_dict[k][i]])), axis=0)
            # end if
    # end for

    # Replace lists with combined results
    tmp_results_dict[k] = combinedXcorrResults
    tmp_windows_count_dict[k] = combinedWindowCountResults

    print(combinedXcorrResults)
    print(combinedXcorrResults.shape)
    print(combinedWindowCountResults)
    print(combinedWindowCountResults.shape)
# end for



print(tmp_results_dict)
print(tmp_windows_count_dict)


# Save Results
for i, k in enumerate(tmp_results_dict.keys()):
    fn = os.path.join(outputDir, '%s.nc'%(k))
    print('\t Writing: %s'%(fn))

    root_grp = Dataset(fn, 'w', format='NETCDF4')
    root_grp.description = 'Cross-correlation results for station-pair: %s' % (k)

    # print(tmp_results_dict[k].shape)
    # Dimensions
    root_grp.createDimension('time', tmp_results_dict[k].shape[0])
    root_grp.createDimension('lag', tmp_results_dict[k].shape[1])

    # Variables
    time = root_grp.createVariable('time', 'f4', ('time',))
    lag = root_grp.createVariable('lag', 'f4', ('lag',))
    nsw = root_grp.createVariable('NumStackedWindows', 'f4', ('time',))
    xc = root_grp.createVariable('xcorr', 'f4', ('time', 'lag',))

    # Populate variables
    time[:] = np.arange(tmp_results_dict[k].shape[0])
    lag[:] = x
    nsw[:] = tmp_windows_count_dict[k]
    xc[:, :] = tmp_results_dict[k][::-1, :]  # Flipping rows
    root_grp.close()


    for day_array in tmp_results_dict[k][::-1, :]:
        print(day_array)
        print(day_array.shape)


        print(x.shape)
        print(tmp_results_dict[k][::-1, :].shape)



    # full_array = np.column_stack(tmp_results_dict[k][::-1, :])



    data_type = "XcorData"
    data_path = k.split(".")[1].split("_")[0] + "/" + k.replace(".", "_")

    # new_ds.add_auxiliary_data(data=tmp_results_dict[k][::-1, :], data_type=data_type, path=data_path, parameters={})

    data = np.random.normal(size=(3,100))

    # print(data)

    # data = np.concatenate(tmp_results_dict[k][::-1, :], np.array(1))

    # print(data.shape)
    new_ds.add_auxiliary_data(data=tmp_results_dict[k][::-1, :], data_type=data_type, path=data_path, parameters={})


del new_ds