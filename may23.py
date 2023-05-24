import numpy as np
import scipy
from matplotlib import pyplot as plt
import obspy
import os
import glob
import pandas as pd


def load_sample_stream():
    return obspy.read('c:/Users/lworthington/Downloads/8J.009.001.991.sgy', unpack_trace_headers=True)

def compare_arrs(*arrs):
    fig, axs = plt.subplots(len(arrs), 1, sharex=True)
    for i,arr in enumerate(arrs):
        axs[i].plot(arr)
    plt.show()

def my_agc(tr, window_len):
    l = int(window_len)
    tr.data /= scipy.ndimage.maximum_filter1d(np.abs(tr.data), l, origin=-l//2, mode='nearest')

# Lindsay's code:
# AGC filter functions
# Use Numba library to compile function in C for better performance
from numba import jit
@jit(nopython=True)
# AGC filter takes [data] as array
# and [agc_len] as filter length in samples
def lindsays_agc(data,agc_len):
    data_out = data
    for i in np.arange(0,len(data_out)):
        s = int(i)
        e = int(i+agc_len)
        # Dealing with very end of trace
        if np.max(np.abs(data_out[s:e])) > 0:
            data_out[s] /= np.max(np.abs(data_out[s:e]))
    return data_out

def test():
    a = st[0].copy()
    a.filter('bandpass', freqmin=2, freqmax=20)
    b = a.copy()
    c = a.copy()
    agc_len = 2.0 * st[0].stats.sampling_rate
    my_agc(b, agc_len)
    c.data = lindsays_agc(c.data, agc_len)
    compare_arrs(st[0], a, b, c, b.data - c.data)

def process_stream(ref_st, window_len_s = 2.0):
    st = ref_st.copy()
    for i in range(len(st)):
        st[i].filter('bandpass', freqmin=3, freqmax=20, zerophase=True)
        my_agc(st[i], window_len_s * st[i].stats.sampling_rate)
    return st

def process_stream_lindsay(ref_st):
    st = ref_st.copy()
    for i in range(len(st)):
        st[i].filter('bandpass', freqmin=3, freqmax=20, zerophase=True)
        st[i].data = lindsays_agc(st[i].data, 2.0 * st[i].stats.sampling_rate)
    return st

# Timing: For one sample stream process_stream took ~0.5 s and process_stream_lindsay took ~14 s.

def process_stream_2(ref_st):
    st = ref_st.copy()
    for i in range(len(st)):
        st[i].data /= np.max(np.abs(st[i].data))
    return st

def plot_stream(st, show=True):
    red_vel = 6.0
    yy = np.arange(len(st[0])) / st[0].stats.sampling_rate
    min_offset = float('inf')
    max_offset = float('-inf')
    for i in range(len(st)):
        offset = st[i].stats.segy.trace_header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group / 1000
        min_offset = min(offset, min_offset)
        max_offset = max(offset, max_offset)
        red_vel_shift = np.abs(offset) / red_vel
        x = st[i].data + offset
        y = yy - red_vel_shift
        z = np.zeros(len(x)) + offset
        plt.fill_betweenx(y, x, z, where=x > z, facecolor='k', zorder=1)
        plt.plot(x, y, 'k', linewidth=0.1)
    plt.xlim(max_offset + 2, min_offset - 2)
    plt.ylim(0, 8)
    if show:
        plt.show()

def load_shot(shotno):
    shotno = str(shotno)
    shot_dir = os.path.join('shots', shotno)
    stream = obspy.core.stream.Stream()
    for path in glob.glob(os.path.join(shot_dir, '*.sgy')):
        stream += obspy.read(path, unpack_trace_headers=True)
    return stream

def load_station_df():
    return pd.read_csv('station_lat_lon_elev.csv', index_col=1)

stat_df = load_station_df()

def stat_lat_lon(stat_code):
    info = stat_df.loc[stat_code]
    return info.lat, info.lon

def load_shot_df():
    return pd.read_excel('PASSCAL_SHOT.xlsx', names=['lineno', 'shotno', 'lat', 'lon', 'unknown', 'time', 'linename'], index_col=1, header=None, usecols='A:G')

shot_df = load_shot_df()

def shot_lat_lon(shotno):
    info = shot_df.loc[shotno]
    return info.lat, -info.lon  # There seems to be a problem with the shot locations. For now just fix it here.

def source_receiver_offset(tr):
    # Ignore elevation of receiver.
    stat_code = tr.stats.segy.trace_header.trace_number_within_the_original_field_record
    shotno = tr.stats.segy.trace_header.energy_source_point_number
    dist_m, _, _ = obspy.geodetics.base.gps2dist_azimuth(*stat_lat_lon(stat_code), *shot_lat_lon(shotno))
    return dist_m

def fix_offsets(old_st):
    st = obspy.core.Stream()
    for i in range(len(old_st)):
        try:
            tr = old_st[i].copy()
            tr.stats.segy.trace_header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group = source_receiver_offset(tr)
            st += tr
        except:
            pass
    return st

def plot_all_shots():
    for shotno in os.listdir('shots'):
        print(f'loading shot {shotno}')
        raw_stream = load_shot(shotno)
        print('...processing')
        processed_stream = process_stream(fix_offsets(raw_stream))
        print('...ploting')
        plt.figure(figsize=(12,9))
        plot_stream(processed_stream, show=False)
        fig_path = os.path.join('shot_plots', f'shot_{shotno}.png')
        plt.savefig(fig_path, dpi=300)
