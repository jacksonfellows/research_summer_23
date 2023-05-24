import numpy as np
import scipy
from matplotlib import pyplot as plt
import obspy


def load_sample_stream():
    return obspy.read('c:/Users/lworthington/Downloads/8J.009.001.991.sgy')

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

def plot_stream(st):
    plt.gca().invert_xaxis()
    red_vel = 6.0
    yy = np.arange(len(st[0])) / st[0].stats.sampling_rate
    for i in range(len(st)):
        offset = st[i].stats.segy.trace_header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group / 1000
        red_vel_shift = np.abs(offset) / red_vel
        x = st[i].data + offset
        y = yy - red_vel_shift
        z = np.zeros(len(x)) + offset
        plt.fill_betweenx(y, x, z, where=x > z, facecolor='k', zorder=1)
        plt.plot(x, y, 'k', linewidth=0.1)
    plt.ylim(0, 8)
    plt.show()
