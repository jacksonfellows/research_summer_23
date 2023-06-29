# TODO: clean this up
# Probably should store the sampling rate in the BinnedTraces object.

import numpy as np
import obspy
import obspy.signal.trigger
import obspy.signal.filter
from matplotlib import pyplot as plt
import pickle

import utils
import binned


def load_binned_traces(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def bin_all_shots_sta_lta(shot_nos):
    min_offset, max_offset = utils.calc_min_max_offsets_km(shot_nos)
    min_offset, max_offset = np.floor(min_offset), np.ceil(max_offset)
    print(f"min_offset = {min_offset}, max_offset = {max_offset}")
    bt = binned.BinnedTraces(min_offset, max_offset, 0.25, 60, 500)
    for shotno in shot_nos:
        print(f"loading shot {shotno}")
        st = utils.load_shot(shotno)
        print(f"applying bandpass to shot {shotno}")
        utils.bandpass_stream_inplace(st)
        print(f"applying sta-lta to shot {shotno}")
        for t in st:
            t.data = obspy.signal.trigger.classic_sta_lta(t.data, 0.05 * 500, 5.0 * 500)
            t.data /= np.abs(t.data).max()  # normalize
            if np.isnan(t.data).any():
                print(f"nan! - skipping trace {t}")
                st.remove(t)
        print(f"adding shot {shotno} to binned traces")
        for t in st:
            offset = utils.source_receiver_offset(t) / 1e3
            bt.add_trace(t.data, offset)
    return bt


def bin_all_shots_envelope_sta_lta(shot_nos):
    min_offset, max_offset = utils.calc_min_max_offsets_km(shot_nos)
    min_offset, max_offset = np.floor(min_offset), np.ceil(max_offset)
    print(f"min_offset = {min_offset}, max_offset = {max_offset}")
    bt = binned.BinnedTraces(min_offset, max_offset, 0.25, 60, 500)
    for shotno in shot_nos:
        print(f"loading shot {shotno}")
        st = utils.load_shot(shotno)
        print(f"applying bandpass to shot {shotno}")
        utils.bandpass_stream_inplace(st)
        print(f"applying envelope and sta-lta to shot {shotno}")
        for t in st:
            t.data = obspy.signal.filter.envelope(t.data)
            t.data = obspy.signal.trigger.classic_sta_lta(t.data, 0.05 * 500, 5.0 * 500)
            t.data /= np.abs(t.data).max()  # normalize
            if np.isnan(t.data).any():
                print(f"nan! - skipping trace {t}")
                st.remove(t)
        print(f"adding shot {shotno} to binned traces")
        for t in st:
            offset = utils.source_receiver_offset(t) / 1e3
            bt.add_trace(t.data, offset)
    return bt


def plot_trace(t):
    x = np.arange(t.shape[0])
    plt.plot(x, t, "k", linewidth=0.1)
    y0 = np.zeros(t.shape)
    plt.fill_between(x, t, y0, where=t > 0, facecolor="k")
    plt.show()


def sta_lta_and_normalize(t, sta_samples, lta_samples):
    t_ = obspy.signal.trigger.classic_sta_lta(t, sta_samples, lta_samples)
    t_ /= t_.max()
    return t_
