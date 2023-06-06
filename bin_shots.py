# TODO: clean this up
# Probably should store the sampling rate in the BinnedTraces object.

import numpy as np
import obspy
import obspy.signal.trigger
from matplotlib import pyplot as plt
import copy
import pickle
import itertools

import utils

class BinnedTraces():
    def round_to_bin(self, x):
        return np.rint(x / self.bin_size) * self.bin_size

    def __init__(self, min_offset, max_offset, bin_size, trace_len):
        self.bin_size = bin_size
        self.offsets = np.arange(self.round_to_bin(min_offset), self.round_to_bin(max_offset) + bin_size, bin_size)
        self.counts = np.zeros(len(self.offsets), dtype=int)
        self.binned = np.zeros((len(self.offsets), trace_len))

    def add_stream(self, st):
        for tr in st:
            offset = utils.source_receiver_offset(tr) / 1000
            j = int((self.round_to_bin(offset) - self.offsets[0]) / self.bin_size)
            self.counts[j] += 1
            self.binned[j] += tr.data

    def plot(self, sampling_rate, red_vel, scale=1.0, show=True, div_by_counts=True, ylim=None):
        fig, axs = plt.subplots(2, 1, sharex=True, height_ratios=[0.8,0.2])
        axs[0].set_xlim(self.offsets.max() + 2, self.offsets.min() - 2)
        if ylim is not None:
            axs[0].set_ylim(*ylim)
        else:
            axs[0].set_ylim(-2, 8)
        axs[0].set_ylabel(f'Reduced time w/ v={red_vel:0.1f} km/s (s)')
        axs[1].set_xlabel('Offset (km)')
        axs[1].set_ylabel('# of traces')
        if div_by_counts:
            x = (scale * self.bin_size * self.binned / self.counts[:,np.newaxis] + self.offsets[:,np.newaxis])
        else:
            x = (scale * self.bin_size * self.binned + self.offsets[:,np.newaxis])
        y = np.subtract.outer(np.arange(self.binned.shape[1]) / sampling_rate, self.offsets / red_vel)
        z = np.zeros(x.shape) + self.offsets[:,np.newaxis]
        where = x > z
        for i in range(len(self.offsets)):
            try:
                axs[0].plot(x[i], y[:,i], 'k', linewidth=0.1)
                axs[0].fill_betweenx(y[:,i], x[i], z[i], where=where[i], facecolor='k')
            except:
                pass
        axs[1].bar(self.offsets, self.counts)
        plt.tight_layout()
        if show:
            plt.show()

    def save(self, path):
        with open(path, 'wb') as f:
            pickle.dump(self, f)

def load_binned_traces(path):
    with open(path, 'rb') as f:
        return pickle.load(f)

def min_dist_to_square(lat, lon, min_lat, max_lat, min_lon, max_lon):
    assert not (min_lat < lat < max_lat and min_lon < lon < max_lon) # assume lat, lon outside of square
    corners = list(itertools.product((min_lat, max_lat), (min_lon, max_lon)))
    sides = []
    if min_lat < lat < max_lat:
        sides += [(lat, min_lon), (lat, max_lon)]
    if min_lon < lon < max_lon:
        sides += [(min_lat, lon), (max_lat, lon)]
    return min(obspy.geodetics.gps2dist_azimuth(lat, lon, lat_, lon_)[0] for lat_, lon_ in corners + sides)

def max_dist_to_square(lat, lon, min_lat, max_lat, min_lon, max_lon):
    assert not (min_lat < lat < max_lat and min_lon < lon < max_lon) # assume lat, lon outside of square
    corners = list(itertools.product((min_lat, max_lat), (min_lon, max_lon)))
    return max(obspy.geodetics.gps2dist_azimuth(lat, lon, lat_, lon_)[0] for lat_, lon_ in corners)

def calc_min_max_offsets_km(shot_nos):
    min_lat, max_lat, min_lon, max_lon = float('inf'), float('-inf'), float('inf'), float('-inf')
    for stat in utils.stat_df.index:
        lat, lon = utils.stat_lat_lon(stat)
        min_lat, max_lat, min_lon, max_lon = min(lat, min_lat), max(lat, max_lat), min(lon, min_lon), max(lon, max_lon)
    min_offset, max_offset = float('inf'), float('-inf')
    for shot in shot_nos:
        try:
            lat, lon = utils.shot_lat_lon(shot)
            try:
                min_offset = min(min_offset, min_dist_to_square(lat, lon, min_lat, max_lat, min_lon, max_lon))
                max_offset = max(max_offset, max_dist_to_square(lat, lon, min_lat, max_lat, min_lon, max_lon))
            except ValueError:
                pass
        except KeyError:
            pass
    return min_offset / 1e3, max_offset / 1e3

def bin_all_shots_sta_lta(shot_nos):
    min_offset, max_offset = calc_min_max_offsets_km(shot_nos)
    min_offset, max_offset = np.floor(min_offset), np.ceil(max_offset)
    print(f'min_offset = {min_offset}, max_offset = {max_offset}')
    bt = BinnedTraces(min_offset, max_offset, 0.25, 60 * 500)
    for shotno in shot_nos:
        print(f'loading shot {shotno}')
        st = utils.load_shot(shotno)
        print(f'applying bandpass to shot {shotno}')
        utils.bandpass_stream_inplace(st)
        print(f'applying sta-lta to shot {shotno}')
        for t in st:
            t.data = obspy.signal.trigger.classic_sta_lta(t.data, 0.05 * 500, 5.0 * 500)
            t.data /= np.abs(t.data).max() # normalize
            assert not np.isnan(t.data).any()
        print(f'adding shot {shotno} to binned traces')
        bt.add_stream(st)
    return bt

def rebin(old_bt, new_bin_size):
    assert(new_bin_size > old_bt.bin_size)
    min_offset, max_offset = old_bt.offsets[0], old_bt.offsets[-1]
    new_bt = BinnedTraces(min_offset, max_offset, new_bin_size, old_bt.binned.shape[1])
    for i in range(len(old_bt.offsets)):
        j = int((new_bt.round_to_bin(old_bt.offsets[i]) - new_bt.offsets[0]) / new_bt.bin_size)
        new_bt.counts[j] += old_bt.counts[i]
        new_bt.binned[j] += old_bt.binned[i]
    return new_bt

def apply_to_binned(bt, func):
    new_bt = copy.deepcopy(bt)
    for i in range(new_bt.binned.shape[0]):
        new_bt.binned[i] = func(new_bt.binned[i])
    return new_bt

def plot_trace(t):
    x = np.arange(t.shape[0])
    plt.plot(x, t, 'k', linewidth=0.1)
    y0 = np.zeros(t.shape)
    plt.fill_between(x, t, y0, where=t>0, facecolor='k')
    plt.show()

def sta_lta_and_normalize(t, sta_samples, lta_samples):
    t_ = obspy.signal.trigger.classic_sta_lta(t, sta_samples, lta_samples)
    t_ /= t_.max()
    return t_
