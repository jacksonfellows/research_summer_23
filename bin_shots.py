import numpy as np
import obspy
from matplotlib import pyplot as plt

import utils

class BinnedTraces():
    def round_to_bin(self, x):
        return np.round(x / self.bin_size) * self.bin_size

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

    def plot(self, sampling_rate, red_vel, scale=1.0, show=True):
        fig, axs = plt.subplots(2, 1, sharex=True, height_ratios=[0.8,0.2])
        axs[0].set_xlim(self.offsets.max() + 2, self.offsets.min() - 2)
        axs[0].set_ylim(0, 8)
        axs[0].set_ylabel(f'Reduced time w/ v={red_vel:0.1f} km/s (s)')
        axs[1].set_xlabel('Offset (km)')
        axs[1].set_ylabel('# of traces')
        x = (scale * self.bin_size * self.binned / self.counts[:,np.newaxis] + self.offsets[:,np.newaxis])
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

def bin_all_shots():
    bt = BinnedTraces(28, 240, 0.1, 60 * 500)
    for shotno in range(991, 1392):
        print(f'loading shot {shotno}')
        st = utils.load_shot(shotno)
        print(f'processing shot {shotno}')
        utils.process_stream_inplace(st)
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
