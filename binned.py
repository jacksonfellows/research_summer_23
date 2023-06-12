import numpy as np
from matplotlib import pyplot as plt


class BinnedTraces:
    """Store binned traces along with counts and some metadata."""

    def round_to_bin(self, offset):
        return np.rint(offset / self.bin_size) * self.bin_size

    def __init__(self, min_offset, max_offset, bin_size, sample_len, sampling_rate):
        self.bin_size = bin_size
        self.sampling_rate = sampling_rate
        self.offsets = np.arange(
            self.round_to_bin(min_offset),
            self.round_to_bin(max_offset) + bin_size,
            bin_size,
        )
        self.counts = np.zeros(len(self.offsets), dtype=int)
        self.binned = np.zeros((len(self.offsets), sample_len * sampling_rate))

    def get_bin_index(self, offset):
        """Get the bin index for offset."""
        return int((self.round_to_bin(offset) - self.offsets[0]) / self.bin_size)

    def add_trace(self, trace, offset):
        """Adds np.ndarray trace with offset to the binned traces."""
        j = self.get_bin_index(offset)
        self.counts[j] += 1
        self.binned[j] += trace

    def plot(
        self,
        red_vel=6.0,
        scale=1.0,
        show=True,
        ylim=(-2, 8),
    ):
        """Plot the binned traces with counts."""
        fig, axs = plt.subplots(2, 1, sharex=True, height_ratios=[0.8, 0.2])
        axs[0].set_xlim(self.offsets.max() + 2, self.offsets.min() - 2)
        axs[0].set_ylim(*ylim)
        axs[0].set_ylabel(f"Reduced time w/ v={red_vel:0.1f} km/s (s)")
        axs[1].set_xlabel("Offset (km)")
        axs[1].set_ylabel("# of traces")
        x = (
            scale * self.bin_size * self.binned / self.counts[:, np.newaxis]
            + self.offsets[:, np.newaxis]
        )
        y = np.subtract.outer(
            np.arange(self.binned.shape[1]) / self.sampling_rate, self.offsets / red_vel
        )
        z = np.zeros(x.shape) + self.offsets[:, np.newaxis]
        where = x > z
        for i in range(len(self.offsets)):
            try:
                axs[0].plot(x[i], y[:, i], "k", linewidth=0.1)
                axs[0].fill_betweenx(y[:, i], x[i], z[i], where=where[i], facecolor="k")
            except:
                pass
        axs[1].bar(self.offsets, self.counts)
        plt.tight_layout()
        if show:
            plt.show()
