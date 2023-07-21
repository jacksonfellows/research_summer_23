import numpy as np
from matplotlib import pyplot as plt


class UnbinnedTraces:
    """Like BinnedTraces but not binned."""

    def __init__(self, min_offset, max_offset, sampling_rate):
        self.sampling_rate = sampling_rate
        # Keep two parallel lists of offset,trace
        self.offsets = []
        self.traces = []

    def add_trace(self, trace, offset):
        self.offsets.append(offset)
        self.traces.append(trace)

    def plot(self, red_vel=6.0, scale=1.0, show=True, ylim=(-2, 8)):
        fig, ax = plt.subplots()
        ax.set_xlim(max(self.offsets) + 2, min(self.offsets) - 2)
        ax.set_ylim(ylim)
        ax.set_xlabel("Offset (km)")
        ax.set_ylabel(f"Reduced time w/ v={red_vel:0.1f} km/s (s)")
        for offset, trace in zip(self.offsets, self.traces):
            x = scale * trace + offset
            y = np.arange(x.shape[0]) / self.sampling_rate - offset / red_vel
            z = np.full(x.shape, offset)
            ax.plot(x, y, "k", linewidth=0.1)
            ax.fill_betweenx(y, x, z, where=x > z, facecolor="k")
        if show:
            plt.show()
        else:
            return fig, [ax]

    def round_to_bin(self, offset):
        # Really round to trace.
        return min(self.offsets, key=lambda x: abs(offset - x))
