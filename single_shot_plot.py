import numpy as np
from matplotlib import pyplot as plt
import utils
import obspy
import obspy.signal.trigger
import os


def plot_stream(st, title, fname):
    print(f"plotting stream to {fname}")
    plt.figure(figsize=(8, 8))
    plt.title(title)
    red_vel = 6.0
    yy = np.arange(len(st[0])) / st[0].stats.sampling_rate
    min_offset = float("inf")
    max_offset = float("-inf")
    for tr in st:
        offset = utils.source_receiver_offset(tr) / 1e3
        min_offset = min(offset, min_offset)
        max_offset = max(offset, max_offset)
        red_vel_shift = np.abs(offset) / red_vel
        x = tr.data + offset
        y = yy - red_vel_shift
        z = np.zeros(len(x)) + offset
        plt.fill_betweenx(y, x, z, where=x > z, facecolor="k", zorder=1)
        plt.plot(x, y, "k", linewidth=0.1)
    plt.xlabel("Offset (km)")
    plt.ylabel("Reduced Travel Time (w/ v=6.0 km/s)")
    plt.xlim(max_offset + 2, min_offset - 2)
    plt.ylim(-8, 8)
    plt.savefig(fname, dpi=250)


def plot_shot(shotno):
    basepath = os.path.join("single_shot_plots", str(shotno))
    st = utils.load_shot(shotno)
    utils.bandpass_stream_inplace(st)
    plot_stream(st, f"Shot {shotno} - Bandpass 3-20 Hz", f"{basepath}_bandpass.png")
    st_ = st.copy()
    for tr in st_:
        utils.agc(tr, 2 * 500)
    plot_stream(st_, f"Shot {shotno} - AGC 2s Window", f"{basepath}_agc.png")
    for tr in st:
        tr.data = obspy.signal.trigger.classic_sta_lta(tr.data, 0.05 * 500, 5.0 * 500)
        tr.data /= np.abs(tr.data).max()  # normalize
    plot_stream(st, f"Shot {shotno} - STA 0.05s LTA 5.0s", f"{basepath}_sta_lta.png")
