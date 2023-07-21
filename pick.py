import code
import datetime
import os
import tempfile
from copy import deepcopy

import numpy as np
import obspy.signal.trigger
import pandas as pd
from matplotlib import pyplot as plt

import binned
import broadband
import megashot
import profile_info
import rayfile
import utils
import wrappers


def threshold_picker(bt, threshold):
    return pd.DataFrame(
        {
            "offset": bt.offsets,
            "tt": np.argmax(bt.binned > threshold, axis=1) / bt.sampling_rate,
        }
    )


def threshold_range_picker(bt, threshold, min_red_tt, max_red_tt, red_vel=6.0):
    yy = np.subtract.outer(
        np.arange(bt.binned.shape[1]) / bt.sampling_rate, bt.offsets / red_vel
    ).T
    in_bounds = (min_red_tt < yy) & (yy < max_red_tt)
    return pd.DataFrame(
        {
            "offset": bt.offsets,
            "tt": np.argmax((bt.binned > threshold) & in_bounds, axis=1)
            / bt.sampling_rate,
        }
    )


def max_picker(bt):
    return pd.DataFrame(
        {"offset": bt.offsets, "tt": np.argmax(bt.binned, axis=1) / bt.sampling_rate}
    )


def sta_lta_binned(old_bt, st_s, lt_s):
    # a little silly
    bt = deepcopy(old_bt)
    # now safe to modify
    for row in range(bt.binned.shape[0]):
        i = np.argmax(bt.binned[row] > 0)
        bt.binned[row, i:] = obspy.signal.trigger.classic_sta_lta(
            bt.binned[row, i:], st_s * bt.sampling_rate, lt_s * bt.sampling_rate
        )
    return bt


def pick(
    bt,
    picks_df=None,
    pick_error=0.1,
    mode="matrix",
    red_vel=6.0,
    ylim=(-2, 8),
    vmin=None,
    vmax=None,
    scale=1,
):
    if mode == "matrix":
        fig, axs = bt.plot_mat(
            show=False, red_vel=red_vel, ylim=ylim, vmin=vmin, vmax=vmax
        )
    elif mode == "squiggle":
        fig, axs = bt.plot(show=False, red_vel=red_vel, ylim=ylim, scale=scale)
    else:
        assert False

    ax_overlay = axs[0]._make_twin_axes(sharex=axs[0], sharey=axs[0])
    ax_overlay.patch.set_visible(False)

    xs = []
    ys = []
    ts = []
    errors = []

    def redraw():
        ax_overlay.clear()
        ax_overlay.errorbar(
            xs,
            ys,
            yerr=errors,
            capsize=5,
            elinewidth=2,
            markeredgewidth=2,
            color="red",
            linestyle="none",
            marker="x",
        )
        fig.canvas.update()

    if picks_df is not None:
        red_tts = picks_df.tt - picks_df.offset / red_vel
        in_bounds = (ylim[0] < red_tts) & (red_tts < ylim[1])
        print(
            f"discarding {in_bounds.shape[0] - np.count_nonzero(in_bounds)} picks out of the visible bounds of the plot"
        )
        xs = list(picks_df.offset[in_bounds])
        ys = list(red_tts[in_bounds])
        ts = list(picks_df.tt[in_bounds])
        errors = list(picks_df.error[in_bounds])
        redraw()

    def onclick_callback(event):
        offset = bt.round_to_bin(event.xdata)
        red_time = event.ydata
        time = red_time + offset / red_vel
        if event.key == "shift":
            if offset in xs:
                i = xs.index(offset)
                xs.pop(i)
                ys.pop(i)
                ts.pop(i)
                errors.pop(i)
            xs.append(offset)
            ys.append(red_time)
            ts.append(time)
            errors.append(pick_error)
            print(f"picked {offset:0.2f} km, {time:0.2f} s, err={pick_error:0.2f}")
            redraw()
        elif event.key == "control":
            i = xs.index(offset)
            xs.pop(i)
            ys.pop(i)
            time = ts.pop(i)
            errors.pop(i)
            print(f"deleted pick {offset:0.2f} km, {time:0.2f} s")
            redraw()

    cid = fig.canvas.mpl_connect("button_press_event", onclick_callback)

    plt.show()

    return pd.DataFrame({"offset": xs, "tt": ts, "error": errors})


def fake_station_number(x_km):
    """Create a fake station number from offset along a profile."""
    # Supports 1-meter precision and up to 999.999 km offset.
    return int(1e6 + 1e3 * x_km)


def fake_station_offset_km(fake_station_number):
    """Get the offset for a fake station number."""
    return 1e-3 * (fake_station_number - 1e6)


def make_pick_file(profile, shotno, phase, picks):
    shot_x_km = profile.project_shots([shotno]).iloc[0].x
    return pd.DataFrame(
        {
            "station": [
                # round to nearest 250 m
                fake_station_number(x)
                for x in np.rint((shot_x_km - picks.offset) / 0.25) * 0.25
            ],
            "shot": shotno,
            "phase": phase,
            "offset": picks.offset,
            "tt": picks.tt,
            "error": picks.error,
        }
    )


def save_picks(picks, filename):
    picks.to_csv(
        filename,
        columns=["station", "shot", "phase", "offset", "tt", "error"],
        sep=" ",
        header=False,
        index=False,
    )


def load_picks(filename):
    return pd.read_csv(
        filename,
        names=["station", "shot", "phase", "offset", "tt", "error"],
        sep=" ",
        header=None,
    )


def pick_megashot(profile, shotno, shots_per_side, min_v, max_v, stacking="mean"):
    pick_file = os.path.join("picks", f"megashot_{shotno}_{shots_per_side}_01")

    picks = None
    if os.path.exists(pick_file):
        print(f"Pick file {pick_file} already exists - loading existing picks!")
        picks = load_picks(pick_file)

    bt = megashot.megashot_all_nodes(shotno, shots_per_side, min_v, max_v, stacking)

    pick_error = 0.1

    def reshoot(min_v, max_v):
        nonlocal bt
        bt = megashot.megashot_all_nodes(shotno, shots_per_side, min_v, max_v, stacking)

    def guess_picks(min_red_tt, max_red_tt, st_s=0.002, lt_s=0.2, threshold=4.0):
        nonlocal picks
        picks = threshold_range_picker(
            sta_lta_binned(bt, st_s, lt_s), threshold, min_red_tt, max_red_tt
        )

    def _save_picks(phase):
        print(f"saving {len(picks)} picks to {pick_file} with {phase=}")
        save_picks(make_pick_file(profile, shotno, phase, picks), pick_file)

    def plot_squiggle(scale=2, **kwargs):
        nonlocal picks
        picks = pick(
            bt, picks, pick_error=pick_error, mode="squiggle", scale=scale, **kwargs
        )

    def plot_matrix(**kwargs):
        nonlocal picks
        picks = pick(bt, picks, pick_error=pick_error, mode="matrix", **kwargs)

    def plot_raw(ylim=(0, 60), **kwargs):
        bt.plot_raw_mat(ylim, **kwargs)

    code.interact(
        local={
            "guess": guess_picks,
            "plot": plot_squiggle,
            "plot_mat": plot_matrix,
            "save": _save_picks,
            "reshoot": reshoot,
            "plot_raw": plot_raw,
        }
    )


def broadband_station_code_to_number(station_code):
    return {
        "KD01": 1,
    }[station_code]


def make_broadband_pick_file(station_code, phase, picks):
    # Need to recover shot numbers.
    shot_offsets = {
        1e-3
        * obspy.geodetics.gps2dist_azimuth(
            *utils.broadband_lat_lon(station_code), *utils.shot_lat_lon(shotno)
        )[0]: shotno
        for shotno in utils.shots_for_line(1)
    }
    shots = [shot_offsets[offset] for offset in picks.offset]
    return pd.DataFrame(
        {
            "station": broadband_station_code_to_number(station_code),
            "shot": shots,
            "phase": phase,
            "offset": picks.offset,
            "tt": picks.tt,
            "error": picks.error,
        }
    )


def pick_broadband(station_code):
    # Assume profile 1

    pick_file = os.path.join("picks", f"broadband_{station_code}_01")

    picks = None
    if os.path.exists(pick_file):
        print(f"Pick file {pick_file} already exists - loading existing picks!")
        picks = load_picks(pick_file)

    traces = broadband.make_unbinned_traces(station_code, 1)

    pick_error = 0.04

    def plot_squiggle(**kwargs):
        nonlocal picks
        picks = pick(traces, picks, pick_error=pick_error, mode="squiggle", **kwargs)

    def _save_picks(phase):
        print(f"saving {len(picks)} picks to {pick_file} with {phase=}")
        save_picks(make_broadband_pick_file(station_code, phase, picks), pick_file)

    code.interact(
        local={
            "plot": plot_squiggle,
            "save": _save_picks,
        }
    )
