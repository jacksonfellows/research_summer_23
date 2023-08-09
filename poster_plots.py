from functools import cache

import numpy as np
from matplotlib import pyplot as plt

import utils
from bin_shots import bin_all_shots_envelope_sta_lta
from inversion import plot_vm_and_rays
from megashot import megashot, megashot_all_nodes
from pickfile import load_picks
from profile_info import profile_1
from rayfile import Rayfile
from tt_curves import merge_pick_calc
from unbinned import UnbinnedTraces

# Set default font size.
plt.rcParams.update({"font.size": 24})

DPI = 300

tt_shots = [1000, 1050, 1287, 1350]


def megashot_comparison_plot(shotno, shots_per_side, min_v, max_v):
    single_bt = get_megashot_bt(shotno, 0, min_v, max_v)
    stacked_bt = get_megashot_bt(shotno, shots_per_side, min_v, max_v)
    red_vel = 6.0
    fig, (single_ax, stacked_ax, counts_ax) = plt.subplots(
        nrows=3, sharex=True, height_ratios=(0.5, 0.5, 0.1), figsize=(8.5, 8)
    )
    ylim = (0, 6)
    single_bt.plot_mat(
        red_vel=red_vel,
        vmin=0,
        show=False,
        fig=fig,
        axs=[single_ax, counts_ax],
        ylim=ylim,
    )
    stacked_bt.plot_mat(
        red_vel=red_vel, fig=fig, axs=[stacked_ax, counts_ax], show=False, ylim=ylim
    )
    single_ax.set_title("Before Stacking (1 shot)")
    stacked_ax.set_title(f"After Stacking ({shots_per_side} shots/side)")
    plt.tight_layout()
    plt.savefig("figures/megashot_comparison.png", dpi=DPI, bbox_inches="tight")


def megashot_single_plot(shotno, min_v, max_v):
    fig, (shots_ax, stacked_ax) = megashot(shotno, 9001, 1, min_v, max_v, show=False)
    shots_ax.set_xlim(6, 9)
    shots_ax.set_title("Individual")
    stacked_ax.set_title("Stacked")
    shots_ax.set_ylim(-0.1, 1.1)
    stacked_ax.set_ylim(-0.1, 1.1)
    plt.tight_layout()
    fig.set_size_inches(8, 3.5)
    shots_ax.legend(
        loc="lower right",
        fontsize=20,
        columnspacing=0.2,
        handletextpad=0.2,
        borderpad=0.2,
        labelspacing=0.2,
    )
    plt.savefig("figures/megashot_single.png", dpi=DPI, bbox_inches="tight")


def model_plots():
    fig, (model_ax, invert_ax) = plt.subplots(nrows=2, figsize=(13, 10), sharex=True)
    cbar_ax = fig.add_axes((0.92, 0.2, 0.05, 0.6))
    plot_vm_and_rays(
        "poster_inversion/vm_000",
        "poster_inversion/rayfile_000",
        show=False,
        ax=model_ax,
        cbar_ax=cbar_ax,
        cbar_orientation="vertical",
    )
    plot_vm_and_rays(
        "poster_inversion/vm_008",
        "poster_inversion/rayfile_008",
        "poster_inversion/dws_vel_009",
        show=False,
        ax=invert_ax,
        cbar_ax=cbar_ax,
        cbar_orientation="vertical",
    )

    tt_shot_locs = profile_1.project_shots(tt_shots)

    for ax in [model_ax, invert_ax]:
        ax.plot(
            tt_shot_locs.x, tt_shot_locs.z, "*", zorder=102, color="gold", markersize=10
        )

    model_ax.set_title("Initial Model")
    invert_ax.set_title("Inverted Model")
    model_ax.set_xlabel(None)
    plt.savefig("figures/model.png", dpi=DPI, bbox_inches="tight")


@cache
def get_megashot_bt(shotno, shots_per_side, min_v, max_v):
    return megashot_all_nodes(shotno, shots_per_side, min_v, max_v, "median")


def plot_tt_curve(
    shotno,
    shots_per_side,
    initial_rayfile,
    final_rayfile,
    min_v,
    max_v,
    red_vel=6.0,
    ylim=(-2, 8),
    show=True,
    fig=None,
    ax=None,
):
    shot_pick = load_picks(f"picks/megashot_{shotno}_{shots_per_side}_01")
    all_tts = []
    for rayfile in [initial_rayfile, final_rayfile]:
        shot_calc = Rayfile.load(rayfile).calculated_tts()
        shot_calc = shot_calc[shot_calc.shot == shotno]
        all_tts.append(merge_pick_calc(shot_pick, shot_calc))

    initial_tts, final_tts = all_tts

    if fig is None or ax is None:
        fig, ax = plt.subplots()

    bt = get_megashot_bt(shotno, shots_per_side, min_v, max_v)
    bt.plot_mat(red_vel=red_vel, ylim=ylim, fig=fig, axs=[ax, None], show=False)

    ax.set_xlabel("Offset (km)")

    ax.errorbar(
        initial_tts.offset,
        initial_tts.tt - initial_tts.offset / red_vel,
        yerr=initial_tts.error,
        capsize=5,
        elinewidth=2,
        markeredgewidth=2,
        color="red",
        linestyle="none",
        marker="x",
        label="Picked",
    )

    ax.plot(
        initial_tts.offset,
        initial_tts.tt_calc - initial_tts.offset / red_vel,
        "x",
        markeredgewidth=2,
        color="purple",
        label="Initial",
    )

    ax.plot(
        final_tts.offset,
        final_tts.tt_calc - final_tts.offset / red_vel,
        "x",
        markeredgewidth=2,
        color="blue",
        label="Inverted",
        zorder=100,
    )

    plt.legend()

    if show:
        plt.show()


plt.set_cmap("viridis_r")


def plot_tt_curves():
    shots = tt_shots
    fig, axs = plt.subplots(nrows=2, ncols=2, sharey="row", figsize=(13, 7))
    axs = axs.flatten()
    for i, shotno in enumerate(shots):
        min_v = 5000 if shotno < 1200 else 7000
        max_v = 7000 if shotno < 1200 else 9000
        ylim = (-1, 5) if shotno < 1200 else (-2, 4)
        plot_tt_curve(
            shotno,
            5,
            "poster_inversion/rayfile_000",
            "poster_inversion/rayfile_008",
            min_v=min_v,
            max_v=max_v,
            ylim=ylim,
            show=False,
            fig=fig,
            ax=axs[i],
        )
        axs[i].set_title(f"Shot {shotno}")

    axs[0].set_xlabel(None)
    axs[1].set_xlabel(None)

    axs[1].set_ylabel(None)
    axs[3].set_ylabel(None)

    axs[0].set_yticks([-1, 1, 3, 5])
    axs[2].set_yticks([-2, 0, 2, 4])

    plt.legend(
        loc="upper left",
        fontsize=20,
        columnspacing=0.2,
        handletextpad=0.2,
        borderpad=0.2,
        labelspacing=0.2,
    )

    plt.tight_layout()
    plt.savefig("figures/tt_curves.png", dpi=DPI, bbox_inches="tight")


@cache
def get_line_binned(lineno):
    return bin_all_shots_envelope_sta_lta(utils.shots_for_line(lineno))


def plot_binned_1_2():
    l1 = get_line_binned(1)
    l2 = get_line_binned(2)

    fig, axs = plt.subplots(
        2, 2, sharex="col", height_ratios=[0.9, 0.1], figsize=(17, 4)
    )

    l1.plot_mat(
        fig=fig,
        axs=[axs[0, 0], axs[1, 0]],
        vmin=0.07,
        vmax=0.2,
        ylim=(-4, 6),
        show=False,
    )

    l2.plot_mat(
        fig=fig,
        axs=[axs[0, 1], axs[1, 1]],
        vmin=0.07,
        vmax=0.25,
        ylim=(-8, 4),
        show=False,
    )

    axs[0, 0].set_yticks([-4, -2, 0, 2, 4, 6])
    axs[0, 1].set_yticks([-8, -6, -4, -2, 0, 2, 4])

    axs[0, 0].set_title("Line 1")
    axs[0, 1].set_title("Line 2")

    plt.savefig("figures/l1_l2_binned.png", dpi=DPI, bbox_inches="tight")


def plot_raw():
    shotno = 1000
    st = utils.load_shot(shotno)
    utils.bandpass_stream_inplace(st)
    offsets = [1e-3 * utils.source_receiver_offset(tr) for tr in st]
    traces = UnbinnedTraces(min(offsets), max(offsets), st[0].stats.sampling_rate)
    for offset, tr in zip(offsets, st):
        tr.data = tr.data / np.max(np.abs(tr.data))
        traces.add_trace(tr.data, offset)
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.set_yticks([-1, 1, 3, 5])
    traces.plot(show=False, scale=0.5, ylim=(-1, 5), fig=fig, ax=ax)
    plt.title(f"Shot {shotno}")
    plt.savefig(f"figures/raw_shot_{shotno}", dpi=DPI, bbox_inches="tight")
    plt.close()

    shotno = 1350
    st = utils.load_shot(shotno)
    utils.bandpass_stream_inplace(st)
    offsets = [1e-3 * utils.source_receiver_offset(tr) for tr in st]
    traces = UnbinnedTraces(min(offsets), max(offsets), st[0].stats.sampling_rate)
    for offset, tr in zip(offsets, st):
        tr.data = tr.data / np.max(np.abs(tr.data))
        traces.add_trace(tr.data, offset)
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.set_yticks([-2, 0, 2, 4])
    traces.plot(show=False, scale=0.5, ylim=(-2, 4), fig=fig, ax=ax)
    plt.title(f"Shot {shotno}")
    plt.savefig(f"figures/raw_shot_{shotno}", dpi=DPI, bbox_inches="tight")
    plt.close()
