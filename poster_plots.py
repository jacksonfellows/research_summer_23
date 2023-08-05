from functools import cache

from matplotlib import pyplot as plt

from inversion import plot_vm_and_rays
from megashot import megashot, megashot_all_nodes
from pickfile import load_picks
from rayfile import Rayfile
from tt_curves import merge_pick_calc

# Set default font size.
plt.rcParams.update({"font.size": 24})

DPI = 300


def megashot_comparison_plot(shotno, shots_per_side, min_v, max_v):
    single_bt = megashot_all_nodes(shotno, 0, min_v, max_v, "median")
    stacked_bt = megashot_all_nodes(shotno, shots_per_side, min_v, max_v, "median")
    red_vel = 6.0
    fig, (single_ax, stacked_ax, counts_ax) = plt.subplots(
        nrows=3, sharex=True, height_ratios=(0.4, 0.4, 0.2), figsize=(8.5, 12)
    )
    single_bt.plot_mat(
        red_vel=red_vel, vmin=0, show=False, fig=fig, axs=[single_ax, counts_ax]
    )
    _, axs = stacked_bt.plot_mat(
        red_vel=red_vel, fig=fig, axs=[stacked_ax, counts_ax], show=False
    )
    single_ax.set_title("Before Stacking (1 shot)")
    stacked_ax.set_title(f"After Stacking ({shots_per_side} shots/side)")
    plt.tight_layout()
    plt.savefig("figures/megashot_comparison.png", dpi=DPI)


def megashot_single_plot(shotno, min_v, max_v):
    fig, (shots_ax, stacked_ax) = megashot(shotno, 9001, 1, min_v, max_v, show=False)
    shots_ax.set_xlim(6, 9)
    shots_ax.set_title("Individual")
    stacked_ax.set_title("Stacked")
    shots_ax.set_ylim(-0.1, 1.1)
    stacked_ax.set_ylim(-0.1, 1.1)
    plt.tight_layout()
    fig.set_size_inches(8, 6)
    plt.savefig("figures/megashot_single.png", dpi=DPI)


def model_plots():
    fig, (model_ax, invert_ax) = plt.subplots(nrows=2, figsize=(15, 11.5), sharex=True)
    cbar_ax = fig.add_axes((0.2, -0.05, 0.6, 0.05))
    plot_vm_and_rays(
        "poster_inversion/vm_000",
        "poster_inversion/rayfile_000",
        show=False,
        ax=model_ax,
        cbar_ax=cbar_ax,
    )
    plot_vm_and_rays(
        "poster_inversion/vm_008",
        "poster_inversion/rayfile_008",
        "poster_inversion/dws_vel_009",
        show=False,
        ax=invert_ax,
        cbar_ax=cbar_ax,
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
    shots = [1000, 1050, 1287, 1350]
    fig, axs = plt.subplots(nrows=2, ncols=2, sharey="row", figsize=(13, 9))
    axs = axs.flatten()
    for i, shotno in enumerate(shots):
        min_v = 5000 if shotno < 1200 else 7000
        max_v = 7000 if shotno < 1200 else 9000
        ylim = (-2, 6)  # (0, 4) if shotno < 1200 else (-2, 2)
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

    axs[1].set_ylabel(None)
    axs[3].set_ylabel(None)

    plt.tight_layout()
    plt.savefig("figures/tt_curves.png", dpi=DPI)
