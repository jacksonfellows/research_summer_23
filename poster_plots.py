from matplotlib import pyplot as plt

from megashot import megashot, megashot_all_nodes

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
