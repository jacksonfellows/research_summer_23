import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from pickfile import load_picks
from rayfile import Rayfile


def merge_pick_calc(pick_df, calc_df):
    tts = pd.merge(
        pick_df,
        calc_df,
        on=("station", "shot", "phase"),
    )
    # Joining on float columns doesn't seem to work.
    assert np.allclose(tts.tt_x, tts.tt_y) and np.allclose(tts.error_x, tts.error_y)
    tts = tts.rename(columns={"tt_x": "tt", "error_x": "error"})
    tts = tts.drop(columns=["tt_y", "error_y"])
    return tts


def plot_tt_curve_megashot(shotno, shots_per_side, rayfile, red_vel=6.0, ylim=(-2, 8)):
    shot_pick = load_picks(f"picks/megashot_{shotno}_{shots_per_side}_01")
    shot_calc = Rayfile.load(rayfile).calculated_tts()
    shot_calc = shot_calc[shot_calc.shot == shotno]
    tts = merge_pick_calc(shot_pick, shot_calc)
    plt.xlim(tts.offset.max() + 2, tts.offset.min() - 2)
    plt.ylim(ylim)
    plt.xlabel("Offset (km)")
    plt.ylabel(f"Reduced time w/ v={red_vel:0.1f} km/s (s)")
    plt.errorbar(
        tts.offset,
        tts.tt - tts.offset / red_vel,
        yerr=tts.error,
        capsize=5,
        elinewidth=2,
        markeredgewidth=2,
        color="red",
        linestyle="none",
        marker="x",
    )
    plt.plot(tts.offset, tts.tt_calc - tts.offset / red_vel, "bx")
    plt.show()
