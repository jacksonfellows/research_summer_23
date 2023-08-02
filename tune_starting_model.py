import itertools
import os
import tempfile
from math import ceil, sqrt

import pandas as pd
from matplotlib import pyplot as plt

from inversion_multi import trace_picks_multi
from pickfile import load_picks
from profile_info import profile_1
from tt_curves import plot_tt_curve_megashot
from velocity_model import build_vm

shots_to_check = [1000, 1050, 1287, 1350]


def check_vm(vm, show=True):
    # Make tt curves for each shot to check.
    with tempfile.TemporaryDirectory() as tmp:
        vm_path = os.path.join(tmp, "vm")
        vm.dump(vm_path)
        shots_per_side = 5
        all_picks = pd.concat(
            [
                load_picks(f"picks/megashot_{shotno}_{shots_per_side}_01")
                for shotno in shots_to_check
            ]
        )
        rayfile_path = os.path.join(tmp, "rayfile")
        trace_picks_multi(profile_1, vm_path, all_picks, rayfile_path, drp=0.5)
        N = len(shots_to_check)
        n = ceil(sqrt(N))
        fig, axs = plt.subplots(nrows=n, ncols=n, sharey=True, figsize=(8, 8))
        axs = axs.flatten()
        for i, shotno in enumerate(shots_to_check):
            plot_tt_curve_megashot(
                shotno, shots_per_side, rayfile_path, show=False, ax=axs[i]
            )
            axs[i].set_title(f"Shot {shotno}")
    plt.tight_layout()
    if show:
        plt.show()


def try_vms():
    start_sed_vs = [2.5, 3.0, 3.5, 4.0]
    final_sed_vs = [4.5, 5.0, 5.5, 6.0]
    sed_thicknesses = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    final_crust_vs = [6.5, 7.0, 7.5]
    crust_thicknesses = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
    for (
        start_sed_v,
        final_sed_v,
        sed_thickness,
        final_crust_v,
        crust_thickness,
    ) in itertools.product(
        start_sed_vs, final_sed_vs, sed_thicknesses, final_crust_vs, crust_thicknesses
    ):
        fname = "-".join(
            str(n)
            for n in [
                start_sed_v,
                final_sed_v,
                sed_thickness,
                final_crust_v,
                crust_thickness,
            ]
        ).replace(".", "_")
        print(f"building {fname}")
        vm = build_vm(
            profile_1.start_lat_lon,
            profile_1.end_lat_lon,
            profile_1.x1,
            profile_1.x2,
            profile_1.z1,
            profile_1.z2,
            10,
            64,
            start_sed_v=start_sed_v,
            final_sed_v=final_sed_v,
            sed_thickness=sed_thickness,
            start_crust_v=final_sed_v,
            final_crust_v=final_crust_v,
            crust_thickness=crust_thickness,
        )
        print(f"checking {fname}")
        check_vm(vm, show=False)
        plt.savefig(f"starting_model_tt_plots/{fname}.png", dpi=250)
