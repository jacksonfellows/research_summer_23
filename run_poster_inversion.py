# Good to document what commands I run.

from inversion import plot_vm_and_rays, run_inversion
from pickfile import load_pick_list
from profile_info import profile_1
from velocity_model import build_vm

vm_test = build_vm(
    profile_1.start_lat_lon,
    profile_1.end_lat_lon,
    profile_1.x1,
    profile_1.x2,
    profile_1.z1,
    profile_1.z2,
    10,
    64,
    start_sed_v=3.0,
    final_sed_v=5.5,
    sed_thickness=1,
    start_crust_v=5.5,
    final_crust_v=6.5,
    crust_thickness=15,
)


def run():
    run_inversion(
        profile_1,
        vm_test,
        load_pick_list("good_inversion_picks"),
        "poster_inversion",
        10,
        drp=0.5,
        ch2n=4,
        ch2n_factor=1 / 1.3,
        reg0=1,
        reg1=1,
        reg2=1,
        xreg=2,
        zreg=4,
        nxc=vm_test.nx // 4,
        nzc=vm_test.nz // 4,
    )
