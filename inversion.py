import datetime
import os
import tempfile

import numpy as np
from matplotlib import pyplot as plt

import pick
import pickfile
import profile_info
import rayfile
import velocity_model
import wrappers
from dws import load_dws_mask
from inversion_multi import trace_picks_multi


def plot_vm_and_rays(
    vm_path, rayfile_path, dws_path=None, show=True, ax=None, cbar_ax=None
):
    vm = velocity_model.VMTOMO_VM.load(vm_path)
    rays = rayfile.Rayfile.load(rayfile_path).rays()
    if ax is None:
        ax = plt.subplot()
    vm.plot(show=False, ax=ax, cbar_ax=cbar_ax)
    ax.set_ylim(40, -2)
    for x, z in rays:
        ax.plot(x, z, color="k", linewidth=0.1)
    if dws_path is not None:
        x, z, d = load_dws_mask(dws_path)
        d_i = np.zeros((*d.shape, 4))
        d_i[d == 0] = np.array([0, 0, 0, 0])
        d_i[d == 1] = np.array([1, 1, 1, 1])
        ax.pcolormesh(x, z, d_i, zorder=100)
        # plot boundaries again
        bx = np.linspace(vm.x1, vm.x2, vm.nx)
        for boundary_i in range(vm.nr):
            b = vm.zrf[boundary_i * vm.nx : (boundary_i + 1) * vm.nx]
            ax.plot(bx, b, color="k", zorder=101)
    if show:
        plt.show()


def invert(
    vm_path,
    rayfile_path,
    new_vm_path,
    dws_vel_path,
    dws_zrf_path,
    mask_path,
    itop,
    ibot,
    nxc=125,
    nzc=80,
    ch2n=12,
    xreg=15.0,
    zreg=8.0,
    crf=1.0,
    reg0=1,
    reg1=1,
    reg2=1,
    dz1=0.1,
):
    wrappers.vm_tomo(
        vm_file=vm_path,
        rayfile=rayfile_path,
        # would be nice if these updated automatically
        itop=itop,
        ibot=ibot,
        # taken from tomo.csh
        cmax=0.08,
        nxc=nxc,
        nzc=nzc,
        dz1=dz1,
        vscal=6.5,
        zscal=15.0,
        xreg=xreg,
        zreg=zreg,
        asr=xreg / zreg,
        reg0=reg0,
        reg1=reg1,
        reg2=reg2,
        crf=crf,
        ch2n=ch2n,
        vpmin=2.0,
        dwsv=dws_vel_path,
        dwsz=dws_zrf_path,
        mask_file=mask_path,
        vmnew=new_vm_path,
    )


def run_inversion(
    profile,
    initial_vm,
    picks,
    inversion_dir,
    n_iters,
    drp=0.1,
    ch2n=16,
    ch2n_factor=0.5,
    **invert_kwargs,
):
    # Keeps all the generated rayfiles and velocity models in the inversion_dir.
    n_iter = 0

    # Get started with initial vm.
    vm_path = os.path.join(inversion_dir, f"vm_{n_iter:03}")
    if not os.path.exists(vm_path):
        print("dumping initial velocity model")
        initial_vm.dump(vm_path)
    else:
        print("initial velocity model already exists - skipping dumping")

    while n_iter < n_iters:
        vm_path = os.path.join(inversion_dir, f"vm_{n_iter:03}")
        rayfile_path = os.path.join(inversion_dir, f"rayfile_{n_iter:03}")
        if not os.path.exists(rayfile_path):
            print(f"raytracing for iteration {n_iter:03}")
            trace_picks_multi(profile, vm_path, picks, rayfile_path, drp=drp)
        else:
            print(f"rayfile {rayfile_path} already exists - skipping raytracing")
        new_vm_path = os.path.join(inversion_dir, f"vm_{n_iter+1:03}")
        if not os.path.exists(new_vm_path):
            print(f"inverting for iteration {n_iter:03} w/ {ch2n=}")
            invert(
                vm_path,
                rayfile_path,
                new_vm_path,
                os.path.join(inversion_dir, f"dws_vel_{n_iter+1:03}"),
                os.path.join(inversion_dir, f"dws_zrf_{n_iter+1:03}"),
                os.path.join(inversion_dir, f"mask_{n_iter+1:03}"),
                itop=picks.phase.min(),
                ibot=picks.phase.max(),
                ch2n=ch2n,
                **invert_kwargs,
            )
            if (
                not os.path.exists(os.path.join(inversion_dir, f"vm_{n_iter+1:03}"))
                or os.path.getsize(os.path.join(inversion_dir, f"vm_{n_iter+1:03}"))
                == 0
            ):
                raise RuntimeError(
                    "Inversion failed! (did not produce a new velocity model)"
                )
            ch2n = max(1, ch2n_factor * ch2n)
            print(f"updated {ch2n=}")
        else:
            print(f"velocity model {new_vm_path} already exists - skipping inversion")
        n_iter += 1
