from matplotlib import pyplot as plt
import tempfile
import os
import datetime

import wrappers
import pick
import profile_info
import velocity_model
import rayfile


def plot_vm_and_rays(vm_path, rayfile_path):
    vm = velocity_model.VMTOMO_VM.load(vm_path)
    rays = rayfile.load(rayfile_path)
    vm.plot(show=False)
    for x, z in rays:
        plt.plot(x, z)
    plt.show()


def trace_picks(profile, vm_path, picks, rayfile_path, drp=0.1):
    with tempfile.TemporaryDirectory() as tmp:
        pick_path = os.path.join(tmp, "p")
        shot_path = os.path.join(tmp, "s")

        all_shots = profile.project_shots(list(picks.shot.unique()))

        new = 1
        for station in list(picks.station.unique()):
            picks_for_station = picks[picks.station == station]
            picks_for_station.to_csv(
                pick_path,
                columns=["station", "shot", "phase", "offset", "tt", "error"],
                sep=" ",
                header=False,
                index=False,
            )

            shots_for_pick = all_shots[
                all_shots.shotno.isin(picks_for_station.shot.unique())
            ]
            shots_for_pick.to_csv(
                shot_path,
                sep=" ",
                columns=["shotno", "x", "y", "z"],
                header=False,
                index=False,
            )

            assert len(str(station)) == 7 and str(station)[0] == "1"
            station_x = pick.fake_station_offset_km(station)
            station_z = 0

            # Use seconds in current day - seconds since epoch
            # overflows 32-bit (?) int in vm_trace :(.
            t = datetime.datetime.now()
            isec = int(
                datetime.timedelta(
                    hours=t.hour, minutes=t.minute, seconds=t.second
                ).total_seconds()
            )

            print(
                f"tracing {len(picks_for_station)} picks for station {station} ({isec=})"
            )

            wrappers.vm_trace(
                vm_file=vm_path,
                vtop=0.8,  # default in raytrace_2D.csh
                drp=drp,
                isec=isec,
                itop=picks_for_station.phase.min(),
                ibot=picks_for_station.phase.max(),
                ins=station,
                xin=station_x,
                zin=station_z,
                pfile=pick_path,
                npx=len(picks_for_station),
                sfile=shot_path,
                nsh=len(shots_for_pick),
                new=new,
                rayfile=rayfile_path,
            )
            new = 0


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
):
    with tempfile.TemporaryDirectory() as tmp:
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
            dz1=0.6,
            vscal=6.5,
            zscal=15.0,
            xreg=15.0,
            zreg=8.0,
            asr=2.5,  # I don't get this since it is not the ratio of xreg and zreg
            reg0=1.51,
            reg1=1.81,
            reg2=1.51,
            crf=0.6,
            ch2n=ch2n,
            vpmin=2.0,
            dwsv=dws_vel_path,
            dwsz=dws_zrf_path,
            mask_file=mask_path,
            vmnew=new_vm_path,
        )


def run_inversion(
    profile, initial_vm, picks, inversion_dir, n_iters, drp=0.1, **invert_kargs
):
    # Keeps all the generated rayfiles and velocity models in the inversion_dir.
    n_iter = 0

    # Get started with initial vm.
    vm_path = os.path.join(inversion_dir, f"vm_{n_iter:03}")
    if not os.path.exists(vm_path):
        print(f"dumping initial velocity model")
        initial_vm.dump(vm_path)
    else:
        print(f"initial velocity model already exists - skipping dumping")

    while n_iter < n_iters:
        vm_path = os.path.join(inversion_dir, f"vm_{n_iter:03}")
        rayfile_path = os.path.join(inversion_dir, f"rayfile_{n_iter:03}")
        if not os.path.exists(rayfile_path):
            print(f"raytracing for iteration {n_iter:03}")
            trace_picks(profile, vm_path, picks, rayfile_path, drp=drp)
        else:
            print(f"rayfile {rayfile_path} already exists - skipping raytracing")
        new_vm_path = os.path.join(inversion_dir, f"vm_{n_iter+1:03}")
        if not os.path.exists(new_vm_path):
            print(f"inverting for iteration {n_iter:03}")
            invert(
                vm_path,
                rayfile_path,
                new_vm_path,
                os.path.join(inversion_dir, f"dws_vel_{n_iter+1:03}"),
                os.path.join(inversion_dir, f"dws_zrf_{n_iter+1:03}"),
                os.path.join(inversion_dir, f"mask_{n_iter+1:03}"),
                itop=picks.phase.min(),
                ibot=picks.phase.max(),
                **invert_kargs,
            )
            if (
                not os.path.exists(os.path.join(inversion_dir, f"vm_{n_iter+1:03}"))
                or os.path.getsize(os.path.join(inversion_dir, f"vm_{n_iter+1:03}"))
                == 0
            ):
                raise RuntimeError(
                    "Inversion failed! (did not produce a new velocity model)"
                )
        else:
            print(f"velocity model {new_vm_path} already exists - skipping inversion")
        n_iter += 1