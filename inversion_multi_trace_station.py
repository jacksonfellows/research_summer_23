import datetime
import os

import pick
import rayfile
import wrappers


def trace_station(trace_args):
    # Silly
    station = trace_args["station"]
    tmp = trace_args["tmp"]
    picks = trace_args["picks"]
    all_shots = trace_args["all_shots"]
    vm_path = trace_args["vm_path"]
    drp = trace_args["drp"]

    # Need unique paths for each station
    pick_path = os.path.join(tmp, f"p_{station}")
    shot_path = os.path.join(tmp, f"s_{station}")
    rayfile_path = os.path.join(tmp, f"rayfile_{station}")

    picks_for_station = picks[picks.station == station]
    picks_for_station.to_csv(
        pick_path,
        columns=["station", "shot", "phase", "offset", "tt", "error"],
        sep=" ",
        header=False,
        index=False,
    )

    shots_for_pick = all_shots[all_shots.shotno.isin(picks_for_station.shot.unique())]
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

    print(f"tracing {len(picks_for_station)} picks for station {station} ({isec=})")

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
        new=1,
        rayfile=rayfile_path,
    )

    return rayfile_path
