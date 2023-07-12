import tempfile, multiprocessing

import profile_info
import pick
import rayfile

# Need to have the function to map over defined in a separate file.
from inversion_multi_trace_station import trace_station


def trace_picks_multi(profile, vm_path, picks, rayfile_path, drp=0.1):
    with tempfile.TemporaryDirectory() as tmp:
        all_shots = profile.project_shots(list(picks.shot.unique()))

        stations = picks.station.unique()

        with multiprocessing.Pool() as pool:
            rf_paths = pool.map(
                trace_station,
                [
                    {
                        "station": station,
                        "tmp": tmp,
                        "picks": picks,
                        "all_shots": all_shots,
                        "vm_path": vm_path,
                        "drp": drp,
                    }
                    for station in stations
                ],
            )
            print(f"combining {len(rf_paths)} rayfiles into 1")
            rayfile.cat_rayfiles(rf_paths, rayfile_path)
