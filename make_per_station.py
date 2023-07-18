import collections
import os

import obspy

import utils


def gather_per_station(shotnos):
    stat_sts = collections.defaultdict(lambda: obspy.Stream())
    for shotno in shotnos:
        print(f"loading shot {shotno}")
        shot_st = utils.load_shot(shotno)
        for t in shot_st:
            stat = (
                t.stats.segy.trace_header.trace_number_within_the_original_field_record
            )
            stat_sts[stat] += t
    return stat_sts


def make_line_per_station(lineno):
    dir = f"line_{lineno}_per_station"
    per_stat_dict = gather_per_station(utils.shots_for_line(lineno))
    for stat, st in per_stat_dict.items():
        path = os.path.join(dir, f"{stat}.segy")
        print(f"writing station {stat} to {path}")
        st.write(path)
