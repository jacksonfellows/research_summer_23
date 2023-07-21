import numpy as np
import obspy

import unbinned
import utils


def make_unbinned_traces(station_code, lineno, freqmin=3, freqmax=15):
    st = utils.load_broadband(station_code, lineno)
    st.filter("bandpass", freqmin=freqmin, freqmax=freqmax)
    offsets = [
        1e-3
        * obspy.geodetics.gps2dist_azimuth(
            *utils.broadband_lat_lon(station_code), *utils.shot_lat_lon(tr.stats.shotno)
        )[0]
        for tr in st
    ]
    traces = unbinned.UnbinnedTraces(
        min(offsets), max(offsets), st[0].stats.sampling_rate
    )
    for offset, tr in zip(offsets, st):
        tr_norm = tr.data / np.max(np.abs(tr.data))
        traces.add_trace(tr_norm, offset)
    return traces
