import numpy as np
import obspy
import obspy.signal.filter
import obspy.signal.trigger

import binned
import utils


def quake_node_offset_km(quake_ns, node_code):
    row = utils._earthquake_df[utils._earthquake_df.ns == quake_ns].iloc[0]
    lat, lon, depth_km = row.origin_latitude, row.origin_longitude, row.origin_depth_km
    surface_offset_km = (
        1e-3
        * obspy.geodetics.gps2dist_azimuth(
            *utils.quake_lat_lon(quake_ns), *utils.node_lat_lon(node_code)
        )[0]
    )
    return np.sqrt(surface_offset_km**2 + depth_km**2)


def bin_quake(quake_ns, component="Z"):
    st = utils.load_quake(quake_ns, component)
    sampling_rate = st[0].stats.sampling_rate
    utils.remove_bad_traces_inplace(st)
    st.filter("bandpass", freqmin=3.0, freqmax=20.0, zerophase=True)
    utils.sta_lta_inplace(st, 0.01, 1.0)
    offsets = [
        quake_node_offset_km(
            quake_ns,
            t.stats.segy.trace_header.trace_number_within_the_original_field_record,
        )
        for t in st
    ]
    bt = binned.BinnedTraces(
        min(offsets),
        max(offsets),
        0.1,
        60,
        int(sampling_rate),
    )
    for t, offset in zip(st.traces, offsets):
        # silly
        if t.data.shape[0] == 30001:
            bt.add_trace(t.data[:-1], offset)
        else:
            bt.add_trace(t.data, offset)
    return bt


def bin_quake_bandpass(quake_ns, component="Z"):
    st = utils.load_quake(quake_ns, component)
    sampling_rate = st[0].stats.sampling_rate
    utils.remove_bad_traces_inplace(st)
    st.filter("bandpass", freqmin=0.5, freqmax=2.0, zerophase=True)
    # utils.sta_lta_inplace(st, 0.01, 1.0)
    offsets = [
        quake_node_offset_km(
            quake_ns,
            t.stats.segy.trace_header.trace_number_within_the_original_field_record,
        )
        for t in st
    ]
    bt = binned.BinnedTraces(
        min(offsets),
        max(offsets),
        0.1,
        60,
        int(sampling_rate),
    )
    for t, offset in zip(st.traces, offsets):
        # silly
        if t.data.shape[0] == 30001:
            bt.add_trace(t.data[:-1], offset)
        else:
            bt.add_trace(t.data, offset)
    return bt
