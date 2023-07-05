import numpy as np
import obspy
import obspy.signal.filter
import obspy.signal.trigger

import utils
import binned


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


def bin_quake(quake_ns):
    st = utils.load_quake(quake_ns)
    sampling_rate = st[0].stats.sampling_rate
    utils.remove_bad_traces_inplace(st)
    st.filter("bandpass", freqmin=3.0, freqmax=20.0, zerophase=True)
    for t in st:
        t.data = obspy.signal.filter.envelope(t.data)
        t.data = obspy.signal.trigger.classic_sta_lta(
            t.data, 0.01 * sampling_rate, 1.0 * sampling_rate
        )
        t.data /= np.abs(t.data).max()  # normalize
        if np.isnan(t.data).any():
            print(f"nan! - skipping trace {t}")
            st.remove(t)
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
        bt.add_trace(t.data[:-1], offset)
    return bt
