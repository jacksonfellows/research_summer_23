import glob
import itertools
import os

import numpy as np
import obspy
import pandas as pd
import scipy


def load_node(station_code, lineno):
    stream = obspy.read(
        os.path.join(f"line_{lineno}_per_station", f"{station_code}.segy"),
        unpack_trace_headers=True,
    )
    return stream


def load_broadband(station_code, lineno):
    # Assumes broadband shots have already been trimmed to so that the
    # start time of each trace is accurate.
    stream = obspy.read(
        os.path.join(f"line_{lineno}_per_broadband", f"{station_code}.pickle")
    )
    return stream


def load_earthquake_df():
    df = pd.read_csv("grace_earthquake_events.csv", sep="\t")
    ns = [obspy.UTCDateTime(row.origin_time).ns for _, row in df.iterrows()]
    df["ns"] = ns
    return df


# locations
_node_df = pd.read_csv("./location_data/node_lat_lon_elev.csv")
_shot_df = pd.read_csv("./location_data/shots_table.csv")
_broadband_df = pd.read_csv("./location_data/broadband_lat_lon_elev.csv")
_volcanoes_df = pd.read_csv(
    "./location_data/volcanoes-2023-06-15_14-58-44_-0600_cleaned.csv"
)
_earthquake_df = load_earthquake_df()


def load_dir(segy_dir, real_starttime, len_s):
    stream = obspy.core.stream.Stream()
    for path in glob.glob(os.path.join(segy_dir, "*.sgy")):
        stream += obspy.read(path, unpack_trace_headers=True)
    # The starttime of the stream is truncated to the second. We need
    # to get the real time of the shot (with >ms precision) and cut
    # the stream so it actually starts at the shot time.
    real_endtime = real_starttime + len_s
    stream.trim(
        real_starttime, real_endtime, pad=True, fill_value=0, nearest_sample=False
    )
    return stream


def load_shot(shotno):
    shot_dir = os.path.join("shots", str(shotno))
    real_starttime = obspy.UTCDateTime(
        _shot_df[_shot_df.shotno == int(shotno)].iloc[0].time
    )
    return load_dir(shot_dir, real_starttime, 60)


def load_quake(quake_ns, component="Z"):
    quake_dir = os.path.join(f"earthquakes_{component}", str(quake_ns))
    real_starttime = obspy.UTCDateTime(
        _earthquake_df[_earthquake_df.ns == int(quake_ns)].iloc[0].origin_time
    )
    return load_dir(quake_dir, real_starttime, 60)


def node_lat_lon(stat_code):
    info = _node_df[_node_df.code == stat_code].iloc[0]
    return info.lat, info.lon


def quake_lat_lon(quake_ns):
    row = _earthquake_df[_earthquake_df.ns == quake_ns].iloc[0]
    return row.origin_latitude, row.origin_longitude


def shots_for_line(lineno):
    return list(_shot_df[_shot_df.lineno == lineno].shotno)


def shot_lat_lon(shotno):
    info = _shot_df[_shot_df.shotno == shotno].iloc[0]
    return (info.lat, info.lon)


def source_receiver_offset(tr):
    # Ignore elevation of receiver.
    stat_code = tr.stats.segy.trace_header.trace_number_within_the_original_field_record
    shotno = tr.stats.segy.trace_header.energy_source_point_number
    dist_m, _, _ = obspy.geodetics.base.gps2dist_azimuth(
        *node_lat_lon(stat_code), *shot_lat_lon(shotno)
    )
    return dist_m


def broadband_lat_lon(station_code):
    row = _broadband_df[_broadband_df.code == station_code].iloc[0]
    return (row.lat, row.lon)


def agc(tr, window_len):
    l = int(window_len)
    tr.data /= scipy.ndimage.maximum_filter1d(
        np.abs(tr.data), l, origin=-l // 2, mode="nearest"
    )


def process_stream_inplace(st):
    for tr in st:
        # Some instruments are messed up.
        if (
            tr.stats.segy.trace_header.trace_number_within_the_original_field_record
            in {3025, 4016}
        ):
            print(
                f"removing trace {tr} with station number {tr.stats.segy.trace_header.trace_number_within_the_original_field_record}"
            )
            st.remove(tr)
        else:
            tr.filter("bandpass", freqmin=3, freqmax=20, zerophase=True)
            agc(tr, 2.0 * tr.stats.sampling_rate)


def bandpass_stream_inplace(st):
    for tr in st:
        # Some instruments are messed up.
        if (
            tr.stats.segy.trace_header.trace_number_within_the_original_field_record
            in {3025, 4016}
        ):
            print(
                f"removing trace {tr} with station number {tr.stats.segy.trace_header.trace_number_within_the_original_field_record}"
            )
            st.remove(tr)
        else:
            tr.filter("bandpass", freqmin=3, freqmax=20, zerophase=True)


def remove_bad_traces_inplace(st):
    for tr in st:
        # Some instruments are messed up.
        if (
            tr.stats.segy.trace_header.trace_number_within_the_original_field_record
            in {3025, 4016}
        ):
            print(
                f"removing trace {tr} with station number {tr.stats.segy.trace_header.trace_number_within_the_original_field_record}"
            )
            st.remove(tr)


def sta_lta_inplace(st, short_window_s, long_window_s):
    sampling_rate = st[0].stats.sampling_rate
    for t in st:
        t.data = obspy.signal.filter.envelope(t.data)
        t.data = obspy.signal.trigger.classic_sta_lta(
            t.data, short_window_s * sampling_rate, long_window_s * sampling_rate
        )
        t.data /= np.abs(t.data).max()  # normalize
        if np.isnan(t.data).any():
            print(f"nan! - skipping trace {t}")
            st.remove(t)


def min_dist_to_square(lat, lon, min_lat, max_lat, min_lon, max_lon):
    assert not (
        min_lat < lat < max_lat and min_lon < lon < max_lon
    )  # assume lat, lon outside of square
    corners = list(itertools.product((min_lat, max_lat), (min_lon, max_lon)))
    sides = []
    if min_lat < lat < max_lat:
        sides += [(lat, min_lon), (lat, max_lon)]
    if min_lon < lon < max_lon:
        sides += [(min_lat, lon), (max_lat, lon)]
    return min(
        obspy.geodetics.gps2dist_azimuth(lat, lon, lat_, lon_)[0]
        for lat_, lon_ in corners + sides
    )


def max_dist_to_square(lat, lon, min_lat, max_lat, min_lon, max_lon):
    assert not (
        min_lat < lat < max_lat and min_lon < lon < max_lon
    )  # assume lat, lon outside of square
    corners = list(itertools.product((min_lat, max_lat), (min_lon, max_lon)))
    return max(
        obspy.geodetics.gps2dist_azimuth(lat, lon, lat_, lon_)[0]
        for lat_, lon_ in corners
    )


def calc_min_max_offsets_km(shot_nos):
    min_lat, max_lat, min_lon, max_lon = (
        float("inf"),
        float("-inf"),
        float("inf"),
        float("-inf"),
    )
    for _, row in _node_df.iterrows():
        lat, lon = row.lat, row.lon
        min_lat, max_lat, min_lon, max_lon = (
            min(lat, min_lat),
            max(lat, max_lat),
            min(lon, min_lon),
            max(lon, max_lon),
        )
    min_offset, max_offset = float("inf"), float("-inf")
    for shot in shot_nos:
        try:
            lat, lon = shot_lat_lon(shot)
            try:
                min_offset = min(
                    min_offset,
                    min_dist_to_square(lat, lon, min_lat, max_lat, min_lon, max_lon),
                )
                max_offset = max(
                    max_offset,
                    max_dist_to_square(lat, lon, min_lat, max_lat, min_lon, max_lon),
                )
            except ValueError:
                pass
        except KeyError:
            pass
    return min_offset / 1e3, max_offset / 1e3
