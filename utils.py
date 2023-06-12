import numpy as np
import scipy
import obspy
import pandas as pd
import os
import glob
import itertools


def load_shot(shotno):
    shotno = str(shotno)
    shot_dir = os.path.join("shots", shotno)
    stream = obspy.core.stream.Stream()
    for path in glob.glob(os.path.join(shot_dir, "*.sgy")):
        stream += obspy.read(path, unpack_trace_headers=True)
    return stream


def load_station_df():
    return pd.read_csv("station_lat_lon_elev.csv", index_col=1)


stat_df = load_station_df()


def stat_lat_lon(stat_code):
    info = stat_df.loc[stat_code]
    return info.lat, info.lon


def load_shot_df():
    return pd.read_excel(
        "PASSCAL_SHOT.xlsx",
        names=["lineno", "shotno", "lat", "lon", "unknown", "time", "linename"],
        index_col=1,
        header=None,
        usecols="A:G",
    )


shot_df = load_shot_df()


def shots_for_line(lineno):
    return list(shot_df[shot_df.lineno == lineno].index)


def shot_lat_lon(shotno):
    info = shot_df.loc[shotno]
    if len(info) != 6:
        info = info.iloc[0]
    return (
        info.lat,
        -info.lon,
    )  # There seems to be a problem with the shot locations. For now just fix it here.


def source_receiver_offset(tr):
    # Ignore elevation of receiver.
    stat_code = tr.stats.segy.trace_header.trace_number_within_the_original_field_record
    shotno = tr.stats.segy.trace_header.energy_source_point_number
    dist_m, _, _ = obspy.geodetics.base.gps2dist_azimuth(
        *stat_lat_lon(stat_code), *shot_lat_lon(shotno)
    )
    return dist_m


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
    for stat in stat_df.index:
        lat, lon = stat_lat_lon(stat)
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
