import numpy as np
import obspy
import obspy.signal.filter
import obspy.signal.trigger
from obspy.clients.fdsn.client import Client

import utils

client = Client("IRIS")

# Keep client.get_waveforms_bulk in mind.


def download_broadband_for_shot(station_code, shotno, window=(0, 60)):
    shot_t = obspy.UTCDateTime(
        utils._shot_df[utils._shot_df.shotno == shotno].iloc[0].time
    )
    t1 = shot_t + window[0]
    t2 = shot_t + window[1]
    get_params = {
        "network": "XO",
        "station": station_code,
        "location": "*",
        "channel": "HHZ",  # high period, high gain seismometer, vertical component
        "starttime": t1,
        "endtime": t2,
    }
    st = client.get_waveforms(**get_params)
    inv = client.get_stations(**get_params, level="response")
    st.remove_response(inventory=inv, output="VEL")
    return st[0]  # Assume 1 trace.


def bandpass(tr):
    tr_ = tr.copy()
    tr_.filter("bandpass", freqmin=3, freqmax=20, zerophase=True)
    return tr_


def sta_lta(tr):
    tr_ = tr.copy()
    tr_.data = obspy.signal.trigger.classic_sta_lta(tr_.data, 0.05 * 500, 5.0 * 500)
    tr_.data /= np.abs(tr_.data).max()  # normalize
    return tr_


def plot_with_shots(tr, shotnos):
    plt.figure()
    x = tr.times("utcdatetime")
    plt.plot(x, tr.data)
    for shotno in shotnos:
        t = obspy.UTCDateTime(
            utils._shot_df[utils._shot_df.shotno == shotno].iloc[0].time
        )
        if tr.stats.starttime < t < tr.stats.endtime:
            print(f"shot {shotno} at {t}")
            plt.axvline(t, color="red")
    plt.show()
