import os

import obspy
from obspy.clients.fdsn.client import Client

import utils


def download_broadband_for_line(station_code, lineno, window=(0, 60)):
    client = Client("IRIS")
    st_all = obspy.Stream()
    inv = None
    for _, shot_row in utils._shot_df[utils._shot_df.lineno == lineno].iterrows():
        shot_t = obspy.UTCDateTime(shot_row.time)
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
        st = None
        while st is None:
            try:
                print(f"trying to download shot {shot_row.shotno}")
                st = client.get_waveforms(**get_params)
            except:
                pass
        if inv is None:
            inv = client.get_stations(**get_params, level="response")
        st.remove_response(inventory=inv, output="VEL")
        print(f"trimming stream")
        st.trim(t1, t2, pad=True, fill_value=0, nearest_sample=False)
        # Save shotno with trace.
        st[0].stats.shotno = shot_row.shotno
        st_all.append(st[0])
    # Need to save as a pickle to keep shotno in stats.
    path = os.path.join(f"line_{lineno}_per_broadband", f"{station_code}.pickle")
    print(f"writing station {station_code} to {path}")
    st_all.write(path)
