import pandas as pd
import obspy
from obspy.clients.fdsn.client import Client
import urllib
import zipfile
import os

import utils

client = Client("IRIS")


def download_earthquake(origin_time, component="Z"):
    t1 = obspy.UTCDateTime(origin_time)
    t2 = t1 + 60
    quake_ns = str(t1.ns)
    target_dir = os.path.join(f"earthquakes_{component}", quake_ns)
    if not os.path.exists(target_dir):
        print(f"downloading earthquake {quake_ns}")
        download_url = f"http://service.iris.edu/ph5ws/dataselect/1/query?reqtype=FDSN&format=segy1&net=8J&sta=*&cha=DP{component}&starttime={t1}&endtime={t2}&length=60&nodata=404"
        path, _ = urllib.request.urlretrieve(download_url)
        with zipfile.ZipFile(path, "r") as zip_ref:
            print(f"extracting files from {path} into {target_dir}")
            zip_ref.extractall(target_dir)
        urllib.request.urlcleanup()
    else:
        print(f"skipping earthquake {quake_ns}")
