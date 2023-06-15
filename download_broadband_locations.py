ofrom obspy import UTCDateTime
from obspy.clients.fdsn.client import Client


def download_broadband_locations(filename):
    client = Client("IRIS")
    starttime = UTCDateTime("2019-05-01")
    endtime = UTCDateTime("2019-05-01")
    inventory = client.get_stations(network="XO", starttime=starttime, endtime=endtime)
    with open(filename, 'w') as f:
        f.write("code,lat,lon,elev_m\n")
        for s in inventory[0]:
            f.write(f"{s.code},{s.latitude},{s.longitude},{s.elevation}\n")
