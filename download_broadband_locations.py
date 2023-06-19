from obspy import UTCDateTime
from obspy.clients.fdsn.client import Client


def dump_inventory(inventory, filename):
    with open(filename, "w") as f:
        f.write("code,lat,lon,elev_m\n")
        for s in inventory[0]:
            f.write(f"{s.code},{s.latitude},{s.longitude},{s.elevation}\n")


def download_broadband_locations(filename):
    client = Client("IRIS")
    starttime = UTCDateTime("2019-05-01")
    endtime = UTCDateTime("2019-05-01")
    inventory = client.get_stations(network="XO", starttime=starttime, endtime=endtime)
    dump_inventory(inventory, filename)


def download_avo_locations(filename):
    client = Client("IRIS")
    inventory = client.get_stations(network="AV")
    dump_inventory(inventory, filename)
