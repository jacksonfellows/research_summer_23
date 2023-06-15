import pygmt
import os
import geopandas

import utils


figs_dir = "figures"


def load_usgs_faults_for_region(region):
    return geopandas.read_file(
        "./Qfaults_GIS/SHP/Qfaults_US_Database.dbf",
        bbox=(region[0], region[2], region[1], region[3]),
    )


def make_overview_map():
    fig = pygmt.Figure()
    fig.coast(
        projection="C-153/57.5/12c",
        region="US.AK+R5",
        shorelines=["1/0.5p", "2/0.5p"],
        frame="afg",
        land="grey",
    )
    fig.savefig(os.path.join(figs_dir, "overview_map.pdf"))


def make_aacse_map():
    fig = pygmt.Figure()
    region = (-165, -147, 51, 60)
    fig.coast(
        projection="M12c",
        region=region,
        shorelines=("1/0.5p", "2/0.5p"),
        frame="a2f",
        land="grey",
    )

    # bathymetry & topography
    grid = pygmt.datasets.load_earth_relief(resolution="15s", region=region)
    fig.grdimage(grid=grid, cmap="geo")

    # faults
    fig.plot(load_usgs_faults_for_region(region))

    # volcanoes
    fig.plot(
        x=utils._volcanoes_df.lon,
        y=utils._volcanoes_df.lat,
        style="t0.2c",
        fill="black",
    )

    # nodal stations
    fig.plot(
        x=utils._node_df.lon,
        y=utils._node_df.lat,
        style="c0.03c",
        fill="red",
    )

    # broadband stations
    fig.plot(
        x=utils._broadband_df.lon,
        y=utils._broadband_df.lat,
        style="s0.2c",
        fill="orange",
        pen="black",
    )

    fig.savefig(os.path.join(figs_dir, "aacse_map.pdf"))
