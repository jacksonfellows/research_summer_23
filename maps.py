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


def load_rupture_zones():
    df = geopandas.read_file("./traced_rupture_zones_2.gpkg")
    return df.sort_values(by=["rupture_name"])


rupture_zones = load_rupture_zones()

plate_boundaries = geopandas.read_file("./tectonicplates-master/PB2002_boundaries.shp")


def make_overview_map():
    fig = pygmt.Figure()
    region = "US.AK+R10"

    # basemap
    center_lon = aacse_map_region[0] + (aacse_map_region[1] - aacse_map_region[0]) / 2
    center_lat = aacse_map_region[2] + (aacse_map_region[3] - aacse_map_region[2]) / 2
    fig.basemap(
        projection=f"C{center_lon}/{center_lat}/12c",
        region=region,
        frame=["x+10f", "y+5f"],
    )

    # bathymetry & topography
    grid = pygmt.datasets.load_earth_relief(resolution="01m", region=region)
    fig.grdimage(grid=grid, cmap="geo")

    # plate boundaries
    fig.plot(
        plate_boundaries,
        pen="0.5p",
    )

    # AACSE map location
    fig.plot(
        data=[
            [aacse_map_region[0], aacse_map_region[2]],
            [aacse_map_region[0], aacse_map_region[3]],
            [aacse_map_region[1], aacse_map_region[3]],
            [aacse_map_region[1], aacse_map_region[2]],
            [aacse_map_region[0], aacse_map_region[2]],
        ],
        pen="1p,red",
        straight_line="mp",
    )

    fig.savefig(os.path.join(figs_dir, "overview_map.pdf"))


aacse_map_region = (-165, -147, 51, 60)


def make_aacse_map():
    fig = pygmt.Figure()

    # set area
    fig.basemap(projection="M12c", region=aacse_map_region, frame="a2f")

    # bathymetry & topography
    grid = pygmt.datasets.load_earth_relief(resolution="15s", region=aacse_map_region)
    fig.grdimage(grid=grid, cmap="geo")

    # faults
    fig.plot(
        load_usgs_faults_for_region(aacse_map_region), pen="0.02c", label="Mapped Fault"
    )

    # rupture zones
    zone_colors = {  # Just picked these colors randomly, should definitely change them.
        "1938": "purple",
        "1946": "red",
        "1948": "blue",
        "1957": "green",
        "1964": "pink",
        "2020": "yellow",
        "2021": "orange",
    }
    zone_magnitudes = {
        "1938": "M 8.2",
        "1946": "M 8.6",
        "1948": "M 7.5",
        "1957": "M 9.1",
        "1964": "M 9.2",
        "2020": "Mw 7.8",
        "2021": "Mw 8.2",
    }
    for i, row in rupture_zones.iterrows():
        fig.plot(
            rupture_zones.loc[[i]],  # a little silly
            pen="0.02c,4_4:4p",
            fill=f"{zone_colors[row.rupture_name]}@50",
            label=f"{row.rupture_name} {zone_magnitudes[row.rupture_name]}",
        )

    # shots
    fig.plot(
        x=utils._shot_df.lon,
        y=utils._shot_df.lat,
        style="c0.02c",
        fill="purple",
        label="Airgun Shot+S0.15c",
    )

    # volcanoes
    fig.plot(
        x=utils._volcanoes_df.lon,
        y=utils._volcanoes_df.lat,
        style="t0.2c",
        fill="black",
        label="Volcano",
    )

    # nodal stations
    fig.plot(
        x=utils._node_df.lon,
        y=utils._node_df.lat,
        style="c0.03c",
        fill="red",
        label="Nodal Station+S0.15c",
    )

    # broadband stations
    fig.plot(
        x=utils._broadband_df.lon,
        y=utils._broadband_df.lat,
        style="s0.2c",
        fill="orange",
        pen="black",
        label="Broadband Station",
    )

    # bathymetry/elevation scale
    with pygmt.config(FONT="8"):
        fig.colorbar(
            cmap="geo",
            truncate=[-8000, 3000],
            frame=["a3000f1000", "x+lBathymetry/Elevation", "y+lm"],
            position="x0.5c/8c+w3.5c+h",
        )

    # legend
    with pygmt.config(FONT="8"):
        fig.legend(position="JBR+jBR+o0.2c", box="+gwhite+p1p")

    # scale
    with pygmt.config(FONT="8"):
        fig.basemap(map_scale="JBL+jBL+o0.5c/0.75+w100")

    fig.savefig(os.path.join(figs_dir, "aacse_map.pdf"))
