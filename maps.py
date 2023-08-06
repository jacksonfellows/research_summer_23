import os

import geopandas
import pandas as pd
import pygmt

import profile_info
import utils

figs_dir = "figures"


def load_usgs_faults_for_region(region):
    return geopandas.read_file(
        "./map_data/Qfaults_GIS/SHP/Qfaults_US_Database.dbf",
        bbox=(region[0], region[2], region[1], region[3]),
    )


def load_rupture_zones():
    df = geopandas.read_file("./map_data/traced_rupture_zones_2.gpkg")
    return df.sort_values(by=["rupture_name"])


rupture_zones = load_rupture_zones()
plate_boundaries = geopandas.read_file(
    "./map_data/tectonicplates-master/PB2002_boundaries.shp"
)
aleut_lines = pd.read_csv("./map_data/aleut_lines.csv")
edge_lines = geopandas.read_file("./map_data/traced_edge_transect.gpkg")
avo_stations = pd.read_csv("./map_data/avo_lat_lon_elev.csv")

model_region = (-155, -150.75, 54.5, 58.5)
aacse_map_region = (-165, -147, 51, 60)


pen_width = "0.015i"
marker_size = "0.15i"


def make_aacse_map():
    fig = pygmt.Figure()

    # set area
    fig.basemap(projection="M8i", region=aacse_map_region, frame="a2f")

    # bathymetry & topography
    grid = pygmt.datasets.load_earth_relief(resolution="15s", region=aacse_map_region)
    pygmt.makecpt(cmap="geo", series=(-8000, 8000))
    fig.grdimage(grid=grid, cmap=True)

    # faults
    fig.plot(
        load_usgs_faults_for_region(aacse_map_region),
        pen=pen_width,
        label="Mapped Fault",
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
            pen=f"{pen_width},4_4:4p",
            fill=f"{zone_colors[row.rupture_name]}@50",
            label=f"{row.rupture_name} {zone_magnitudes[row.rupture_name]}",
        )

    # EDGE lines
    fig.plot(edge_lines, pen=f"{pen_width},gold", label="EDGE Line")

    # ALEUT lines
    for i, row in aleut_lines.iterrows():
        fig.plot(
            x=(row.start_longitude, row.stop_longitude),
            y=(row.start_latitude, row.stop_latitude),
            pen=f"{pen_width},green",
            label="ALEUT Line" if i == 0 else None,
        )

    # shots
    fig.plot(
        x=utils._shot_df.lon,
        y=utils._shot_df.lat,
        style=f"c{pen_width}",
        fill="red",
        label="Airgun Shot+S0.15c",
    )

    # volcanoes
    fig.plot(
        x=utils._volcanoes_df.lon,
        y=utils._volcanoes_df.lat,
        style=f"t{marker_size}",
        fill="black",
        label="Volcano",
    )

    # AVO stations
    fig.plot(
        x=avo_stations.lon,
        y=avo_stations.lat,
        style=f"i{marker_size}",
        pen="black",
        fill="turquoise",
        label="AVO Station",
    )

    # nodal stations
    fig.plot(
        x=utils._node_df.lon,
        y=utils._node_df.lat,
        style=f"c{pen_width}",
        fill="purple",
        label="Nodal Station+S0.15c",
    )

    # broadband stations
    fig.plot(
        x=utils._broadband_df.lon,
        y=utils._broadband_df.lat,
        style=f"s{marker_size}",
        fill="orange",
        pen="black",
        label="Broadband Station",
    )

    # model map
    fig.plot(
        data=[
            [model_region[0], model_region[2]],
            [model_region[0], model_region[3]],
            [model_region[1], model_region[3]],
            [model_region[1], model_region[2]],
            [model_region[0], model_region[2]],
        ],
        pen=f"{pen_width},black",
        straight_line="mp",
    )

    # overview inset
    center_lon = aacse_map_region[0] + (aacse_map_region[1] - aacse_map_region[0]) / 2
    center_lat = aacse_map_region[2] + (aacse_map_region[3] - aacse_map_region[2]) / 2
    inset_region = "US.AK+R10"
    with fig.inset(
        position="jTL+o0.2c",
        box="+gwhite+p1p",
        region=inset_region,
        projection=f"C{center_lon}/{center_lat}/2i",
    ):
        # bathymetry & topography
        inset_grid = pygmt.datasets.load_earth_relief(
            resolution="10m", region=inset_region
        )
        fig.grdimage(grid=inset_grid, cmap=True)

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
            pen=f"{pen_width},black",
            straight_line="mp",
        )

    # bathymetry/elevation scale
    with pygmt.config(FONT="12"):
        fig.colorbar(
            cmap="geo",
            truncate=[-8000, 3000],
            frame=["a3000f1000", "x+lBathymetry/Elevation", "y+lm"],
            position="x0.25i/5i+w2i+h",
        )

    # legend
    with pygmt.config(FONT="12"):
        fig.legend(position="JBR+jBR+o0.2c", box="+gwhite+p1p")

    # scale
    with pygmt.config(FONT="12"):
        fig.basemap(map_scale="JBL+jBL+o0.5c/0.75+w100")

    fig.savefig(os.path.join(figs_dir, "aacse_map.pdf"))


def make_model_map():
    fig = pygmt.Figure()

    region = model_region

    # set area
    fig.basemap(projection="M6i", region=region, frame="a1f")

    # bathymetry & topography
    grid = pygmt.datasets.load_earth_relief(resolution="15s", region=region)
    pygmt.makecpt(cmap="geo", series=(-8000, 8000))
    fig.grdimage(grid=grid, cmap=True)

    # faults
    fig.plot(
        load_usgs_faults_for_region(aacse_map_region), pen="0.02c", label="Mapped Fault"
    )

    # nodal stations
    fig.plot(
        x=utils._node_df.lon,
        y=utils._node_df.lat,
        style=f"c{pen_width}",
        fill="purple",
        label="Nodal Station+S0.15c",
    )

    # shots
    fig.plot(
        x=utils._shot_df.lon,
        y=utils._shot_df.lat,
        style=f"c{pen_width}",
        fill="red",
        label="Airgun Shot+S0.15c",
    )

    # broadband stations
    fig.plot(
        x=utils._broadband_df.lon,
        y=utils._broadband_df.lat,
        style=f"s{marker_size}",
        fill="orange",
        pen="black",
        label="Broadband Station",
    )

    fig.text(
        x=utils._broadband_df.lon,
        y=utils._broadband_df.lat,
        text=utils._broadband_df.code,
        font="12p",
    )

    # profile line
    fig.plot(
        x=(
            profile_info.profile_1.start_lat_lon[1],
            profile_info.profile_1.end_lat_lon[1],
        ),
        y=(
            profile_info.profile_1.start_lat_lon[0],
            profile_info.profile_1.end_lat_lon[0],
        ),
        pen=f"{pen_width},gold,4_4:4p",
        label="Kodiak Profile",
    )

    # legend
    with pygmt.config(FONT="12"):
        fig.legend(position="JBL+jBL+o0.2c", box="+gwhite+p1p")

    # scale
    with pygmt.config(FONT="12"):
        fig.basemap(map_scale="JBR+jBR+o0.5c/0.75+w50")

    fig.savefig(os.path.join(figs_dir, "model_map.pdf"))
