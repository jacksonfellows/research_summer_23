import pygmt
import os
import geopandas
import pandas as pd

import utils
import profile_info


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
    pygmt.makecpt(cmap="geo", series=(-8000, 8000))
    fig.grdimage(grid=grid, cmap=True)

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
    pygmt.makecpt(cmap="geo", series=(-8000, 8000))
    fig.grdimage(grid=grid, cmap=True)

    # slab contours
    fig.plot(
        "./map_data/alu_slab2_dep_02.23.18_contours.in",
        pen="0.5p,pink",
        label="Slab 20 km",
    )

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

    # EDGE lines
    fig.plot(edge_lines, pen="0.02c,gold", label="EDGE Line")

    # ALEUT lines
    for i, row in aleut_lines.iterrows():
        fig.plot(
            x=(row.start_longitude, row.stop_longitude),
            y=(row.start_latitude, row.stop_latitude),
            pen="0.02c,green",
            label="ALEUT Line" if i == 0 else None,
        )

    # shots
    fig.plot(
        x=utils._shot_df.lon,
        y=utils._shot_df.lat,
        style="c0.02c",
        fill="red",
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

    # AVO stations
    fig.plot(
        x=avo_stations.lon,
        y=avo_stations.lat,
        style="i0.2c",
        pen="black",
        fill="turquoise",
        label="AVO Station",
    )

    # nodal stations
    fig.plot(
        x=utils._node_df.lon,
        y=utils._node_df.lat,
        style="c0.03c",
        fill="purple",
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
    with pygmt.config(FONT="6"):
        fig.legend(position="JBR+jBR+o0.2c", box="+gwhite+p1p")

    # scale
    with pygmt.config(FONT="8"):
        fig.basemap(map_scale="JBL+jBL+o0.5c/0.75+w100")

    fig.savefig(os.path.join(figs_dir, "aacse_map.pdf"))


def make_aacse_mag_map():
    fig = pygmt.Figure()

    # set area
    fig.basemap(projection="M12c", region=aacse_map_region, frame="a2f")

    # magnetic
    grid = pygmt.datasets.load_earth_magnetic_anomaly(
        resolution="03m",
        region=aacse_map_region,
        data_source="wdmam",
    )
    pygmt.makecpt(cmap="mag", series=(-500, 500))
    fig.grdimage(grid=grid, cmap=True)

    # shorelines
    fig.coast(region=aacse_map_region, shorelines="1/0.5p")

    # faults
    fig.plot(
        load_usgs_faults_for_region(aacse_map_region), pen="0.02c", label="Mapped Fault"
    )

    # rupture zones
    fig.plot(
        rupture_zones,
        pen="0.02c,4_4:4p",
        label="Historical Rupture Zone",
    )

    # magnetic scale
    with pygmt.config(FONT="8"):
        fig.colorbar(
            cmap=True,
            frame=["a250f125", "x+lMagnetic Anomaly", "y+lnT"],
            position="x0.5c/8c+w3.5c+h",
            box="+gwhite+p1p",
        )

    # scale
    with pygmt.config(FONT="8"):
        fig.basemap(map_scale="JBL+jBL+o0.5c/0.75+w100")

    fig.savefig(os.path.join(figs_dir, "aacse_mag_map.pdf"))


def make_aacse_grav_map():
    fig = pygmt.Figure()

    # set area
    fig.basemap(projection="M12c", region=aacse_map_region, frame="a2f")

    # magnetic
    grid = pygmt.datasets.load_earth_free_air_anomaly(
        resolution="01m",
        region=aacse_map_region,
    )
    pygmt.makecpt(cmap="polar", series=(-275, 275))
    fig.grdimage(grid=grid, cmap=True)

    # shorelines
    fig.coast(region=aacse_map_region, shorelines="1/0.5p")

    # faults
    fig.plot(
        load_usgs_faults_for_region(aacse_map_region), pen="0.02c", label="Mapped Fault"
    )

    # rupture zones
    fig.plot(
        rupture_zones,
        pen="0.02c,4_4:4p",
        label="Historical Rupture Zone",
    )

    # gravity scale
    with pygmt.config(FONT="8"):
        fig.colorbar(
            cmap=True,
            frame=["a100f50", "x+lFree-Air Anomaly", "y+lmGal"],
            position="x0.5c/8c+w3.5c+h",
        )

    # scale
    with pygmt.config(FONT="8"):
        fig.basemap(map_scale="JBL+jBL+o0.5c/0.75+w100")

    fig.savefig(os.path.join(figs_dir, "aacse_grav_map.pdf"))


def make_model_map():
    fig = pygmt.Figure()
    region = (-155, -150.75, 54.5, 58.5)

    # set area
    fig.basemap(projection="M12c", region=region, frame="a1f")

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
        style="c0.06c",
        fill="purple",
        label="Nodal Station+S0.15c",
    )

    # shots
    fig.plot(
        x=utils._shot_df.lon,
        y=utils._shot_df.lat,
        style="c0.06c",
        fill="red",
        label="Airgun Shot+S0.15c",
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
    fig.text(
        x=utils._broadband_df.lon,
        y=utils._broadband_df.lat,
        text=utils._broadband_df.code,
    )

    # profile lines
    for profile, color in (
        (profile_info.profile_1, "saddlebrown"),
        (profile_info.profile_2, "limegreen"),
    ):
        fig.plot(
            x=(profile["start_lat_lon"][1], profile["end_lat_lon"][1]),
            y=(profile["start_lat_lon"][0], profile["end_lat_lon"][0]),
            pen=f"0.06c,{color},4_4:4p",
            label=profile["name"],
        )

    # legend
    with pygmt.config(FONT="10"):
        fig.legend(position="JBL+jBL+o0.2c", box="+gwhite+p1p")

    fig.savefig(os.path.join(figs_dir, "model_map.pdf"))


def make_quake_map():
    fig = pygmt.Figure()
    region = (-155.5, -150.5, 56, 59)

    # set area
    fig.basemap(projection="M12c", region=region, frame="a1f")

    # bathymetry & topography
    grid = pygmt.datasets.load_earth_relief(resolution="15s", region=region)
    pygmt.makecpt(cmap="geo", series=(-8000, 8000))
    fig.grdimage(grid=grid, cmap=True)

    # nodal stations
    fig.plot(
        x=utils._node_df.lon,
        y=utils._node_df.lat,
        style="c0.06c",
        fill="purple",
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
    fig.text(
        x=utils._broadband_df.lon,
        y=utils._broadband_df.lat,
        text=utils._broadband_df.code,
    )

    # earthquakes
    for i, row in utils._earthquake_df.iterrows():
        fig.plot(
            x=row.origin_longitude,
            y=row.origin_latitude,
            style="x0.4c",
            pen="0.05c,red",
        )
        fig.text(x=row.origin_longitude, y=row.origin_latitude, text=f"{i}")
        if i >= 10:
            break

    # legend
    with pygmt.config(FONT="10"):
        fig.legend(position="JBL+jBL+o0.2c", box="+gwhite+p1p")

    fig.savefig(os.path.join(figs_dir, "quake_map.pdf"))
