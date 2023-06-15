import pygmt
import os

import utils


figs_dir = "figures"


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
    region = [-165, -147, 51, 60]
    fig.coast(
        projection="M12c",
        region=region,
        shorelines=["1/0.5p", "2/0.5p"],
        frame="a2f",
        land="grey",
    )
    grid = pygmt.datasets.load_earth_relief(resolution="15s", region=region)
    fig.grdimage(grid=grid, cmap="geo")
    fig.savefig(os.path.join(figs_dir, "aacse_map.pdf"))
