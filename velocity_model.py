import numpy as np
import pygmt
from matplotlib import pyplot as plt
from geographiclib.geodesic import Geodesic

import utils

trev = lambda x: tuple(reversed(x))


def make_model_space():
    # Let's make the model space by drawing a line from the first shot
    # of line 1 to the last shot of line 1. Note that GMT expects
    # lon,lat ordering (i.e. x,y).
    line1 = utils.shots_for_line(1)
    first_shot, last_shot = line1[0], line1[-1]
    # pointing towards Kodiak
    line_start = trev(utils.shot_lat_lon(last_shot))
    line_end = trev(utils.shot_lat_lon(first_shot))
    node_locs = np.column_stack(
        (utils._node_df.lon, utils._node_df.lat, utils._node_df.elev_m)
    )
    print(line_start, line_end)
    projected = pygmt.project(
        data=node_locs, center=line_start, endpoint=line_end, unit=True
    )
    print(
        f"distance from node to farthest shot: closest={projected[3].min():0.2f} km, farthest={projected[3].max():0.2f} km"
    )
    round_to_nearest_km = 50
    max_dist_km = (
        np.ceil(projected[3].max() / round_to_nearest_km) * round_to_nearest_km
    )
    print(f"chosen round max distance of {max_dist_km}")
    # geographiclib wants lat,lon
    start_point_dict = Geodesic.WGS84.InverseLine(
        *trev(line_start), *trev(line_end)
    ).Position(max_dist_km * 1e3)
    # Round to same precision as shot location.
    start_point = (
        np.round(start_point_dict["lat2"], decimals=6),
        np.round(start_point_dict["lon2"], decimals=6),
    )
    end_point = trev(line_start)
    print(f"the overall line will be from {start_point} to {end_point}")
    print(
        f"sanity check (should be {max_dist_km}): {Geodesic.WGS84.InverseLine(*start_point, *end_point).s13 * 1e-3}"
    )


# Want to have these written out in the source. Obtained by running
# make_model_space. The line points from Kodiak out towards the nodes.
model_start_lat_lon = (57.953815, -152.717426)
model_end_lat_lon = (55.910548, -151.01544)
