import numpy as np
import pygmt
from matplotlib import pyplot as plt
from geographiclib.geodesic import Geodesic
import pygmt
import scipy
import functools

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


class VelocityModel:
    def __init__(
        self, start_lat_lon, end_lat_lon, length_km, min_depth_km, max_depth_km
    ):
        """
        Create a velocity model along a profile from start_lat_lon to
        end_lat_lon. The profile should be length_km km long. It will
        store depths from min_depth_km to max_depth_km.
        """
        self.start_lat_lon = start_lat_lon
        self.end_lat_lon = end_lat_lon
        self.length_km = length_km
        self.min_depth_km = min_depth_km
        self.max_depth_km = max_depth_km

        # Store a mapping of layer number -> (xx, yy, vv). xx and yy
        # should span the entire model grid. Layers numbers start at 1
        # and increase with depth.
        self.layer_velocities = {}

        # Store a mapping of boundary number -> (x, z). x should span
        # the length of the model. Boundary numbers start at 1 and
        # increase with depth. E.g. the boundary between layers 1 and
        # 2 is boundary number 1.
        self.boundary_depths = {}

    def plot(self):
        """
        Plots the velocity model.
        """
        plt.xlim(0, self.length_km)
        plt.ylim(self.max_depth_km, self.min_depth_km)
        plt.set_cmap("gist_rainbow")

        first_layer = min(self.layer_velocities.keys())
        for layer_num in sorted(self.layer_velocities.keys()):
            xx, zz, vv = self.layer_velocities[layer_num]
            x1, z1 = (
                (xx[:, 0], np.full(xx.shape[0], self.min_depth_km))
                if layer_num == first_layer
                else self.boundary_depths[layer_num - 1]
            )
            z1_interp = scipy.interpolate.interp1d(x1, z1)(xx[:, 0])
            vv_masked = np.ma.masked_array(vv, zz < z1_interp[:, np.newaxis])
            plt.pcolormesh(xx, zz, vv_masked, vmin=0, vmax=10)
        plt.xlabel("x (km)")
        plt.ylabel("z (km)")
        plt.colorbar(label="v (km/s)")
        plt.show()

    def plot_layer(self, layer_num):
        """
        Plot a single layer.
        """
        plt.xlim(0, self.length_km)
        plt.ylim(self.max_depth_km, self.min_depth_km)
        plt.set_cmap("gist_rainbow")

        xx, zz, vv = self.layer_velocities[layer_num]
        plt.pcolormesh(xx, zz, vv, vmin=0, vmax=10)
        # plt.plot(xx, zz, marker="x", color="k", linestyle="none")

        plt.xlabel("x (km)")
        plt.ylabel("z (km)")
        plt.colorbar(label="v (km/s)")
        plt.show()

    def add_layer_f(self, layer_num, v_f, n_x_samples, n_z_samples):
        """
        Adds a layer with velocity v_f(x, z) to the velocity model.
        Takes evenly spaced n_x_samples along the length and
        n_z_samples along the depth. Positions are assumed to be in km
        and velocity is assumed to be in km/s.
        """
        x = np.linspace(0, self.length_km, n_x_samples)
        z = np.linspace(self.min_depth_km, self.max_depth_km, n_z_samples)
        xx, zz = np.meshgrid(x, z, indexing="ij")
        self.layer_velocities[layer_num] = (xx, zz, v_f(xx, zz))

    def add_boundary_f(self, boundary_num, z_f, n_x_samples):
        """
        Adds a boundary with depth z_f(x) to the velocity model. Takes
        evenly spaced n_x_samples along the length. Positions are
        assumed to be in km.
        """
        x = np.linspace(0, self.length_km, n_x_samples)
        z = z_f(x)
        assert x.shape == z.shape == (n_x_samples,)
        self.boundary_depths[boundary_num] = (x, z)

    def add_boundary(self, boundary_num, x, z):
        assert x.shape == z.shape
        self.boundary_depths[boundary_num] = (x, z)


def create_boundary_from_grid(grid, start_lat_lon, end_lat_lon, length_km, n_x_samples):
    """
    Create a boundary from a grid.
    """
    sample_interval_km = length_km / (n_x_samples - 1)
    df = pygmt.grdtrack(
        # +i sets the sampling interval
        # +d tells it to output distances along profile
        profile=f"{start_lat_lon[1]}/{start_lat_lon[0]}/{end_lat_lon[1]}/{end_lat_lon[0]}+i{sample_interval_km}+d",
        grid=grid,
    )

    x = df[2].to_numpy()
    z = df[3].to_numpy()
    assert x.shape == z.shape == (n_x_samples,)
    return x, z


@functools.cache
def make_elev_boundary(region, start_lat_lon, end_lat_lon, length_km, n_x_samples):
    elev_grid = pygmt.datasets.load_earth_relief(resolution="01s", region=region)
    elev_grid.data = (
        -1e-3 * elev_grid.data
    )  # convert from m of elevation to km of depth
    return create_boundary_from_grid(
        elev_grid, start_lat_lon, end_lat_lon, length_km, n_x_samples
    )


model_start_lat_lon = (57.953815, -152.717426)
model_end_lat_lon = (55.910548, -151.01544)


def make_v_model():
    # Start, end, and distance obtained by running make_model_space.
    # The line points from Kodiak out towards the nodes.
    vm = VelocityModel(
        start_lat_lon=model_start_lat_lon,
        end_lat_lon=model_end_lat_lon,
        length_km=250,
        min_depth_km=0,
        max_depth_km=60,
    )

    x_samples, z_samples = 10 * vm.length_km, 10 * (vm.max_depth_km - vm.min_depth_km)

    # Layer 1 is water with a constant v=1.5 km/s.
    vm.add_layer_f(1, lambda x, z: 1.5 + x * 0 + z * 0, x_samples, z_samples)

    # Layer 3 is the mantle. I'll start by giving it a constant v=8.0 km/s.
    vm.add_layer_f(3, lambda x, z: 8.0 + x * 0 + z * 0, x_samples, z_samples)

    # Boundary 1 is just elevation/bathymetry.
    grid_region = (
        vm.start_lat_lon[1],
        vm.end_lat_lon[1],
        vm.end_lat_lon[0],
        vm.start_lat_lon[0],
    )
    elev_x, elev_z = make_elev_boundary(
        grid_region,
        vm.start_lat_lon,
        vm.end_lat_lon,
        vm.length_km,
        1 * vm.length_km,
    )

    vm.add_boundary(1, elev_x, elev_z)

    # Layer 2 is the crust.
    vm.add_layer_f(
        2,
        lambda x, z: 2
        + np.clip(
            1
            * np.sqrt(
                (elev_x[:, np.newaxis] - x) ** 2 + (elev_z[:, np.newaxis] - z) ** 2
            ),
            0,
            4,
        ),
        1 * vm.length_km,
        z_samples,
    )

    # Boundary 2 is the Moho (slab depth). I'll start by setting it to 10 km.
    vm.add_boundary_f(2, lambda x: 10 + 0 * x, vm.length_km // 5)

    return vm
