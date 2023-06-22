import numpy as np
import struct
from matplotlib import pyplot as plt
import pygmt
import functools
from obspy.clients.fdsn.client import Client


import utils
import profile_info


class VMTOMO_VM:
    def __init__(self, nx, nz, nr, x1, x2, z1, z2, zrf, idr, vel):
        self.nx = nx
        self.nz = nz
        self.nr = nr
        self.x1 = x1
        self.x2 = x2
        self.z1 = z1
        self.z2 = z2
        assert zrf.shape == (nx * nr,)
        self.zrf = zrf
        assert idr.shape == (nx * nr,)
        self.idr = idr
        assert vel.shape == (nx * nz * (nr + 1),)
        self.vel = vel

    def plot(self, show=True):
        x = np.linspace(self.x1, self.x2, self.nx)
        z = np.linspace(self.z1, self.z2, self.nz)
        xx, zz = np.meshgrid(x, z)
        velocity_cmap = "gist_rainbow"
        vert_exag = 2
        contour_levels = np.linspace(0, 10, 21)
        contour_levels_label = np.linspace(0, 10, 11)
        im = None  # Needed for the colorbar.
        # Plot the layers.
        for layer_i in range(self.nr + 1):
            b = (
                self.z1
                if layer_i == 0
                else self.zrf[(layer_i - 1) * self.nx : layer_i * self.nx]
            )
            vv = (
                self.vel[
                    layer_i * self.nx * self.nz : (layer_i + 1) * self.nx * self.nz
                ]
                .reshape((self.nx, self.nz))
                .T
            )
            vv_masked = np.ma.masked_array(vv, zz < b)
            im = plt.imshow(
                vv_masked,
                vmin=0,
                vmax=10,
                extent=(self.x1, self.x2, self.z2, self.z1),
                cmap=velocity_cmap,
                aspect=vert_exag,
            )
            CS = plt.contour(
                xx,
                zz,
                vv_masked,
                levels=contour_levels,
                colors="k",
                linestyles="dashed",
                linewidths=0.5,
            )
            plt.clabel(CS, inline=1, fontsize=8, levels=contour_levels_label)
        # Plot the boundaries.
        for boundary_i in range(self.nr):
            b = self.zrf[boundary_i * self.nx : (boundary_i + 1) * self.nx]
            plt.plot(x, b, color="k")
        # Add title and axis labels
        plt.title(f"Vertical Exaggeration = {vert_exag:0.1f}x")
        plt.xlabel("Profile Distance (km)")
        plt.ylabel("Depth (km)")
        # Add scale bar
        plt.colorbar(im, location="bottom", shrink=0.7, label="Velocity (km/s)")
        if show:
            plt.show()

    # load and dump assume little-endian with 4-byte integers and
    # 4-byte floats.

    @classmethod
    def load(cls, filename):
        with open(filename, "rb") as f:
            nx = struct.unpack("<i", f.read(4))[0]
            nz = struct.unpack("<i", f.read(4))[0]
            nr = struct.unpack("<i", f.read(4))[0]
            x1 = struct.unpack("<f", f.read(4))[0]
            x2 = struct.unpack("<f", f.read(4))[0]
            z1 = struct.unpack("<f", f.read(4))[0]
            z2 = struct.unpack("<f", f.read(4))[0]
            zrf_n_bytes = nx * nr * 4
            zrf = np.frombuffer(f.read(zrf_n_bytes), dtype="<f")
            idr_n_bytes = zrf_n_bytes
            idr = np.frombuffer(f.read(idr_n_bytes), dtype="<i")
            nv = nx * nz
            vel_n_bytes = nv * (nr + 1) * 4
            vel = np.frombuffer(f.read(vel_n_bytes), dtype="<f")
            return cls(nx, nz, nr, x1, x2, z1, z2, zrf, idr, vel)

    def dump(self, filename):
        with open(filename, "wb") as f:
            f.write(struct.pack("<i", self.nx))
            f.write(struct.pack("<i", self.nz))
            f.write(struct.pack("<i", self.nr))
            f.write(struct.pack("<f", self.x1))
            f.write(struct.pack("<f", self.x2))
            f.write(struct.pack("<f", self.z1))
            f.write(struct.pack("<f", self.z2))
            f.write(self.zrf.tobytes())
            f.write(self.idr.tobytes())
            f.write(self.vel.tobytes())


def create_boundary_from_grid(grid, start_lat_lon, end_lat_lon, length_km, n_x_samples):
    """
    Create a boundary from a grid.
    """
    sample_interval_km = length_km / (n_x_samples + 1)
    df = pygmt.grdtrack(
        # +i sets the sampling interval
        # +d tells it to output distances along profile
        profile=f"{start_lat_lon[1]}/{start_lat_lon[0]}/{end_lat_lon[1]}/{end_lat_lon[0]}+i{sample_interval_km}+d",
        grid=grid,
    )

    z = df[3].to_numpy()
    assert z.shape == (n_x_samples,)
    return z


@functools.cache
def make_elev_boundary(region, start_lat_lon, end_lat_lon, length_km, n_x_samples):
    elev_grid = pygmt.datasets.load_earth_relief(resolution="01s", region=region)
    elev_grid.data = (
        -1e-3 * elev_grid.data
    )  # convert from m of elevation to km of depth
    return create_boundary_from_grid(
        elev_grid, start_lat_lon, end_lat_lon, length_km, n_x_samples
    )


@functools.cache
def make_slab_boundary(region, start_lat_lon, end_lat_lon, length_km, n_x_samples):
    sample_interval_km = length_km / (n_x_samples + 1)
    df = pygmt.grdtrack(
        grid="./alu_slab2_dep_02.23.18.grd",
        profile=f"{start_lat_lon[1]}/{start_lat_lon[0]}/{end_lat_lon[1]}/{end_lat_lon[0]}+l{length_km}k+i{sample_interval_km}+d",
    )
    z = df[3].to_numpy()
    z = -z  # convert to depth
    assert z.shape == (n_x_samples,)
    return z


def build_vm(start_lat_lon, end_lat_lon, nx, nz, nr, x1, x2, z1, z2):
    # Model location.
    region = (
        start_lat_lon[1],
        end_lat_lon[1],
        end_lat_lon[0],
        start_lat_lon[0],
    )

    # Velocity model info.
    zrf = np.zeros((nr, nx), dtype="f")
    idr = np.zeros((nr, nx), dtype="i")
    vel = np.zeros((nr + 1, nx, nz), dtype="f")

    # Might need this.
    x = np.linspace(x1, x2, nx)
    z = np.linspace(z1, z2, nz)
    xx, zz = np.meshgrid(x, z)

    # Set layer ids. Not really sure what these are for.
    for i in range(nr):
        idr[i] = i + 1

    # Set the boundary locations.

    # Boundary 1: Sea level
    zrf[0] = 0

    # Boundary 2: Elevation/bathymetry
    elev_z = make_elev_boundary(region, start_lat_lon, end_lat_lon, x2 - x1, nx)
    zrf[1] = elev_z

    # Boundary 3: Moho
    slab_z = make_slab_boundary(region, start_lat_lon, end_lat_lon, x2 - x1, nx)
    moho_z = slab_z + 8
    undefined = np.isnan(moho_z)
    x_rest = x[undefined]
    last_z = moho_z[np.logical_not(undefined)][-1]
    final_z = 15
    final_x = x[-1]
    moho_z[undefined] = last_z + (final_z - last_z) / (final_x - x_rest[0]) * (
        x_rest - x_rest[0]
    )
    zrf[2] = moho_z

    # Set the layer velocities.

    # Layer 1: Air
    vel[0] = 0.3

    # Layer 2: Water
    vel[1] = 1.5

    # Layer 3: Crust
    start_sed_v = 2
    final_sed_v = 5
    sed_thickness = 2
    start_crust_v = 5
    final_crust_v = 7.1  # Add 0.1 so we see the contour.
    crust_thickness = 8
    depth = zz - elev_z
    vel[2] = np.clip(
        np.where(
            depth < 0,
            start_sed_v,  # above the surface
            np.where(
                depth < sed_thickness,
                # in sediment
                start_sed_v + (final_sed_v - start_sed_v) / sed_thickness * depth,
                # in crust
                start_crust_v
                + (final_crust_v - start_crust_v)
                / crust_thickness
                * (depth - sed_thickness),
            ),
        ),
        start_sed_v,
        final_crust_v,
    ).T

    # Layer 4: Mantle
    vel[3] = 8

    return VMTOMO_VM(
        nx, nz, nr, x1, x2, z1, z2, zrf.flatten(), idr.flatten(), vel.flatten()
    )


trev = lambda x: tuple(reversed(x))


@functools.cache
def find_earthquakes(min_lat, max_lat, min_lon, max_lon, min_mag=0):
    """
    Finds earthquake events within the given bounds. Returns a 2d
    array with lon, lat, z, mag as columns. z will be depth in km.
    """
    client = Client("IRIS")
    cat = client.get_events(
        minlatitude=min_lat,
        maxlatitude=max_lat,
        minlongitude=min_lon,
        maxlongitude=max_lon,
        minmagnitude=min_mag,
    )
    locs = np.zeros((len(cat.events), 4))
    for i, event in enumerate(cat.events):
        origin = event.preferred_origin()
        if (
            origin.longitude is None
            or origin.latitude is None
            or origin.depth is None
            or event.preferred_magnitude().mag is None
        ):
            locs[i] = np.nan
        else:
            locs[i][0] = origin.longitude
            locs[i][1] = origin.latitude
            locs[i][2] = 1e-3 * origin.depth
            locs[i][3] = event.preferred_magnitude().mag
    locs_masked = np.ma.masked_array(locs, np.isnan(locs))
    return locs_masked


@functools.cache
def find_earthquakes_along_profile(start_lat_lon, end_lat_lon, min_mag=0):
    min_lat = min(start_lat_lon[0], end_lat_lon[0])
    max_lat = max(start_lat_lon[0], end_lat_lon[0])
    min_lon = min(start_lat_lon[1], end_lat_lon[1])
    max_lon = max(start_lat_lon[1], end_lat_lon[1])
    return find_earthquakes(
        min_lat - 0.5, max_lat + 0.5, min_lon - 0.5, max_lon + 0.5, min_mag=min_mag
    )


class Profile:  # Could have a better name.
    """
    A study profile. Stores the associated velocity model and
    metadata.
    """

    def __init__(
        self, name, start_lat_lon, end_lat_lon, vm, node_codes=None, shot_nos=None
    ):
        self.name = name
        self.start_lat_lon = start_lat_lon
        self.end_lat_lon = end_lat_lon
        self.vm = vm
        # Set some helper variables from vm
        self.x1, self.x2, self.z1, self.z2 = vm.x1, vm.x2, vm.z1, vm.z2
        # Nodes
        if node_codes is not None:
            self.node_x, self.node_z = self.project_nodes(node_codes)
        else:
            self.node_x, self.node_z = np.zeros(0), np.zeros(0)
        # Shots
        if shot_nos is not None:
            self.shot_x, self.shot_z = self.project_shots(shot_nos)
        else:
            self.shot_x, self.shot_z = np.zeros(0), np.zeros(0)

    def project_nodes(self, node_codes):
        """
        Helper function to project nodes to profile. Returns arrays x, z in
        the model space.
        """
        df = utils._node_df[utils._node_df.code.isin(node_codes)]
        node_locs = np.column_stack((df.lon, df.lat, df.elev_m))
        projected = pygmt.project(
            data=node_locs,
            center=trev(self.start_lat_lon),
            endpoint=trev(self.end_lat_lon),
            unit=True,
        )
        x = projected[3].to_numpy()
        z = projected[2].to_numpy()
        z = -1e-3 * z  # Convert z from m of elevation to km of depth.
        return x, z

    def project_shots(self, shot_nos):
        """
        Helper function to project shots to profile. Returns arrays x, z in
        the model space.
        """
        df = utils._shot_df[utils._shot_df.shotno.isin(shot_nos)]
        node_locs = np.column_stack((df.lon, df.lat))
        projected = pygmt.project(
            data=node_locs,
            center=trev(self.start_lat_lon),
            endpoint=trev(self.end_lat_lon),
            unit=True,
        )
        x = projected[2].to_numpy()
        z = np.zeros(x.shape)  # All shots happen at sea level.
        return x, z

    def project_earthquakes(self, earthquake_locs, max_dist_from_profile_km):
        """
        Projects an array of earthquake lon, lat, z, mag to profile.
        Returns arrays x, z, mag in the model space.
        """
        projected = pygmt.project(
            data=earthquake_locs,
            center=trev(self.start_lat_lon),
            endpoint=trev(self.end_lat_lon),
            unit=True,
        )
        x = projected[4]
        dist_from_profile = projected[5]
        z = projected[2].to_numpy()
        mag = projected[3].to_numpy()
        assert x.shape == z.shape == mag.shape == dist_from_profile.shape
        # Only return events that are in the profile
        valid = (
            (self.x1 < x)
            & (x < self.x2)
            & (self.z1 < z)
            & (z < self.z2)
            & (np.abs(dist_from_profile) < max_dist_from_profile_km)
        )
        return np.column_stack((x[valid], z[valid], mag[valid]))

    def plot(self, projected_earthquakes=None):
        plt.suptitle(self.name)
        self.vm.plot(show=False)
        # Earthquakes. I think it makes sense to not store them with
        # the profile since we'll only need them for plotting.
        if projected_earthquakes is not None:
            quake_x, quake_z, quake_mag = (
                projected_earthquakes[:, 0],
                projected_earthquakes[:, 1],
                projected_earthquakes[:, 2],
            )
            quake_s = (10 * quake_mag / quake_mag.max()) ** 2
            plt.scatter(
                quake_x,
                quake_z,
                quake_s,
                facecolors="none",
                edgecolors="white",
                label="Earthquake",
            )
        # Nodes
        plt.plot(
            self.node_x,
            self.node_z,
            "o",
            color="purple",
            label="Nodal Station",
            markersize=4.0,
        )
        # Shots
        plt.plot(
            self.shot_x,
            self.shot_z,
            "x",
            color="blue",
            label="Airgun Shot",
            markersize=4.0,
        )
        plt.legend()
        plt.show()
