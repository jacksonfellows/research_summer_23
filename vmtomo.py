import numpy as np
import struct
from matplotlib import pyplot as plt
import pygmt
import functools


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
        plt.title(f"Kodiak Cross-Section (Vertical Exaggeration = {vert_exag:0.1f}x)")
        plt.xlabel("Profile Distance (km)")
        plt.ylabel("Depth (km)")
        # Add scale bar
        plt.colorbar(im, location="bottom", shrink=0.7, label="Velocity (km/s)")
        if show:
            plt.show()

    @classmethod
    def load(cls, filename):
        with open(filename, "rb") as f:
            # Assumes little-endian with 4-byte integers and 4-byte
            # floats.
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
        grid="/Users/jackson/Downloads/Slab2_AComprehe/alu_slab2_dep_02.23.18.grd",
        profile=f"{start_lat_lon[1]}/{start_lat_lon[0]}/{end_lat_lon[1]}/{end_lat_lon[0]}+l{length_km}k+i{sample_interval_km}+d",
    )
    z = df[3].to_numpy()
    z = -z  # convert to depth
    assert z.shape == (n_x_samples,)
    return z


def build_vm():
    # Model location.
    start_lat_lon = (57.953815, -152.717426)
    end_lat_lon = (55.910548, -151.01544)
    region = (
        start_lat_lon[1],
        end_lat_lon[1],
        end_lat_lon[0],
        start_lat_lon[0],
    )

    # Velocity model info.
    x1, x2 = 0, 250
    z1, z2 = -2, 60
    nx = 4 * (x2 - x1)
    nz = 8 * (z2 - z1)
    nr = 3  # Number of boundaries.
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
