import numpy as np
import struct
from matplotlib import pyplot as plt


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


def load_vm_from_file(filename):
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
        return VMTOMO_VM(nx, nz, nr, x1, x2, z1, z2, zrf, idr, vel)
