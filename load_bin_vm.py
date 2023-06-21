import numpy as np
import struct
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


class VMTOMO_VM:
    def __init__(self, filename):
        with open(filename, "rb") as f:
            # Assumes little-endian with 4-byte integers and 4-byte
            # floats.
            self.nx = struct.unpack("<i", f.read(4))[0]
            self.nz = struct.unpack("<i", f.read(4))[0]
            self.nr = struct.unpack("<i", f.read(4))[0]
            self.x1 = struct.unpack("<f", f.read(4))[0]
            self.x2 = struct.unpack("<f", f.read(4))[0]
            self.z1 = struct.unpack("<f", f.read(4))[0]
            self.z2 = struct.unpack("<f", f.read(4))[0]
            zrf_n_bytes = self.nx * self.nr * 4
            self.zrf = np.frombuffer(f.read(zrf_n_bytes), dtype="<f")
            idr_n_bytes = zrf_n_bytes
            self.idr = np.frombuffer(f.read(idr_n_bytes), dtype="<i")
            nv = self.nx * self.nz
            vel_n_bytes = nv * (self.nr + 1) * 4
            self.vel = np.frombuffer(f.read(vel_n_bytes), dtype="<f")

    def plot(self):
        x = np.linspace(self.x1, self.x2, self.nx)
        z = np.linspace(self.z1, self.z2, self.nz)
        xx, zz = np.meshgrid(x, z)
        velocity_cmap = "gist_rainbow"
        vert_exag = 2
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
            plt.imshow(
                vv_masked,
                vmin=0,
                vmax=10,
                extent=(self.x1, self.x2, self.z2, self.z1),
                cmap=velocity_cmap,
                aspect=vert_exag,
            )
        # Plot the boundaries.
        for boundary_i in range(self.nr):
            b = self.zrf[boundary_i * self.nx : (boundary_i + 1) * self.nx]
            plt.plot(x, b, color="k")
        # Add title and axis labels
        plt.title(f"Kodiak Cross-Section (Vertical Exaggeration = {vert_exag:0.1f}x)")
        plt.xlabel("Profile Distance (km)")
        plt.ylabel("Depth (km)")
        # Add scale bar
        plt.colorbar(location="bottom", shrink=0.7, label="Velocity (km/s)")
        plt.show()
