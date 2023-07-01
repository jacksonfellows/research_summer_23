import numpy as np
import struct


def load(filename):
    with open(filename, "rb") as f:
        rays = []
        num = struct.unpack("<i", f.read(4))[0]  # number of rays
        for _ in range(num):
            inr = struct.unpack("<i", f.read(4))[0]  # instrument ID
            npr = struct.unpack("<i", f.read(4))[0]  # number of picks
            nrr = struct.unpack("<i", f.read(4))[0]  # number of ray points
            # Not sure what these are for.
            isk = np.frombuffer(f.read(4 * npr), dtype="<i")  # Shot number?
            fas = np.frombuffer(f.read(4 * npr), dtype="<i")
            ttp = np.frombuffer(f.read(4 * npr), dtype="<f")  # Picked travel time?
            etp = np.frombuffer(f.read(4 * npr), dtype="<f")  # Error in pick?
            tca = np.frombuffer(f.read(4 * npr), dtype="<f")  # Calculated travel time?
            len_ = np.frombuffer(f.read(4 * npr), dtype="<i")
            xry = np.frombuffer(f.read(4 * nrr), "<f")
            zry = np.frombuffer(f.read(4 * nrr), "<f")
            start = 0
            for l in len_.flat:
                rays.append((xry[start : start + l], zry[start : start + l]))
                start += l
        # is there still stuff in the ray file?
        return rays
