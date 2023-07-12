import numpy as np
import struct
from dataclasses import dataclass


@dataclass
class RayfileChunk:
    inr: int
    npr: int
    nrr: int
    # I don't think dataclass actually checks the types.
    isk: "int[npr]"
    fas: "int[npr]"
    ttp: "float[npr]"
    etp: "float[npr]"
    tca: "float[npr]"
    len_: "int[npr]"
    xry: "float[nrr]"
    zry: "float[nrr]"


class Rayfile:
    def __init__(self, num, chunks):
        self.num = num
        self.chunks = chunks

    @classmethod
    def load(cls, filename):
        with open(filename, "rb") as f:
            num = struct.unpack("<i", f.read(4))[0]  # number of rays
            # Each chunk could have a different size so I can't just have a big array.
            chunks = []
            for _ in range(num):
                inr = struct.unpack("<i", f.read(4))[0]  # instrument ID
                npr = struct.unpack("<i", f.read(4))[0]  # number of picks
                nrr = struct.unpack("<i", f.read(4))[0]  # number of ray points
                # Not sure what these are for.
                # Shot number?
                isk = np.frombuffer(f.read(4 * npr), dtype="<i")
                fas = np.frombuffer(f.read(4 * npr), dtype="<i")
                # Picked travel time?
                ttp = np.frombuffer(f.read(4 * npr), dtype="<f")
                # Error in pick?
                etp = np.frombuffer(f.read(4 * npr), dtype="<f")
                # Calculated travel time?
                tca = np.frombuffer(f.read(4 * npr), dtype="<f")
                len_ = np.frombuffer(f.read(4 * npr), dtype="<i")
                xry = np.frombuffer(f.read(4 * nrr), "<f")
                zry = np.frombuffer(f.read(4 * nrr), "<f")
                chunks.append(
                    RayfileChunk(inr, npr, nrr, isk, fas, ttp, etp, tca, len_, xry, zry)
                )

        return cls(num, chunks)

    def rays(self):
        rays = []
        chunk: RayfileChunk
        for chunk in self.chunks:
            start = 0
            for l in chunk.len_.flat:
                rays.append(
                    (chunk.xry[start : start + l], chunk.zry[start : start + l])
                )
                start += l
        return rays
