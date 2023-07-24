import struct
from dataclasses import dataclass

import numpy as np
import pandas as pd


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
    def __init__(self, num_chunks, chunks):
        self.num_chunks = num_chunks
        self.chunks = chunks

    @classmethod
    def load(cls, filename):
        with open(filename, "rb") as f:
            num_chunks = struct.unpack("<i", f.read(4))[0]  # Number of chunks
            # Each chunk could have a different size so I can't just have a big array.
            chunks = []
            for _ in range(num_chunks):
                inr = struct.unpack("<i", f.read(4))[0]  # Instrument ID
                npr = struct.unpack("<i", f.read(4))[0]  # Number of picks
                nrr = struct.unpack("<i", f.read(4))[0]  # Number of ray points
                # Shot number
                isk = np.frombuffer(f.read(4 * npr), dtype="<i")
                # Phase
                fas = np.frombuffer(f.read(4 * npr), dtype="<i")
                # Picked travel time
                ttp = np.frombuffer(f.read(4 * npr), dtype="<f")
                # Error in pick
                etp = np.frombuffer(f.read(4 * npr), dtype="<f")
                # Calculated travel time
                tca = np.frombuffer(f.read(4 * npr), dtype="<f")
                len_ = np.frombuffer(f.read(4 * npr), dtype="<i")
                xry = np.frombuffer(f.read(4 * nrr), "<f")
                zry = np.frombuffer(f.read(4 * nrr), "<f")
                chunks.append(
                    RayfileChunk(inr, npr, nrr, isk, fas, ttp, etp, tca, len_, xry, zry)
                )

        return cls(num_chunks, chunks)

    def dump(self, filename):
        with open(filename, "wb") as f:
            f.write(struct.pack("<i", self.num_chunks))
            for chunk in self.chunks:
                f.write(struct.pack("<i", chunk.inr))
                f.write(struct.pack("<i", chunk.npr))
                f.write(struct.pack("<i", chunk.nrr))
                f.write(chunk.isk.tobytes())
                f.write(chunk.fas.tobytes())
                f.write(chunk.ttp.tobytes())
                f.write(chunk.etp.tobytes())
                f.write(chunk.tca.tobytes())
                f.write(chunk.len_.tobytes())
                f.write(chunk.xry.tobytes())
                f.write(chunk.zry.tobytes())

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

    def calculated_tts(self):
        # Goal is to return enough info to draw travel-time curves.
        stations = []
        shots = []
        phases = []
        tt_picks = []
        tt_calcs = []
        errors = []
        for chunk in self.chunks:
            stations.extend(chunk.npr * [chunk.inr])
            shots.extend(chunk.isk)
            phases.extend(chunk.fas)
            tt_picks.extend(chunk.ttp)
            tt_calcs.extend(chunk.tca)
            errors.extend(chunk.etp)
        return pd.DataFrame(
            {
                "station": stations,
                "shot": shots,
                "phase": phases,
                "tt": tt_picks,
                "tt_calc": tt_calcs,
                "error": errors,
            }
        )


def cat_rayfiles(input_paths, output_path):
    num_chunks = 0
    chunks = []
    for path in input_paths:
        rf = Rayfile.load(path)
        num_chunks += rf.num_chunks
        chunks.extend(rf.chunks)
    Rayfile(num_chunks, chunks).dump(output_path)
