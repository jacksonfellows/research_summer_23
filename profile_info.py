from dataclasses import dataclass
import utils
import numpy as np
import pandas as pd
import pygmt


@dataclass
class Profile:
    """Place to keep the metadata for each profile"""

    start_lat_lon: tuple
    end_lat_lon: tuple
    x1: int
    x2: int
    z1: int
    z2: int
    name: str

    def project_shots(self, shot_nos):
        """
        Helper function to project shots to profile.
        """
        df = utils._shot_df[utils._shot_df.shotno.isin(shot_nos)]
        shot_locs = np.column_stack((df.lon, df.lat))
        projected = pygmt.project(
            data=shot_locs,
            center=(self.start_lat_lon[1], self.start_lat_lon[0]),
            endpoint=(self.end_lat_lon[1], self.end_lat_lon[0]),
            unit=True,
        )
        shot_depth_km = 1e-3 * 12
        return pd.DataFrame(
            {
                # Need the to_numpy() because otherwise pandas tries
                # to create an index based on the index of the
                # shot_df.
                "shotno": df.shotno.to_numpy(),
                "x": projected[2],
                "y": projected[3],
                "z": np.full(len(projected[2]), shot_depth_km),
            }
        )


profile_1 = Profile(
    start_lat_lon=(57.953815, -152.717426),
    end_lat_lon=(55.910548, -151.01544),
    x1=0,
    x2=250,
    z1=-2,
    z2=60,
    name="Kodiak Profile 1",
)

profile_2 = Profile(
    start_lat_lon=(57.94745, -154.210595),
    end_lat_lon=(55.107529, -151.77111),
    x1=0,
    x2=350,
    z1=-2,
    z2=84,  # Same ratio of length/depth as profile 1.
    name="Kodiak Profile 2",
)
