from dataclasses import dataclass


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
