# Place to put the metadata for each profile.

profile_1 = {
    "start_lat_lon": (57.953815, -152.717426),
    "end_lat_lon": (55.910548, -151.01544),
    "nx": 1000,
    "nz": 496,
    "nr": 3,
    "x1": 0,
    "x2": 250,
    "z1": -2,
    "z2": 60,
    "name": "Kodiak Profile 1",
}

profile_2 = {
    "start_lat_lon": (57.94745, -154.210595),
    "end_lat_lon": (55.107529, -151.77111),
    "nx": 4 * 350,
    "nz": 8 * (82 - -2),
    "nr": 3,
    "x1": 0,
    "x2": 350,
    "z1": -2,
    "z2": 84,  # Same ratio of length/depth as profile 1.
    "name": "Kodiak Profile 2",
}
