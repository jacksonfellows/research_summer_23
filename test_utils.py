import utils

# Let's make sure I don't mess up stat_lat_lon before I rewrite stuff.
# I created these test cases by taking a random sample of stations.

nodal_station_test_stat_lat_lon = (
    (3087, (57.7368, -152.5144)),
    (2025, (57.61771, -152.356901)),
    (2025, (57.61771, -152.356901)),
    (9062, (57.507258, -152.450571)),
    (3056, (57.7826, -152.5714)),
    (9093, (57.582677, -152.456459)),
    (9004, (57.435117, -152.347248)),
    (9126, (57.643565, -152.476271)),
    (2008, (57.596825, -152.422067)),
    (9055, (57.488592, -152.452586)),
)


def test_stat_lat_lon():
    for stat, lat_lon in nodal_station_test_stat_lat_lon:
        assert utils.stat_lat_lon(stat) == lat_lon


def test_calc_min_max_offsets_km():
    # I made these test cases by looking at a few per-shot plots I had
    # made previously.
    min_offset, max_offset = utils.calc_min_max_offsets_km([991])
    assert min_offset < 30
    assert 70 < max_offset

    min_offset, max_offset = utils.calc_min_max_offsets_km([1300])
    assert min_offset < 160
    assert 200 < max_offset

    min_offset, max_offset = utils.calc_min_max_offsets_km([3497])
    assert min_offset < 100
    assert 130 < max_offset
