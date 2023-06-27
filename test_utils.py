import obspy

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
        assert utils.node_lat_lon(stat) == lat_lon


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


# Again, good to have some sanity checks before I rewrite stuff.


def test_shots_for_line():
    # Manually created these by looking at the shot spreadsheet.
    l1 = utils.shots_for_line(1)
    assert l1[0] == 991 and l1[-1] == 1391
    l15 = utils.shots_for_line(15)
    assert l15[0] == 16175 and l15[-1] == 17685


def test_shot_lat_lon():
    # Manually created these by looking at the shot spreadsheet. Note
    # that in the spreadsheet the longitudes are positive when they
    # should be negative.
    assert utils.shot_lat_lon(991) == (57.220942, -152.08304)
    assert utils.shot_lat_lon(5079) == (56.288621, -153.34762)
    assert utils.shot_lat_lon(5323) == (55.496843, -152.67804)
    assert utils.shot_lat_lon(7093) == (55.077562, -152.92654)
    assert utils.shot_lat_lon(9305) == (55.214477, -153.63619)
    assert utils.shot_lat_lon(11324) == (55.554297, -154.53278)
    assert utils.shot_lat_lon(13205) == (55.963331, -155.36349)
    assert utils.shot_lat_lon(15256) == (55.058355, -155.29734)
    assert utils.shot_lat_lon(23167) == (54.201261, -156.90926)
    assert utils.shot_lat_lon(22084) == (53.684722, -156.42038)


def test_shot_starttime():
    for shotno in [1000, 1100, 1200, 1300]:
        st = utils.load_shot(shotno)
        real_starttime = obspy.UTCDateTime(
            utils._shot_df[utils._shot_df.shotno == shotno].iloc[0].time
        )
        time_diff_s = abs(st[0].stats.starttime - real_starttime)
        # Assert time is the same within the error of our sampling rate.
        assert time_diff_s < (1 / st[0].stats.sampling_rate)
