import numpy as np
import obspy
import obspy.signal.trigger
import obspy.signal.filter
from matplotlib import pyplot as plt
import scipy

import utils
import binned


def load_trace(per_shot_st, node_code):
    for tr in per_shot_st:
        if (
            tr.stats.segy.trace_header.trace_number_within_the_original_field_record
            == node_code
        ):
            return tr


def process_trace(tr):
    tr_ = tr.copy()
    # bandpass
    tr_.filter("bandpass", freqmin=3, freqmax=20, zerophase=True)
    # sta / lta
    sta_s = 0.05
    lta_s = 5.0
    tr_.data = obspy.signal.trigger.classic_sta_lta(
        tr_.data,
        int(sta_s * tr_.stats.sampling_rate),
        int(lta_s * tr_.stats.sampling_rate),
    )
    tr_.data /= np.abs(tr_.data).max()
    # trim off what the sta/lta messes up
    # tr_.trim(tr_.stats.starttime + lta_s, tr_.stats.endtime)
    return tr_


def correlate_traces(t1, t2):
    assert (
        t1.data.shape == t2.data.shape
        and t1.stats.sampling_rate == t2.stats.sampling_rate
    )
    offset1_m = utils.source_receiver_offset(t1)
    offset2_m = utils.source_receiver_offset(t2)
    diff_m = offset2_m - offset1_m
    print(f"t2 has {diff_m:0.2f} m more offset then t1")
    min_v = 1000  # m/s
    most_lag = -diff_m / min_v * t1.stats.sampling_rate
    max_lag = max(most_lag, 0)
    min_lag = min(most_lag, 0)
    print(min_lag, max_lag)
    correlation = scipy.signal.correlate(t1.data, t2.data)
    lags = scipy.signal.correlation_lags(t1.data.shape[0], t2.data.shape[0])
    correlation_masked = np.ma.masked_array(
        correlation, (lags < min_lag) | (lags > max_lag)
    )
    lag = lags[np.argmax(correlation_masked)]
    print(lag)
    if lag > 0:
        t1d = t1.data[lag:]
        t2d = t2.data
    else:
        t1d = t1.data
        t2d = t2.data[-lag:]
    fig, axs = plt.subplots(2, 1, sharex=True)
    axs[0].plot(t1.data, label="t1")
    axs[0].plot(t2.data, label="t2")
    axs[1].plot(t1d)
    axs[1].plot(t2d)
    axs[0].legend()
    plt.show()


def megashot(center_shot, node_code, shots_per_side):
    t0 = process_trace(load_trace(utils.load_shot(center_shot), node_code))

    t_stacked = np.zeros(t0.data.shape)

    shots_to_stack = tuple(
        range(center_shot - shots_per_side, center_shot + shots_per_side + 1)
    )

    for other_shot in shots_to_stack:
        t_ = process_trace(load_trace(utils.load_shot(other_shot), node_code))
        offset1_m = utils.source_receiver_offset(t0)
        offset2_m = utils.source_receiver_offset(t_)
        diff_m = offset2_m - offset1_m
        min_v = 1500  # m/s
        max_v = 8000  # m/s
        least_lag = diff_m / max_v * t0.stats.sampling_rate
        most_lag = diff_m / min_v * t0.stats.sampling_rate
        max_lag = max(most_lag, least_lag)
        min_lag = min(most_lag, least_lag)
        correlation = scipy.signal.correlate(t0.data, t_.data)
        lags = scipy.signal.correlation_lags(t0.data.shape[0], t_.data.shape[0])
        correlation_masked = np.ma.masked_array(
            correlation, (lags < min_lag) | (lags > max_lag)
        )
        lag = lags[np.argmax(correlation_masked)]
        lag_s = lag / t0.stats.sampling_rate
        print(f"lag of {lag_s} s for change in offset of {diff_m}")
        if lag == 0:
            t_stacked += t_
        elif lag > 0:
            # shift t_ to the left
            t_stacked[:-lag] += t_[lag:]
        else:
            # shift t_ to the right
            t_stacked[-lag:] += t_[:lag]

    fig, axs = plt.subplots(2, 1, sharex=True)
    for other_shot in shots_to_stack:
        axs[0].plot(
            t0.times(),
            process_trace(load_trace(utils.load_shot(other_shot), node_code)).data,
            label=f"{other_shot}",
        )
    axs[0].legend()
    axs[1].plot(t0.times(), t_stacked / len(shots_to_stack))
    plt.show()


def megashot_all_nodes(center_shot, shots_per_side, min_v, max_v):
    # First load all the shots.
    shotnos = range(center_shot - shots_per_side, center_shot + shots_per_side + 1)
    shot_sts = {}
    for shotno in shotnos:
        print(f"loading shot {shotno}")
        st = utils.load_shot(shotno)
        print(f"applying bandpass to shot {shotno}")
        utils.bandpass_stream_inplace(st)
        print(f"applying sta-lta to shot {shotno}")
        for t in st:
            t.data = obspy.signal.trigger.classic_sta_lta(
                t.data, 0.05 * t.stats.sampling_rate, 5.0 * t.stats.sampling_rate
            )
            t.data /= np.abs(t.data).max()  # normalize
            if np.isnan(t.data).any():
                print(f"nan! - skipping trace {t}")
                st.remove(t)
        # Sort by node.
        print(f"sorting traces by node")
        st.traces.sort(
            key=lambda t: t.stats.segy.trace_header.trace_number_within_the_original_field_record
        )
        print(f"removing duplicate traces")
        last_n = None
        for t in st:
            n = t.stats.segy.trace_header.trace_number_within_the_original_field_record
            if n == last_n:
                print(f"removing duplicate trace {t}")
                st.remove(t)
            last_n = n
        shot_sts[shotno] = st

    min_offset, max_offset = utils.calc_min_max_offsets_km(shotnos)
    sample_len = 60  # s
    sampling_rate = int(shot_sts[center_shot][0].stats.sampling_rate)
    bt = binned.BinnedTraces(min_offset, max_offset, 0.25, sample_len, sampling_rate)

    # Go through each station and merge those shots.
    for ts in zip(*[st.traces for st in shot_sts.values()]):
        t0 = [
            t
            for t in ts
            if t.stats.segy.trace_header.energy_source_point_number == center_shot
        ][0]
        assert all(
            t.stats.segy.trace_header.trace_number_within_the_original_field_record
            == t0.stats.segy.trace_header.trace_number_within_the_original_field_record
            for t in ts
        )
        t_stacked = np.zeros(t0.data.shape)

        for i, t_ in enumerate(ts):
            offset1_m = utils.source_receiver_offset(t0)
            offset2_m = utils.source_receiver_offset(t_)
            diff_m = offset2_m - offset1_m
            least_lag = diff_m / max_v * t0.stats.sampling_rate
            most_lag = diff_m / min_v * t0.stats.sampling_rate
            max_lag = max(most_lag, least_lag)
            min_lag = min(most_lag, least_lag)
            correlation = scipy.signal.correlate(t0.data, t_.data)
            lags = scipy.signal.correlation_lags(t0.data.shape[0], t_.data.shape[0])
            correlation_masked = np.ma.masked_array(
                correlation, (lags < min_lag) | (lags > max_lag)
            )
            lag = lags[np.argmax(correlation_masked)]
            if lag == 0:
                t_stacked += t_
            elif lag > 0:
                # shift t_ to the left
                t_stacked[:-lag] += t_[lag:]
            else:
                # shift t_ to the right
                t_stacked[-lag:] += t_[:lag]

        bt.add_trace(t_stacked / len(shotnos), 1e-3 * utils.source_receiver_offset(t0))

    return bt
