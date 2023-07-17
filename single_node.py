import numpy as np

import utils
import unbinned


def make_unbinned_traces(node_code, lineno):
    st = utils.load_node(node_code, lineno)
    st.filter("bandpass", freqmin=3, freqmax=20)
    offsets = [1e-3 * utils.source_receiver_offset(tr) for tr in st]
    traces = unbinned.UnbinnedTraces(
        min(offsets), max(offsets), st[0].stats.sampling_rate
    )
    for offset, tr in zip(offsets, st.traces):
        tr_norm = tr.data / np.max(np.abs(tr.data))
        traces.add_trace(tr_norm, offset)
    return traces
