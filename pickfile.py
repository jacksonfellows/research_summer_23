import pandas as pd

PICKFILE_COLUMNS = ["station", "shot", "phase", "offset", "tt", "error"]


def save_picks(picks, filename):
    picks.to_csv(
        filename,
        columns=PICKFILE_COLUMNS,
        sep=" ",
        header=False,
        index=False,
        float_format="%.3f",
    )


def load_picks(filename):
    return pd.read_csv(
        filename,
        names=PICKFILE_COLUMNS,
        sep=" ",
        header=None,
    )


def load_pick_list(filename):
    """Load picks from a list of pick files in filename."""
    pick_dfs = []
    with open(filename, "r") as f:
        for pickfile in f.readlines():
            pick_dfs.append(load_picks(f"picks/{pickfile.strip()}"))
    return pd.concat(pick_dfs)
