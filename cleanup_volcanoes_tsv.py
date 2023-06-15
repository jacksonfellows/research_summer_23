# Helper to cleanup the .tsv returned by the NOAA service (https://www.ngdc.noaa.gov/hazel/view/hazards/volcano/loc-search).

import pandas as pd
import os


def cleanup_tsv(filename):
    df = pd.read_csv(filename, sep="\t")
    df = df.rename(
        columns={"Volcano Name": "name", "Latitude": "lat", "Longitude": "lon"}
    )
    df.to_csv(
        f"{os.path.splitext(filename)[0]}_cleaned.csv",
        columns=["name", "lat", "lon"],
        index=False,
    )
