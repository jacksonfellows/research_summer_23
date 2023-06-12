import multiprocessing
from download_shots import download_shot

# multiprocessing can't handle running a function defined in the same module so I have to make a new module and import download_shot.


def download_shots_multi(shots):
    with multiprocessing.Pool(8) as p:
        p.map(download_shot, shots, chunksize=1)


def download_shots_multi_retry(shots):
    try:
        download_shots_multi(shots)
    except:
        download_shots_multi_retry(shots)
