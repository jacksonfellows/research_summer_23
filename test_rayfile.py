import numpy as np
import pickle

import rayfile


def test_load_rayfile():
    rays_loaded = rayfile.load("test_data/test_rayfile")
    with open("test_data/test_rayfile_rays.pickle", "rb") as f:
        rays_pickled = pickle.load(f)
    assert all(
        np.all(x1 == x2) and np.all(z1 == z2)
        for (x1, z1), (x2, z2) in zip(rays_loaded, rays_pickled)
    )
