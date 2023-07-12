import numpy as np
import pickle
import tempfile
import os

import rayfile


def test_load_rayfile():
    rays_loaded = rayfile.Rayfile.load("test_data/test_rayfile").rays()
    with open("test_data/test_rayfile_rays.pickle", "rb") as f:
        rays_pickled = pickle.load(f)
    assert all(
        np.all(x1 == x2) and np.all(z1 == z2)
        for (x1, z1), (x2, z2) in zip(rays_loaded, rays_pickled)
    )


def test_roundtrip():
    rayfile_path = "test_data/test_rayfile"
    rf = rayfile.Rayfile.load(rayfile_path)
    with tempfile.TemporaryDirectory() as tmp:
        rayfile_roundtrip_path = os.path.join(tmp, "r")
        rf.dump(rayfile_roundtrip_path)
        with open(rayfile_path, "rb") as f1:
            with open(rayfile_roundtrip_path, "rb") as f2:
                assert f1.read() == f2.read()
