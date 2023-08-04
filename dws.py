import numpy as np
from matplotlib import pyplot as plt


def load_dws_mask(path):
    a = np.loadtxt(path)
    x, z, d = a[:, 0], a[:, 1], a[:, 3]
    # hacky
    n0s = len(x[x == 0])
    new_shape = (len(x) // n0s, n0s)
    x = x.reshape(new_shape)
    z = z.reshape(new_shape)
    d = d.reshape(new_shape)
    d[d == -0.5] = 1
    return x, z, d
