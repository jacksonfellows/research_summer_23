import numpy as np
from matplotlib import pyplot as plt

# 2d example


def plot_source_stations(source, stations):
    plt.xlim(0, 10)
    plt.ylim(0, 10)
    plt.scatter(source[0], source[1])
    plt.scatter(stations[:, 0], stations[:, 1])
    plt.show()


def simulate_arrivals(source, stations, v):
    return (1 / v) * np.linalg.norm(source[:2] - stations, axis=1) + source[2]


def calc_G(source, stations, v):
    d = source[:2] - stations
    t = np.sum(d**2, axis=1) ** (-1 / 2)
    G = np.column_stack(
        ((d[:, 0] / v) * t, (d[:, 1] / v) * t, np.ones(t.shape))  # x  # y  # t
    )
    return G


def invert(source_0, stations, arrivals, v):
    G = calc_G(source_0, stations, v)
    arrivals_0 = simulate_arrivals(source_0, stations, v)
    arrivals_delta = arrivals - arrivals_0
    source_delta, residuals, _, _ = np.linalg.lstsq(G, arrivals_delta)
    source_1 = source_0 + source_delta
    return source_1


def demo_inversion(source, stations, arrivals, v, N):
    source_guess = np.random.normal(0, 0.01, 3) + np.concatenate(
        (stations[arrivals.argmin()], np.array([arrivals.min()]))
    )
    guesses = np.zeros((N + 1, 3))
    for n in range(N):
        guesses[n] = source_guess
        source_guess = invert(source_guess, stations, arrivals, v)
    guesses[N] = source_guess
    plt.plot(guesses[:, 0], guesses[:, 1], "-x")
    plt.plot(source[0], source[1], "rx")
    plt.plot(stations[:, 0], stations[:, 1], "bo")
    for i in range(stations.shape[0]):
        plt.annotate(f"{arrivals[i]:0.2f}", (stations[i, 0], stations[i, 1]))
    plt.legend()
    plt.show()


def demo_with_random_data():
    plt.xlim(0, 10)
    plt.ylim(0, 10)
    source = np.concatenate(
        (10 * np.random.random(2), np.zeros(1))
    )  # earthquake source (x,y,t)
    stations = 10 * np.random.random((7, 2))  # station locations
    v = 1
    arrivals = np.random.normal(0, 0.1, stations.shape[0]) + simulate_arrivals(
        source, stations, v
    )
    demo_inversion(source, stations, arrivals, v, 10)
