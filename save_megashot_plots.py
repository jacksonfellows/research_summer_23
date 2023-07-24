import os

from matplotlib import pyplot as plt

import megashot


def make_plots():
    plt.set_cmap("viridis_r")
    for shotno in range(996, 1392 - 5):
        path = f"megashot_plots/{shotno}.png"
        if os.path.exists(path):
            print(f"skipping shot {shotno}")
        else:
            print(f"plotting shot {shotno}")
            try:
                if shotno < 1250:
                    bt = megashot.megashot_all_nodes(shotno, 5, 5000, 7000, "median")
                else:
                    bt = megashot.megashot_all_nodes(shotno, 5, 7000, 9000, "median")
                bt.plot_mat(show=False)
                print(f"saving shot {shotno}")
                plt.savefig(path)
            except:
                pass
