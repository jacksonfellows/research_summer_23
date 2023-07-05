from matplotlib import pyplot as plt

import binned


def pick(bt, mode="matrix", red_vel=6.0, ylim=(-2, 8), vmin=None, vmax=None, scale=1):
    if mode == "matrix":
        fig, axs = bt.plot_mat(
            show=False, red_vel=red_vel, ylim=ylim, vmin=vmin, vmax=None
        )
    elif mode == "squiggle":
        fig, axs = bt.plot(show=False, red_vel=red_vel, ylim=ylim, scale=scale)
    else:
        assert False

    xs = []
    ys = []
    ts = []
    (line,) = axs[0].plot([], [], "x", color="red")

    def onclick_callback(event):
        offset = bt.round_to_bin(event.xdata)
        red_time = event.ydata
        time = red_time + offset / red_vel
        if event.key == "shift":
            print(f"picking {offset:0.2f} km, {time:0.2f} s")
            xs.append(offset)
            ys.append(red_time)
            ts.append(time)
            line.set_data(xs, ys)
            axs[0].draw_artist(line)
            axs[0].figure.canvas.update()
        elif event.key == "control":
            i = xs.index(offset)
            xs.pop(i)
            ys.pop(i)
            time = ts.pop(i)
            print(f"deleting pick {offset:0.2f} km, {time:0.2f} s")
            line.set_data(xs, ys)
            axs[0].draw_artist(line)
            axs[0].figure.canvas.update()

    cid = fig.canvas.mpl_connect("button_press_event", onclick_callback)

    plt.show()

    print(xs, ts)
