import pygmt

import utils

def plot_lines(linenos):
    fig = pygmt.Figure()
    region = [-158, -150, 54, 59]
    grid = pygmt.datasets.load_earth_relief(resolution = '15s', region=region)
    fig.grdimage(grid=grid, projection='M15c', cmap='geo', frame=True)
    pygmt.makecpt(cmap='categorical', series=(min(linenos), max(linenos) + 1, 1))
    print('plotting nodes')
    fig.plot(x=utils.stat_df.lon, y=utils.stat_df.lat, style='c0.06c', fill='red')
    for lineno in linenos:
        print(f'plotting line {lineno}')
        shot_nos = utils.shots_for_line(lineno)
        df = utils.shot_df.loc[shot_nos]
        fig.plot(x=-df.lon, y=df.lat, cmap=True, zvalue=lineno, pen='thick,+z', label=lineno)
    fig.legend()
    fig.show()
