#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 15:55:49 2023

Make a plot of LFEs

@author: amt
"""

import pandas as pd
import numpy as np
import pygmt
import datetime

# load data
locs=np.load("lfe_locations.npy")  
df = pd.DataFrame(locs, columns = ['id','lat','lon','dep'])

# load slab model
fm=pd.read_csv('mccrory_model_2012.dat',names=["lon", "lat", "depth"])

# geographic region
region=[np.min(df['lon'])-0.1, np.max(df['lon'])+0.1, np.min(df['lat'])-0.15, np.max(df['lat'])+0.1]

# Load sample grid (3 arc-minutes global relief) in target area
grid = pygmt.datasets.load_earth_relief(resolution="03s", region=region)

# plot map
fig = pygmt.Figure()
fig.basemap(region=region, 
            projection="M15c", frame=True)

# plot bathy
watergrid = pygmt.grdclip(grid, above=[0, np.nan])
watercpt=pygmt.makecpt(
    cmap="/Users/amt/Documents/GMT_scripts/colorpalettes/bath_111.cpt",
    series='-700,0,10',
    continuous=True)
fig.grdimage(grid=watergrid, cmap=watercpt, nan_transparent=True)

# plot topo
landgrid = pygmt.grdclip(grid, below=[0, np.nan])
#dlandgrid = pygmt.grdgradient(grid=landgrid, radiance=[270, 30], normalize='e0.4')
landcpt=pygmt.makecpt(
    cmap="/Users/amt/Documents/GMT_scripts/colorpalettes/gray.cpt",
    series='0,1000,10',
    continuous=True)
#fig.grdimage(grid=landgrid, shading=dlandgrid, cmap=landcpt, nan_transparent=True)
fig.grdimage(grid=landgrid, cmap=landcpt, nan_transparent=True)

# # plot shorelines
# fig.coast(shorelines="1p,103/103/103,solid",
#           area_thresh=0)

# make colormap
pygmt.makecpt(cmap="plasma", series=[28,47])

fig.plot(
    x=df['lon'],
    y=df['lat'],
    fill=df['dep'],
    cmap=True,
    style="c0.3c",
    pen="black",
)

# plot slab contours
for d in np.arange(-50,-20,5):
    tmp=fm[fm['depth']==d]
    fig.plot(
        x=tmp['lon'],
        y=tmp['lat'],
        pen="2p,103/103/103",
        transparency=50,
        )
  
# colorbar if you want it
fig.colorbar(position="n0.7/0.1c+w5c/0.5c+h", frame="af+lDepth (km)")

with fig.inset(
    position="jBL+o0.1c",
    box="+gwhite+p2p",
    region=[-126, -121, 46, 51.5],
    projection="M3c",
):
    # Highlight the Japan area in "lightbrown"
    # and draw its outline with a pen of "0.2p".
    fig.coast(
        dcw="US,CA+p0.2p",
        borders="2/thin",
        shorelines="thin",
        land="lightyellow",
        water="lightblue",
        area_thresh=10000,
    )
    # Plot a rectangle ("r") in the inset map to show the area of the main
    # figure. "+s" means that the first two columns are the longitude and
    # latitude of the bottom left corner of the rectangle, and the last two
    # columns the longitude and latitude of the uppper right corner.
    rectangle = [np.min(df['lon'])-0.1, np.max(df['lon'])+0.1, np.min(df['lat'])-0.1, np.max(df['lat'])+0.1]
    fig.plot(data=rectangle, style="r+s", pen="2p,blue")

fig.show()