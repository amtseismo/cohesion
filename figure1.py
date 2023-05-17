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
import matplotlib

pygmt.config(COLOR_MODEL="rgb")

# load data
locs=np.load("lfe_locations.npy")  
df = pd.DataFrame(locs, columns = ['id','lat','lon','dep'])

# load stations
stas=pd.read_csv('svi_stations.dat', delim_whitespace=True)

# load slab model
fm=pd.read_csv('mccrory_model_2012.dat',names=["lon", "lat", "depth"])

# geographic region
region=[np.min(df['lon'])-0.1, np.max(df['lon'])+0.1, np.min(df['lat'])-0.125, np.max(df['lat'])+0.1]

# Load sample grid (3 arc-minutes global relief) in target area
grid = pygmt.datasets.load_earth_relief(resolution="03s", region=region)

# plot map
fig = pygmt.Figure()
fig.basemap(region=region, 
            projection="M15c", 
            frame=True)

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

# plot stations
fig.plot(
    x=stas['LON'],
    y=stas['LAT'],
    fill='50/50/50',
    style="i0.45c",
    pen="black",
)

# plot lfes
fig.plot(
    x=df['lon'],
    y=df['lat'],
    fill=df['dep'],
    cmap=True,
    style="c0.35c",
    pen="0/0/0",
)

lfeid=254

# Plotting text annotations for LFE
fig.text(text=lfeid, 
         x=df[df['id']==lfeid]['lon'], 
         y=df[df['id']==lfeid]['lat'],
         font="18p,Helvetica-Bold",
         offset='0/0.5')

# Plotting text annotations for stations
fig.text(text=['PGC','LZB','TWKB','SSIB','SILB'], 
         x=[stas[stas['STATION']=='PGC']['LON'].values[0], 
            stas[stas['STATION']=='LZB']['LON'].values[0], 
            stas[stas['STATION']=='TWKB']['LON'].values[0],
            stas[stas['STATION']=='SSIB']['LON'].values[0],
            stas[stas['STATION']=='SILB']['LON'].values[0]], 
         y=[stas[stas['STATION']=='PGC']['LAT'].values[0], 
            stas[stas['STATION']=='LZB']['LAT'].values[0], 
            stas[stas['STATION']=='TWKB']['LAT'].values[0],
            stas[stas['STATION']=='SSIB']['LAT'].values[0],
            stas[stas['STATION']=='SILB']['LAT'].values[0]],
         font="17p,Helvetica-Bold",
         offset='0/0.5')

fig.text(text='A', 
         x=-124.36, 
         y=48.96,
         font="36p,Helvetica",
         offset='0/0')

# plot slab contours
for d in np.arange(-50,-20,5):
    tmp=fm[fm['depth']==d]
    fig.plot(
        x=tmp['lon'],
        y=tmp['lat'],
        pen="2p,103/103/103",
        transparency=50,
        )

# ADD SCALEBAR
pygmt.config(MAP_SCALE_HEIGHT="10p", FONT_ANNOT_PRIMARY=14, MAP_TICK_PEN_PRIMARY="1p")
fig.basemap(map_scale="g-124/49+w20k")
  
# colorbar if you want it
fig.colorbar(position="n0.6/0.09c+w5c/0.5c+h", frame="af+lDepth (km)")

# get outline color from colormap
cmap = matplotlib.colormaps['plasma']
rgba=cmap(0.25,bytes=True)

# plot distance along strike bar
meanlat=48.510
meanlon=-123.733
fig.plot(
    x=[meanlon],
    y=[meanlat],
    style="v0.6c+e",
    direction=[[132],[7]],
    pen="2p",
    fill="black")
fig.text(
    text='+40 km',
    x=-124.25,
    y=48.8,
    font="14p,Helvetica-Bold,50/50/50")
fig.plot(
    x=[meanlon],
    y=[meanlat],
    style="v0.6c+e",
    direction=[[-48],[7]],
    pen="2p",
    fill="black")
fig.text(
    text='-40 km',
    x=-123.23,
    y=48.23,
    font="14p,Helvetica-Bold,50/50/50")

# Make inset
fig.shift_origin(xshift="-1.5c",yshift="-1.5c")
fig.coast(
    projection="G-124/48/4/5c", 
    dcw="US,CA+p1p,100/100/100",
    resolution="h",
    region="g", 
    frame='a',
    land="200/200/200",
    area_thresh=100,
    water="143/212/255") 

inset_region=[-126, -121, 47, 50.5]
fig.plot(data=[[region[0], region[2], region[1], region[3]]], 
          style="r+s", 
          region=inset_region, 
          pen="2p,"+str(rgba[0])+"/"+str(rgba[1])+"/"+str(rgba[2])
          )

# Plotting text annotations US
fig.text(text='US', 
         x=-121, 
         y=46,
         font="18p,Helvetica-Bold,50/50/50",
         offset='0/0.5')

# Plotting text annotations US
fig.text(text='CA', 
         x=-122, 
         y=49.7,
         font="18p,Helvetica-Bold,50/50/50",
         offset='0/0.5')

# pygmt.config(MAP_FRAME_TYPE="plain")
# fig.shift_origin(xshift="0.5c",yshift="0.5c")
# inset_proj="M4c"
# inset_region=[-126, -121, 47, 50.5]
# fig.coast(
#         region=inset_region,
#         dcw="US,CA+p0.5p,150/150/150",
#         projection=inset_proj,
#         land="white",
#         water="148/216/255",
#         area_thresh=10000,
#         frame=["WSne"]
#     )

# # ADD COASTLINES TO INSET
# fig.plot(data=[[region[0], region[2], region[1], region[3]]], 
#           style="r+s", 
#           projection=inset_proj, 
#           region=inset_region, 
#           pen="2p,"+str(rgba[0])+"/"+str(rgba[1])+"/"+str(rgba[2])
#           )

# fig.plot(data=[[inset_region[0], inset_region[2], inset_region[1], inset_region[3]-0.05]], 
#          style="r+s", 
#          projection=inset_proj, 
#          region=inset_region, 
#          pen="4p,150/150/150"
#          )

fig.show()
fig.savefig('f1pA.png',transparent=True, dpi=300)
fig.savefig('f1pA.eps',dpi=300)