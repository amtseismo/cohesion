#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 17:23:15 2022

Download SVI elevation

@author: amt
"""

import pygmt

# reg=[-124.7,-122.8,48,49.4]
reg=[-125,-122,47,50]
grid = pygmt.datasets.load_earth_relief(resolution='15s', region=reg)
grid.to_netcdf('svi_topo.nc')