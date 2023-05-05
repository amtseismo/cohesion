#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 08:45:15 2022

Moment rates from LFEs

@author: amt
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as dates
import datetime
from obspy import geodetics
import utm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import matplotlib.patches as patches

def mag2mom(mag):
    mom=10**(1.5*(mag+10.7))/(1e7)
    return mom

def mom2mag(mom):
    mag=np.log10(mom*(1e7))/1.5-10.7
    return mag


# for 2003 slow slip event
start=datetime.datetime(2012,9,6)
finish=datetime.datetime(2012,9,22)

# read slab model
fm=pd.read_csv('~/Documents/sse_nucleation/mccrory_model_2012.dat',names=["lon", "lat", "depth"])
fm['depth']*=-1

# read tremor data
tdf = pd.read_csv('tremor_events-2009-08-06T00_00_00-2022-06-08T23_59_59.csv',header=0,delim_whitespace=(True))
tdf=tdf[(tdf['lat']>48.119166666666665-0.1) & (tdf['lat']<48.930166666666665+0.1)
        & (tdf['lon']>-124.33683333333333-0.1) & (tdf['lon']<-123.10216666666666+0.1)]
dateform=tdf['starttime1']+' '+tdf['starttime2']
datenums=[datetime.datetime.strptime(d,"%Y-%m-%d %H:%M:%S") for d in dateform]
tdf['datenums']=datenums
tdf['datetime']=pd.to_datetime(datenums)
tdf=tdf[(tdf['datetime']>start) & (tdf['datetime']<finish)]

# plot LFE family locations
locs=np.load("lfe_locations.npy")   
# plt.scatter(locs[:,2], locs[:,1], c=locs[:,3], s=40, alpha=0.5)
# for ii in range(len(locs)):
#     plt.text(locs[ii,2], locs[ii,1], str(int(locs[ii,0])))
# plt.colorbar()
contour=fm[fm['depth']==35]
contour=contour[(contour['lat']>48.475) & (contour['lat']<48.55)]
# plt.plot(contour['lon'].values,contour['lat'].values,'ko',markersize=12)


# get local strike
_,az1,az2=geodetics.base.gps2dist_azimuth(contour.iloc[0]['lat'],contour.iloc[0]['lon'],contour.iloc[1]['lat'],contour.iloc[1]['lon'])
meanlat=np.mean(locs[:,1])
meanlon=np.mean(locs[:,2])
# plt.plot(meanlon,meanlat,'ro',markersize=12)

coords=utm.from_latlon(locs[:,1],locs[:,2])
tcoords=utm.from_latlon(tdf['lat'].values,tdf['lon'].values)
meanx=np.mean(coords[0])
meany=np.mean(coords[1])
# make rotation tensor
theta = np.radians(az2)
c, s = np.cos(theta), np.sin(theta)
R = np.array(((c, -s), (s, c)))
rot=np.matmul(R,np.vstack((coords[0]-np.mean(coords[0]),coords[1]-np.mean(coords[1]))))
trot=np.matmul(R,np.vstack((tcoords[0]-np.mean(coords[0]),tcoords[1]-np.mean(coords[1]))))
# plt.figure()
# plt.scatter(rot[0], rot[1], c=locs[:,3], s=40, alpha=0.5)
# for ii in range(len(locs)):
#     plt.text(rot[0][ii], rot[1][ii], str(int(locs[ii,0])))
# plt.colorbar()
# plt.axis('equal')

lfes=pd.DataFrame()
lfes['id']=locs[:,0]
lfes['lat']=locs[:,1]
lfes['lon']=locs[:,2]
lfes['dep']=locs[:,3]
lfes['utmx']=coords[0]
lfes['utmy']=coords[1]
lfes['rotx']=rot[0]/1000
lfes['roty']=rot[1]/1000
lfes=lfes.sort_values(by=['roty'])

# add new coords to tremor dataframe
tdf['utmx']=tcoords[0]
tdf['utmy']=tcoords[1]
tdf['rotx']=trot[0]/1000
tdf['roty']=trot[1]/1000

# Load stats
stats=pd.read_csv('total_mag_detect_0000_cull.txt', delim_whitespace=True, names=('famid','ymd','hour','secofhour','mag','no_phases'))

# get things needed for datetime
year=np.floor((stats['ymd'].values+20000000)/10000)
month=stats['ymd'].astype(str).str[-4:-2].astype(np.int64)
day=stats['ymd'].astype(str).str[-2:].astype(np.int64)
minute=np.floor(stats['secofhour'].values/60)
second=stats['secofhour'].values-60*minute

# make new df and fill it
df = pd.DataFrame() #columns = ['year','month','day','hour','minute','second'])
df['year']=year
df['month']=month
df['day']=day
df['hour']=stats['hour']
df['minute']=minute
df['second']=second
df['datetime']=pd.to_datetime(df)
plt_dates = dates.date2num(df['datetime'].values)
stats['dates']= plt_dates
stats['moment']=mag2mom(stats['mag'].values)
stats['datetime']=df['datetime']

# make the mask 
td=finish-start
td_hours=td/datetime.timedelta(hours=1)
atd=np.arange(0,td_hours,4)
asd=np.arange(-68,78,4)
mask=np.zeros((len(atd),len(asd)))
for ii, hours in enumerate(atd):
    for jj, ydist in enumerate(asd):
        tmp=tdf[((tdf['datetime']-start)/np.timedelta64(1, 'h') > hours) &
                ((tdf['datetime']-start)/np.timedelta64(1, 'h') < hours+4) &
                (tdf['roty'] > ydist) &
                (tdf['roty'] < ydist+4)]
        if len(tmp) > 0:
            mask[ii,jj]=1
plt.imshow(mask)
        
        
'''
fig = plt.figure(tight_layout=True, figsize=(15,10))
gs = gridspec.GridSpec(11, 25, hspace=0.1, wspace=0.1)
ax = fig.add_subplot(gs[1:,:24])
ax1 = fig.add_subplot(gs[0,:24])
ax2 = fig.add_subplot(gs[1:,24])
ax.set_xlabel('Date',fontsize=16)
ax.set_ylabel('Along strike distance (km)',fontsize=16)
mr=np.zeros(int(td_hours))
mrt = np.arange(start, finish, datetime.timedelta(hours=1)).astype(datetime.datetime)
for ii in range(int(td_hours)):
    tmp=stats[(stats['datetime']>=start+datetime.timedelta(hours=ii)) 
              & (stats['datetime']<start+datetime.timedelta(hours=ii+1))]  
    mr[ii]=np.sum(tmp['moment'].values)
    print(len(tmp))
for ii, fid in enumerate(lfes['id'].values):
    # get dates
    tmp=stats[(stats['famid']==fid) & (stats['datetime']>=start)
              & (stats['datetime']<finish)]
    if len(tmp) > 0:
        dates = [pd.to_datetime(d) for d in tmp['datetime']]
        ec=ax.scatter(dates,lfes.iloc[ii]['roty']*1/1000*np.ones(len(tmp)),
                   s=40,c=tmp['mag'].values,vmin=np.min(tmp['mag'].values),
                   vmax=np.max(tmp['mag'].values), cmap='viridis_r',
                   edgecolors='k',linewidth=0.2)

ax.set_xlim((start,finish))
ax.tick_params(axis='both', which='major', labelsize=14)
ax.tick_params(axis='both', which='minor', labelsize=14)
ax.set_ylim((-40,40))

# rect = patches.Rectangle((start+datetime.timedelta(days=10.6), 0), 
#                           datetime.timedelta(days=8), 40, 
#                           linewidth=1, facecolor='white',alpha=0.8)
# ax.add_patch(rect)
# rect = patches.Rectangle((start, -40), 
#                           datetime.timedelta(days=5.5), 35, 
#                           linewidth=1, facecolor='white',alpha=0.8)
# ax.add_patch(rect)

ax1.fill_between(mrt,mr,color=(.2,.2,.2), alpha=0.5)
ind=np.where(mr==np.max(mr))[0][0]
ax1.arrow(mrt[ind],mr[ind],0.00,0,width=0.1,color='k')
ax1.text(mrt[ind]+datetime.timedelta(days=0.3),0.9*mr[ind],
         "Peak Hourly Moment Rate: "+"{:.2E}".format(mr[ind]),fontsize=14)
ax1.set_xlim((start,finish))
ax1.axis('off')
ax1.text(start,0.6*np.max(mr),'Hourly Moment Rate',fontsize=16)
ax1.set_xticklabels([]) 
ax1.set_yticklabels([]) 

cbar = plt.colorbar(ec, cax=ax2)
cbar.set_label('Magnitude', fontsize=14)
cbar.ax.tick_params(labelsize=14)

'''