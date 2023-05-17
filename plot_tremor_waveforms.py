#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 20:48:02 2022

Plot tremor

@author: amt
"""

import pandas as pd
from obspy import Stream, Trace, UTCDateTime
from datetime import datetime, timedelta
from obspy.core import read
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
plt.rcParams['font.family'] = 'Helvetica'
import matplotlib.gridspec as gridspec

# setup plot
fig = plt.figure(figsize=(9,7), tight_layout=True)
gs = gridspec.GridSpec(7, 5)
ax = fig.add_subplot(gs[:5,:])
ax1 = fig.add_subplot(gs[5:,0])
ax2 = fig.add_subplot(gs[5:,1])
ax3 = fig.add_subplot(gs[5:,2])
ax4 = fig.add_subplot(gs[5:,3])
ax5 = fig.add_subplot(gs[5:,4])

# # color
cmap = matplotlib.cm.get_cmap('plasma')
rgba = cmap((36.75-28)/(47-28))

files=['20050913.CN.LZB..BHN.mseed','20050913.PO.SSIB..HHN.mseed',
       '20050913.CN.PGC..BHN.mseed','20050913.PO.SILB..HHN.mseed','20050913.PO.TWKB..HHN.mseed']

t1=UTCDateTime(2005,9,13)

st=Stream()
pre_filt = (0.9, 1, 10, 15.0)
for file in files:
    tr=read(file)
    tr.merge()
    # inv = read_inventory("RESP."+tr[0].stats.station+"."+tr[0].stats.network+".."+tr[0].stats.channel)
    # tr.remove_response(inventory=inv, output='VEL', pre_filt=pre_filt)
    st+=tr
st.detrend()

# filter the data
st.filter('bandpass',freqmin=1,freqmax=10)

# decimate
st=st.resample(5)

# trim to hour
st.trim(UTCDateTime(2005,9,13,3,0),UTCDateTime(2005,9,13,4,0))

# trim to lfe time
lfetimes=[122.275,155.075,220.275,238.000,265.900,306.675,316.275,361.525,
          394.350,435.750,478.400,487.975,493.750,532.800,593.675,611.975,
          883.925,916.975,1147.100,1188.050,1234.650,1282.900,1334.825,
          1385.275,1438.700,1489.025,1592.300,1657.825,1983.950,2088.450,
          2128.300,2443.425,2491.025,2734.775,2800.500,2898.200,3224.825,
          3305.075,3474.100,3498.675,3503.850]

# label figure
ax.text(900,-0.5,'B',weight='bold', fontsize=32)

# plot tremor waveforms
fac=20
for ii in range(len(st)):
    rgba = cmap(ii/(len(st)))
    ax.plot(st[ii].times('relative'),st[ii].data/np.median(fac*np.abs(st[ii].data))+ii,color=(0.6,0.6,0.6))
    for jj in range(len(lfetimes)):
        #lfest=st.copy().trim(UTCDateTime(lfetimes[ind])-10,UTCDateTime(lfetimes[ind])+10)
        start=int(lfetimes[jj]*st[ii].stats.sampling_rate)-int(4*st[ii].stats.sampling_rate)
        print(start)
        finish=int(lfetimes[jj]*st[ii].stats.sampling_rate)+int(4*st[ii].stats.sampling_rate)
        ax.plot(st[ii].times('relative')[start:finish],st[ii].data[start:finish]/np.median(fac*np.abs(st[ii].data))+ii,color=rgba)
    text=ax.text(st[ii].times('relative')[0]+950,ii-0.3,st[ii].stats.station,fontsize=12, weight='bold',color=rgba)
    text.set_path_effects([matplotlib.patheffects.Stroke(linewidth=1.5, foreground=(0.8,0.8,0.8)),
                       matplotlib.patheffects.Normal()])
ax.plot([1200,2400],[-0.59,-0.59],'k',linewidth=2,color=(0.25,0.25,0.25))
ax.text(1720,-0.65,'20 minutes',color=(0.25,0.25,0.25),fontsize=14)
ax.set_xlim((898,2700))
ax.set_ylim((-0.6,len(st)-0.46))
ax.invert_yaxis()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.annotate('', xy=(1488, 0.1), xytext=(1320, 2.5),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color=(0.25,0.25,0.25)), fontsize=14, color=(0.25,0.25,0.25))
ax.annotate('', xy=(1488, 1.1), xytext=(1320, 2.5),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color=(0.25,0.25,0.25)), fontsize=14)
ax.annotate('', xy=(1488, 2.1), xytext=(1320, 2.5),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color=(0.25,0.25,0.25)), fontsize=14)
ax.annotate('', xy=(1488, 3.1), xytext=(1320, 2.5),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color=(0.25,0.25,0.25)), fontsize=14)
ax.annotate('', xy=(1488, 4.1), xytext=(1320, 2.5),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color=(0.25,0.25,0.25)), fontsize=14)
ax.text(1130,2.55,'detections',color=(0.25,0.25,0.25),fontsize=14)

s1=UTCDateTime(2005,9,13,3,15)
ax.plot([901,901,920],[4,4.43,4.43],color=(0.25,0.25,0.25), linewidth=3)
ax.text(928,4.48,str(s1.year)+'/'+str(s1.month)+'/'+str(s1.day)
        +' '+str(s1.hour).zfill(2)+':'+str(s1.minute).zfill(2)+':'+str(s1.second).zfill(2),color=(0.25,0.25,0.25), fontsize=12)


start=UTCDateTime(2005,9,13,3,0)
# ax.set_xlabel('Time relative to '+str(start.year)+'/'+str(start.month)+'/'+str(start.day)+' '
#               +str(start.hour).zfill(2)+':'+str(start.minute).zfill(2)+':'+str(start.second).zfill(2)+
#               ' (s)',fontsize=14)

# label figure
ax1.text(0,1.1,'C',weight='bold', fontsize=32)

# plot first stack
with open('LZB_Nstack.npy', 'rb') as f:
    clip = np.load(f)
t=np.arange(-11,19.01,0.01)
ax1.plot(t,clip,color=cmap(0), linewidth=2, path_effects=[matplotlib.patheffects.Stroke(linewidth=1, foreground=(0.8,0.8,0.8)), 
                        matplotlib.patheffects.Normal()])
# ax1.plot(t,clip,color=rgba,
#          path_effects=[matplotlib.patheffects.Stroke(linewidth=1, foreground='black'), 
#                        matplotlib.patheffects.Normal()], alpha=0.5)
ax1.set_xlim((0,8))
ax1.get_xaxis().set_visible(False)
ax1.get_yaxis().set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.set_ylim((-1.,1.2))
text=ax1.text(0.5,0.6,st[0].stats.station,fontsize=12, weight='bold',color=cmap(0))
text.set_path_effects([matplotlib.patheffects.Stroke(linewidth=1.5, foreground=(0.8,0.8,0.8)),
                   matplotlib.patheffects.Normal()])

# plot second stack
with open('SSIB_Nstack.npy', 'rb') as f:
    clip = np.load(f)
ax2.plot(t,clip,color=cmap(0.2), linewidth=2, path_effects=[matplotlib.patheffects.Stroke(linewidth=1, foreground=(0.8,0.8,0.8)), 
                        matplotlib.patheffects.Normal()])
# ax1.plot(t,clip,color=rgba,
#          path_effects=[matplotlib.patheffects.Stroke(linewidth=1, foreground='black'), 
#                        matplotlib.patheffects.Normal()], alpha=0.5)
ax2.set_xlim((0,8))
ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.set_ylim((-1.,1.2))
text=ax2.text(0.5,0.6,st[1].stats.station,fontsize=12, weight='bold',color=cmap(0.2))
text.set_path_effects([matplotlib.patheffects.Stroke(linewidth=1.5, foreground=(0.8,0.8,0.8)),
                   matplotlib.patheffects.Normal()])

# plot third stack
with open('PGC_Nstack.npy', 'rb') as f:
    clip = np.load(f)
ax3.plot(t,clip,color=cmap(0.4), linewidth=2, path_effects=[matplotlib.patheffects.Stroke(linewidth=1, foreground=(0.8,0.8,0.8)), 
                        matplotlib.patheffects.Normal()])
ax3.plot([0,8],[1.1,1.1],'k',linewidth=2,color=(0.25,0.25,0.25))
ax3.text(2.,1.2,'8 seconds',color=(0.25,0.25,0.25),fontsize=14)
ax3.set_xlim((0,8))
ax3.get_xaxis().set_visible(False)
ax3.get_yaxis().set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.set_ylim((-1.,1.2))
text=ax3.text(0.5,0.6,st[2].stats.station,fontsize=12, weight='bold',color=cmap(0.4))
text.set_path_effects([matplotlib.patheffects.Stroke(linewidth=1.5, foreground=(0.8,0.8,0.8)),
                   matplotlib.patheffects.Normal()])

# plot fourth stack
with open('SILB_Nstack.npy', 'rb') as f:
    clip = np.load(f)
ax4.plot(t,clip,color=cmap(0.6), linewidth=2, path_effects=[matplotlib.patheffects.Stroke(linewidth=1, foreground=(0.8,0.8,0.8)), 
                        matplotlib.patheffects.Normal()])
# ax1.plot(t,clip,color=rgba,
#          path_effects=[matplotlib.patheffects.Stroke(linewidth=1, foreground='black'), 
#                        matplotlib.patheffects.Normal()], alpha=0.5)
ax4.set_xlim((0,8))
ax4.get_xaxis().set_visible(False)
ax4.get_yaxis().set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.spines['bottom'].set_visible(False)
ax4.spines['left'].set_visible(False)
ax4.set_ylim((-1.,1.2))
text=ax4.text(0.5,0.6,st[3].stats.station,fontsize=12, weight='bold',color=cmap(0.6))
text.set_path_effects([matplotlib.patheffects.Stroke(linewidth=1.5, foreground=(0.8,0.8,0.8)),
                   matplotlib.patheffects.Normal()])

# plot final stack
with open('TWKB_Nstack.npy', 'rb') as f:
    clip = np.load(f)
ax5.plot(t,clip,color=cmap(0.8), linewidth=2, path_effects=[matplotlib.patheffects.Stroke(linewidth=1, foreground=(0.8,0.8,0.8)), 
                        matplotlib.patheffects.Normal()])
# ax1.plot(t,clip,color=rgba,
#          path_effects=[matplotlib.patheffects.Stroke(linewidth=1, foreground='black'), 
#                        matplotlib.patheffects.Normal()], alpha=0.5)
ax5.set_xlim((0,8))
ax5.get_xaxis().set_visible(False)
ax5.get_yaxis().set_visible(False)
ax5.spines['top'].set_visible(False)
ax5.spines['right'].set_visible(False)
ax5.spines['bottom'].set_visible(False)
ax5.spines['left'].set_visible(False)
ax5.set_ylim((-1.,1.2))
text=ax5.text(0.5,0.6,st[4].stats.station,fontsize=12, weight='bold',color=cmap(0.8))
text.set_path_effects([matplotlib.patheffects.Stroke(linewidth=1.5, foreground=(0.8,0.8,0.8)),
                   matplotlib.patheffects.Normal()])

plt.savefig('f2pBC.pdf', dpi=300)