import numpy as np
import pylab as pl
import matplotlib.patches as patches
# header
fo = open('hist_gxy.txt')
for k in range(11): fo.readline()
nx = int( fo.readline().split()[1] )
ny = int( fo.readline().split()[1] )
# main
x,y,p = [], [], []
for kk in range(ny):
    xraw = []
    yraw = []
    praw = []
    for k in range(nx):
        a,b,_,_,c = (float(v) for v in fo.readline().split())
        xraw += [a]
        yraw += [b]
        praw += [c]
    x += [xraw]
    y += [yraw]
    p += [praw]
x,y,p = (np.array(v) for v in [x,y,p])
# plot
fig, ax = pl.subplots()
ax.set_aspect('equal')
patch = patches.Circle((0,0), radius=min(x.max(),y.max()), transform=ax.transData)
ax.pcolor(x,y,p,clip_path=patch)
pl.show()
