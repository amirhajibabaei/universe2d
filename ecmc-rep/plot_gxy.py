import numpy as np
import pylab as pl
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
print(list(v.shape for v in (x,y,p)))
pl.axes().set_aspect('equal')
pl.pcolor(x,y,p)
pl.show()
