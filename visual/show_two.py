import numpy as np
import pylab as pl
import matplotlib.patches as patches
import sys

inp = sys.argv[1]
out = sys.argv[2]
#
I, x, y = np.loadtxt(inp,float,skiprows=5,unpack=True)
I2, x2, y2 = np.loadtxt(out,float,skiprows=5,unpack=True)

with open(inp) as ff:
   ff.readline()
   ff.readline()
   ff.readline()
   nop = int( ff.readline().split()[0] )
   lx, ly = (float(v) for v in ff.readline().split() )

#
try:
   xlim, ylim = int(sys.argv[3]), int(sys.argv[4])
except IndexError:
   xlim, ylim = lx, ly
reg = np.logical_and(x<xlim,y<ylim)
reg2 = np.logical_and(x2<xlim,y2<ylim)
   
#
fig = pl.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.scatter(x[reg],y[reg],color='b',label='i')
ax.scatter(x2[reg2],y2[reg2],color='r',label='f')
ax.set_aspect('equal')

#
if '-n' in sys.argv:
   for k in range(len(I)):
    if reg[k]:
      ax.text(x[k]+0.1,y[k]+0.1,int(I[k]),color='b')
    if reg2[k]:
      ax.text(x2[k]+0.1,y2[k]-0.1,int(I2[k]),color='r')


ax.set_xlim(0,xlim)
ax.set_ylim(0,ylim)
ax.set_xticks(np.arange(0,xlim,1.0))
ax.set_yticks(np.arange(0,ylim,1.0))
ax.grid()
ax.legend()


pl.show()
