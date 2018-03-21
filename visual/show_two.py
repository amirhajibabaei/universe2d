import numpy as np
import pylab as pl
import matplotlib.patches as patches
import sys

inp = sys.argv[1]
out = sys.argv[2]
#
I, x, y = np.loadtxt(inp,float,skiprows=5,unpack=True)
I2, x2, y2 = np.loadtxt(out,float,skiprows=5,unpack=True)

with open("fin.txt") as ff:
   ff.readline()
   ff.readline()
   ff.readline()
   nop = int( ff.readline().split()[0] )
   lx, ly = (float(v) for v in ff.readline().split() )

   
#
fig = pl.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.scatter(x,y,color='b',label='i')
ax.scatter(x2,y2,color='r',label='f')
ax.set_aspect('equal')

for k in range(len(I)):
      ax.text(x[k]+0.1,y[k]+0.1,int(I[k]),color='b')
      ax.text(x2[k]+0.1,y2[k]-0.1,int(I2[k]),color='r')


ax.set_xlim(0,lx)
ax.set_ylim(0,ly)
ax.set_xticks(np.arange(0,lx,1.0))
ax.set_yticks(np.arange(0,ly,1.0))
ax.grid()
ax.legend()


pl.show()
