import numpy as np
import sys
import os
if not os.path.isfile(sys.argv[1]):
    print "# "+sys.argv[1]+" doesn't exist!"
    sys.exit()
t,u,p,o1,o2 = np.loadtxt(sys.argv[1],float,unpack=True)
o = o1**2+o2**2
ndat = len(t)
nskip = ndat/10
avgs = [ s[nskip:].mean() for s in [u,p,o] ]
print t[-1], "{0} {1} {2}".format(*avgs)
