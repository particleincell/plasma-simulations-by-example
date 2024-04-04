# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 10:21:42 2018

@author: lbrieda
"""

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab as pl
from matplotlib import cm
from matplotlib.colors import ListedColormap

vth = 1
ni = 21
nj = 21
f = np.zeros((ni,nj))

u,v = np.meshgrid(np.linspace(-5*vth,5*vth,ni),np.linspace(-5*vth,5*vth,nj))
f = np.exp(-((u-1)**2+v**2)/vth**2)

fig = pl.figure(figsize=(8,6))
ax = fig.gca(projection='3d')
ax._axis3don = False

#transform for plotting
f = np.power(f,1/4);

# transition to transparency on first 10 entries
cmap = pl.cm.binary
my_cmap = cmap(np.arange(cmap.N))
print(cmap.N)
my_cmap[0:10,-1] = np.linspace(0,1,10);
my_cmap = ListedColormap(my_cmap)

surf = ax.plot_surface(u, v, f-0.01, cmap=my_cmap,
                       linewidth=0, antialiased=True)
z = np.zeros_like(u);
ax.plot_wireframe(u, v, z, color='888')
ax.plot(u[10,:],v[10,:],zs=f[10,:],color='k',LineWidth=4)
ax.plot(u[:,10],v[:,10],zs=f[:,10],color='k',LineWidth=4)
ax.text(0,-6*vth,0,"u",color='222',size='20');
ax.text(-6*vth,0,0,"v",color='222',size='20');
ax.set_zlim3d(0, 1.2) 
pl.savefig('vdf.png')

