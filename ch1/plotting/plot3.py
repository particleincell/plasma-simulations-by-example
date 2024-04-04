import pylab as pl
import numpy as np
import csv

with open('results.csv', 'r') as f:
  reader = csv.reader(f)
  entries = list(reader)

#convert data to float, skip first line which contains headers
data = [[float(c) for c in r] for r in entries[1:]]

#transpose into a numpy array 
data = np.transpose(np.array(data))

#compute analytical solution
x = data[0]
rho = data[-1]
e_th = rho/(8.85418782e-12)*(x-x[-1]/2)

#create a new figure
fig, ax1 = pl.subplots(figsize=(7,4))
#ax2 = ax1.twinx()

#plot column 0 (x) vs. column 1 (phi)
ax1.plot(x,e_th,lineWidth=2,color='black',label="analytical")
ax1.plot(x,data[2],lineStyle='-.',color='black', marker='o',
         markerFaceColor='gray',markerSize=8,label="first order")
ax1.plot(x,data[3],lineStyle='-.',color='black', marker='v',
         markerFaceColor='white',markerSize=6,label="second order")

#ax2.plot(x,rho,'--',lineWidth=2,color='gray')

ax1.set_xlabel('x (m)')
ax1.set_ylabel('electric field (V/m)')
#ax2.set_ylabel('rho (C/m^3)', color='gray')
ax1.legend();
pl.show()
pl.grid()


