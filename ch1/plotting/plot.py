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

#create a new figure
fig, ax1 = pl.subplots(figsize=(7,4))
ax2 = ax1.twinx()

#plot column 0 (x) vs. column 1 (phi)
ax1.plot(data[0],data[1],lineWidth=2,color='black',
        marker='o',markerFaceColor='white',markerSize=8)

ax2.plot(data[0],data[3],'--',lineWidth=2,color='gray')

ax1.set_xlabel('x (m)')
ax1.set_ylabel('phi (V)')
ax2.set_ylabel('rho (C/m^3)', color='gray')
pl.show()


