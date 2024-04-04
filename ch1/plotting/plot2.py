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
phi_th = rho/(2*8.85418782e-12)*x*(x[-1]-x)

#create a new figure
fig, ax1 = pl.subplots(figsize=(7,4))
#ax2 = ax1.twinx()

#plot column 0 (x) vs. column 1 (phi)
ax1.plot(x,phi_th,lineWidth=2,color='black',label="analytical")
ax1.plot(x,data[1],lineStyle='None',color='black', marker='o',
         markerFaceColor='gray',markerSize=8,label="direct")
ax1.plot(x,data[2],lineStyle='-.',color='black', marker='v',
         markerFaceColor='white',markerSize=6,label="GS 20")
ax1.plot(x,data[3],lineStyle='-.',color='black', marker='x',
         markerFaceColor='white',markerSize=6,label="GS 60")
ax1.plot(x,data[4],lineStyle='-.',color='black', marker='p',
         markerFaceColor='white',markerSize=6,label="GS 400")

#ax2.plot(x,rho,'--',lineWidth=2,color='gray')

ax1.set_xlabel('x (m)')
ax1.set_ylabel('phi (V)')
#ax2.set_ylabel('rho (C/m^3)', color='gray')
ax1.legend();
pl.show()
pl.grid()


