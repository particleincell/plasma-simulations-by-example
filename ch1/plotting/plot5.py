import pylab as pl
import numpy as np
import csv

pl.rcdefaults()
pl.rc('text', usetex=False)
#pl.rc('font', family='sans-serif')

with open('trace.csv', 'r') as f:
  reader = csv.reader(f)
  entries = list(reader)

#convert data to float, skip first line which contains headers
data = [[float(c) for c in r] for r in entries[1:]]

#transpose into a numpy array 
data = np.transpose(np.array(data))


with open('trace2.csv', 'r') as f:
  reader = csv.reader(f)
  entries = list(reader)

#convert data to float, skip first line which contains headers
data2 = [[float(c) for c in r] for r in entries[1:]]

#transpose into a numpy array 
data2 = np.transpose(np.array(data2))

#theory
t = data[0]
E = -100
QE = 1.602176565e-19		# C, electron charge
EPS_0 = 8.85418782e-12  # F/m, vacuum permittivity
ME = 9.10938215e-31	
v_th = -QE/ME*E*t
x_th = -0.5*QE/ME*E*t**2



x1 = np.array(data[1])
v1 = np.array(data[2])
x2 = np.array(data[3])
v2 = np.array(data[4])
x3 = np.array(data[5])
v3 = np.array(data[6])

x1s = np.array(data2[1])
v1s = np.array(data2[2])
x2s = np.array(data2[3])
v2s = np.array(data2[4])


#create a new figure
pl.close('all')
fig, ax1 = pl.subplots(figsize=(7,4))
#ax2 = ax1.twinx()

#plot column 0 (x) vs. column 1 (phi)
#ln1 = ax1.plot(x2,t,lineWidth=1,lineStyle="-",color='888',label="forward")
#ln2 = ax1.plot(x,t,lineWidth=2,color='black',label="central")
#ln3 = ax1.plot(x3,t,lineWidth=1,lineStyle="-.",color='444',label="backward")
ax1.plot(x_th,v_th,'o',markerFaceColor='1',lineWidth=1,color='0',label="theory")
ax1.plot(x1,v1,lineWidth=2,color='0.3',label="forward")
ax1.plot(x2,v2,lineWidth=2,color='0.7',label="backward")
ax1.plot(x3,v3,lineWidth=2,color='0',label="central")

ax1.plot(x1s,v1s,'--',lineWidth=2,color='0.3',label="forward, half dt")
ax1.plot(x2s,v2s,'--',lineWidth=2,color='0.7',label="backward, half dt")

#ax1.plot(t,x3,lineWidth=1,color='blue',label="pos")

#ln2 = ax2.plot(v2,t,'--',lineWidth=2,color='gray',label="vel")



#pl.xticks(rotation=45)

ax1.set_ylabel('veloctiy (m/s)')
ax1.set_xlabel('position (m)')
#ax2.set_xlabel('vel (m/s)', color='gray')
#pl.xticks(np.arange(t[0],t[-1],(t[-1]-t[0])/2) )
#lns = ln1+ln2+ln3
#labs = [l.get_label() for l in lns]
ax1.legend(loc='lower right')
pl.tight_layout()

ax1.grid(axis='both')
pl.show()


