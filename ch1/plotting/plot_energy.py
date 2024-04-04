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

t = data[0]*1e6
x = np.array(data[1])
v = np.array(data[2])
ke = np.array(data[3])
pe = np.array(data[4])

#create a new figure
pl.close('all')
fig, ax1 = pl.subplots(figsize=(7,4))
#ax2 = ax1.twinx()

#plot column 0 (x) vs. column 1 (phi)
#ln1 = ax1.plot(x2,t,lineWidth=1,lineStyle="-",color='888',label="forward")
#ln2 = ax1.plot(x,t,lineWidth=2,color='black',label="central")
#ln3 = ax1.plot(x3,t,lineWidth=1,lineStyle="-.",color='444',label="backward")
ax1.plot(t,pe,':',lineWidth=2,color='0.3',label="potential energy")
ax1.plot(t,ke,'--',lineWidth=2,color='0',label="kinetic energy")
ax1.plot(t,ke+pe,'-',lineWidth=2,color='0.6',label="total energy")


#ax1.plot(t,x3,lineWidth=1,color='blue',label="pos")

#ln2 = ax2.plot(v2,t,'--',lineWidth=2,color='gray',label="vel")



#pl.xticks(rotation=45)

ax1.set_xlabel('time (microsecond)')
ax1.set_ylabel('energy (eV)')
#ax2.set_xlabel('vel (m/s)', color='gray')
#pl.xticks(np.arange(t[0],t[-1],(t[-1]-t[0])/2) )
#lns = ln1+ln2+ln3
#labs = [l.get_label() for l in lns]
ax1.legend(loc='lower right')
pl.tight_layout()

ax1.grid(axis='both')
pl.show()


