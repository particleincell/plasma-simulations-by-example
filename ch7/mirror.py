#simulation of magnetic mirror

import numpy as np
import pylab as pl


def evalB(z0):
    z = z0/1e-5
    B = np.exp(-2*(z-1)**2) + np.exp(-2*(z+1)**2)
    return B

def evalGradB(z0):
    z = z0/1e-5
    gB = np.exp(-2*(z-1)**2)*(-4*(z-1)) + np.exp(-2*(z+1)**2)*(-4*(z+1))
    return gB*1e5  #because of unit conversion


vz=1e1

q = -1.602e-19
m = 9.19e-31	#electron mass

#electric field
E = np.asarray([0,0,0])

#initial position and velocity
max_it=20000
x=np.zeros((max_it,3))
v=np.zeros((max_it,3))

Bmag = np.zeros(max_it)
x[0,:] = [0,0,0]
v[0,:] = [1e4,0,1.6e4]
dt0 = 1e-12

#compute some values
Rm = (1+np.exp(-8))/(2*np.exp(-2))	#Bm/B0
theta = np.arcsin(np.sqrt(1/Rm))*180/np.pi
vperp2 = v[0,0]**2+v[0,1]**2
vpar2 = v[0,2]**2
print("theta: %.3g,   1/Rm: %.3g, vp2/v02: %.3g"%(theta,1/Rm, vperp2/(vperp2+vpar2)))

B0 = evalB(0);
vtot2 = vperp2+vpar2
print("Particle mirrors at Bm: %.3g, vtot:%.3g"%(B0*vtot2/vperp2, np.sqrt(vtot2)))
#magnetic moment, should be conserved
mu0 = vperp2/evalB(x[0,2])
mu = 0.5*m*mu0
Bmag[0] = vperp2/mu0	

for it in range (max_it-1):

	#evaluate B at particle position
	z = x[it,2]
	B = np.asarray([0,0,evalB(z)])
	gradB = np.asarray([0,0,evalGradB(z)])
	
	#force
	F=q*E - mu*gradB
    	
     #velocity rewind for leapfrog
	if (it==0):
         dt = 0.5*dt0
	else:
         dt = dt0
         
	#v minus
	v_minus = v[it,:] + (F/m)*0.5*dt
	
	#v prime
	t = q/m*B*0.5*dt	#t vector
	v_prime = v_minus + np.cross(v_minus,t)
	
	#v plus
	s = 2*t/(1+np.dot(t,t))#    #s vector
	v_plus = v_minus + np.cross(v_prime,s)
	
	#v n+1/2
	v[it+1,:] = v_plus + (F/m)*0.5*dt
	
	#position update
	x[it+1,:] = x[it,:] + v[it+1,:]*dt0
	
	#save bmag for plotting
	Bmag[it+1] = B[2]
	
	#screen output
	if (np.mod(it,500)==0):
		#perpendicular and tangential velocities
		vperp2 = v[it,0]**2 + v[it,1]**2
		vpar2 = v[it,2]**2
		
		print("it:%d, z:%.2e, v_par:%.2e, v_perp:%.2e, v_tot:%.4e"%
		(it,x[it,2],np.sqrt(vperp2),np.sqrt(vpar2),np.sqrt(vperp2+vpar2)))
	
pl.figure(1)
pl.plot (x[:,0],x[:,1])
pl.title("x-y view")
pl.axis('equal')

pl.figure(2)
pl.plot (np.arange(0,max_it),x[:,2])
pl.title("it vs. z")

pl.figure(3)
pl.plot (np.arange(0,max_it),(v[:,0]**2 + v[:,1]**2)/B[2])
pl.title("it vs. mu0")

pl.figure(4)
pl.plot(x[:,2],Bmag)
pl.title("z vs Bmag")

