# -*- coding: utf-8 -*-
"""
Demo magnetic field solver: uniformly magnetized sphere
Follows Jackson 5.10
"""

import numpy
import pylab as pl
import math

#sphere properties
a = 0.4
M0 = 1
mu0=1

#mesh properties
nr = 41
nz = 41
dr = 0.1
dz = 0.1
z0 = -0.5*(nz-1)*dz

#--- theoretical solution -----------------
phiM_th = numpy.zeros([nz,nr])

#evaluate
for i in range (nz):
    for j in range (nr):
        pos = [z0+i*dz, j*dr]
        R = math.sqrt(pos[0]**2+pos[1]**2)
        if R==0: continue
        cos_theta = pos[0]/R
        phiM_th[i][j] = (1.0/3.0)*M0*a*a*(min(R,a)/max(R,a)**2)*cos_theta

#---- numerical solution ------------------
#set fixed M in the sphere
M = numpy.zeros([nz,nr,2])
for i in range (nz):
    for j in range (nr):
        pos = [z0+i*dz, j*dr]
        R = math.sqrt(pos[0]**2+pos[1]**2)
        if (R<=a):
            M[i][j] = [M0,0]

#compute div(M)
divM = numpy.zeros([nz,nr])
for i in range (nz):
    for j in range (nr):
        jp = (j+1) if j<(nr-2) else (nr-1)
        jm = (j-1) if j>0 else 0
        rp = jp*dr
        rm = jm*dr
        r = j*dr;
        divM[i][j] = -(1/r)*((rp*M[i][jp][1]-rm*M[i][jm][1])/(rp-rm)) if r>0 else 0
        
        #z component
        ip = (i+1) if i<(nz-2) else (nz-1)
        im = (i-1) if i>0 else 0
        zp = ip*dz
        zm = im*dz
        divM[i][j] -= (M[ip][j][0]-M[im][j][0])/(zp-zm)
        
#Poisson solver
phiM = numpy.zeros_like(divM)
dz2 = dz*dz
dr2 = dr*dr

#set radia
r = numpy.zeros_like(phiM)
for i in range(nz):
    for j in range(nr):
        r[i][j] = j*dr
    
for it in range (200):
    #regular form inside        
    phiM[1:-1,1:-1] = (divM[1:-1,1:-1] + 
                   (phiM[1:-1,2:]+phiM[1:-1,:-2])/dr2 +
                   (phiM[1:-1,0:-2]-phiM[1:-1,2:])/(2*dr*r[1:-1,1:-1]) +
                   (phiM[2:,1:-1] + phiM[:-2,1:-1])/dz2) / (2/dr2 + 2/dz2)
        
    #boundaries
    phiM[0] = 0   #left 
    phiM[-1] = 0  #right
    phiM[:,-1] = 0 #top
    phiM[:,0] = phiM[:,1] #bottom
        

#compute H=-grad(phiM)
Hz = numpy.zeros_like(phiM)
Hr = numpy.zeros_like(phiM)

#central difference on internal nodes
Hz[1:-1] = (phiM[0:nz-2]-phiM[2:nz+1])/(2*dz)
Hr[:,1:-1] = (phiM[:,0:nr-2]-phiM[:,2:nr+1])/(2*dr)
    
#one sided difference on boundaries
Hz[0,:] = (phiM[0,:]-phiM[1,:])/dz
Hz[-1,:] = (phiM[-2,:]-phiM[-1,:])/dz
Hr[:,0] = (phiM[:,0]-phiM[:,1])/dr
Hr[:,-1] = (phiM[:,-2]-phiM[:,-1])/dr

#B=mu0*(H+M)
mu0 = 1
Bz = mu0*(Hz+M[:,:,0])
Br = mu0*(Hr+M[:,:,1])

#plotting
pl.close("all")
fig = pl.figure(num=None, figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')
sub = (pl.subplot(221),pl.subplot(222),pl.subplot(223),pl.subplot(224))

sub[0].contourf(numpy.transpose(phiM_th),8,alpha=.75,linewidth=1,cmap='jet',vmin=-1./4.,vmax=1./4.)
sub[0].set_aspect('equal', adjustable='box')       
sub[0].set_title("phi_M theory")    
sub[1].contourf(numpy.transpose(phiM),8,alpha=.75,linewidth=1,cmap='jet',vmin=-1./4.,vmax=1./4.)
sub[1].set_aspect('equal', adjustable='box')       
sub[1].set_title("phi_M numerical")    
Z, R = pl.meshgrid( pl.arange(0,nz,1),pl.arange(0,nr,1) )
sub[2].streamplot(Z,R,numpy.transpose(Bz),numpy.transpose(Br),density=2)
sub[2].set_aspect('equal', adjustable='box')       
sub[2].set_title("B")    
sub[3].streamplot(Z,R,numpy.transpose(Hz),numpy.transpose(Hr),density=2)
sub[3].set_aspect('equal', adjustable='box')   
sub[3].set_title("H")    
         
