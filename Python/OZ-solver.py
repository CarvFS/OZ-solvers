"""
Program which solves the Ornstein-Zernike (OZ) equation using the Picard iteration method 
Author: Felipe Silva Carvalho
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

global dk,dr
pi=math.pi

########################## Defining functions ##########################
# Fourier transform
def FT(r,q,f):
    global dk,dr
    Ch=[]
    for k in q:
        kernel=(4.0*pi/k)*np.sin(k*r)*f*r*dr
        integral=np.sum(kernel,axis=0)
        Ch.append(integral)
    Ch=np.array(Ch)
    return Ch

# inverse Fourier transform
def FTinv(r,q,f):
    global dk,dr
    C2=[]
    for k in r:
        kernel=1.0/(k*2.0*pi**2.0)*np.sin(k*q)*f*q*dk
        integral=np.sum(kernel,axis=0)
        C2.append(integral)
    C2=np.array(C2)
    return C2

# defining energy potential function
def pot(r,T):
    Ur=(4.0/T)*((1.0/r) **12 - (1.0/r)**6)
    return Ur

# defining bridge function
def bridge(r):
    ## Bridge function which gives HNC approximation
    Br = np.zeros(np.shape(r))
    return Br

# Picard iteration process to solve OZ equation
def oz(r,q,rho,fgamma,V):
    # Calculating C given a bridge function
    C=np.exp(fgamma - V + bridge(r)) - fgamma - 1.0
    # Calculating Fourier transform of C, C_hat
    C_hat=FT(r,q,C)
    # Acquring Fourier transform of gamma function, g_hat, from OZ equation
    g_hat=(rho*np.square(C_hat))/(1.0-rho*C_hat)
    # Calculating inverse Fourier transform of g_hat to acquire new value for gamma function
    g_novo=FTinv(r,q,g_hat)
    # defining mixing parameter (0 < alfa <= 1)
    alfa=1.0
    g_novo = alfa*g_novo + (1-alfa)*fgamma
    return g_novo

def excess_chemical_potential(r, g, C, rho, T):
    int_kernel = (0.5*(g-1)*(g-1) - C - 0.5*(g-1)*C) * r*r
    I = np.trapz(int_kernel, r, axis = 0)
    print(r.shape, g.shape, C.shape)
    # mu = T * rho * 4 * pi * I # This gives mu/epsilon
    mu = rho * 4 * pi * I # This gives beta * mu
    return mu

def d_pot(r,T):
    dUdr=(4.0/T)*(-12.0*(1.0/r)**13 + 6.0*(1.0/r)**7)
    return dUdr

def pressure_v(r, rho, T, g):
    int_kernel = r*r*r*g*d_pot(r,T)
    I = np.trapz(int_kernel, r, axis = 0)
    p_v = 1 - (rho/6) * 4 * pi * I # This gives beta*p_v/rho
    return p_v

def isotherm_comp(r, rho, g):
    int_kernel = r*r * (g - 1)
    I = np.trapz(int_kernel, r, axis = 0)
    chi_T = 1 + rho * 4 * pi * I # This gives rho*chi_T/beta
    return chi_T

########################## All functions defined ##########################


Rmax=15.0 # maximum value for variable r
N=800.0 # number of points to be generated

# generating r and q coordinates

imax=N-1.0
dr=Rmax/N
dk=pi/Rmax

r=[]
for i in range(int(imax)):
    valor=(i+1)*dr
    r.append(valor)
r=np.reshape(r, (len(r),1))

q=[]
for i in range(int(imax)):
    valor=(i+1)*dk
    q.append(valor)
q=np.reshape(q, (len(q),1))

rho=0.9 # define liquid density
T=1.5 # define temperature

# calculating potential
V=pot(r,T)

# initial estimative for gamma function
fgammai=np.zeros(np.shape(r))

tol=10.0
count=0
# loop for solve OZ with the desired tolerance
while tol > 1e-5:
    fgamma=oz(r,q,rho,fgammai,V)
    tol=np.sqrt(np.sum(np.square(fgammai-fgamma)))
    fgammai=fgamma
    count=count+1
    print(tol)

data=pd.read_csv('gammaHNCtest.txt') #Reading expected result for testing accuracy

# Plotting results

plot1=plt.figure(1)
plt.plot(r,fgammai,'-k',data['r'],data['gr'],'--r')
plt.legend(('Result','Expected'))
plt.xlim([0, 5])
plt.ylim([-2, 21])
plt.xlabel('r/$\sigma$',fontsize=16)
plt.ylabel('$\gamma$(r)',fontsize=16)
plt.tick_params(labelsize=16)
plt.tight_layout()
#plt.savefig('gamma.eps', format='eps', dpi=600)

# Calculating g(r) and C(r)
Cr = np.exp(fgammai - V + bridge(r)) - fgammai - 1.0
gr = fgammai + Cr + 1.0

plot2=plt.figure(2)
plt.plot(r,Cr,'-k')
plt.xlim([0, 5])
plt.ylim([-22, 2])
plt.xlabel('r/$\sigma$',fontsize=16)
plt.ylabel('C(r)',fontsize=16)
plt.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig('Cr.png')
plt.savefig('Cr.pdf')
plt.show()

plot3=plt.figure(3)
plt.plot(r,gr-1,'-k')
plt.xlim([0, 5])
plt.ylim([-0.5, 3])
plt.xlabel('r/$\sigma$',fontsize=16)
plt.ylabel('h(r)',fontsize=16)
plt.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig('hr.png')
plt.savefig('hr.pdf')
plt.show()

## calculate excess chemical potential
print(f'beta * mu = {excess_chemical_potential(r,gr,Cr,rho,T)}')
print(f'beta * p_v / rho = {pressure_v(r, rho, T, gr)}')
print(f'rho * chi_T / beta = {isotherm_comp(r, rho, gr)}')