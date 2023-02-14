"""
Using MDIIS to solve OZ equation for a simple LJ fluid
Author: Felipe Silva Carvalho
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

pi = math.pi

def main():
    global dr,dq,R

    # loading expected result acquired using Picard iteration method
    data = pd.read_csv("hrHNC.txt")

    # generating coordinates r and q

    N = 1000 # number of grid points
    rmax = 15.0 # maximum value of r
    dr = rmax/N # increment of r
    dq = pi/rmax # increment of q

    r=np.array([])
    q=np.array([])

    for i in range(N-1):
        r = np.append(r,[(i+1)*dr],axis = 0)
        q = np.append(q,[(i+1)*dq],axis = 0)

    r = np.reshape(r,[len(r),1])
    q = np.reshape(q,[len(r),1])

    # defining density and temperature - set up to match results contained in hrHNC.txt
    rho = 0.8
    T = 2.0

    # defining initial guess for c(r)
    
    V = pot(r,T)

    # defining parameters for MDIIS method
    m = 1 # related to the matrix size (m x m)
    n = 1 # number of iterations
    eta = 0.5 # damping parameter


    # creating vectors and matrices to store information
    sto_h = np.empty((len(r),1)) # store h(r) calculated from OZ equation
    sto_h2 = np.empty((len(r),1)) # store h(r) calculated from closure relation
    R = np.empty((len(r),1)) # store residual
    s = np.empty((m,m)) # store S matrix elements
    c = np.zeros((len(r),1)) # store c(r) values

    # initiating iterative procedure
    check = 1 # parameter to help fixing S matrix maximum size
    while n <= 200:
        # calculating Fourier Tranform (FT) of c(r)
        ch = FT(r,q,c[:,[m-1]])

        # calculating FT of h(r)
        hh = ch/(1.0 - rho*ch)

        # calculating first h(r)
        h = invFT(r,q,hh)
        
        sto_h[:,[m-1]] = h # store h(r)
        
        # calculating second h(r) from closure relation
        h2 = np.exp(-V + h - c[:,[m-1]]) - 1.0
        
        sto_h2[:,[m-1]] = h2 # store h(r)
        
        # calculating residual
        R[:,[m-1]] = h2 - h

        # setting stop criterium
        if np.sqrt(np.dot(R[:,m-1],R[:,m-1])) < 1e-5:
            break
        print("it ", n, "error = ", np.sqrt(np.dot(R[:,m-1],R[:,m-1]))) # print iteration and residual norm

        # calculating dot products and constructing S matrix
        for i in range(m):
            for j in range(m):
                s[i][j] =   np.dot(R[:,i],R[:,j])

        # adding -1's and 0 to S matrix while not reaching maximum dimensions
        if check < 2:
            s = np.append(s,-np.ones((m,1)),axis=1)
            s = np.append(s,-np.ones((1,m+1)),axis=0)
            s[m][m] = 0
        
        # building vector b
        b = np.zeros((m+1,1))
        b[m] = -1

        # acquiring coefficients
        cs = np.linalg.solve(s, b)
        
        # calculatiing R*
        r_star = np.matmul(R,cs[0:m])
            
        # calculating new value for c(r)
        c_new = np.matmul(c,cs[0:m]) + eta*r_star
        
        # increment S matrix dimension
        m = m + 1
        
        # check if dimension reached the maximum value disired, redefining check value and discarding old values
        if m > 10:
            check = 2
            m=10
            for i in range(m-1):
                c[:,i] = c[:,i+1]
            c[:,[m-1]] = c_new
            for i in range(m-1):
                R[:,i] = R[:,i+1]
        else: # if did not reach maximum value, increment R and c matrices
            R = np.append(R,np.ones((len(r),1)),axis=1)
            c = np.append(c,c_new,axis=1)

        # append new values calculated
        sto_h = np.append(sto_h,np.ones((len(r),1)),axis=1)
        sto_h2 = np.append(sto_h2,np.ones((len(r),1)),axis=1)

        # increment iteration value
        n = n + 1

        
    #plot final result along with the one acquired prevously
    plt.plot(data["r"],data["hr"],"-k",r,h,"--r")  
    plt.xlim([0, 6])
    plt.ylim([-1.1, 1.6])
    plt.show()  
        
##################### Defining functions for calculating potential, direct and inverse Fourier transforms

def pot(r,T):
    Ur=(4.0/T)*((1.0/r)**12 - (1.0/r)**6) # reduced coordinates
    return Ur

def FT(r,q,f):
    func = []
    for k in q:
        kernel = (4.0*pi/k)*np.sin(k*r)*dr*r*f
        integral=np.sum(kernel,0)
        func.append(integral)
    return np.array(func)

def invFT(r,q,f):
    func = []
    for k in r:
        kernel = (1.0/(k*2*pi*pi))*dq*np.sin(q*k)*q*f
        integral=np.sum(kernel,0)
        func.append(integral)
    return np.array(func)

# call main function to run code
main()
