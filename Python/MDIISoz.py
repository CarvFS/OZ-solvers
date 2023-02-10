"""
Using MDIIS to solve OZ equation for a simple LJ fluid
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

    #print(data)

    # generating coordinates r and k

    N = 1000 # number of grid points
    rmax = 15.0 # maximum value of r
    dr = rmax/N # increment of r
    dq = pi/rmax # increment of q

    r=np.array([])
    q=np.array([])

    for i in range(N-1):
        value_r = np.array([(i+1)*dr])
        value_q = np.array([(i+1)*dq])
        r = np.append(r,value_r,axis = 0)
        q = np.append(q,value_q,axis = 0)

    r = np.reshape(r,[len(r),1])
    q = np.reshape(q,[len(r),1])

    # defining density and temperature
    rho = 0.8
    T = 2.0

    # defining initial guess for c(r)
    
    V = pot(r,T)

    # defining parameters for MDIIS method
    m = 1
    n = 1
    eta = 0.5

    sto_h = np.empty((len(r),1))
    sto_h2 = np.empty((len(r),1))
    R = np.empty((len(r),1))
    s = np.empty((m,m))
    c = np.zeros((len(r),1))


    check = 1
    while n <= 200:
        # calculating Fourier Tranform (FT) of c(r)
        ch = FT(r,q,c[:,[m-1]])

        # calculating FT of h(r)
        hh = ch/(1.0 - rho*ch)

        # calculating first h(r)
        h = invFT(r,q,hh)
        
        sto_h[:,[m-1]] = h
        
        #print(sto_h)
        #input("press any bottom to continue")
        
        
        # calcuçating second h(r) from closure relation
        h2 = np.exp(-V + h - c[:,[m-1]]) - 1.0
        
        sto_h2[:,[m-1]] = h2
        
        #print(sto_h2)
        #input("press any bottom to continue")
        

        # calculating res
        R[:,[m-1]] = h2 - h
    
        if np.sqrt(np.dot(R[:,m-1],R[:,m-1])) < 1e-5:
            break
        print("it ", n, "error = ", np.sqrt(np.dot(R[:,m-1],R[:,m-1])))
        #input("press any key")

        # calculating dot products
        for i in range(m):
            for j in range(m):
                s[i][j] =   np.dot(R[:,i],R[:,j])

        #print(s)
        if check < 2:
            s = np.append(s,-np.ones((m,1)),axis=1)
            s = np.append(s,-np.ones((1,m+1)),axis=0)
            s[m][m] = 0
        
        

        b = np.zeros((m+1,1))
        b[m] = -1
        #print(s)
        #input("press any key")
        

        cs = np.linalg.solve(s, b)

        #print(s.shape, R.shape, c.shape, cs[0:m].shape)

        #print(cs[0:m])
        
        #print(R.shape,cs[0:m].shape)
        #input()
        
        r_star = np.matmul(R,cs[0:m])
            
        #print(r_star.shape)
        c_new = np.matmul(c,cs[0:m]) + eta*r_star
        
        #plt.plot(r,h,"--r")
        #plt.pause(.01)
        
        m = m + 1
        
        if m > 10:
            check = 2
            m=10
            for i in range(m-1):
                c[:,i] = c[:,i+1]
            c[:,[m-1]] = c_new
            for i in range(m-1):
                R[:,i] = R[:,i+1]
        else:
            R = np.append(R,np.ones((len(r),1)),axis=1)
            c = np.append(c,c_new,axis=1)

        
        sto_h = np.append(sto_h,np.ones((len(r),1)),axis=1)
        sto_h2 = np.append(sto_h2,np.ones((len(r),1)),axis=1)

        n = n + 1

        

    plt.plot(data["r"],data["hr"],"-k",r,h,"--r")  
    plt.show()  
        

def pot(r,T):
    Ur=(4.0/T)*((1.0/r)**12 - (1.0/r)**6)
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


main()
