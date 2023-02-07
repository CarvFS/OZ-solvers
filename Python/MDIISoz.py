import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

pi = math.pi

def main():
    global dr,dq

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
    c=np.array([])
    for i in range(N-1):
        c_ini = np.array([0.0])
        #c_ini = np.exp(-r[i])
        c = np.append(c,c_ini)

    c=np.reshape(c,[len(c),1])
    
    V = pot(r,T)

    # defining parameters for MDIIS method
    m = 1
    n = 1
    eta = 0.005

    sto_h = np.empty((len(r),1))
    sto_h2 = np.empty((len(r),1))
    R = np.empty((len(r),1))
    s = np.empty((m,m))

    while n <= 1:
        # calculating Fourier Tranform (FT) of c(r)
        ch = FT(r,q,c)

        # calculating FT of h(r)
        hh = ch/(1.0 - rho*ch)

        # calculating first h(r)
        h = invFT(r,q,hh)
        sto_h = np.append(sto_h,h,axis=1) ## TENTAR MODIFICAR PRA NÃO TER UMA COLUNA VAZIA... TALVEZ USAR RESIZE...

        #print(sto_h[:,[m]])
        
        # calcuçating second h(r) from closure relation
        h2 = np.exp(-V + h - c) - 1.0
        sto_h2 = np.append(sto_h2,h2,axis=1) ## TENTAR MODIFICAR PRA NÃO TER UMA COLUNA VAZIA... TALVEZ USAR RESIZE...

        #print(sto_h2.shape)

        # calculating res
        R = np.append(R,sto_h2[:,[m]] - sto_h[:,[m]],axis=1) ## TENTAR MODIFICAR PRA NÃO TER UMA COLUNA VAZIA... TALVEZ USAR RESIZE...
        #print(R)

        # calculating dot products
        for i in range(m):
            for j in range(m):
                s[i][j] =   np.dot(R[:,i+1],R[:,j+1])

        s = np.resize(s,(m+1,m+1))

        for i in range(m):
            s[i][m] = -1.0
            s[m][i] = -1.0
        
        s[m][m] = 0

        b = np.zeros((m+1,1))
        b[m] = -1

        #print(s)
        #print(b)

        cs = np.linalg.solve(s, b)

        #print(cs[0:m])
        #################################### CONTINUAR DAQUI!
        r_star =np.matmul(c,cs[0:m]) + eta* 
        new_c = np.matmul(c,cs[0:m]) + eta*

        print(new_c)

        n = n + 1
        m = m + 1


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
