from math import pi
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, interpolate

x0 = 0; sigma0 = 0.1; dx = 0.005; k0=0;dt = 0.001
L = 2; x = np.arange(-L,L,dx); N = len(x); V0 = 10**6;
j1=complex(0,1)

def Gaussian(x,x0,sigma0,k0):
    L = 2;dx = 0.005; x = np.arange(-L,L,dx); N = len(x); j1=complex(0,1)
    RePsi0 = np.zeros(N)
    ImPsi0 = np.zeros(N)
    Psi0 = np.zeros(N)
    RePsi0 = (2*pi*sigma0**2)**(-0.25)*np.exp(-(x-x0)**2/(4*sigma0**2))*np.cos(k0*x)
    ImPsi0 = (2*pi*sigma0**2)**(-0.25)*np.exp(-(x-x0)**2/(4*sigma0**2))*np.sin(k0*x)
    Psi0 = RePsi0+j1*ImPsi0;
    plt.subplot(1,2,1)
    plt.plot(x,RePsi0);
    plt.subplot(1,2,2)
    plt.plot(x,ImPsi0)
    return Psi0

#Psi0 = Gaussian(x,x0,sigma0,k0)

def Potential(V0,x,L):
    N = len(x)
    V = np.zeros(N)
#    V[0] = V0; V[N-1] = V0
    plt.subplot(1,1,1)
    plt.plot(x,V)
    return V
#V = Potential(V0,x,N)
#
def Hmatrix(V,dx,dt,N):
    f = V*dt/2;
    g = dt/(4*dx**2);
    H = np.diag(f+2*g)+np.diag(-g*np.ones(N-1),-1)+np.diag(-g*np.ones(N-1),1);
    return H
#H = Hmatrix(V,dx,dt,N)
#
# 
def deltax(x,Psi):
    f = abs(np.conj(Psi[0:N])*Psi[0:N]);
    meanx = integrate.simpson(x,f*x);
#    meanx = interpolate.interp1d(x, f*x, kind='cubic')
#    print(meanx)
    meanxsqr = abs(integrate.simpson(x,f*x**2));
#    meanxsqr = interpolate.interp1d(x, f*x**2, kind='cubic')
#    print(meanxsqr)
    w = round(np.sqrt(meanxsqr - meanx**2),5);
    return  meanx, w
               
def CrankNicholsonTDSE(k0,n):
    x0 = 0; sigma0 = 0.1; dx = 0.005; #k0=0;
    L = 2; x = np.arange(-L,L,dx); N = len(x); V0 = 10**6;
    j1=complex(0,1)
    psi0 = Gaussian(x,x0,sigma0,k0);
    V = Potential(V0,x,L);
    dt = 0.001; 
    t = dt;
    I = np.eye(N);
    H = Hmatrix(V,dx,dt,N);
#    Psi = []#
    Psi = np.zeros(N,dtype='complex_')
    i = 1
#    Psi.append(psi0)#
    Psi[0:N]= psi0[0:N];
#    Psi[0:N-1][i] = psi0[0:N-1][i-1];
    #print(shape(Psi))
    i = 2; 
    p = 1;
    j = 1; 
    M = np.linalg.inv(I + j1*H) @ (I - j1*H);
    #print(shape(M))
    while j <= n:
#        Psi.append(M@Psi[j-1])
        Psi[0:N] = M@Psi[0:N]#M.dot(Psi[0:N-1][i]);
        if j%20 == 1:
            plt.subplot(3,3,p)
            plt.plot(x[0:N],np.conj(Psi[0:N])*Psi[0:N])#,rect=[-L, 0, L ,1.2])
#            plt.plot(x,np.conj(Psi[j-1])*Psi[j-1])
            plt.ylim(0,4)
            plt.title('For j={}'.format(j))
#            f = abs(np.conj(Psi[0:N])*Psi[0:N]);
            meanx,delx = deltax(x,Psi[0:N]);
            print("<x> =",meanx,"width =", delx)
            p = p + 1;
        i = i+1;  
        t = t+dt; 
        j = j+1;
    return Psi
    
Psi = CrankNicholsonTDSE(0,180)
