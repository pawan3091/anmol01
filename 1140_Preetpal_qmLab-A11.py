import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate
from scipy import optimize,stats
from scipy.optimize import fsolve
import pandas as pd
from scipy.integrate import solve_ivp
import array

def a_k(xi,xf,N,B):
    array=[];psi=[]
    x=np.linspace(xi,xf,N)
    k=np.linspace(-10,10,N)
    for i in x:
        if abs(i)<=B:
            psi.append(1/np.sqrt(2*B))
        else:
            psi.append(0)
    for i in k:
        array_=(1/np.sqrt(2*np.pi))*integrate.simps(psi*np.exp(complex(0,-1)*x*i),x)
        array.append(array_)
    u_norm=array/np.sqrt(integrate.simps(np.power(array,2),x))
    plt.scatter(x,**2)
    plt.xlabel("x")
    plt.ylabel("$Ψ(x)^2$")
    plt.grid()
    plt.title("Probability Density for Momentum of Particle at t=0")
    plt.show()
    return u_norm, x

def wv_fn(xi,xf,t,N,B):
    k=np.linspace(-10,10,N)
    u_norm,x = a_k(-3,3,N,1)
    for j in t:
        wv_fn=[]        
        for i in x:
            wv_fn_=(1/np.sqrt(2*np.pi))*integrate.simps(u_norm*np.exp(complex(0,-1)*(i*k-(k**2)*j)),k)
            wv_fn.append(wv_fn_)
        wv_norm=wv_fn/np.sqrt(integrate.simps(np.power(wv_fn,2),k))
        plt.plot(k,wv_norm**2,label=f't={np.round(j,2)}')
        plt.grid()
        plt.xlabel("a(k)")
        plt.ylabel("$Ψ(a(k))^2$")
        plt.title("Probability Density at different τ for |x|<b/2")
        plt.legend()
    plt.show()

'''-----------------------------------------Q(a)---------------------------------------------------'''
# xi=-3;xf=3;N=1000;B=1
# t=np.arange(0,0.1,0.1)
# wv_fn(xi,xf,t,N,B/2)

'''-------------------------------------------Q(b)-----------------------------------------------'''
# xi=-3;xf=3;N=1000;B=1
# t=np.arange(0,2.1,0.1)
# wv_fn(xi,xf,t,N,B/2)

def a_k_gaussian(xi,xf,N,A,a):
    array=[];psi=[]
    x=np.linspace(xi,xf,N)
    k=np.linspace(-10,10,N)
    for i in x:
        psi.append(np.power(2*a/np.pi,1/4))
    for i in k:
        array_=(1/np.sqrt(2*np.pi))*integrate.simps(psi*np.exp(-a*(x**2))*np.exp(complex(0,-1)*i*x),x)
        array.append(array_)
    u_norm=array/np.sqrt(integrate.simps(np.power(array,2),x))
    # plt.scatter(x,u_norm**2)
    # plt.xlabel("x")
    # plt.ylabel("$Ψ(x)^2$")
    # plt.grid()
    # plt.title("Probability Density for Momentum of Particle at t=0")
    # plt.show()
    return u_norm, x

def wv_fn_gaussian(xi,xf,t,N,A,a):
    k=np.linspace(-10,10,N)
    u_norm,x = a_k_gaussian(xi,xf,N,A,a)
    for j in t:
        wv_fn=[]        
        for i in x:
            wv_fn_=(1/np.sqrt(2*np.pi))*integrate.simps(u_norm*np.exp(complex(0,-1)*(i*k-(k**2)*j)),k)
            wv_fn.append(wv_fn_)
        wv_norm=wv_fn/np.sqrt(integrate.simps(np.power(wv_fn,2),k))
        plt.plot(k,wv_norm**2,label=f't={np.round(j,2)}')
        plt.grid()
        plt.xlabel("a(k)")
        plt.ylabel("$Ψ(a(k))^2$")
        plt.title("Probability Density of Gaussian Wave Packet at different τ ")
        plt.legend()
    plt.show()
xi=-3;xf=3;N=1000;B=1;a=1
t=np.arange(0,2.1,0.1)
wv_fn_gaussian(xi,xf,t,N,B,a)