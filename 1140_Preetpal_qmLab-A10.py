import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate
from scipy import optimize,stats
from scipy.optimize import fsolve
import pandas as pd
from scipy.integrate import solve_ivp
import array
def numerov(x,E_min, E_max):    
    c_i =[];u=[]
    h = x[1]-x[0]
    Alpha = 2/x+((l)*(l+1))/(np.power(x,2))-1
    ddx_12 = (h**2)/12
    for i in range(0,len(x)):
        c_i_ = 1 + np.multiply(ddx_12,Alpha[i])
        c_i.append(c_i_)
    u_0 = 0
    u_1=h
    u.append(u_0);u.append(u_1)
    for i in range(2,len(x)):
        u_ = (1/c_i[i])*(((12-10*c_i[i-1])*u[i-1])-c_i[i-2]*u[i-2])
        u.append(u_)
    return u,  x

def numerov_(x,E_min, E_max):
    c_i =[];u=[]
    h = -(x[1]-x[0])
    Alpha = 2*(((-x**2)/2)+(E_min+E_max)/2)
    ddx_12 = (h**2)/12

    for i in range(0,len(x)):
        c_i_ = 1 + np.multiply(ddx_12,Alpha[i])
        c_i.append(c_i_)
    u_0 = 0
    u_1=h
    u.append(u_0);u.append(u_1)
    for i in range(2,len(x)):
        u_ = (1/c_i[i])*(((12-10*c_i[i-1])*u[i-1])-c_i[i-2]*u[i-2])
        u.append(u_)
    return u,  x


def e_range(u,E_min,E_max) :
    I = []
    E = (E_min+E_max)/2
    for i in range(len(u)):
        if (u[i-1]*u[i]) < 0:
           I.append(i)
    N_node = len(I)

    if N_node > 0:
       E_max = E
    else:
       E_min = E
    return len(I),E_min,E_max
def matching(p,p1,x,x1):
    rescale_fac = p[-1]/p1[-1]
    p2=[]
    for i in p1:
        p2.append(i*rescale_fac)
    print(type(x))
    plt.scatter(x,p)
    print(len(p))
    plt.scatter(x1,p2)
    print(len(p2))
    plt.show()    
    p = p[:-1]
    p2=p2[:-1]
    x = x[:-1]
    x1=x1[:-1]
    p_final=p+p2
    print(p_final)
    print("final p",len(p_final))
    t=np.append(x,x1)
    x = t.tolist()
    print(x)
    print(p_final)
    print("x",len(x))
    plt.scatter(x,p_final,label="Combined")
    plt.legend()
    plt.show()
    return p,p2

def E(E_min,E_max,tol):
    for i in range(1000):        
        p=numerov(x,E_min, E_max  )
        p1=numerov_(x1,E_min, E_max)
        U_, U_back_ = matching(p,p1,x,x1)
        U = U_[:-1]
        U_back = U_back_[::-1]
        for i in U_back:
            U.append(i)
        u_norm=U/np.sqrt(integrate.simps(np.power(U,2),np.linspace(xi_1,xf_2,len(U))))
        plt.scatter(np.linspace(xi_1,xf_2,len(U)),u_norm,label="Final with E")
        plt.legend()
        plt.show()
        I ,E_min_new,E_max_new = e_range(u_norm,E_min,E_max)
        if abs(E_max_new - E_min_new)<tol:
            break
        else:
           E_min = E_min_new
           E_max = E_max_new
    return E_min_new,E_max_new,U,u_norm

'''-----------------------------------------------------------------------------------------'''
def cl_trn_pts(xi_1,xf_1,N,n):
    x = np.linspace(xi_1,xf_1,N+1)
    index=0;xf_=0
    for i in x:
        index+=1
        if round(i,2) == n**2:
            xf_ = i
            break
        x_for,x_back=x[:index+1],x[index:]
    return xf_,index,x_for,x_back[::-1]

# def phi(E):
#     E_min_new,E_max_new,U,parity=E(n_node,E_min,E_max,tol)
#     array = parity(n_node,E,E_max,tol)
#     c_i=[]
#     h = x[1]-x[0]
#     Alpha = 2*(((-x**2)/2)+E)
#     ddx_12 = (h**2)/12
#     for i in range(0,len(x)):
#         c_i_ = 1 + np.multiply(ddx_12,Alpha[i])
#         c_i.append(c_i_)
#     t= int(len(u_norm)/2)
#     p=len(x)
#     G = (1/h)*(u_norm[t+p+1]+u_norm[t+p-1]-((12*c_i[-1])-10)*u_norm[t+p])
#     return abs(G)



xi_1=0.1;xf_2=20;N=200;n=1;l=0;E_min=0;E_max=1/n**2;tol=0.4
xf_,index,x,x1 = cl_trn_pts(xi_1,xf_2,N,n)
# print("Classical Turning Point",xf_)
# print("X forward",x)
# print("xbackward",x1)
# xf_1=xf_;xi_2=xf_1_;xf_2=xf_
p=numerov(x,E_min, E_max)[0]
p1=numerov(x1,E_min,E_max)[0]
#p1,p2=matching(p,p1,x,x1)    
E(E_min,E_max,tol)

# num_eig_val=[];n=[];anal_eig_val=[]
# #E_min_new,E_max_new,U=E(n_node,E_min,E_max,tol)
# p,x=numerov(x,E_min,E_max)
# plt.plot(x,p)
# plt.show()
#     num_eig_val.append(num_eig_val_)
#     n.append(i)
#     anal_eig_val.append(i+0.5)
# print("Table for Eigen Values for xmax = 10")
# data = {
#     "N":n,
#     "Numerical Eigen Value":num_eig_val,
#     "Analytical Eigen Value":anal_eig_val,
# }
# print(pd.DataFrame(data))