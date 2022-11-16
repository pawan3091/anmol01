import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import pandas as pd
from scipy.special import assoc_laguerre
from scipy.optimize import fsolve
def diag_mat(xi,xf,N,l):
    X = np.linspace(xi,xf,N+2)
    x=X[1:-1]
    h = x[1]-x[0]
    a,v=np.zeros((len(x),len(x))),np.zeros((len(x),len(x)))
    for i in range(len(x)):
        for j in range(len(x)):
            if i==j:
                a[i][i]=2/h**2
                v[i][i]= ((x[i]**2)) +((2/3)*(const)*(x[i]**3))
            elif i==j+1:
                a[i][j]=-1/h**2
            elif i == j-1:
                a[i][j] = -1/h**2
    A=(a+v)
    eig = eigh(A) 
    return eig,x

def graph(x,y,label,xlabel,ylabel,title):
    plt.plot(x,y,label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()
    plt.legend()

N=500;const_=[0];no_eigen_values=10;conversion_factor=(197.3/2)*np.sqrt(100/940)
for i in range(0,4,1):
    const_.append(10**-i)
k=100;m=940;h_cut=6.5821*(10**-22);r0=np.power((h_cut**2)/m*k,1/4)
for const in const_:
    U_a=[];anal_eig_val=[];state=[]
    for i in range(0,no_eigen_values,1):
        state.append(i+1)
        anal_eig_val.append(((2*i+1)-(1/8)*((const)**2)*(15*np.power((2*i+1),2)+7))*conversion_factor)
        xi=-5;xf=5
        U_,x=diag_mat(xi,xf,N,0)
        U=U_[1][:,i]
        U_a=(U_[0][:no_eigen_values])*conversion_factor
        u_norm=U/np.sqrt(integrate.simps(np.power(U,2),x))
        # graph(x,np.power(u_norm,1),f'n={i+1}',"$\u03BE$","$(u_r(\u03BE))$",f'alpha={const}')
        # graph(x,np.power(u_norm,2),f'n={i+1}',"$\u03BE$","$(u_r(\u03BE)^2)$",f'alpha={const}')
        graph(np.arange(1,no_eigen_values+1,1),U_a,None,"n","Eigen Value","Eigen Value vs n")        
    plt.show()

#---------------------Q-a(ii)--------------------------#
    print("Eigen Values for different alpha=",const)
    data ={
        "n":state,
        "Numerical Eigen Values ": U_a,
        "Analytical Eigen Values ":anal_eig_val ,
    }
    print(pd.DataFrame(data))