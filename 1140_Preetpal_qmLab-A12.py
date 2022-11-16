import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import pandas as pd
from scipy.special import assoc_laguerre
from scipy.optimize import fsolve
def diag_mat(xi,xf,N,l,ratio):
    X = np.linspace(xi,xf,N+2)
    x=X[1:-1]
    h = x[1]-x[0]
    a,v=np.zeros((len(x),len(x))),np.zeros((len(x),len(x)))
    for i in range(len(x)):
        for j in range(len(x)):
            if i==j:
                a[i][i]=2/h**2
                v[i][i]=(-2/x[i])*np.exp(np.divide(-x[i],ratio))+((l)*(l+1))/(2*(x[i]**2))
            elif i==j+1:
                a[i][j]=-1/h**2
            elif i == j-1:
                a[i][j] = -1/h**2
    A=(a+v)
    eig = eigh(A) 
    return eig,x
def V(x,ratio):
    v_x=(-2/x)*np.exp(np.divide(-x,ratio))
    v_coulomb = (-2/x)
    return v_x, v_coulomb
def graph(x,y,label,xlabel,ylabel,title):
    plt.scatter(x,y,label=f'ratio={j}')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()
    plt.legend()

N=500
for i in range(0,1,1):
    xi=0;xf=20
    U_a=[];ratio_=[];v_x=[];v_coulomb=[]
    ratio_=[2,5,10,20,100]
    for j in ratio_:
        ratio=j
        U_,x=diag_mat(xi,xf,N,0,j)
        U=U_[1][:,i]
        U_a.append(U_[0][:1][0])
        u_norm=U/np.sqrt(integrate.simps(np.power(U,2),x))
        v_x, v_coulomb=V(x,j)
        # graph(x,v_x,j,"$\u03BE$","$V(\u03BE)$","Screened Coulomb Potential")    
        # graph(x,v_coulomb,j,"$\u03BE$","$V_c$","Coulomb Potential")    
        # graph(x,np.power(u_norm,1),j,"x","$(u_r(\u03BE))$","Radial Wavefunction for n=1,l=0")
        # graph(x,np.power(u_norm,2),j,"x","$(u_r(\u03BE))^2$","Radial Probability for n=1,l=0")
# plt.show()

#---------------------Q-a(ii)--------------------------#
p=[]
for i in range(1,len(ratio_)+1):
    p.append(-1/1**2)
data ={
    "ratio":ratio_,
    "Numerical Eigen Values ": U_a,
    # "Analytical Eigen Values ":p ,
}
print("Bound Energy State Eigen Value for n=0")
print(pd.DataFrame(data))
# plt.scatter(ratio_,U_a)
# plt.grid()
# plt.xlabel("ratio")
# plt.ylabel("Eigen Value")
# plt.show()