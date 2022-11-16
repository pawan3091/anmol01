import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import pandas as pd
def V():
    return 0
def analytical(x,n):
        if (i%2)== 0:
            return 2*(np.cos(n*np.pi*x))**2
        else:
            return 2*(np.sin(n*np.pi*x))**2
def analytical_1(x,n):
        if (i%2)== 0:
            return np.sqrt(2)*(np.cos(n*np.pi*x))
        else:
            return np.sqrt(2)*(np.sin(n*np.pi*x))

def sub_plot(ax,a,b,d,title,x_label,y_label,key=1):
    ax.scatter(a,b,label="Numerical Value",marker="*",color="red")
    if key == 1:
        ax.plot(a,d,label="Inbuilt Solution")
        ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.legend()
    ax.grid(True)

def diag_mat(n,xi,xf):
    h = (xf-xi)/(N+1)
    a,v,e=np.zeros((n,n)),np.zeros((n,n)),np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i==j:
                a[i][i]=-2
                a[i][i-1]=1
                a[i-1][i]=1
                v[i][i]=0
            else:
                a[i][j]=0
                v[i][j]=0
                e[i][i]=1
    A=a*(-1/(h**2))+v
    eig = eigh(A) 
    return eig

n=2;N=2
xi=-1/2;xf=1/2
x=np.linspace(xi,xf,N)
U=diag_mat(n,xi,xf)[0]
data={
    "Eigen Value": diag_mat(n,xi,xf)[0],
    "Eigen Vector 1":diag_mat(n,xi,xf)[1][0],
    "Eigen Vector 2":diag_mat(n,xi,xf)[1][1]
}
print(pd.DataFrame(data))


n=100;N=100
xi=-1/2;xf=1/2
x=np.linspace(xi,xf,N)
print("First 11 Eigen Values")
data={
    "Eigen Value": diag_mat(n,xi,xf)[0][:11],
}
print(pd.DataFrame(data))


fig2, ((axx1, axx2), (axx3, axx4)) = plt.subplots(2,2)
fig2.suptitle("Probability Density")
dict = {'0': axx1,'1': axx2,'2': axx3,'3': axx4}
for i in range(0,4,1):
    U=diag_mat(n,xi,xf)[1][:,i]
    u_norm=U/np.sqrt(integrate.simps(np.power(U,2),x))
    sub_plot(dict[str(int(i))],x,np.power(u_norm,2),analytical(x,i+1),f'N={i+1}',"\u03BE","$u(\u03BE)^2$")
plt.tight_layout()
plt.show()

fig2, ((axx1, axx2), (axx3, axx4)) = plt.subplots(2,2)
fig2.suptitle("Wavefunction")
dict = {'0': axx1,'1': axx2,'2': axx3,'3': axx4}
for i in range(0,4,1):
    U=diag_mat(n,xi,xf)[1][:,i]
    u_norm=U/np.sqrt(integrate.simps(np.power(U,2),x))
    sub_plot(dict[str(int(i))],x,np.power(u_norm,1),analytical_1(x,i+1),f'N={i+1}',"\u03BE","$u(\u03BE)^2$",key=0)
plt.tight_layout()
plt.show()


