import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import pandas as pd
from scipy.special import assoc_laguerre
def diag_mat(xi,xf,N,l):
    X = np.linspace(xi,xf,N+2)
    x=X[1:-1]
    h = x[1]-x[0]
    a,v=np.zeros((len(x),len(x))),np.zeros((len(x),len(x)))
    for i in range(len(x)):
        for j in range(len(x)):
            if i==j:
                a[i][i]=2/h**2
                v[i][i]=(-2/x[i])+((l)*(l+1))/(x[i]**2)
            elif i==j+1:
                a[i][j]=-1/h**2
            elif i == j-1:
                a[i][j] = -1/h**2
    A=(a+v)
    eig = eigh(A) 
    return eig,x

def Analytic(x,n,l):
    return ((2*x/n)**(l)*assoc_laguerre(2*x/n,n-l-1,2*l+1))/(np.exp(x/n))*x

m=1.67*10**(-27)
def Veff(x, l):
    Vef = (l * (l + 1) / (x ** 2)) - (2 / x)
    V = - (2 / x)

    return Vef, V

def plot(i, l, power):
    H, x = diag_mat(0.01, 30, 1000, l)
    u = H[1][:, i]
    c = integrate.simps(u ** 2, x)
    N = u / np.sqrt(c)
    plt.plot(x, N ** power,label=f'l,n={l,i+1}')
    plt.title("Radial Probability Density for different n and l")
    plt.xlabel("x")
    plt.ylabel("$(u_r(\u03BE))^2$")
    plt.grid()
    plt.legend()

    # plt.show()

N=1000
xi=10**-14;xf=30

#A-i
# for i in range(0, 4):
#     H, x = diag_mat(0.1, 50, N, i)
#     Vef, V = Veff(x, i)
#     plt.scatter(x, V,label=f'V,l={i}')
#     # plt.plot(x, Vef,label=f'V_eff,l={i}')
#     plt.xlabel("x")
#     plt.ylabel("V_eff")
#     plt.legend()
# plt.show()

# for i in range(0, 4):
#     H, x = diag_mat(0.1, 50, N, i)
#     Vef, V = Veff(x, i)
#     plt.plot(x, V,label=f'V,l={i}')
#     plt.xlabel("x")
#     plt.ylabel("V")
#     plt.legend()
# plt.show()

#---------------------------------Q1_ii----------------------------------------#
for i in range(0,4,1):
    U_,x=diag_mat(xi,xf,N,0)
    U=U_[1][:,i]
    u_norm=U/np.sqrt(integrate.simps(np.power(U,2),x))
    U_anal = Analytic(x,i,0)
    u_anal_normalised=U_anal/np.sqrt(integrate.simps(np.power(U_anal,2),x))
    plt.scatter(x,np.power(u_norm,1),label=f'Numerical, n={i+1}')
    plt.plot(x,np.power(u_anal_normalised,1),label=f'Analytical, n={i+1}')
    plt.xlabel("x")
    plt.ylabel("$(u_r(\u03BE))$")
    plt.title("Radial Wavefunction for l=0")
plt.legend()
plt.show()

# #---------------------Q-a(ii)--------------------------#
# print("First 10 Energy Eigen Values for l=0 and r_max=30")
# p=[]
# for i in range(1,11):
#     p.append(-1/i**2)
# data ={
#     "Numerical Eigen Values ": U_[0][:10],
#     "Analytical Eigen Values ":p 
# }
# print(pd.DataFrame(data))

# #---------------------------------Qb_i----------------------------------------#
# for i in range(0,4,1):
#     U_,x=diag_mat(xi,xf,N,1)
#     U=U_[1][:,i]
#     u_norm=U/np.sqrt(integrate.simps(np.power(U,2),x))
#     U_anal = Analytic(x,i,1)
#     u_anal_normalised=U_anal/np.sqrt(integrate.simps(np.power(U_anal,2),x))
#     plt.scatter(x,np.power(u_norm,1),label=f'Numerical, n={i+1}')
#     plt.plot(x,np.power(u_anal_normalised,1),label=f'Analytical, n={i+1}')
# plt.legend()
# plt.show()

# print("First 10 Energy Eigen Values for l=1 and r_max=30")
# p=[]
# for i in range(1,11):
#     p.append(-1/(i+1)**2)
# data ={
#     "Numerical Eigen Values ": U_[0][:10],
#     "Analytical Eigen Values ":p 
# }
# print(pd.DataFrame(data))

#---------------------------------Qb_i----------------------------------------#
# for i in range(0,4,1):
#     U_,x=diag_mat(xi,xf,N,2)
#     U=U_[1][:,i]
#     u_norm=U/np.sqrt(integrate.simps(np.power(U,2),x))
#     U_anal = Analytic(x,i,2)
#     u_anal_normalised=U_anal/np.sqrt(integrate.simps(np.power(U_anal,2),x))
#     plt.scatter(x,np.power(u_norm,1),label=f'Numerical, n={i+1}')
#     plt.plot(x,np.power(u_anal_normalised,1),label=f'Analytical, n={i+1}')
# plt.legend()
# plt.show()

# print("First 10 Energy Eigen Values for l=2 and r_max=30")
# p=[]
# for i in range(1,11):
#     p.append(-1/(i+2)**2)
# data ={
#     "Numerical Eigen Values ": U_[0][:10],
#     "Analytical Eigen Values ":p 
# }
# print(pd.DataFrame(data))

# C

# '''for n=1:'''
# plot(0, 0, 2)
# plt.show()
# '''for n=2:'''
# plot(1, 0, 2)
# plot(1, 1, 2)
# plt.show()
#n=3
# plot(2, 0, 2)
# plot(2, 1, 2)
# plot(2, 2, 2)
# plt.show()


