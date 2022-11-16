import math
import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import scipy.integrate as integrate
import pandas as pd
def alpha(n,x,del_e):
    return 2*(n+(1/2)+del_e-(x**2)/2)
def func_(x,x_vec):
    ans_vec = np.zeros((2))
    ans_vec[0] = x_vec[1]
    ans_vec[1] = -2*(n+(1/2)+(1e-2)-(x**2)/2)*x_vec[0]
    return ans_vec
def sub_plot(ax,a,b,d,title,x_label,y_label):
    ax.scatter(a,b,label="Numerical Value",marker="*",color="red")
    ax.plot(a,d,label="Inbuilt Solution")
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend()
    ax.grid(True)

def numerov(x_min, x_max, u_0, u_prime,n, N,del_e):
    c_i =[];u=[]
    x = np.linspace(x_min,x_max,N+1)
    Alpha = alpha(n,x,del_e)
    h = x[1]-x[0]
    u_1 = 1 + ((h**2)/math.factorial(2)) + 3*((h**4)/math.factorial(4))
    u.append(u_0);u.append(u_1)
    ddx_12 = (h**2)/12
    for i in range(0,N+1):
        c_i_ = 1 + np.multiply(ddx_12,Alpha[i])
        c_i.append(c_i_)
    for i in range(2,N+1):
        u_ = (1/c_i[i])*(((12-10*c_i[i-1])*u[i-1])-c_i[i-2]*u[i-2])
        u.append(u_)
    u_norm=u/np.sqrt(integrate.simps(np.power(u,2),x))
    sol = solve_ivp(func_, [x_min,x_max], [u_0,u_prime],dense_output=True)
    inbuilt = sol.sol(x)
    sol_ = inbuilt[0]/np.sqrt(integrate.simps(np.power(inbuilt[0],2),x))
    return x, u_norm, sol_, Alpha
'''--------------------------------------------------N=2----------------------------------------------'''
def combine(list1, list2):
    list1_ = np.delete(list1, len(list1)-1)
    result_list = []
    result_list = list(list1_) 
    for item in list2:
        result_list.append(item)
    return result_list

def parity(x_min, x_max, u_0, u_prime,n, N,del_e):
    p = numerov(x_min, x_max, u_0, u_prime,n, N,del_e)
    t = p[1][::-1]
    t_1 = p[2][::-1]
    array = [];array_1=[]
    if n%2 == 0:
        array_1 = combine(t_1, p[2])        
        array = combine(t, p[1])
    elif n%2 != 0:
        array_1 = combine((np.multiply(-1,t_1)),p[2])
        array = combine((np.multiply(-1,t)),p[1])
    #print("u_norm",array)
    # plt.scatter(np.linspace(-x_max,x_max,2*N+1),array_1)
    # plt.plot(np.linspace(-x_max,x_max,2*N+1),array)
    # plt.show()
    return array, array_1, p[3]
# '''----------------------------u vs x Plotting--------------------------------------'''
# x_min = 0; n = 0; x_max = n+1;N=50
# if n%2 ==0:
#     u_0 = 1
#     u_prime=0
# else:
#     u_0 = 0
#     u_prime = 1
# fig2, ((axx1, axx2), (axx3, axx4)) = plt.subplots(2,2)
# fig2.suptitle(f'N={n}')
# dict = {'0': axx1,'1': axx2,'2': axx3,'3': axx4}
# for i in range(2,10,2):
#     par_1, par_2 = parity(x_min, x_max, u_0, u_prime,n, N,10**-i)
#     if n ==2 or n==3:
#         sub_plot(dict[str(int((i-2)/2))],np.linspace(-x_max,x_max,2*N +1),np.multiply(-1,par_1),np.multiply(-1,par_2),f'del_e={10**-i}',"\u03BE","u(\u03BE)")
#     else:
#         sub_plot(dict[str(int((i-2)/2))],np.linspace(-x_max,x_max,2*N +1),par_1,par_2,f'del_e={10**-i}',"\u03BE","u(\u03BE)")
# plt.tight_layout()
# plt.show()

'''----------------------------Probability Density vs x Plotting--------------------------------------'''
# x_min = 0 ;N=100;del_e=1e-6
# fig2, ((axx1, axx2), (axx3, axx4)) = plt.subplots(2,2)
# fig2.suptitle("Probability Density")
# dict = {'0': axx1,'1': axx2,'2': axx3,'3': axx4}
# for n in range(0,4,1):
#     x_max = n+1
#     if n%2 ==0:
#         u_0 = 1
#         u_prime=0
#     else:
#         u_0 = 0
#         u_prime = 1
#     par_1, par_2 = parity(x_min, x_max, u_0, u_prime,n, N,del_e)
#     if n ==2 or n==3:
#         sub_plot(dict[str(int(n))],np.linspace(-x_max,x_max,2*N +1),np.power(np.multiply(-1,par_1),2),np.power(np.multiply(-1,par_2),2),f'N={n}',"\u03BE","$u(\u03BE)^2$")
#     else:
#         sub_plot(dict[str(int(n))],np.linspace(-x_max,x_max,2*N +1),np.power(par_1,2),np.power(par_2,2),f'N={n}',"\u03BE","$u(\u03BE)^2$")
# plt.tight_layout()
# plt.show()
'''----------------------------------------Q-C----------------------------------------------'''
x_min = 0; n = 0; x_max = n+1;N=10;del_e=1e-6;e = 1.6e-19; h_cut = 1.0545e-34;umega=5.5e14
if n%2 ==0:
    u_0 = 1
    u_prime=0
else:
    u_0 = 0
    u_prime = 1
par_1, par_2, alpha_ = parity(x_min, x_max, u_0, u_prime,n, N,del_e)
print("                            Energy in eV at N=0")
data = {
    "Approximated Eigen Values":np.multiply(alpha_,(h_cut*umega)/e),
    "Analytical Eigen Values":np.multiply(alpha_-del_e,(h_cut*umega)/e)
}
print(pd.DataFrame(data))
'''------------------------------------------------------Q_d-------------------------------------------'''
# x_min = 0; n = 0; x_max = n+4;N=100;del_e=1e-6
# if n%2 ==0:
#     u_0 = 1
#     u_prime=0
# else:
#     u_0 = 0
#     u_prime = 1
# x_test = np.linspace(x_min,x_max,N+1)
# alfa = alpha(n,x_test,del_e)
# # par_1, par_2, alpha_, x_ = parity(x_min, x_max, u_0, u_prime,n, N,del_e)
# j=0;prob_density=[]
# for i in alfa:
#     j+=1
#     if i<0:
#         print("Probability of finding electron in the classically forbidden region when it is in the ground state")
#         # print(len(np.power(par_1[j-1:],2)))
#         # print(integrate.simps(np.power(par_1[j-1:],2),x_[j-1:]))
#         print(len(alfa[j-1:]))
#         print((alfa[j-1:]))
#         p_ = numerov(x_test[j-1], x_max, u_0, u_prime,n, len(x_test[j-1:])-1,del_e)
#         print(len(x_test[j-1:]))
#         print(len(p_[1]))
#         print("probability",integrate.simps(np.power(p_[1],2),x_test[j-1:]))
#         # par_1, par_2, alpha_, x_ = parity(x_test[j-1], x_max, u_0, u_prime,n, len(x_test[j-1:]),del_e)
# #       print(np.power(par_1,2))
#         # t = p_[1][::-1]
#         # print(integrate.simps(np.power(t,2),np.linspace(-x_test[j-1],-x_max,len(x_test[j-1:])))+integrate.simps(np.power(p_[1],2),np.linspace(x_test[j-1],x_max,len(x_test[j-1:]))))
        
#         # print(integrate.simps(np.power(t,2),np.linspace(-x_test[j-1],-x_max,len(x_test[j-1:]))))
#         break
    