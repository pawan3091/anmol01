import math
import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pandas as pd
def alpha(x):
    return -(1+x**2)
def func_(x,x_vec):
    ans_vec = np.zeros((2))
    ans_vec[0] = x_vec[1]
    ans_vec[1] = (1+x**2)*x_vec[0]
    return ans_vec

def sub_plot(ax,a,b,d,title):
    ax.scatter(a,b,label="Numerical Value")
    ax.plot(a,d,label="Inbuilt Solution")
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("u")
    ax.legend()
    ax.grid(True)

def numerov(x_min, x_max, u_0, u_prime, N):
    c_i =[];u=[]
    x = np.linspace(x_min,x_max,N+1)
    Alpha = alpha(x)
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
    sol = solve_ivp(func_, [x_min,x_max], [u_0,u_prime],dense_output=True)
    inbuilt = sol.sol(x)
    return x, u, inbuilt[0],c_i
'''--------------------------------------------------N=2----------------------------------------------'''
p = numerov(0,1,1,0,2)
data = {
    "x":p[0],
    "u_num":p[1],
    "u_inbuilt": p[2],
    "E = u_inbuilt - u_num": abs(p[2]-p[1])
}
print(pd.DataFrame(data))

'''--------------------------------------------------N=4----------------------------------------------'''
p = numerov(0,1,1,0,4)
data = {
    "x":p[0],
    "u_num":p[1],
    "u_inbuilt": p[2],
    "E = u_inbuilt - u_num": abs(p[2]-p[1])
}
print(pd.DataFrame(data))
'''-----------------------------------------------------------------------------------------------------'''
n = np.arange(1,7,1)
# fig2 = plt.figure(figsize=(12,12))
# fig, (ax1, ax2) = plt.subplots(2)
fig2, ((axx1, axx2), (axx3, axx4), (axx5, axx6)) = plt.subplots(3,2)
fig2.suptitle('u(x) vs x for N intervals')
dict = {'0': axx1,'1': axx2,'2': axx3,'3': axx4,'4': axx5,'5': axx6}
for i in range(1,7):
    p = numerov(0,1,1,0,2**i)
    sub_plot(dict[str(i-1)],p[0],p[1],p[2],f'N={2**i}')
plt.tight_layout()
plt.show()
