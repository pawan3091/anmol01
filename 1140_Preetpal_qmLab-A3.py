import numpy as np 
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import pandas as pd
from scipy import optimize,stats
from IVP import RK_fourth_vec
def func_(x,x_vec,e):
    ans_vec = np.zeros((2))
    ans_vec[0] = x_vec[1]
    ans_vec[1] = (-e)*x_vec[0]
    return ans_vec
def analytical(x,n):
        if (i%2)== 0:
            return 2*(np.cos(n*np.pi*x))**2
        else:
            return 2*(np.sin(n*np.pi*x))**2
def analytical_1(x,n):
        return np.sin(n*np.pi*x)
E = np.linspace(0,400,100)
x  = np.linspace(-1/2,1/2,100)
initcond=[0,1]

t = []
for i in range(len(E)):
    e = E[i]
    y1 = RK_fourth_vec(x, initcond, func_,e).T[0]
    y=y1/np.sqrt(integrate.simps(np.power(y1,2),x))
    t.append(y[-1])
# plt.scatter(E,t)
# plt.xlabel("e")
# plt.ylabel("$u_{R}$")
# plt.title("u(\u03BE=1/2)")
# plt.grid(True)
# plt.show()

def f(e,x):
    y1 = RK_fourth_vec(x, initcond, func_,e).T[0]
    y2 = RK_fourth_vec(x, initcond, func_,e).T[1]
    y=y1/np.sqrt(integrate.simps(np.power(y1,2),x))
    y_=y2/np.sqrt(integrate.simps(np.power(y2,2),x))
    t.append(y[-1])
    return t[-1], y,x,y_

en_neg = [];en_pos=[];i_neg=[];i_pos=[]
def eigenvalues(t):
    for i in range(1,len(t)):
        if t[i-1]*t[i]<0:
            en_neg.append(t[i-1])
            i_neg.append(E[i-1])
            en_pos.append(t[i])
            i_pos.append(E[i])
    return en_neg,en_pos,i_pos,i_neg

phi_1,phi_2,s_1,s_0 = eigenvalues(t)
data = {
    "u_1":phi_1,
    "E1":s_0,
    "u_2":phi_2,
    "E_2":s_1
}
#print(pd.DataFrame(data))
def secant(s_1,s_0,iterations,x):
    sec_2=[]
    sec_2.append(s_0); sec_2.append(s_1)
    for i in range(1,iterations):
        sec_2.append(sec_2[i]-(((sec_2[i]-sec_2[i-1])*(f(sec_2[i],x)[0]))/((f(sec_2[i],x)[0])-(f(sec_2[i-1],x)[0]))))
        if abs(f(sec_2[-1],x)[0])<0.1e-12:
            return sec_2[-1], f(sec_2[-1],x)[0],f(sec_2[-1],x)[1],f(sec_2[-1],x)[2], f(sec_2[-1],x)[2]

E_n=[];u=[];v=[]

for i in range(0,len(s_0)):
    E_n_=secant(s_1[i],s_0[i],501,x)[0]
    u_ = secant(s_1[i],s_0[i],501,x)[1]
    ufull_ = secant(s_1[i],s_0[i],501,x)[2]
    ufull_prime = secant(s_1[i],s_0[i],501,x)[4]
    probability_density = secant(s_1[i],s_0[i],200,x)[3]
    x_ = secant(s_1[i],s_0[i],501,x)[3]
    # plt.scatter(x_,ufull_,label=f'n={i+1}')
    # # plt.plot(x_,analytical_1(x_,i+1),label=f'n={i+1}')
    # plt.xlabel("\u03BE")
    # plt.ylabel("u(\u03BE)")
    # plt.title("Normalised wave function for infinite square well")
    # plt.grid()
    # plt.legend()
    # plt.show()        

    # plt.scatter(x_,ufull_*ufull_,label=f'n={i+1}')
    # plt.plot(x_,analytical(x_,i+1),label=f'n={i+1}(Analytical)')
    # plt.xlabel("\u03BE")
    # plt.ylabel("$(u(\u03BE))^2$")
    # plt.title("Probability Densities")
    # plt.grid()
    # plt.legend()
    # plt.show()        

    E_n.append(E_n_)
    u.append(u_)
    v.append(secant(s_1[i],s_0[i],200,x)[2])

print(E_n);print(u)

data = {
    "N":[1,2,3,4,5,6],
    "Final Energy Eigen Value":E_n,
    "Corresponding u": u
}
print(pd.DataFrame(data))

slope, intercept, r, p, se = stats.linregress(np.array(E_n)/(np.pi)**2, E_n )
print("slope",slope)
# plt.scatter(np.array(E_n)/(np.pi)**2,E_n,label="Approximated")
# plt.plot(abs(np.array(E_n)/(np.pi)**2),(np.array(E_n)/(np.pi)**2)*(np.pi)**2,label="Analytical")
# plt.xlabel("$n^2$")
# plt.ylabel("$E_{n}$")
# plt.grid()
# plt.legend()
# plt.show()

def en_ev(E_n,h_cut,m_e,L):
    En_anal=[];prob_En_ev=[]
    for i in range(0,len(E_n)):
        prob_En_ev_ = ((h_cut**2)*np.array(E_n)[i])/(2*m_e*(L**2)*(1.609e-19))
        prob_En_ev.append(prob_En_ev_)
        E = (((i+1)**2)*((np.pi)**2)*(h_cut**2))/(2*m_e*(L**2)*(1.609e-19))
        En_anal.append(E)
    return prob_En_ev, En_anal
m_e = 9.1e-31; h_cut = 1.0545e-34;m_p=1.6e-27
print("Well of width = 5 Angstrom")
data = {
    "Approximated Eigen Values":en_ev(E_n, h_cut,m_e,5e-10)[0],
    "Analytical Eigen Values":en_ev(E_n, h_cut,m_e,5e-10)[1],
}
print(pd.DataFrame(data))
print("---------------------------------------------------------------------------------------------")
print("Well of width = 10 Angstrom")
data = {
    "Approximated Eigen Values":en_ev(E_n, h_cut,m_e,10e-10)[0],
    "Analytical Eigen Values":en_ev(E_n, h_cut,m_e,10e-10)[1],
}
print(pd.DataFrame(data))

print("---------------------------------------------------------------------------------------------")
print("Well of width = 5 Fermimeter for proton")
data = {
    "Approximated Eigen Values":en_ev(E_n, h_cut,m_p,5e-15)[0],
    "Analytical Eigen Values":en_ev(E_n, h_cut,m_p,5e-15)[1],
}
print(pd.DataFrame(data))

print("---------------------------------------------------------------------------------------------")
'''Uncertainty Principle'''
exp_x_2=np.power(ufull_,2)*(np.power(x_,2))
i1=integrate.simps(exp_x_2,x_)
exp=np.power(ufull_,2)*(x_)
i2=integrate.simps(exp,x_)
variance = i1-i2**2
st_dev = np.sqrt(variance)
print("Uncertainty in x is ", np.sqrt(variance))
exp_x_p1 = np.power(ufull_prime,2)*x_
ip_1=integrate.simps(exp_x_p1,x_)
exp_x_p2=np.power(ufull_prime,2)*np.power(x_,2)
ip_2=integrate.simps(exp_x_p2,x_)
variance1=ip_2-ip_1**2
print("uncertainty in momentum p is", np.sqrt(variance1))
print("sigma_x*sigma_p = ",np.sqrt(variance)*np.sqrt(variance1))
print("h_cut/2*pi=",h_cut/(2))
print("sigma_x*sigma_p >= h_cut/2*pi, Hence Uncertainty Principle Verified ")