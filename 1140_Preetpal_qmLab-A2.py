import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import scipy.integrate as integrate

def My_RK4(Y0,func,xi,xf,n):
    x = np.linspace(xi,xf,n)
    h = (xf-xi)/n
    Y = np.zeros((n,len(Y0)))
    Y[0,:] = Y0
    for i in range(n-1):
        k1 = h*func(x[i],Y[i,:])
        k2 = h*func(x[i]+h*0.5,Y[i,:]+k1*0.5)
        k3 = h*func(x[i]+h*0.5,Y[i,:]+k2*0.5)
        k4 = h*func(x[i]+h,Y[i,:]+k3)
        Y[i+1,:] = Y[i,:]+(k1+2*k2+2*k3+k4)/6
    return Y
    
def func_(x,x_vec):
    ans_vec = np.zeros((2))
    ans_vec[0] = x_vec[1]
    ans_vec[1] = (-n**2)*(8)*x_vec[0]
    return ans_vec

def analytical(x,n):
    u_ana = []
    for i in range(n):
        if (i%2)== 0:
            u = np.sqrt(2)*np.cos(n*np.pi*x)
            u_ana.append(u)
        else:
            u = np.sqrt(2)*np.sin(n*np.pi*x)
            u_ana.append(u)
    return u

x  = np.linspace(-1/2,1/2,100)
for n in range(1,5):
    y1  = My_RK4([0,1],func_,0,1,100).T[0]
    y2  = My_RK4([0,1],func_,0,1,100).T[1]
    #normalized
    y=y1/np.sqrt(integrate.simps(y1**2,x))
    plt.rcParams["figure.figsize"] = (15,10)
    plt.plot(x,y,label ="Normalized E =8,u'=1 ")
    plt.plot(x,y1,label ="Not Normalized E =8,u'=1 ")
    
    y1  = My_RK4([0,1.5],func_,0,1,100).T[0]
    y2  = My_RK4([0,1.5],func_,0,1,100).T[1]
    #normalized
    y=y1/np.sqrt(integrate.simps(y1**2,x))
    plt.rcParams["figure.figsize"] = (15,10)
    plt.plot(x,y,label = "Normalized  E =8,u'=1.5")
    plt.plot(x,y1,label = "Not Normalized  E =8,u'=1.5")
    
    #Analytical solution
    u = analytical(x,n)
    plt.plot(x,y,marker = ".",label = "Analytical solution")
    plt.xlabel("\u03BE")
    plt.ylabel("u(\u03BE)")
    plt.title(label=f'n={n}')
    plt.grid()
    plt.legend()
    plt.show()

def func_(x,x_vec):
    ans_vec = np.zeros((2))
    ans_vec[0] = x_vec[1]
    ans_vec[1] = (-n**2)*(11)*x_vec[0]
    return ans_vec

def analytical(x,n):
    u_ana = []
    for i in range(n):
        if (i%2)== 0:
            u = np.sqrt(2)*np.cos(n*np.pi*x)
            u_ana.append(u)
        else:
            u = np.sqrt(2)*np.sin(n*np.pi*x)
            u_ana.append(u)
    return u

x  = np.linspace(-1/2,1/2,100)

for n in range(1,5):
    y1  = My_RK4([0,1],func_,0,1,100).T[0]
    y2  = My_RK4([0,1],func_,0,1,100).T[1]
    #normalized
    y=y1/np.sqrt(integrate.simps(y1**2,x))
    plt.rcParams["figure.figsize"] = (15,10)
    plt.plot(x,y,label ="Normalized E =11,u'=1 ")
    plt.plot(x,y1,label ="Not Normalized E =11,u'=1 ")
    
    y1  = My_RK4([0,1.5],func_,0,1,100).T[0]
    y2  = My_RK4([0,1.5],func_,0,1,100).T[1]
    #normalized
    y=y1/np.sqrt(integrate.simps(y1**2,x))
    plt.rcParams["figure.figsize"] = (15,10)
    plt.plot(x,y,label ="Normalized E =11,u'=1 ")
    plt.plot(x,y1,label ="Not Normalized E =11,u'=1 ")
    
    #Analytical solution
    u = analytical(x,n)
    plt.plot(x,y,marker = ".",label = "Analytical solution")
    plt.grid()
    plt.xlabel("\u03BE")
    plt.ylabel("u(\u03BE)")
    plt.title(label=f'n={n}')
    plt.legend()
    plt.show()


def func_(x,x_vec):
    ans_vec = np.zeros((2))
    ans_vec[0] = x_vec[1]
    ans_vec[1] = (-n**2)*(E**2)*x_vec[0]
    return ans_vec

x  = np.linspace(-1/2,1/2,100)
epsilonvalue = np.linspace(np.pi*0.9,1.1*np.pi,10)

for n in range(1,5):
    for i in range(len(epsilonvalue)):
        E = epsilonvalue[i]
        y1  = My_RK4([0,1],func_,0,1,100).T[0]
        y2  = My_RK4([0,1],func_,0,1,100).T[1]
        plt.rcParams["figure.figsize"] = (15,10)
        plt.plot(x,y1,label=f'$e$ ={epsilonvalue[i]**2}')
    plt.title(label=f'n={n}')
    plt.xlabel("\u03BE")
    plt.ylabel("u(\u03BE)")
    plt.grid()
    plt.legend()
    plt.show()
