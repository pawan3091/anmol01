import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

m_e = 9.1e-34; q_e  = 1.6e-19; epsilon=8.854187817e-12; h_cut = 1.0545e-34
def f_cln(n):
    return (m_e*(q_e)**4)/(32*((np.pi)**3)*(epsilon**2)*(h_cut**3)*(n**3))
    
def f_qn(n):
    return (m_e*((q_e)**4)*((2*n)-1))/(64*((np.pi)**3)*(epsilon**2)*(h_cut**3)*(n**2)*((n-1)**2))

def E(n):
    return -13.6/(n**2)

t = [];rel_dif =[];f_quantum=[];f_classical = []
p = np.arange(0.5,200,0.5)
for i in 10**p:
    del_f = (f_qn(i)-f_cln(i))
    t.append(i );rel_dif.append(np.abs(del_f/f_qn(i))*100);f_quantum.append(f_qn(i));f_classical.append(f_cln(i))
    if np.abs(del_f/f_qn(i)) < 10**-5:
        break;
table = {
    "n" : t,
    "f_qn" : f_quantum,
    "f_cl" : f_classical,
    "Relative Difference(%)" : rel_dif,
}
print(pd.DataFrame(table))
plt.plot(np.log10(t ),rel_dif, marker="*")
plt.xlabel("ln(n)")
plt.ylabel("Relative Difference(%)")
plt.title("Relative Difference(%) vs ln(n)")
plt.grid(True)
plt.legend()
plt.show()
#----------------------------Energy Level Diagram ---------------------------------#
n = [];E_n=[]
for i in range(1,11):
    t = np.ones(10)*(E(i))
    E_n.append(t)
    n.append(i)
for i in range(1,9):
    plt.plot(n,E_n[i],label=f'n={i}')
    plt.ylabel("En(in eV)")
    plt.xlabel("n")
    plt.title("Energy Level Diagram")
    plt.legend()
plt.show()

