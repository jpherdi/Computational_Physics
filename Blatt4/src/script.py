import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
# Set fontsize larger for plots
matplotlib.rcParams.update({'font.size': 20})

# Exercise 1
for i in ["",10,20,50]:
    # Generate data from files
    a = np.genfromtxt("output/data"+str(i)+".txt", unpack=True)
    a.reshape(512,512)
    # Plotting of the apes with k=512, 10, 20, 50
    plt.figure(figsize=(12,12))
    plt.imshow(a, cmap="gray")
    plt.tight_layout()
    plt.savefig("output/Bild"+str(i)+".pdf")

# Exercise 2, generate data from file
t_random, t_LU, t_solve = np.genfromtxt("output/times.txt", unpack=True)

# Fitting just for fun, to get an idea for the order of the exponent (because of log-log we fitted with a linear function)
num = np.arange(1,len(t_random)+1)
num_new = np.arange(12, 2000)
param_LU, err_LU = curve_fit(lambda x, a, N: a+x*N, np.log(num), np.log(t_LU))
param_random, err_random = curve_fit(lambda x, a, N: a+x*N, np.log(num), np.log(t_random))
param_solve, err_solve = curve_fit(lambda x, a, N: a+x*N, np.log(num), np.log(t_solve))

# Plotting for three different times
plt.figure(figsize=(12,8))
plt.plot(num, t_random, "x", label=r"$t_\mathrm{random}$", color="C0")
plt.plot(num_new, np.exp(param_random[0]+np.log(num_new)*param_random[1]), label=r"$N^{"+f"{param_random[1]:.2f}"+r"}$- Random", color="C0")
plt.plot(num, t_LU, "x", label=r"$t_\mathrm{LU}$", color="C1")
plt.plot(num_new, np.exp(param_LU[0]+np.log(num_new)*param_LU[1]), label=r"$N^{"+f"{param_LU[1]:.2f}"+r"}$- LU Zerlegung", color="C1")
plt.plot(num, t_solve, "x", label=r"$t_\mathrm{solve}$", color="C2")
plt.plot(num_new, np.exp(param_solve[0]+np.log(num_new)*param_solve[1]), label=r"$N^{"+f"{param_solve[1]:.2f}"+r"}$- Solve", color="C2")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"Dimension $N$")
plt.ylabel(r"Zeit $t\, / \,\mathrm{s}$")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("output/times.pdf")