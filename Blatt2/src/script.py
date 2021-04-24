import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
# Set fontsize larger for latex plots
matplotlib.rcParams.update({'font.size': 20})

# Generate data from file
for i in ["",10,20,50]:
    a = np.genfromtxt("bin/data"+str(i)+".txt", unpack=True)
    a.reshape(512,512)
    # Plotting
    plt.figure(figsize=(12,12))
    plt.imshow(a, cmap="gray")
    plt.tight_layout()
    plt.savefig("bin/Bild"+str(i)+".pdf")

t_random, t_LU, t_solve = np.genfromtxt("bin/times.txt", unpack=True)

num = np.arange(1,len(t_random)+1)

plt.figure(figsize=(12,8))
plt.plot(num, t_random, "x", label=r"$t_\mathrm{random}$")
plt.plot(num, t_LU, "x", label=r"$t_\mathrm{LU}$")
plt.plot(num, t_solve, "x", label=r"$t_\mathrm{solve}$")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"Size of Vector $N$")
plt.ylabel(r"time $t / \mathrm{s}$")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("bin/times.pdf")