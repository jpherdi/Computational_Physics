import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
# Set fontsize larger for latex plots
matplotlib.rcParams.update({'font.size': 20})

# Generate data from file
x, y = np.genfromtxt("bin/python_Aufgabe2.txt", unpack=True)
m, n = x[-1], y[-1]


# Plotting
plt.figure(figsize=(12,7))
plt.grid()
plt.xlabel("x")
plt.ylabel("y")
x_new = np.linspace(min(x)-x[:-1].std()/2, max(x)+x[:-1].std()/2)
plt.plot(x[:-1], y[:-1], "x", mew=2., alpha=2, label="Datenpunkte")
plt.plot(x_new, m*x_new+n, "-", linewidth=3, label="Ausgleichsgerade")
plt.legend()
plt.tight_layout()
plt.savefig("bin/figure.pdf", dpi=1200)
