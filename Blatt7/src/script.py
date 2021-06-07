import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Set fontsize larger for plots
matplotlib.rcParams.update({"font.size": 20})

a_x, a_y = np.genfromtxt("output/output_a.txt", unpack=True)
b_x, b_y = np.genfromtxt("output/output_b.txt", unpack=True)
c_x, c_y = np.genfromtxt("output/output_c.txt", unpack=True)
# Plotting for three different times

x_min, x_max = -2, 2
y_min, y_max = -2, 2

x = np.linspace(x_min, x_max, 1000)
y = np.linspace(y_min, y_max, 1000)
x, y = np.meshgrid(x, y)
z = (1 - x) ** 2 + 100 * (y - x ** 2) ** 2

fig = plt.figure(figsize=(15, 15))
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
contour_ = plt.contourf(x, y, z, vmin=0, vmax=4000)
plt.plot(a_x, a_y, "x", label=r"Inverse Hesse-Matrix", color="C0")
plt.plot(b_x, b_y, "x", label=r"Diag. Inv. Hesse-Matrix", color="C1")
plt.plot(c_x, c_y, "x", label=r"Einheitsmatrix", color="C2")
plt.plot(1, 1, ".r", mew=10, label="Minimum")

plt.xlabel(r"$x_0$")
plt.ylabel(r"$x_1$")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
fig.colorbar(contour_)
plt.savefig("output/bfgs.pdf")
