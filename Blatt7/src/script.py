import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Set fontsize larger for plots
matplotlib.rcParams.update({"font.size": 20})

# Plotting for three different times
plt.figure(figsize=(12, 8))
plt.plot(num, t_random, "x", label=r"$t_\mathrm{random}$", color="C0")
plt.xlabel(r"Dimension $N$")
plt.ylabel(r"Zeit $t\, / \,\mathrm{s}$")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("output/times.pdf")
