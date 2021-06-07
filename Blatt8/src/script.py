import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Set fontsize larger for plots
matplotlib.rcParams.update({"font.size": 20})


fig = plt.figure(figsize=(15, 15))
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)
plt.plot(1, 1, ".r", mew=10, label="Minimum")

plt.xlabel(r"$x_0$")
plt.ylabel(r"$x_1$")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("output/bfgs.pdf")
