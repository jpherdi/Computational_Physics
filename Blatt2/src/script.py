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
    plt.savefig("bin/Bild"+str(i)+".pdf"),
