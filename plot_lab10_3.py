import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("bme_unstable_lab10.txt")

t = data[:, 0]
bme = data[:, 1]
ana = data[:, 2]

plt.figure()
plt.suptitle("BME z niestabilnym krokiem")
plt.plot(t, bme, marker='o', linestyle='None', label="BME")
plt.plot(t, ana, label="EXACT")

plt.xlabel("t")
plt.ylabel("y(t)")
plt.legend()
plt.grid(True, which="both")

plt.show()