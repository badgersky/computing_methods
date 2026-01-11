import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("bme_unstable_lab10.txt")
data2 = np.loadtxt("bme_almost_stable_lab10.txt")

t1 = data[:, 0]
bme_u = data[:, 1]
bme_b = data2[:, 1]
t2 = data2[:, 0]
ana = data[:, 2]

plt.figure()
plt.suptitle("BME niestabilne i na granicy stabilności")
plt.plot(t1, bme_u, marker='o', linestyle='None', label="BME niestabilne")
plt.plot(t1, ana, label="EXACT")
plt.plot(t2, bme_b, marker='o', linestyle='None', label="BME na granicy")

plt.xlabel("t")
plt.ylabel("y(t)")
plt.legend()
plt.grid(True, which="both")

plt.show()