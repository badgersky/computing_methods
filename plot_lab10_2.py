import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("results_lab10.txt")

t = data[:, 0]
bme = data[:, 1]
pme = data[:, 2]
pmt = data[:, 3]
ana = data[:, 3]

plt.figure()
plt.suptitle("Porównanie metod z wynikiem analitycznym")
plt.plot(t, bme, marker='o', linestyle='None', label="BME")
plt.plot(t, pme, marker='o', linestyle='None', label="PME")
plt.plot(t, pmt, marker='o', linestyle='None', label="PMT")
plt.plot(t, ana, label="EXACT")

plt.xlabel("t")
plt.ylabel("y(t)")
plt.legend()
plt.grid(True, which="both")

plt.show()