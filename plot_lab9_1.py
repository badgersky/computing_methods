import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("results_lab9.txt")

x  = data[:, 0]
y1 = data[:, 1]
y2 = data[:, 2]
y3 = data[:, 3]

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(6, 8))
fig.suptitle("Porównanie metod numerycznych z rozwiązaniem analitycznym")

axs[0].plot(x, y1)
axs[0].set_title("Metoda strzałów")
axs[0].set_ylabel("y")
axs[0].grid(True)

axs[1].plot(x, y2)
axs[1].set_title("Algorytm Thomasa")
axs[1].set_ylabel("y")
axs[1].grid(True)

axs[2].plot(x, y3)
axs[2].set_title("Rozwiązanie analityczne")
axs[2].set_xlabel("x")
axs[2].set_ylabel("y")
axs[2].grid(True)

plt.tight_layout()
plt.show()