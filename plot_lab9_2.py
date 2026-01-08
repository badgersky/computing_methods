import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("errors_lab9.txt")

h    = data[:, 0]
err1 = data[:, 1]
err2 = data[:, 2]

plt.figure()
plt.suptitle("Błąd metod w zależności od kraku h")
plt.loglog(h, err1, label="Błąd metoda strzałów")
plt.loglog(h, err2, label="Błąd algorytm Thomasa")

plt.xlabel("h")
plt.ylabel("błąd")
plt.legend()
plt.grid(True, which="both")

plt.show()