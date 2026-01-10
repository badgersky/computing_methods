import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("errors_lab10.txt")

h    = data[:, 0]
err1 = data[:, 1]
err2 = data[:, 2]
err3 = data[:, 3]

plt.figure()
plt.suptitle("Błąd metod w zależności od kroku h")
plt.loglog(h, err1, label="Błąd BME")
plt.loglog(h, err2, label="Błąd PME")
plt.loglog(h, err3, label="Błąd PMT")

plt.xlabel("h")
plt.ylabel("błąd")
plt.legend()
plt.grid(True, which="both")

plt.show()