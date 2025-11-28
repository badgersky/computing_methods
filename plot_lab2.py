import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("wyniki_porownania_lab_2.txt")

logx = data[:, 0]
logerr = data[:, 1]

plt.plot(logx, logerr)
plt.xlabel("log10(x)")
plt.ylabel("log10(błąd względny)")
plt.title("Porównanie błędu względnego funkcji")
plt.grid(True)
plt.show()