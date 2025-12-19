import numpy as np
import matplotlib.pyplot as plt

def moving_average(y, window_size=5):
    return np.convolve(y, np.ones(window_size)/window_size, mode='valid')

double_data = np.loadtxt('double_errors.txt')
long_double_data = np.loadtxt('longdouble_errors.txt')

labels = ['progresywna 2', 'wsteczna 2', 'progresywna 3', 'wsteczna 3', 'centralna 3']

plt.figure(figsize=(10, 6))

for i in range(1, 6):
    y_smooth = moving_average(double_data[:, i])
    x_smooth = double_data[:len(y_smooth), 0]
    # y = double_data[:, i]
    # x = double_data[: len(y), 0]
    plt.loglog(x_smooth, y_smooth, label=f'DOUBLE {labels[i-1]}')
    # plt.loglog(x, y, label=f'DOUBLE {labels[i-1]}')

for i in range(1, 6):
    y_smooth = moving_average(long_double_data[:, i])
    x_smooth = long_double_data[:len(y_smooth), 0]
    # y = long_double_data[:, i]
    # x = long_double_data[: len(y), 0]
    plt.loglog(x_smooth, y_smooth, label=f'LONG DOUBLE {labels[i-1]}', linestyle='--')
    # plt.loglog(x, y, label=f'LONG DOUBLE {labels[i-1]}', linestyle='--')

plt.xlabel('Log10(h)')
plt.ylabel('Log10(err)')
plt.title('Błędy bezwzględne funkcji różnicowej do wielkości kroku')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.show()