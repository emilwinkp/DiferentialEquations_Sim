import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Sistema_2DOF.csv", delimiter=",", skiprows=1)

t = data[:, 0]
y = data[:, 1]
y2 = data[:, 2]

plt.plot(t, y2, label="Masa 2")
plt.plot(t, y, label="Masa 1")

plt.xlabel("Tiempo (s)")
plt.ylabel("Desplazamiento (m)")
plt.title("Oscilación de un resorte")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
