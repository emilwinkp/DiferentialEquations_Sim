''' El siguiente codigo define el calculo matemático de un sistema masa resorte de dos grados de livertad. '''

''' ----------------- Autores -> Emil Winkler Partida, David Alfonso Bueno Nachon ------------------------ '''


# Importamos librerias para computacion numerica y graficacion
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Masas del sistema
masa_1 = 40/1000 # Peso del aro intermedio en kilogramos
masa_2 = 30/1000 # Peso de la lampara en kilogramos

# Constantes elasticas de las ligas obtenidas a partir de experimentacion
k1 = 4.5  # Constante elastica de la liga superior en N
k2 = 3.5  # Constante elastica de la liga inferior en N

# Resistencias de amortiguamiento
resistencia_1 = 0.5  # Resistencia de amortiguamiento de la masa 1 en Ns/m
resistencia_2 = 0.3  # Resistencia de amortiguamiento de la

def Derivadas(t, y):
    x1, x2, v1, v2 = y
    dx1 = v1
    dx2 = v2
    dv1 = ( -resistencia_1  *v1 - k1*x1 - k2*(x1 - x2))/masa_1
    dv2 = ( -resistencia_2  *v2 - k2*(x2 - x1))/masa_2
    return np.array([dx1, dx2, dv1, dv2])

#Funcion para el metodo numerico Runge Kutta
def Runge_Kutta_4(f,t,y,h):
    k1 = f(t,y)
    k2 = f(t + h/2, y + h/2 * k1)
    k3 = f(t + h/2, y + h/2 * k2)
    k4 = f(t + h, y + h * k3)
    return y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

# Condiciones iniciales
y0 = np.array([0.01,0.0,0.0,0.0]) # [x1, x2, v1, v2]
t0 = 0.0 
tf = 10.0 
dt = 0.001 

# Vector para almacenar suluciones
N = int((tf - t0)/dt) # Cantidad de iteraciones
valores_t = np.linspace(t0, tf, N)
valores_y = np.zeros((N,4))
valores_y[0] = y0 # Asigno valor inicial

for i in range(1,N):
    valores_y[i] = Runge_Kutta_4(Derivadas, valores_t[i-1], valores_y[i-1], dt)

x1_valores = valores_y[:, 0]
x2_valores = valores_y[:, 1]


# Graficacion 
plt.figure(figsize=(10, 5))
plt.plot(valores_t, x1_valores, label="x1 (masa intermedia)")
plt.plot(valores_t, x2_valores, label="x2 (lampara)")
plt.xlabel("Tiempo (s)")
plt.ylabel("Desplazamiento (m)")
plt.title("Sistema de 2 DOF masa-resorte amortiguado")
plt.grid(True)
plt.legend()
plt.tight_layout()

''' Animacion del sistema masa resorte amortiguado '''

#Por Checar 

# fig, ax = plt.subplots(figsize=(5, 6))
# ax.set_xlim(-0.2, 0.2)
# ax.set_ylim(-0.25, 0.05)
# ax.set_title("Animación del sistema masa-resorte-amortiguado")
# ax.set_xlabel("Eje X")
# ax.set_ylabel("Altura (m)")
# ax.grid()

# # Dibujar elementos iniciales
# support, = ax.plot([0], [0], marker='s', markersize=10, color='brown')  # soporte superior
# line1, = ax.plot([], [], 'k-', lw=2)  # resorte 1
# mass1, = ax.plot([], [], 'ro', markersize=10)  # masa intermedia
# line2, = ax.plot([], [], 'k-', lw=2)  # resorte 2
# mass2, = ax.plot([], [], 'bo', markersize=12)  # botella PET

# # Función de actualización para la animación
# def update(frame):
#     x1 = x1_valores[frame]
#     x2 = x2_valores[frame]
    
#     # Coordenadas verticales
#     y_support = 0
#     y_mass1 = y_support - x1
#     y_mass2 = y_support - x2
    
#     # Actualizar resortes y masas
#     line1.set_data([0, 0], [y_support, y_mass1])
#     mass1.set_data([0], [y_mass1])
    
#     line2.set_data([0, 0], [y_mass1, y_mass2])
#     mass2.set_data([0], [y_mass2])
    
#     return line1, mass1, line2, mass2

# # Crear animación
# ani = FuncAnimation(fig, update, frames=len(valores_t), interval=10, blit=True)

plt.show()

    