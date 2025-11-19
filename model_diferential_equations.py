''' El siguiente codigo define el calculo matemático de un sistema masa resorte de dos grados de livertad. '''

''' ----------------- Autores -> Emil Winkler Partida, David Alfonso Bueno Nachon ------------------------ '''


# Importamos librerias para computacion numerica y graficacion
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Masas del sistema
masa_1 = 200/1000 # Peso del aro intermedio en kilogramos
masa_2 = 60/1000 # Peso de la lampara en kilogramos

# Constantes elasticas de las ligas obtenidas a partir de experimentacion
k1 = 4.5  # Constante elastica de la liga superior en N
k2 = 3.5  # Constante elastica de la liga inferior en N

# Resistencias de amortiguamiento
resistencia_1 = 0.5  # Resistencia de amortiguamiento de la masa 1 en Ns/m
resistencia_2 = 0.3  # Resistencia de amortiguamiento de la masa 2 en Ns/m

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

# Longitudes de los resortes 
L1 = 0.09
L2 = 0.09

# Dimensiones para tabla horizontal (masa 1)
ancho_tabla = 0.18
alto_tabla = 0.03

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-0.4, 0.4)
ax.set_ylim(-0.4, 0.2)
ax.set_title("Animación del sistema masa-resorte-amortiguado")
ax.set_xlabel("Eje X")
ax.set_ylabel("Altura (m)")
ax.grid(True)
ax.set_aspect('equal')

# Dibujar techo
ax.axhline(y=0, color = 'black', linewidth=4)

def dibujar_resorte(x0,y0,x1,y1,width, numero_de_curvas = 5):
    dx = x1 - x0
    dy = y1 - y0
    distancia = np.sqrt(dx**2 + dy**2)

    # Puntos intermedios del resorte
    t = np.linspace(0,1,100)
    oscilacion = width * np.sin(2* np.pi * numero_de_curvas * t)

    # Direccion perpendicular 
    nx = -dy/distancia
    ny = dx/distancia

    # Calcular puntos del resorte
    x_resorte = x0  + dx * t + nx * oscilacion
    y_resorte = y0 + dy * t + ny * oscilacion
    return x_resorte, y_resorte

# Elementos de la animacion 
resorte1_linea, = ax.plot([],[], 'b-', lw = 1.5)
masa1_tabla = plt.Rectangle((-ancho_tabla/2, 0), ancho_tabla, alto_tabla, fc="brown", ec="black", zorder=10)
resorte2_linea, = ax.plot([],[], 'g-', lw = 1.5)
masa2_circulo = plt.Circle((0,0), 0.04, fc = "blue", ec = "black")

# Masas en la grafica
ax.add_patch(masa1_tabla)
ax.add_patch(masa2_circulo)

# Etiquetas de las masas
etiqueta_m1 = ax.text(0, 0, 'Masa 1', ha='center', va='center', color='gray', fontweight='bold', zorder=11)
etiqueta_m2 = ax.text(0, 0, 'Masa 2', ha='center', va='center', color='gray', fontweight='bold', zorder=11)

# Inicializacion
def init():
    resorte1_linea.set_data([],[])
    resorte2_linea.set_data([],[])
    masa1_tabla.set_xy((-ancho_tabla/2, -alto_tabla/2))
    masa2_circulo.center = (0,-L1 -L2)
    etiqueta_m1.set_position((0, -alto_tabla/2))
    etiqueta_m2.set_position((0, -L1 -L2))
    return resorte1_linea, resorte2_linea, masa1_tabla, masa2_circulo, etiqueta_m1, etiqueta_m2

# Actualizacion para animacion 

def actualizar(frame):
    factor_visual = 2
    x1 = x1_valores[frame] * factor_visual
    x2 = x2_valores[frame] * factor_visual

    # Posiciones verticales
    y_techo = 0
    y_masa1 = y_techo - L1 - x1
    y_masa2 = y_masa1 - L2 - (x2 - x1)

    # Actualizar resortes
    resorte1_x, resorte1_y = dibujar_resorte(0, y_techo, 0, y_masa1, 0.025, 8)
    resorte1_linea.set_data(resorte1_x, resorte1_y)
    
    resorte2_x, resorte2_y = dibujar_resorte(0, y_masa1, 0, y_masa2, 0.025, 8)
    resorte2_linea.set_data(resorte2_x, resorte2_y)

    # Actualizar masas 
    masa1_tabla.set_xy((-ancho_tabla/2, y_masa1 - alto_tabla/2))
    masa2_circulo.center = (0, y_masa2)

    etiqueta_m1.set_position((0, y_masa1))
    etiqueta_m2.set_position((0, y_masa2))
    return resorte1_linea, resorte2_linea, masa1_tabla, masa2_circulo, etiqueta_m1, etiqueta_m2

# Animacion 

animacion = FuncAnimation(fig, actualizar, frames= range(0, len(valores_t), 5), init_func=init, blit=True, interval=20, repeat=True)
plt.tight_layout()
plt.show()