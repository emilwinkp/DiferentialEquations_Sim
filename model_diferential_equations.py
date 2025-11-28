''' El siguiente codigo define el calculo matemático de un sistema masa resorte de dos grados de livertad. '''

''' ----------------- Autores -> Emil Winkler Partida, David Alfonso Bueno Nachon ------------------------ '''


# Importamos librerias para computacion numerica y graficacion
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
from matplotlib.animation import FFMpegWriter
import csv

# Masas del sistema
masa_1 = 74/1000 # Peso del aro intermedio en kilogramos
masa_2 = 271/1000 # Peso de la lampara en kilogramos

# Constantes elasticas de las ligas obtenidas a partir de experimentacion
k1 = 80.23  # Constante elastica de la liga superior en N
k2 = 82.36 # Constante elastica de la liga inferior en N

# Resistencias de amortiguamiento
resistencia_1 = 0.02     # Resistencia de amortiguamiento de la masa 1 en Ns/m
resistencia_2 = 0.02  # Resistencia de amortiguamiento de la masa 2 en Ns/m

radio = 0.03

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
y0 = np.array([0.01,0.02,0.0,0.0]) # [x1, x2, v1, v2] Modifican condiciones iniciales dependiendo de lo requerido
t0 = 0.0 
tf = 50.3 # Para poder comparar con el experimento
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

# Transcribir datos obtenidos en un archivo csv

with open("datos_sistema_2DOF.csv", mode="w", newline="") as archivo:
    escritor = csv.writer(archivo)
    
    # Encabezados
    escritor.writerow(["t", "x1", "x2", "v1", "v2"])
    
    # Guardar datos
    for i in range(N):
        escritor.writerow([
            valores_t[i],
            valores_y[i, 0],   # x1
            valores_y[i, 1],   # x2
            valores_y[i, 2],   # v1
            valores_y[i, 3]    # v2
        ])

print("Archivo CSV 'datos_sistema_2DOF.csv' creado.")

# Graficacion 
plt.figure(figsize=(10, 5))
plt.plot(valores_t, x1_valores, label="x1 (masa intermedia)", linewidth=0.6)
plt.plot(valores_t, x2_valores, label="x2 (lampara)",linewidth=1)
plt.xlabel("Tiempo (s)")
plt.ylabel("Desplazamiento (m)")
plt.title("Sistema de 2 DOF masa-resorte amortiguado")
plt.grid(True)
plt.legend()
plt.tight_layout()

''' Animacion del sistema masa resorte amortiguado '''

# Longitudes de los resortes en metros
L1 = 0.064 #En reposo 0.046
L2 = 0.06 #En reposo 0.031

# Dimensiones para tabla horizontal (masa 1) en metros
ancho_tabla = 0.1
alto_tabla = 0.02

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
masa2_circulo = plt.Circle((0,0), radio , fc = "blue", ec = "black")

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
    
    resorte2_x, resorte2_y = dibujar_resorte(0, y_masa1, 0, y_masa2 + radio, 0.025, 8)
    resorte2_linea.set_data(resorte2_x, resorte2_y)

    # Actualizar masas 
    masa1_tabla.set_xy((-ancho_tabla/2, y_masa1 - alto_tabla/2))
    masa2_circulo.center = (0, y_masa2)

    etiqueta_m1.set_position((0, y_masa1))
    etiqueta_m2.set_position((0, y_masa2))
    return resorte1_linea, resorte2_linea, masa1_tabla, masa2_circulo, etiqueta_m1, etiqueta_m2

# Animacion 

writer = FFMpegWriter(fps=30)
tiempo_max_animacion = 30  # segundos
N_anim = int(tiempo_max_animacion / dt)

frames = range(0, N_anim, 5)
animacion = FuncAnimation(fig, actualizar, frames, init_func=init, blit=True, interval=20, repeat=True)
plt.tight_layout()
writer = PillowWriter(fps=30)
#animacion.save("sistema_masa_resorte.gif", writer=writer)

plt.show()


