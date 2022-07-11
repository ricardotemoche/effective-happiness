import rbdl
import numpy as np


# Lectura del modelo del robot a partir de URDF (parsing)
modelo = rbdl.loadModel('../urdf/robot_proyecto.urdf')
# Grados de libertad
ndof = modelo.q_size

# Configuracion articular
q = np.array([0.1, 0.2, 0.3, 0.8, 0.5, 0.6])
# Velocidad articular
dq = np.array([0.05, 0.8, 0.8, 0.6, 0.9, 0.9])
# Aceleracion articular
ddq = np.array([0, 0, 0.04, 0.3, 1.0, 0.5])

# Arrays numpy
zeros = np.zeros(ndof)          # Vector de ceros
tau = np.zeros(ndof)          # Para torque
g = np.zeros(ndof)          # Para la gravedad
c = np.zeros(ndof)          # Para el vector de Coriolis+centrifuga
M = np.zeros([ndof, ndof])  # Para la matriz de inercia
e = np.eye(6)               # Vector identidad

# Torque dada la configuracion del robot
rbdl.InverseDynamics(modelo, q, dq, ddq, tau)

# Parte 1: Calcular vector de gravedad, vector de Coriolis/centrifuga,
# y matriz M usando solamente InverseDynamics

# El vector de gravedad es:
rbdl.InverseDynamics(modelo, q, zeros, zeros, g)
g = np.round(g, 2)
# El vector de fuerza centrifuga y Coriolis es:
rbdl.InverseDynamics(modelo, q, dq, zeros, c)
c = np.round(c-g, 2)
# La matriz de inercia es:
for i in range(6):
    rbdl.InverseDynamics(modelo, q, zeros, e[i, :], M[i, :])
    M[i, :] = M[i, :]-g
M = np.round(M, 2)
Mtrans=np.transpose(M)
print("El vector de gravedad es:\n{} \n".format(g))
print("El vector de fuerza centrifuga y Coriolis es:\n{}\n".format(c))
print("La matriz de inercia es:\n {} \n".format(M))
print("La matriz de inercia transpuesta es:\n {} \n".format(Mtrans))


# Parte 2: Calcular M y los efectos no lineales b usando las funciones
# CompositeRigidBodyAlgorithm y NonlinearEffects. Almacenar los resultados
# en los arreglos llamados M2 y b2
b2 = np.zeros(ndof)          # Para efectos no lineales
M2 = np.zeros([ndof, ndof])  # Para matriz de inercia


rbdl.CompositeRigidBodyAlgorithm(modelo, q, M2)
M2 = np.round(M2, 2)
rbdl.NonlinearEffects(modelo, q, dq, b2)
b2 = np.round(b2, 2)
print("Parte 2 \n")
print("La matriz de inercia es:\n {} \n".format(M2))
print("El vector de efectos no lineales es:\n {}\n".format(b2))

# Parte 3: Verificacion de valores

print("Parte 3\n")
print("\nComprobacion por resta de las matrices de inercia obtenidas por ambos metodos\n")
print(M-M2)
print("\nComprobacion por resta de [c(q,dq) + g(d) - b(q,dq)]:\n")
print(c+g-b2)

# Parte 4: Verificacion de la expresion de la dinamica
print("\nParte 4\n")
tau2 = np.round(M.dot(ddq) + c + g, 2)
tau3 = np.round(M2.dot(ddq) + b2, 2)
print("\nVector de torques obtenidos con la funcion InverseDynamics directamente\n")
print(np.round(tau, 2))
print("\nVector de torques obtenidos con M+c+g obtenidos de la funcion InverseDynamics\n")
print(tau2)
print("\nVector de torques obtenidos con M2 y b2 las funciones CompositeRigidBodyAlgorithm NonlinearEffects\n")
print(tau3)
