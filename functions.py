import numpy as np
import rbdl
from copy import copy


pi = np.pi

class Robot(object):
    def __init__(self, q0, dq0, ndof, dt):
        self.q = q0    # numpy array (ndof x 1)
        self.dq = dq0  # numpy array (ndof x 1)
        self.M = np.zeros([ndof, ndof])
        self.b = np.zeros(ndof)
        self.dt = dt
        self.robot = rbdl.loadModel('../urdf/robot_proyecto.urdf')

    def send_command(self, tau):
        rbdl.CompositeRigidBodyAlgorithm(self.robot, self.q, self.M)
        rbdl.NonlinearEffects(self.robot, self.q, self.dq, self.b)
        ddq = np.linalg.inv(self.M).dot(tau-self.b)
        self.q = self.q + self.dt*self.dq
        self.dq = self.dq + self.dt*ddq

    def read_joint_positions(self):
        return self.q

    def read_joint_velocities(self):
        return self.dq



def dh(d, theta, a, alpha):
    """
    Matriz de transformacion homogenea asociada a los parametros DH.
    Retorna una matriz 4x4
    """
    sth = np.sin(theta)
    cth = np.cos(theta)
    sa = np.sin(alpha)
    ca = np.cos(alpha)
    T = np.array([[cth, -ca*sth,  sa*sth, a*cth],
                  [sth,  ca*cth, -sa*cth, a*sth],
                  [0.0,      sa,      ca,     d],
                  [0.0,     0.0,     0.0,   1.0]])
    return T


def fkine(q):
    """
    Calcular la cinematica directa del robot UR5 dados sus valores articulares. 
    q es un vector numpy de la forma [q1, q2, q3, q4, q5, q6]
    """
    # Matrices DH
    T1 = dh(  q[0],      pi/2, 0.195,    0)
    T2 = dh( 0.148,      q[1], 0.410,    0)
    T3 = dh(-0.148, q[2]+pi/2,      0, pi/2)
    T4 = dh( 0.430,   q[3]+pi,      0, pi/2)
    T5 = dh( 0.148,   q[4]+pi,      0,  pi/2)
    T6 = dh( 0.111, q[5]+pi/2,      0,     0)
    # Efector final con respecto a la base
    T_inercial = np.array([[ -1,       0,       0,     0],
                           [  0,       0,       1,     0],
                           [0.0,       1,       0,     0],
                           [0.0,     0.0,     0.0,   1.0]]) 
    T = T_inercial.dot(T1).dot(T2).dot(T3).dot(T4).dot(T5).dot(T6)
    return T#, T1, T2, T3, T4, T5, T6


def jacobian_ur5(q, delta=0.0001):
    """
    Jacobiano analitico para la posicion. Retorna una matriz de 3x6 y toma como
    entrada el vector de configuracion articular q=[q1, q2, q3, q4, q5, q6]
    """
    # Crear una matriz 3x6
    J = np.zeros((3, 6))
    # Transformacion homogenea inicial (usando q)
    T_inicial = fkine(q)

    # Iteracion para la derivada de cada columna
    for i in range(6):
        # Copiar la configuracion articular inicial
        dq = copy(q)
        # Incrementar la articulacion i-esima usando un delta
        dq[i] = q[i]+delta
        # Transformacion homogenea luego del incremento (q+delta)
        T_final = fkine(dq)
        # Aproximacion del Jacobiano de posicion usando diferencias finitas
        J[:, i] = (T_final[:3, 3]-T_inicial[:3, 3])/delta
    return J


def ikine_ur5(xdes, q0):
    """
    Calcular la cinematica inversa de UR5 numericamente a partir de la configuracion articular inicial de q0. 
    Emplear el metodo de newton
    """
    epsilon = 0.001
    max_iter = 1000
    delta = 0.00001

    q = copy(q0)
    for i in range(max_iter):
        # Main loop
        J = jacobian_ur5(q, delta)
        f = fkine(q)[:3, 3]
        e = xdes-f
        q = q+np.dot(np.linalg.pinv(J), e)

        if(np.linalg.norm(e) < epsilon):
            break

    return q

def jacobian_position(q, delta=0.0001):
    # Alocacion de memoria
    J = np.zeros((3,6))
    # Transformacion homogenea inicial (usando q)
    T = fkine(q)
    # Iteracion para la derivada de cada columna
    for i in range(6):
        # Copiar la configuracion articular inicial (usar este dq para cada
        # incremento en una articulacion)
        dq = copy(q)
        # Incrementar la articulacion i-esima usando un delta
        dq[i] = q[i] + delta
        # Transformacion homogenea luego del incremento (q+dq)
        Tf = fkine(dq)
        # Aproximacion del Jacobiano de posicion usando diferencias finitas
        J[:,i] = (Tf[:3,3] - T[:3,3])/delta
    return J
 
def jacobian_pose(q, delta=0.0001):
    J = np.zeros((6,6))
    # Implementar este Jacobiano aqui
    T = fkine(q)
    x = T[:3,3]
    Q = rot2quat(T[:3,:3])
    for i in range(6):
        dq = copy(q)
        dq[i] = q[i] + delta
        Tf = fkine(dq)
        xf = Tf[:3,3]
        Qf = rot2quat(Tf[:3,:3])
        J[:3,i] = (xf - x)/delta
        J[3:,i] =(Qf - Q)/delta
    return J
 
def rot2quat(R):
    dEpsilon = 1e-6
    quat = 4*[0.,]
    quat[0] = 0.5*np.sqrt(R[0,0]+R[1,1]+R[2,2]+1.0)
    if ( np.fabs(R[0,0]-R[1,1]-R[2,2]+1.0) < dEpsilon ):
        quat[1] = 0.0
    else:
        quat[1] =0.5*np.sign(R[2,1]-R[1,2])*np.sqrt(R[0,0]-R[1,1]-R[2,2]+1.0)
    if ( np.fabs(R[1,1]-R[2,2]-R[0,0]+1.0) < dEpsilon ):
        quat[2] = 0.0
    else:
        quat[2] =0.5*np.sign(R[0,2]-R[2,0])*np.sqrt(R[1,1]-R[2,2]-R[0,0]+1.0)
    if ( np.fabs(R[2,2]-R[0,0]-R[1,1]+1.0) < dEpsilon ):
        quat[3] = 0.0
    else:
        quat[3] =0.5*np.sign(R[1,0]-R[0,1])*np.sqrt(R[2,2]-R[0,0]-R[1,1]+1.0)
    return np.array(quat)
 
def TF2xyzquat(T):
    quat = rot2quat(T[0:3,0:3])
    res = [T[0,3], T[1,3], T[2,3], quat[0], quat[1], quat[2], quat[3]]
    return np.array(res)
 
def skew(w):
    R = np.zeros([3,3])
    R[0,1] = -w[2]; R[0,2] = w[1]
    R[1,0] = w[2]; R[1,2] = -w[0]
    R[2,0] = -w[1]; R[2,1] = w[0]
    return R


def ikine(xdes, q0):
    # Error
    epsilon = 0.001
    # Maximas iteraciones
    max_iter = 1000
    # Delta de la jacobiana
    delta = 0.00001
    # Copia de las articulaciones
    q = copy(q0)
    # Almacenamiento del error
    ee = []
    # Transformacion homogenea (usando q)
    To = fkine(q)
    To = To[0:3, 3]  # vector posicion
    # Resetear cuando se llega a la cantidad maxima de iteraciones
    restart = True
 
    while restart:
        for i in range(max_iter):
            # Hacer el for 1 vez
            restart = False
            # Pseudo-inversa del jacobiano
            J = jacobian_position(q)
            J = np.linalg.pinv(J)
            # Error entre el x deseado y x actual
            e = xdes - To
            # q_k+1
            q = q + np.dot(J, e)
            # Nueva mtransformada homogenea
            To = fkine(q)
            To = To[0:3, 3]  # vector posicion
 
            # Norma del error
            enorm = np.linalg.norm(e)
            ee.append(enorm)    # Almacena los errores
            # Condicion de termino
            if (enorm < epsilon):
                print("Error en la iteracion ", i, ": ", np.round(enorm, 4))
                break
            if (i == max_iter-1 and enorm > epsilon):
                print("Iteracion se repite")
                print("Error en la iteracion ", i, ": ", enorm)
            restart = True
    return q,ee
 
def ik_gradient_sawyer(xdes, q0):
    """
    Calcular la cinematica inversa de UR5 numericamente a partir de la
    configuracion articular inicial de q0.
    Emplear el metodo gradiente
    """
    epsilon = 0.001
    max_iter = 1000
    delta = 0.00001
    alpha = 0.5
    q = copy(q0)
    for i in range(max_iter):
        # Main loop
        J = jacobian_ur5(q,delta)
        f = fkine(q)[:3,3]
        e = xdes - f
        q = q + alpha*np.dot(J.T,e)
        if (np.linalg.norm(e)<epsilon):
            break
    return q
