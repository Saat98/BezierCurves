# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 13:56:43 2018

@author: SaulAlvarez
"""
import numpy as np
import matplotlib as plt

def Derivada(f, t_0):
    h = 1e-6
    
    return (f(t_0 + h) - f(t_0 - h)) / (2 * h)

def Curvatura(f, t_0):
    def e_1(t):
        f_prime_t = Derivada(f, t)
        return f_prime_t / np.linalg.norm(f_prime_t, 2)
    
    return np.linalg.norm(e_1(t_0), 2)
    

class polynomial:
    '''
    Irrelevante por el momento
    '''

    def __init__(self, polinomio = [0], variable = 'x'):
        self.vector = np.array(polinomio)
        self.deg = len(self.vector) - 1
        self.variable = variable
        
    def __add__(self, other):
        propio, otro = self.vector, other.vector
        propio.extend([0 for i in range(max(0, len(otro) - len(propio)))])
        otro.extend(  [0 for i in range(max(0, len(propio) - len(otro)))])
        
        return polynomial(propio + otro, self.variable)

    def __str__(self):
        ret = ''
        
        for i,elem in enumerate(self.vector):
            if elem >= 0:
                ret = '%s + %s%s^%s' %(ret,  elem, self.variable, i)
            else:
                ret = '%s - %s%s^%s' %(ret, -elem, self.variable, i)
                
        return ret
    
    def __multiply__(self, other):
        ret = polynomial(variable = self.variable)
        
        for i,elem in enumerate(other.vector):
            cur = [0 for j in range(i)]
            cur.extend(elem * self.vector)
            ret = ret + polynomial(cur,self.variable)
            
        return ret
    
    def evaluate(self,t):
        ret = 0
        
        for elem in self.vector[::-1]:
            ret = ret*t + elem
            
        return ret
       

def grafica_bezier_plano(puntos, num_divs = 1000):
    '''
    Grafica la curva de Bézier en geometría plana cuyos puntos de referencia
    son puntos, en ese orden.
    El argumento num_divs representa la división del intervalo temporal. Debe
    ser un número natural.
    
    Procurar que el número de divisiones multiplicado por el cuadrado del número
    de puntos de referencia sea pequeño (menor a 1e7)
    '''

    num_puntos = len(puntos)
    divs = [i/num_divs for i in range(num_divs + 1)]
    
    partial_beziers = [[None for i in range(num_puntos)] for j in range(num_puntos)]
    for i in range(num_puntos):
        partial_beziers[i][i] = [puntos[i] for t in divs]
    
    for i in range(1, num_puntos):
        for j in range(num_puntos-i):
            partial_beziers[j][j + i] = [np.array([0, 0]) for t in divs]
            
            for t in range(num_divs + 1):
                partial_beziers[j][j + i][t] = (1 - divs[t])*partial_beziers[j][j + i - 1][t] + \
                divs[t]*partial_beziers[j + 1][j + i][t]
    
    x = [partial_beziers[0][num_puntos - 1][i][0] for i in range(num_divs + 1)]
    y = [partial_beziers[0][num_puntos - 1][i][1] for i in range(num_divs + 1)]
    
    plt.pyplot.plot(x, y)
    
def generar_puntos_naive(func, num_puntos, left = 0, right = 1):
    '''
    Genera un vector con num_puntos+1 puntos, que corresponden a las evaluaciones
    de func en t=left+(right-left)*i/num_puntos, con i=0,1,...,num_puntos.
    '''
    
    return [func(left + (right - left) * i / num_puntos) for i in range(num_puntos + 1)]

def generar_puntos_curvatura(func, num_puntos, left = 0, right = 1):
    '''
    Genera un vector con num_puntos+1 puntos, que corresponden a las evaluaciones
    de func en t=z_i, donde \int_{z_{i-1}}^{z_i} \kappa(func(t))dt es aproximadamente
    constante.
    '''
    
    divs = [(left + (right - left) * i) / 10000 for i in range(10001)]
    kappas = [Curvatura(func,i) for i in divs]
    kappas[0] = 0
    
    kappasac=[kappas[0]]
    for i in range(1,10001):
        kappasac.append(kappas[i] + kappasac[i-1])
    
    total = kappasac[-1]
    corte = total / num_puntos
    tiempos = [left]
    
    for i in range(10001):
        if len(tiempos)*corte <= kappasac[i]:
            tiempos.append(divs[i])
            
    if tiempos[-1]-1e-9 <= right:
        tiempos.append(right)
        
    return [func(t) for t in tiempos]

def grafica_original(func, left = 0, right = 1):
    '''
    Genera dos vectores que se usen para graficar func entre t=left y t=right
    '''
    
    x = [func(left + (right - left) * t / 10000)[0] for t in range(10001)]
    y = [func(left + (right - left) * t / 10000)[1] for t in range(10001)]
    
    plt.pyplot.plot(x, y)

def circulo(t):
    '''
    Da puntos en el circulo unitario
    '''
    
    return np.array([np.cos(2 * t * np.pi), np.sin(2 * t * np.pi)])

def hiperbola(t):
    '''
    Da puntos en la hiperbola xy=1. t debe ser distinto a 0
    '''
    
    return np.array([1 / t, t])

def seno(t):
    '''
    Da puntos en la grafica de y=sin(x)
    '''
    
    return np.array([t, np.sin(t)])

def contraejemplo(t):
    '''
    Da puntos en la grafica de y=x*sin(1/x)
    '''
    
    if t == 0:
        return np.array([0.,0.])
    
    return np.array([t,t * np.sin(1 / t)])

if __name__ == '__main__':
    func = contraejemplo
    left = 0.
    right = 0.1
    num_puntos = 50
    
#    grafica_bezier_plano(generar_puntos_naive(func, num_puntos, left, right))
#    grafica_bezier_plano(generar_puntos_curvatura(func, num_puntos, left, right))
#    grafica_original(func, left, right)