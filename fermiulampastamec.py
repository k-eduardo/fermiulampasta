# -*- coding: utf-8 -*-
#
# Este código reproduce los resultados que reportaron Fermi, Ulam y Pasta
# un poco antes de la muerte de Fermi.
# También intentaré calcular los periodos de oscilación y ver cómo cambian
# estos con el tamaño del sistema.
#

from math import sin,pi
from numpy import *	
from progressbar import ProgressBar
from visual.graph import *
import time
import sys
import datetime
import subprocess
import os

p = ProgressBar('cyan', width=20, block='▣', empty='□')

##################### SISTEMA
# Tamaño del sistema
N = 128
n = N + 2
# Pasos de integración
Ti = 1000000
# Tamaño de la perturbación
epsilon = 0.1
# Masa del sistema
m = 1.
# Constante de los resortes (no se usa en este programa)
k = 1.
##################### SISTEMA



##################### PARÁMETROS DE INTEGRACIÓN
# Paso temporal
h = 0.1
# ti = pasos de integración
ti = 0.
# t = tiempo físico en segundos
t = 0.
# Energía por modo (variable temporal)
E = []
for i in range(0,n):
	E.append(0.)
# Variable de escritura (variable temporal)
salida = ''
# Vector de coordenadas normales
a = []
for i in range(0,n):
	a.append(0.)
#Vector de velocidades normales
av = []
for i in range(0,n):
	av.append(0.)
##################### PARÁMETROS DE INTEGRACIÓN


##################### Cada cuántos pasos de calcula la energía de los modos
modimpr = 1000
##################### Cada cuántos pasos de calcula la energía de los modos

##################### ARCHIVO DE ESCRITURA
archivo = open('.datos128','w')
prueba = open('.prueba128','w')
#prueba2 = open('.prueba128','w')
##################### ARCHIVO DE ESCRITURA


#################### ESTADO INICIAL
# Vector de estado
# (las primeras n coordenadas son posición)
# (de la n+1 a la 2n son las de velocidad)
x = []
for i in range(0,n):
		x.append(sin(2.*pi*(float(i))/float(N+1)))
# Los extremos están fijos
x[0] = 0.
x[n-1] = 0.

# Velocidades iniciales
for i in range(0,n):
	x.append(0.)
x = array(x)
#################### ESTADO INICIAL

#for i in range(0,len(x)):
#	prueba2.write(str(i)+' '+str(x[i])+'\n')
#prueba2.close()

#################### FUNCIÓN DE FUERZA
def F(t,x):
	y = []
	y.append(0.)
	for i in range(1,n - 1):
		auxiliar=x[i + n]
		y.append(auxiliar)
	y.append(0.)
	y.append(0.)
	for k in range(1,n - 1):
		a = (x[k+1]-x[k]) + (x[k-1]-x[k])
		b = epsilon * ((x[k+1] - x[k])**2 - (x[k] - x[k-1])**2)
		c = a+b
		y.append(c)
	y.append(0.)
	fuerza = array(y)
	return fuerza
#################### FUNCIÓN DE FUERZA

#################### RUNGE KUTTA 4
def RK4(t_init,x_init,F,h):
	state = x_init
	t = t_init
	k1 = F(t,state)
	k2 = F(t + .5*h,state + .5*h*k1)
	k3 = F(t + .5*h,state + .5*h*k2)
	k4 = F(t + h,state + h*k3)
	phi = (1/6.)*k1 + (1/3.)*k2 + (1/3.)*k3 + (1/6.)*k4
	state = state+h*phi
	return state
#################### RUNGE KUTTA 4

#################### INTEGRACIÓN
while (ti<Ti):
	################ Barra de progereso
	E2 = time.clock()
	test = E2*float(Ti-ti)/float(ti+1)
	p.render(int((ti+1)*100/Ti),'\nTiempo restante %.1f min.'%(test/60.))
	################ Barra de progereso
	x = RK4(t,x,F,h)
	if (ti%modimpr == 0):
		############ Cálculo de los modos normales
		salida = ''
		for i in range(0,n):
			a[i] = 0.
			av[i] = 0.
			for j in range(1,n-1):
				a[i] += x[j]*sin(float(i+1)*float(j)*pi/float(N))
				av[i] += x[j+n]*sin(float(i+1)*float(j)*pi/float(N))
			E[i] = (av[i]**2)/2. + 2*(a[i]**2)*sin(pi*float(i+1)/float(2*N))**2
		salida = str(t)+' '+str(E[0])
		for i in range(1,n):
			salida += ' '+str(E[i])
		salida += '\n'
		archivo.write(salida)
		############ Cálculo de los modos normales
	################ Se escribe un modo de prueba
#	prueba.write(str(t)+' '+str(x[0])+' '+str(x[2])+' '+str(x[n-1])+'\n')
	################ Se escribe un modo de prueba
	ti += 1 ######## Actualización del tiempo de integración
	t += h  ######## Actualización del tiempo
#################### INTEGRACIÓN
