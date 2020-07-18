import struct
import numpy as np
import pickle
import csv
import os
import csv

from os import listdir
SalidaSimulacion=[]
Dat=[]
for out in listdir("./DatosBruto"):
	A=out
	Dat.append(A)
print('Salidas de simulación encontradas:')
Dat.sort()
print(Dat)
for out in range(len(Dat)):
	SalidaSimulacion.append('DatosBruto/'+Dat[out])
print(SalidaSimulacion)



try:
	os.mkdir('d5004')
except FileExistsError:
	print('La carpeta d5004 ya existe')


from VolI import *
s=len(SalidaSimulacion)
for i in range(s):
	print('Salida en lectura: ' +Dat[i])
	root=Dat[i]
	try:
		Ruta='d5004/'+root
		os.mkdir(Ruta)
		volcado(Ruta,SalidaSimulacion[i])
	except FileExistsError:
		print('Fichero creado. Asegúrese que contiene 3 Búfers y 1 Diccionario con los Parámetros de la simulación')
	
print('Lectura completada con éxito. Ejecute VolII para asegurarse de que no falta contenido.')

