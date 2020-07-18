import numpy as np 
import csv
from os import listdir
from shutil import rmtree
Ruta=[]
for out in listdir("./DatosBruto"):
	Ruta.append(out)
Ruta.sort()
lista=['Data_ibuf.npy', 'Data_mrv.npy' , 'Data_bar.npy']
a=len(lista)
con=1
w=input('Datos que desea mostrar de cada fichero:' '\n  por comodidad introduce enteros bajos (Ej:20)    ')
problema=False 
problema2=False
for j in range(len(Ruta)):
	try:
		w=int(w)
	except ValueError:
		print('No es un número entero')
		print('//////	BREAK 	//////')
		break

	print('||||||||||||||||||||||||||||||||||||||||||||||||||CAMBIO DE SALIDA||||||||||||||||||||||||||||||||||||||||||||||||||')
	print('Salida de la simulación: ' + Ruta[j])

	for i in range((a)):
		print('Se van a mostrar los valores del fichero de datos:   ' + lista[i])
		try:
			F=np.load('d5004/'+Ruta[j]+'/'+lista[i])
			print(F[0:w])
		except: 
			print('NO EXISTE EL FICHERO ' + lista[i] + ' EN LA SALIDA '+ Ruta[j])
			problema=True
	try:	
		Values = {}
		for key, val in csv.reader(open('d5004/'+Ruta[j]+ "/Values.csv")):
	    		Values[key] = val
		print('Valores de los parámetros de la  simulación de interés:	')
		print(Values)
	except:
		print('NO EXISTE EL FICHERO VALUES EN LA SALIDA '+ Ruta[j])
		problema=True
	try:
		SI_Val={}
		for key, val in csv.reader(open('d5004/'+Ruta[j]+'/SI_Val.csv')):
			SI_Val[key]= val
		print('Parámetros con el cambio de unidades:		')
		print(SI_Val)
	except FileNotFoundError:
		print('No se ha creado el fichero SI_Val en ' + Ruta[j])
		from Conversions import *
		Conversion()
		print('Reinicie el programa')
		problema2=True
	if problema2==True:
		print('No se ha creado el fichero SI_Val en ' + Ruta[j] + ' y se ha activado el programa Conversions para crear /SI_Val \n Reinicie el programa de lectura')
		break
	if problema==True:
		print('ERROR' '\n ELIMINADA LA CARPETA ' + Ruta[j] +'  EJECUTE EL PROGRAMA VOLCADO''\n BREAK')
		rmtree('d5004/'+Ruta[j])

		break

