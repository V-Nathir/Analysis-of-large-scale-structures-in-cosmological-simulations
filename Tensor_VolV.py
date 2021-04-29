
from mpl_toolkits.mplot3d import Axes3D
from os import listdir
import pylab as pl
from pylab import *
import sys
import csv
import os
import pickle
import itertools
import cv2
import statistics
import math
import matplotlib.animation as animation
import matplotlib.cm as cm
from astropy import *
import time as time_module
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import Ellipse

#This will not work unless you download the module at the beginning
#of this file and put it in the same directory or the executable directory
#for python on your machine

try:
	os.mkdir('Particle/Tensor')
except:
	print('#INFO: Directorio Particle/Tensor creado')
try:
	os.mkdir('Particle/Tensor/Loads')
except:
	print('#INFO: Directorio Particle/Tensor creado')

dir='Particle/CMs'
Icabecera=['Ixx','Iyy','Izz','Ixy','Ixz','Izy']
webs=listdir('Particle/CMs')

Ruta=[]
for out in listdir("./DatosBruto"):
	Ruta.append(out)
Ruta.sort()
webs.remove('Cabecera.txt')
webs=['sat_nopersistentes_41132.dat_web.txt', 'sat_persistentes_41132.dat_web.txt', 'web2_thres19_ids_500_500_500_5004_3667.dat_web.txt', 'web3_thres19_ids_500_500_500_5004_3667.dat_web.txt']
print(webs)
print('Habrá dos ecuaciones: \n I_ii= sum[ masa*(x_k^2+x_j^2)/r²] |\n| I_ij= sum[ masa*(x_i*x_j)/r²] \n')
select='h'
while select!='max' and select!= 'min' and select!='non':
	select=input('Representar el autovector asociado al mayor autovalor, menor o intermedio: max,min,non---> ')

cm1=[]
cmper=[]
cmwEb=[]
cmwebb=[]
contadorm=0
Satdic={}
for rut in range(len(Ruta)):
	F=np.load('d5004/'+ Ruta[rut]+'/Data_mrv.npy')
	print(f'Salida: {Ruta[rut]} || Creando diccionario con los datos en bruto-> {round((rut+1)*100/int(len(Ruta)))}%')
	Satdic[Ruta[rut]]=np.array(F)
MegaSat={}
substring='sat'
for tn in webs:
	megabuf=[]

	if substring in tn:
		F=open('Particle/'+str(tn)[0:-8]+'/web.txt','r')
		fsat=F.read()
		fsat=np.array(fsat.split())
		for t in range((len(fsat))):
			megabuf.append([])
			megabuf[t]=fsat[t]
		MegaSat[tn]=np.array(fsat)
	else:
		F=open('Particle/Webs/'+str(tn)[0:-8],'r')
		fsat=F.read()
		fsat=np.array((fsat.split()))
		MegaSat[tn]=np.delete(array(fsat),0)
medi='h'
while medi!='yes' and medi!='no':
	medi=input('¿Calcular el Momento Angular mediante la mediana? (yes/no): \n La alternativa:media --> ')

print(MegaSat['web2_thres19_ids_500_500_500_5004_3667.dat_web.txt'])
masatotal={}
elip1=[]
elip2=[]
elip3=[]
elip4=[]
velcmx=0
velcmy=0
velcmz=0
contadorangulo=0

for w in webs:
	W=np.loadtxt('Particle/CMs/'+w)
	#print('Valores del Cm, masa, partículas implicadas y totales\n')
	#print(W)
	posiciones=[]
	if contadorm==0:
		cm1=np.array(W[:,0:3])
	if contadorm==1:
		cmper=np.array(W[:,0:3])
	if contadorm==2:
		cmwEb=np.array(W[:,0:3])
	if contadorm==3:
		cmwebb=np.array(W[:,0:3])
	I=np.zeros((int(len(W)),6))
	elipsoide=[]
	newvelx=[]
	newvely=[]
	newvelz=[]
	estrellabuf=np.zeros(61, dtype=int)
	gasbuf=np.zeros(61, dtype=int)
	for ruta, v in zip(Ruta,range(len(Ruta))):
		print(f'---- \n Calculando Tensor de Inercia para {ruta}: ... {round(100*(v+1)/61)}%')

		Ixx=0
		Iyy=0
		Izz=0
		Ixy=0
		Ixz=0
		Izy=0
		mass=0
		velbufx=0
		velbufy=0
		velbufz=0
		gas=0
		estrellas=0
		otros=0
	
		for i in MegaSat[w]:
			#print(len(MegaSat[w]))
			if contadorm<=1:
				i=int(i)
			if contadorm>=2:
				i=int(i)-1
			if int(round(float(Satdic[ruta][i,7])))==int(1):
				gas=gas+1	
			if int(round(float(Satdic[ruta][i,7])))==int(-1):
				estrellas=estrellas+1
			if int(round(float(Satdic[ruta][i,7])))!= -1 and int(Satdic[ruta][i,7])!= 1:
				otros=otros+1

			velbufx=velbufx+Satdic[ruta][i,0]*Satdic[ruta][i,4]
			velbufy=velbufy+Satdic[ruta][i,0]*Satdic[ruta][i,5]
			velbufz=velbufz+Satdic[ruta][i,0]*Satdic[ruta][i,6]
			mass=mass+Satdic[ruta][i,0]
			x=Satdic[ruta][i,1]-W[v,0]
			y=Satdic[ruta][i,2]-W[v,1]
			z=Satdic[ruta][i,3]-W[v,2]
			Ixx=Satdic[ruta][i,0]*(y**2+z**2)/(x**2+y**2+z**2)+Ixx
			Iyy=Satdic[ruta][i,0]*(x**2+z**2)/(x**2+y**2+z**2)+Iyy
			Izz=Satdic[ruta][i,0]*(x**2+y**2)/(x**2+y**2+z**2)+Izz
			Ixy=-Satdic[ruta][i,0]*(x*y)/(x**2+y**2+z**2)+Ixy
			Ixz=-Satdic[ruta][i,0]*(x*z)/(x**2+y**2+z**2)+Ixz
			Izy=-Satdic[ruta][i,0]*(y*z)/(x**2+y**2+z**2)+Izy
		estrellabuf[v]=int(estrellas)
		gasbuf[v]=int(gas)
		velcm1=velbufx/mass
		velcm2=velbufy/mass
		velcm3=velbufz/mass
		VCM=velcm1,velcm2,velcm3
		print(f'Para la web {w} salida {ruta} hay: \n ESTRELLAS--> {estrellas} GAS--> {gas} OTROS--> {otros} ')
		newvelx.append(velcm1)
		newvely.append(velcm2)
		newvelz.append(velcm3)
		
		I[v,:]=Ixx,Iyy,Izz,Ixy,Ixz,Izy
	newvelx=np.array(newvelx)
	newvely=np.array(newvely)
	newvelz=np.array(newvelz)
	angular2=[]
	print(f'Cantidad de partículas {len(MegaSat[w])}')
	longit=np.round(np.linspace(0,61,61))
	R=plt.figure(figsize=(13.0, 10.0))
	ax=R.add_subplot(111)
	estrellabuf=np.array(estrellabuf)
	gasbuf=np.array(gasbuf)
	ax.bar(longit,estrellabuf)
	plt.xticks(longit,Ruta)
	plt.xticks(rotation=75)
	rects=ax.patches
	labels=[ i for i in estrellabuf]
	for rect,label in zip(rects,labels):
		height=rect.get_height()
		ax.text(rect.get_x() + rect.get_width() / 2, height + 10, label,ha='center', va='bottom',rotation=90)
	ax.set_ylabel('Poblaciones estelares')

	if contadorm==0:
		ax.set_title('Variación de Estrellas -- Sat. No Persistentes')
	if contadorm==1:
		ax.set_title('Variación de Estrellas -- Sat. Persistentes')
	if contadorm==2:
		ax.set_title('Variación de Estrellas -- Web2')
	if contadorm==3:
		ax.set_title('Variación de Estrellas -- Web3')
	plt.savefig('Particle/Tensor/Poblacion_estelar'+w+'.png')
	plt.close()	
	R=plt.figure(figsize=(13.0, 10.0))
	ax=R.add_subplot(111)
	estrellabuf=np.array(estrellabuf)
	gasbuf=np.array(gasbuf)
	ax.bar(longit-0.3,estrellabuf, label='Estrellas')
	ax.bar(longit+0.3,gasbuf,color='black',label='Gas')
	plt.legend(prop={'size': 12},shadow=True)
	comb=np.append(estrellabuf,gasbuf)
	plt.xticks(longit,Ruta)
	plt.xticks(rotation=75)
	rects=ax.patches
	labels=[ i for i in comb]
	for rect,label in zip(rects,labels):
		height=rect.get_height()
		ax.text(rect.get_x() + rect.get_width() / 2, height + 10, label,ha='center', va='bottom',rotation=90)
	ax.set_ylabel('Poblaciones estelares')

	if contadorm==0:
		ax.set_title('Variación de Estrellas & Gas -- Sat. No Persistentes')
	if contadorm==1:
		ax.set_title('Variación de Estrellas & Gas  -- Sat. Persistentes')
	if contadorm==2:
		ax.set_title('Variación de Estrellas & Gas  -- Web2')
	if contadorm==3:
		ax.set_title('Variación de Estrellas & Gas  -- Web3')
	plt.savefig('Particle/Tensor/GAS_Estrellas_'+w+'.png')
	plt.close()	
	R=plt.figure(figsize=(13.0, 10.0))
	ax=R.add_subplot(111)

	ax.bar(longit,gasbuf)

	plt.xticks(longit,Ruta)
	plt.xticks(rotation=75)
	rects=ax.patches
	labels=[ i for i in gasbuf]
	for rect,label in zip(rects,labels):
		height=rect.get_height()
		ax.text(rect.get_x() + rect.get_width() / 2, height + 10, label,ha='center', va='bottom',rotation=90)
	ax.set_ylabel('Poblaciones estelares')
	
	if contadorm==0:
		ax.set_title('Variación de Gas -- Sat. No Persistentes')
	if contadorm==1:
		ax.set_title('Variación de Gas -- Sat. Persistentes')
	if contadorm==2:
		ax.set_title('Variación de Gas -- Web2')
	if contadorm==3:
		ax.set_title('Variación de Gas -- Web3')
	plt.savefig('Particle/Tensor/Poblacion_gas'+w+'.png')
	plt.close()	
	for r in range(len(Ruta)):
		f=0
		angular=[]
		print(f'Calculando Momento angular para las salidas: ... {round(100*r/61)}%')
		for i in MegaSat[w]:
			if contadorm==0 or contadorm==1:
				i=int(i)
			if contadorm==2 or contadorm==3:
				i=int(i)-1

			#Velocidad angular
			X=np.array([[Satdic[Ruta[r]][i,0]*((Satdic[Ruta[r]][i,2]-W[r,1])*(Satdic[Ruta[r]][i,6]-newvelz[r])-(Satdic[Ruta[r]][i,3]-W[r,2])*(Satdic[Ruta[r]][i,5]-newvely[r]))],[Satdic[Ruta[r]][i,0]*(-(Satdic[Ruta[r]][i,1]-W[r,0])*(Satdic[Ruta[r]][i,6]-newvelz[r])+(Satdic[Ruta[r]][i,3]-W[r,2])*(Satdic[Ruta[r]][i,4]-newvelx[r]))],[Satdic[Ruta[r]][i,0]*((Satdic[Ruta[r]][i,1]-W[r,0])*(Satdic[Ruta[r]][i,5]-newvely[r])-(Satdic[Ruta[r]][i,2]-W[r,1])*(Satdic[Ruta[r]][i,4]-newvelx[r]))]])


			angular.append([])
			angular[f]= X[0]/(np.sqrt(X[0]**2+X[1]**2+X[2]**2)),X[1]/(np.sqrt(X[0]**2+X[1]**2+X[2]**2)),X[2]/(np.sqrt(X[0]**2+X[1]**2+X[2]**2))
			posiciones.append(int(r))
			f=f+1
		angular=np.array(angular)
		if medi=='no':
			momentoangularx=np.sum(angular[:,0])/len(angular[:,0])
			momentoangulary=np.sum(angular[:,1])/len(angular[:,1])
			momentoangularz=np.sum(angular[:,2])/len(angular[:,2])
		if medi=='yes':
			momentoangularx=float(statistics.median(angular[:,0]))
			momentoangulary=float(statistics.median(angular[:,1]))
			momentoangularz=float(statistics.median(angular[:,2]))
		
		
		angular2.append([])
		angular2[r]=momentoangularx,momentoangulary,momentoangularz


	angular2=np.array(angular2)

	print(f'Longitud del angular {len(angular2)}')
	masatotal[w]=mass
	print('\n   -----')

	print('\n   -----')
	
	buftensor={}
	bufautovectores={}
	bufautovalores={}
	elipsoide=[]
	trialxi=[]
	ellipticity=[]
	prolatness=[]
	print(f'Longitud de la velocidad angular pariculas {len(angular)}')
	for ruta, n in zip(Ruta,range(len(Ruta))):
		elipsoide.append([])
		trialxi.append([])
		ellipticity.append([])
		prolatness.append([])
		bufI=np.array([[I[n,0],I[n,3],I[n,4]],[I[n,3],I[n,1],I[n,5]],[I[n,4],I[n,5],I[n,2]]])
		#print(f'--------------\n {ruta}')

		#print(f'Tensor de inercia \n {bufI}\n ')

	
		autovalores, autovect = np.linalg.eig(bufI)
		#print('\n Autovectores')
		#print(autovect)
		#print('\n Autovalores del tensor ')
		#print(autovalores)
		lamda=np.sort(autovalores)
		a=np.sqrt(5*(lamda[1]-lamda[0]+lamda[2])/(2*masatotal[w]))

		b=np.sqrt(5*(lamda[2]-lamda[1]+lamda[0])/(2*masatotal[w]))
		c=np.sqrt(5*(lamda[0]-lamda[2]+lamda[1])/(2*masatotal[w]))
		T=(1-b**2/a**2)/(1-c**2/a**2)
		ell=(a**2-c**2)/(a**2+b**2+c**2)
		prol=(a**2+c**2-2*b**2)/(a**2+b**2+c**2)
		trialxi[n]=T
		ellipticity[n]=ell
		prolatness[n]=prol
		elipsoide[n]=a,b,c
		bufautovalores[ruta]=autovalores
		bufautovectores[ruta]=autovect
		buftensor[ruta]=bufI
	elipsoide=np.array(elipsoide)
	trialxi=np.array(trialxi)
	ellipticity=np.array(ellipticity)
	prolatness=np.array(prolatness)
	if contadorm==0:
		elip1=elipsoide
		triaxi1=trialxi
		ellipticity1=ellipticity
		prolatness1=prolatness
	if contadorm==1:
		elip2=elipsoide
		triaxi2=trialxi
		ellipticity2=ellipticity
		prolatness2=prolatness
	if contadorm==2:
		elip3=elipsoide
		triaxi3=trialxi
		ellipticity3=ellipticity
		prolatness3=prolatness
	if contadorm==3:
		elip4=elipsoide
		triaxi4=trialxi
		ellipticity4=ellipticity
		prolatness4=prolatness
	contadorm=contadorm+1
	print('Elipsoide a>b>c')
	print('a         b         c ')
	print(elipsoide)
	fig=plt.figure(figsize=(13.0, 10.0))
	ax=fig.add_subplot(111)
	colors = cm.rainbow(np.linspace(0,1,len(angular2[:,0])))
	angular2=np.array(angular2)
	for hui in range(len(angular2[:,0])):

		c=[colors[hui]]
		
		ax.scatter(angular2[hui,0],angular2[hui,1],marker='o',c=c)
	

	ax.set_xlabel('X')

	ax.set_ylabel('Y')
	plt.legend()
	plt.title('Momento Angular L')
	plt.grid(True)
	plt.savefig('Particle/Tensor/XYvelocidad_angular'+w+'.png')
	plt.show()
	plt.close()
	fig=plt.figure(figsize=(13.0, 10.0))
	ax=fig.add_subplot(111)
	colors = cm.rainbow(np.linspace(0,1,len(angular2[:,0])))
	angular2=np.array(angular2)
	if contadorangulo==0:
		MomAngular1=angular2
	if contadorangulo==1:
		MomAngular2=angular2
	if contadorangulo==2:
		MomAngular3=angular2
	if contadorangulo==3:
		MomAngular4=angular2

	for hui in range(len(angular2[:,0])):

		c=[colors[hui]]
		ax.scatter(angular2[hui,0],angular2[hui,2],marker='o',c=c)

	ax.set_xlabel('X')

	ax.set_ylabel('Z')
	plt.legend(prop={'size': 12},shadow=True,facecolor='grey')
	plt.title('Momento Angular L')
	plt.grid(True)
	plt.savefig('Particle/Tensor/XZvelocidad_angular'+w+'.png')
	plt.show()
	plt.close()
	fig=plt.figure(figsize=(13.0, 10.0))
	ax=plt.axes(projection='3d')
	colors = cm.rainbow(np.linspace(0,1,len(angular2[:,0])))
	angular2=np.array(angular2)
	print(f'Momento Angular.  {angular2}')
	#for hui in range(len(angular2[:,0])):

		#c=[colors[hui]]

	ax.scatter3D(angular2[:,0],angular2[:,1],angular2[:,2],marker='o',c=colors,cmap='Greens',s=35)
	ax.plot3D(angular2[:,0],angular2[:,1],angular2[:,2],'gray')

	ax.set_xlabel('X')

	ax.set_ylabel('Y')

	ax.set_zlabel('Z')

	plt.legend()
	if contadorangulo==0:
		plt.title('Momento Angular L -- Satélites No Persistentes')
	if contadorangulo==1:
		plt.title('Momento Angular L -- Satélites Persistentes')
	if contadorangulo==2:
		plt.title('Momento Angular L -- Web2')
	if contadorangulo==3:
		plt.title('Momento Angular L -- Web3')
	plt.grid(True)
	plt.savefig('Particle/Tensor/momento_angular'+w+'.png')
	plt.show()
	plt.close()
	contadorangulo=contadorangulo+1

	np.savetxt('Particle/Tensor/'+str(w)+'_Elipsoide.txt',elipsoide,fmt='%s')
	with open('Particle/Tensor/'+w+'_TensorInertia.txt','w') as file:
		file.write(str(buftensor))
	file.close()
	with open('Particle/Tensor/'+w+'_Autovectores.txt','w') as file:
		file.write(str(bufautovectores))
	file.close()
	with open('Particle/Tensor/'+w+'_Autovalores.txt','w') as file:
		file.write(str(bufautovalores))
	file.close()

	pickle_obj=open('Particle/Tensor/Loads/'+w+'_TensorInertia.pickle','wb')
	pickle.dump(buftensor,pickle_obj)
	pickle_obj.close()
	pickle_obj=open('Particle/Tensor/Loads/'+w+'_Autovectores.pickle','wb')
	pickle.dump(bufautovectores,pickle_obj)
	pickle_obj.close()
	pickle_obj=open('Particle/Tensor/Loads/'+w+'_Autovalores.pickle','wb')
	pickle.dump(bufautovalores,pickle_obj)
	pickle_obj.close()

fig=plt.figure(figsize=(13.0, 10.0))
ax=plt.axes(projection='3d')
colors = cm.rainbow(np.linspace(0,1,len(angular2[:,0])))
angular2=np.array(angular2)
print(f'Momento Angular.  {angular2}')
print(f' Desviaciones típicas : \n SNP {np.std(MomAngular1,0)} \n SP {np.std(MomAngular2,0)} \n Web2 {np.std(MomAngular3,0)} \n Web3 {np.std(MomAngular4,0)}')
#for hui in range(len(angular2[:,0])):

	#c=[colors[hui]]

ax.scatter3D(MomAngular1[:,0],MomAngular1[:,1],MomAngular1[:,2],marker='o',c=colors,s=15)
ax.plot3D(MomAngular1[:,0],MomAngular1[:,1],MomAngular1[:,2],'red',label='Satélites no persistentes')
ax.scatter3D(MomAngular2[:,0],MomAngular2[:,1],MomAngular2[:,2],marker='o',c=colors,s=15)
ax.plot3D(MomAngular2[:,0],MomAngular2[:,1],MomAngular2[:,2],'green',label='Satélites persistentes')
#ax.scatter3D(MomAngular3[:,0],MomAngular3[:,1],MomAngular3[:,2],marker='o',c=colors,s=15)
ax.plot3D(MomAngular3[:,0],MomAngular3[:,1],MomAngular3[:,2],'blue',label='Web2')
#ax.scatter3D(MomAngular4[:,0],MomAngular4[:,1],MomAngular4[:,2],marker='o',c=colors,s=15)
ax.plot3D(MomAngular4[:,0],MomAngular4[:,1],MomAngular4[:,2],'magenta', label='Web3')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.view_init(elev=35,azim=90)
plt.legend()
plt.title('Momento Angular L De las Cuatro Webs')
plt.grid(True)
plt.savefig('Particle/Tensor/Momento_angular.png')
plt.show()
plt.close()

AR1=[]
AR2=[]
AR3=[]
AR4=[]
DEC1=[]
DEC2=[]
DEC3=[]
DEC4=[]

for bu in range(len(MomAngular1[:,0])):
	x=MomAngular1[bu,0]
	y=MomAngular1[bu,1]
	z=MomAngular1[bu,2]
	if z >0:
		DEC=(np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if z<0:
		DEC=(-np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if y>0 and x>0:
		AR=np.arctan(y/x)

	if y<0 and x>0:
		AR=np.arctan(y/x)

	if x<0:
		AR=np.pi+np.arctan(y/x)
		if AR>np.pi:
			AR=AR-2*np.pi
	if DEC<0:
		if AR<0:
			AR=(AR)+np.pi
		if AR>0:
			AR=(AR)-np.pi
		DEC=np.absolute(DEC)
	AR1.append(AR)
	DEC1.append(DEC)

	x=MomAngular2[bu,0]
	y=MomAngular2[bu,1]
	z=MomAngular2[bu,2]
	if z >0:
		DEC=(np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if z<0:
		DEC=(-np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if y>0 and x>0:
		AR=np.arctan(y/x)

	if y<0 and x>0:
		AR=np.arctan(y/x)

	if x<0:
		AR=np.pi+np.arctan(y/x)
		if AR>np.pi:
			AR=AR-2*np.pi
	if DEC<0:
		if AR<0:
			AR=(AR)+np.pi
		if AR>0:
			AR=(AR)-np.pi
		DEC=np.absolute(DEC)
	AR2.append(AR)
	DEC2.append(DEC)

	x=MomAngular3[bu,0]
	y=MomAngular3[bu,1]
	z=MomAngular3[bu,2]
	if z >0:
		DEC=(np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if z<0:
		DEC=(-np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if y>0 and x>0:
		AR=np.arctan(y/x)

	if y<0 and x>0:
		AR=np.arctan(y/x)

	if x<0:
		AR=np.pi+np.arctan(y/x)
		if AR>np.pi:
			AR=AR-2*np.pi
	if DEC<0:
		if AR<0:
			AR=(AR)+np.pi
		if AR>0:
			AR=(AR)-np.pi
		DEC=np.absolute(DEC)
	AR3.append(AR)
	DEC3.append(DEC)

	x=MomAngular4[bu,0]
	y=MomAngular4[bu,1]
	z=MomAngular4[bu,2]
	if z >0:
		DEC=(np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if z<0:
		DEC=(-np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if y>0 and x>0:
		AR=np.arctan(y/x)

	if y<0 and x>0:
		AR=np.arctan(y/x)

	if x<0:
		AR=np.pi+np.arctan(y/x)
		if AR>np.pi:
			AR=AR-2*np.pi
	if DEC<0:
		if AR<0:
			AR=(AR)+np.pi
		if AR>0:
			AR=(AR)-np.pi
		DEC=np.absolute(DEC)
	AR4.append(AR)
	DEC4.append(DEC)


print(f' Desviaciones típicas : \n SNP {np.std(AR1)} {np.std(DEC1)} \n SP {np.std(AR2)} {np.std(DEC2)} \n Web2 {np.std(AR3)} {np.std(DEC3)} \n Web3 {np.std(AR4)} {np.std(DEC4)}')
ceros=np.zeros(61)
a=input('Aitoff')
colors = cm.rainbow(np.linspace(0, 1, 61))

plt.figure(figsize=(13.0, 10.0))

plt.subplot(221, projection="aitoff")
for hui in range(0,61):
	c=[colors[hui]]
	plt.scatter(AR1[hui],(DEC1[hui]),color=c,marker='o',s=12)
plt.title('Satétiles no persistentes',fontsize=20)
plt.grid(True)
plt.legend()

plt.subplot(222, projection="aitoff")
for hui in range(0,61):
	c=[colors[hui]]
	plt.scatter(AR2[hui],DEC2[hui],color=c,marker='o',s=12)
plt.title('Satélites persistentes',fontsize=20)
plt.grid(True)
plt.legend()


plt.subplot(223, projection="aitoff")
for hui in range(0,61):
	c=[colors[hui]]
	plt.scatter(AR3[hui],DEC3[hui],color=c,marker='o',s=12)

plt.title('Web2',fontsize=20)
plt.grid(True)
plt.legend()

plt.subplot(224, projection="aitoff")
for hui in range(0,61):
	c=[colors[hui]]
	plt.scatter(AR4[hui],DEC4[hui],color=c,marker='o',s=12)

plt.title('Web3',fontsize=20)


plt.grid(True)



plt.savefig('Particle/Tensor/Angular_Aitoff.png')
plt.show()

plt.close()





longit=np.round(np.linspace(0,61,61))
R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)
ax.plot(longit,elip1[:,0],linestyle='--',marker='o',color='red',markersize=4,label='a -- Sat.No persistentes')
ax.plot(longit,elip1[:,1],linestyle='--',marker='o',color='darkred',markersize=4, label='b -- Sat.No persistentes')
ax.plot(longit,elip1[:,2],linestyle='--',marker='o',color='indianred',markersize=4, label='c -- Sat.No persistentes')	
ax.plot(longit,elip2[:,0],linestyle='--',marker='o',color='green',markersize=4,label='a -- Sat. persistentes')
ax.plot(longit,elip2[:,1],linestyle='--',marker='o',color='darkgreen',markersize=4, label='b -- Sat. persistentes')
ax.plot(longit,elip2[:,2],linestyle='--',marker='o',color='lime',markersize=4, label='c -- Sat. persistentes')	

plt.xticks(longit,Ruta)
plt.xticks(rotation=90)
rects=ax.patches
labels=[ i for i in elip1]
for rect,label in zip(rects,labels):
	height=rect.get_height()
	ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
ax.set_ylabel('Valores del semieje')
ax.set_title('Variaciones de los semiejes')
plt.legend(shadow=True)
plt.grid(True)

plt.savefig('Particle/Tensor/SemiejesSAt.png')
plt.show()
plt.close()
R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)
ax.plot(longit,elip3[:,0],linestyle='--',marker='o',color='blue',markersize=4,label='a -- Web2')
ax.plot(longit,elip3[:,1],linestyle='--',marker='o',color='mediumblue',markersize=4, label='b -- Web2')
ax.plot(longit,elip3[:,2],linestyle='--',marker='o',color='navy',markersize=4, label='c -- Web2')	
ax.plot(longit,elip4[:,0],linestyle='--',marker='o',color='magenta',markersize=4,label='a -- Web3')
ax.plot(longit,elip4[:,1],linestyle='--',marker='o',color='darkmagenta',markersize=4, label='b -- Web3')
ax.plot(longit,elip4[:,2],linestyle='--',marker='o',color='mediumvioletred',markersize=4, label='c -- Web3')	

plt.xticks(longit,Ruta)
plt.xticks(rotation=90)
rects=ax.patches
labels=[ i for i in elip1]
for rect,label in zip(rects,labels):
	height=rect.get_height()
	ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
ax.set_ylabel('Valores del semieje')
ax.set_title('Variaciones de los semiejes')
plt.legend(shadow=True)
plt.grid(True)

plt.savefig('Particle/Tensor/SemiejesWEB.png')
plt.show()
plt.close()
R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)

ax.scatter(elip1[:,1]/elip1[:,0],elip1[:,2]/elip1[:,0],marker='o',c='red',s=20,label='Satélites no persistentes')
ax.scatter(elip2[:,1]/elip2[:,0],elip2[:,2]/elip2[:,0],marker='o',c='green',s=20, label='Satélites persistentes')
ax.scatter(elip3[:,1]/elip3[:,0],elip3[:,2]/elip3[:,0],marker='o',c='blue',s=20, label='Web2')	
ax.scatter(elip4[:,1]/elip4[:,0],elip4[:,2]/elip4[:,0],marker='o',c='magenta',s=20,label='Web3')

ax.plot(np.linspace(0,1,61),np.linspace(0,1,61),color='magenta', linestyle='solid', label='Prolate Spheroids')


ax.set_ylabel('Valor  de c/a')
ax.set_xlabel('Valor  de b/a')

ax.set_title('Relación entre b/a y c/a')
plt.legend(shadow=True)
plt.grid(True)

plt.savefig('Particle/Tensor/Semieje_bca.png')
plt.show()
plt.close()
longit=np.round(np.linspace(0,61,61))
R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)
ax.plot(longit,elip1[:,1]/elip1[:,0],marker='o',color='red',markersize=4,label='b/a --Sat.No persistentes')
ax.plot(longit,elip1[:,2]/elip1[:,0],marker='o',color='indianred',markersize=4, label='c/a --Sat.No persistentes')
ax.plot(longit,elip2[:,1]/elip2[:,0],marker='o',color='green',markersize=4, label='b/a --Sat. persistentes')	
ax.plot(longit,elip2[:,2]/elip2[:,0],marker='o',color='darkgreen',markersize=4,label='c/a --Sat. persistentes')

plt.xticks(longit,Ruta)
plt.xticks(rotation=75)
rects=ax.patches
labels=[ i for i in elip1]
for rect,label in zip(rects,labels):
	height=rect.get_height()
	ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
ax.set_ylabel('Valor  de las relaciones')
ax.set_title('Comportamiento de los semiejes')
plt.legend(shadow=True)
plt.grid(True)

plt.savefig('Particle/Tensor/Semiejes_comparativaSAt.png')
plt.show()
plt.close()
R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)
ax.plot(longit,elip3[:,1]/elip3[:,0],linestyle='--',marker='o',color='blue',markersize=4,label='b/a --Web2')
ax.plot(longit,elip3[:,2]/elip3[:,0],linestyle='--',marker='o',color='navy',markersize=4, label='c/a --Web2')
ax.plot(longit,elip4[:,1]/elip4[:,0],linestyle='--',marker='o',color='magenta',markersize=4, label='b/a --Web3')	
ax.plot(longit,elip4[:,2]/elip4[:,0],linestyle='--',marker='o',color='darkmagenta',markersize=4,label='c/a --Web3')

plt.xticks(longit,Ruta)
plt.xticks(rotation=75)
rects=ax.patches
labels=[ i for i in elip1]
for rect,label in zip(rects,labels):
	height=rect.get_height()
	ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
ax.set_ylabel('Valor  de las relaciones')
ax.set_title('Comportamiento de los semiejes')
plt.legend(shadow=True)
plt.grid(True)

plt.savefig('Particle/Tensor/Semiejes_comparativaWebs.png')
plt.show()
plt.close()
R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)
ax.plot(longit,elip1[:,1]/elip1[:,0],linestyle='--',marker='o',color='red',markersize=4,label='Satélites no persistentes')
ax.plot(longit,elip2[:,1]/elip2[:,0],linestyle='--',marker='o',color='green',markersize=4, label='Satélites persistentes')
ax.plot(longit,elip3[:,1]/elip3[:,0],linestyle='--',marker='o',color='blue',markersize=4, label='Web2')	
ax.plot(longit,elip4[:,1]/elip4[:,0],linestyle='--',marker='o',color='magenta',markersize=4,label='Web3')

plt.xticks(longit,Ruta)
plt.xticks(rotation=75)
rects=ax.patches
labels=[ i for i in elip1]
for rect,label in zip(rects,labels):
	height=rect.get_height()
	ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
ax.set_ylabel('Valor  de b/a')
ax.set_title('Relación del los semiejes a y b')
plt.legend(shadow=True)
plt.grid(True)

plt.savefig('Particle/Tensor/Semieje_ba.png')
plt.show()
plt.close()
R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)
ax.plot(longit,elip1[:,2]/elip1[:,0],linestyle='--',marker='o',color='red',markersize=4,label='Satélites no persistentes')
ax.plot(longit,elip2[:,2]/elip2[:,0],linestyle='--',marker='o',color='green',markersize=4, label='Satélites persistentes')
ax.plot(longit,elip3[:,2]/elip3[:,0],linestyle='--',marker='o',color='blue',markersize=4, label='Web2')	
ax.plot(longit,elip4[:,2]/elip4[:,0],linestyle='--',marker='o',color='magenta',markersize=4,label='Web3')

plt.xticks(longit,Ruta)
plt.xticks(rotation=75)
rects=ax.patches
labels=[ i for i in elip1]
for rect,label in zip(rects,labels):
	height=rect.get_height()
	ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
ax.set_ylabel('Valor  c/a')
ax.set_title('Relación del los semiejes a y c')
plt.legend(shadow=True)
plt.grid(True)

plt.savefig('Particle/Tensor/Semieje_ca.png')
plt.show()
plt.close()
R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)
ax.plot(longit,triaxi1,linestyle='--',marker='o',color='red',markersize=4,label='Satélites no persistentes')
ax.plot(longit,triaxi2,linestyle='--',marker='o',color='green',markersize=4, label='Satélites persistentes')
ax.plot(longit,triaxi3,linestyle='--',marker='o',color='blue',markersize=4, label='Web2')	
ax.plot(longit,triaxi4,linestyle='--',marker='o',color='magenta',markersize=4,label='Web3')

plt.xticks(longit,Ruta)
plt.xticks(rotation=75)
rects=ax.patches
labels=[ i for i in triaxi1]
for rect,label in zip(rects,labels):
	height=rect.get_height()
	ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
ax.set_ylabel('Triaxility parameter')
ax.set_title('Triaxility parameter T')
plt.legend(shadow=True)
plt.grid(True)

plt.savefig('Particle/Tensor/T.png')
plt.show()
plt.close()
R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)
ax.plot(longit,ellipticity1,linestyle='--',marker='o',color='red',markersize=4,label='Satélites no persistentes')
ax.plot(longit,ellipticity2,linestyle='--',marker='o',color='green',markersize=4, label='Satélites persistentes')
ax.plot(longit,ellipticity3,linestyle='--',marker='o',color='blue',markersize=4, label='Web2')	
ax.plot(longit,ellipticity4,linestyle='--',marker='o',color='magenta',markersize=4,label='Web3')

plt.xticks(longit,Ruta)
plt.xticks(rotation=75)
rects=ax.patches
labels=[ i for i in ellipticity1]
for rect,label in zip(rects,labels):
	height=rect.get_height()
	ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
ax.set_ylabel('Ellipticity parameter')
ax.set_title('Ellipticity parameter e')
plt.legend(shadow=True)
plt.grid(True)

plt.savefig('Particle/Tensor/Ellipticity.png')
plt.show()
plt.close()
R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)
ax.plot(longit,prolatness1,linestyle='--',marker='o',color='red',markersize=4,label='Satélites no persistentes')
ax.plot(longit,prolatness2,linestyle='--',marker='o',color='green',markersize=4, label='Satélites persistentes')
ax.plot(longit,prolatness3,linestyle='--',marker='o',color='blue',markersize=4, label='Web2')	
ax.plot(longit,prolatness4,linestyle='--',marker='o',color='magenta',markersize=4,label='Web3')

plt.xticks(longit,Ruta)
plt.xticks(rotation=75)
rects=ax.patches
labels=[ i for i in prolatness1]
for rect,label in zip(rects,labels):
	height=rect.get_height()
	ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
ax.set_ylabel('Prolatness parameter')
ax.set_title('Prolatness parameter p')
plt.legend(shadow=True)
plt.grid(True)

plt.savefig('Particle/Tensor/Prolatness.png')
plt.show()
plt.close()
R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)
ax.plot(longit,elip1[:,0],linestyle='--',marker='o',color='red',markersize=4,label='Satélites no persistentes')
ax.plot(longit,elip2[:,0],linestyle='--',marker='o',color='green',markersize=4, label='Satélites persistentes')
ax.plot(longit,elip3[:,0],linestyle='--',marker='o',color='blue',markersize=4, label='Web2')	
ax.plot(longit,elip4[:,0],linestyle='--',marker='o',color='magenta',markersize=4,label='Web3')

plt.xticks(longit,Ruta)
plt.xticks(rotation=75)
rects=ax.patches
labels=[ i for i in elip1]
for rect,label in zip(rects,labels):
	height=rect.get_height()
	ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
ax.set_ylabel('Valor del semieje a')
ax.set_title('Variación del semieje a')
plt.legend(shadow=True)
plt.grid(True)

plt.savefig('Particle/Tensor/Semieje_a.png')
plt.show()
plt.close()
R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)
ax.plot(longit,elip1[:,1],linestyle='--',marker='o',color='red',markersize=4,label='Satélites no persistentes')
ax.plot(longit,elip2[:,1],linestyle='--',marker='o',color='green',markersize=4, label='Satélites persistentes')
ax.plot(longit,elip3[:,1],linestyle='--',marker='o',color='blue',markersize=4, label='Web2')	
ax.plot(longit,elip4[:,1],linestyle='--',marker='o',color='magenta',markersize=4,label='Web3')

plt.xticks(longit,Ruta)
plt.xticks(rotation=75)
rects=ax.patches
labels=[ i for i in elip1]
for rect,label in zip(rects,labels):
	height=rect.get_height()
	ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
ax.set_ylabel('Valor del semieje b')
ax.set_title('Variación del semieje b')
plt.legend(shadow=True)
plt.grid(True)
plt.savefig('Particle/Tensor/Semieje_b.png')
plt.show()
plt.close()
R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)
ax.plot(longit,elip1[:,2],linestyle='--',marker='o',color='red',markersize=4,label='Satélites no persistentes')
ax.plot(longit,elip2[:,2],linestyle='--',marker='o',color='green',markersize=4, label='Satélites persistentes')
ax.plot(longit,elip3[:,2],linestyle='--',marker='o',color='blue',markersize=4, label='Web2')	
ax.plot(longit,elip4[:,2],linestyle='--',marker='o',color='magenta',markersize=4,label='Web3')

plt.xticks(longit,Ruta)
plt.xticks(rotation=75)
rects=ax.patches
labels=[ i for i in elip1]
for rect,label in zip(rects,labels):
	height=rect.get_height()
	ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
ax.set_ylabel('Valor del semieje c')
ax.set_title('Variación del semieje c')
plt.legend(shadow=True)
plt.grid(True)
plt.savefig('Particle/Tensor/Semieje_c.png')
plt.show()
plt.close()
fig=plt.figure(figsize=(30.0, 25.0))
ax=fig.add_subplot(221, projection='3d')

plt.plot(cm1[:,0],cm1[:,1],cm1[:,2],linestyle='--',marker='o',color='red', label='Satélites no persistentes')
ax.scatter3D(cm1[60,0],cm1[60,1],cm1[60,2],marker='*',c='black',s=25)

plt.title('Satéites no persistentes',fontsize=22)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.legend(shadow=True)
ax.view_init(elev=45,azim=250)

plt.grid(True)

ax=fig.add_subplot(222, projection='3d')
plt.plot(cmper[:,0],cmper[:,1],cmper[:,2],linestyle='--',marker='o',color='green',label='Satétiles persistentes')

ax.scatter3D(cmper[60,0],cmper[60,1],cmper[60,2],marker='*',c='black',s=25)

plt.title('Satéites persistentes',fontsize=22)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.grid(True)
plt.legend()
ax.view_init(elev=22,azim=250)

ax=fig.add_subplot(223, projection='3d')

plt.plot(cmwEb[:,0],cmwEb[:,1],cmwEb[:,2],linestyle='--',marker='o',color='blue',label='Web2')
plt.title('Web2',fontsize=22)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.grid(True)
plt.legend()
ax.view_init(elev=20,azim=250)

ax=fig.add_subplot(224, projection='3d')

plt.plot(cmwebb[:,0],cmwebb[:,1],cmwebb[:,2],linestyle='--',marker='o',color='magenta',label='Web3')
plt.title('Web3',fontsize=22)
ax.set_xlabel('X')
plt.legend()
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.grid(True)
ax.view_init(elev=20,azim=250)

plt.savefig('Particle/Tensor/Centros_de_masas.png')
plt.show()
plt.close()
fig=plt.figure(figsize=(30.0, 25.0))
ax=fig.add_subplot(111, projection='3d')
plt.plot(cm1[:,0],cm1[:,1],cm1[:,2],linestyle='--',marker='o',markersize=4,color='red', label='Sat. no persistentes')

plt.plot(cmper[:,0],cmper[:,1],cmper[:,2],linestyle='--',markersize=4,marker='o',color='green',label='Sat. persistentes')

plt.plot(cmwEb[:,0],cmwEb[:,1],cmwEb[:,2],linestyle='--',markersize=4,marker='o',color='blue',label='Web2')

plt.plot(cmwebb[:,0],cmwebb[:,1],cmwebb[:,2],linestyle='--',markersize=4,marker='o',color='magenta',label='Web3')
ax.scatter3D(cmper[60,0],cmper[60,1],cmper[60,2],marker='*',c='black')
ax.scatter3D(cm1[60,0],cm1[60,1],cm1[60,2],marker='*',c='black')
#ax.scatter3D(cmper[37,0],cmper[37,1],cmper[37,2],marker='*',c='grey',s=15, label='Salida 41132')
#ax.scatter3D(cm1[37,0],cm1[37,1],cmper[37,2],marker='*',c='grey',s=15, label='Salida 41132')

ax.legend(prop={'size': 20},shadow=True)
ax.set_xlabel('X')

ax.set_ylabel('Y')
plt.title('Centros de Masas de cada Web',fontsize=25)
ax.set_zlabel('Z')
ax.view_init(elev=20,azim=200)

plt.savefig('Particle/Tensor/All_Centros_de_masas.png')
plt.show()
plt.close()




s=open('Particle/Tensor/Loads/'+webs[0]+'_Autovalores.pickle','rb')
	
A=pickle.load(s)
s.close()
s=open('Particle/Tensor/Loads/'+webs[1]+'_Autovalores.pickle','rb')

B=pickle.load(s)
s.close()
s=open('Particle/Tensor/Loads/'+webs[2]+'_Autovalores.pickle','rb')

C=pickle.load(s)
s.close()
s=open('Particle/Tensor/Loads/'+webs[3]+'_Autovalores.pickle','rb')

D=pickle.load(s)
s.close()

print('#\n VALORES DE LOS AUTOVALORES/ MOMENTOS PRINCIPALES DE INERCIA PARA CADA SAILIDA\n#')
acumulador1=[]
acumulador2=[]
acumulador3=[]
acumulador4=[]

for i in Ruta:
	print(f'Salida de la simulación {i}\n Ixx      Iyy      Izz')
	#print(webs[0])
	print(np.around(A[i],decimals=6))
	if select=='max':
		amx=np.amax(A[i])
	if select=='min':
		amx=np.amin(A[i])
	if select=='non':
		if A[i][0]!=np.amax(A[i]) and A[i][0]!=np.amin(A[i]) :
			amx=A[i][0]
		if A[i][1]!=np.amax(A[i]) and A[i][1]!=np.amin(A[i]) :
			amx=A[i][1]
		if A[i][2]!=np.amax(A[i]) and A[i][2]!=np.amin(A[i]) :
			amx=A[i][2]		
	if amx==A[i][0]:
		acumulador1.append([0])
	if amx==A[i][1]:
		acumulador1.append([1])
	if amx==A[i][2]:
		acumulador1.append([2])
	print(webs[1])
	print(np.around(B[i],decimals=6))
	if select=='max':
		amx=np.amax(B[i])
	if select=='min':
		amx=np.amin(B[i])
	if select=='non':
		if B[i][0]!=np.amax(B[i]) and B[i][0]!=np.amin(B[i]) :
			amx=B[i][0]
		if B[i][1]!=np.amax(B[i]) and B[i][1]!=np.amin(B[i]) :
			amx=B[i][1]
		if B[i][2]!=np.amax(B[i]) and B[i][2]!=np.amin(B[i]) :
			amx=B[i][2]	
	if amx==B[i][0]:
		acumulador2.append([0])
	if amx==B[i][1]:
		acumulador2.append([1])
	if amx==B[i][2]:
		acumulador2.append([2])
	print(webs[2])
	print(np.around(C[i],decimals=6))
	if select=='max':

		amx=np.amax(C[i])
	if select=='min':
		amx=np.amin(C[i])
	if select=='non':
		if C[i][0]!=np.amax(C[i]) and C[i][0]!=np.amin(C[i]) :
			amx=C[i][0]
		if C[i][1]!=np.amax(C[i]) and C[i][1]!=np.amin(C[i]) :
			amx=C[i][1]
		if C[i][2]!=np.amax(C[i]) and C[i][2]!=np.amin(C[i]) :
			amx=C[i][2]

	if amx==C[i][0]:
		acumulador3.append([0])
	if amx==C[i][1]:
		acumulador3.append([1])
	if amx==C[i][2]:
		acumulador3.append([2])
	print(webs[3])
	print(np.around(D[i],decimals=6))
	if select=='max':
		amx=np.amax(D[i])
	if select=='min':
		amx=np.amin(D[i])
	if select=='non':
		if D[i][0]!=np.amax(D[i]) and D[i][0]!=np.amin(D[i]) :
			amx=D[i][0]
		if D[i][1]!=np.amax(D[i]) and D[i][1]!=np.amin(D[i]) :
			amx=D[i][1]
		if D[i][2]!=np.amax(D[i]) and D[i][2]!=np.amin(D[i]) :
			amx=D[i][2]
	if amx==D[i][0]:
		acumulador4.append([0])
	if amx==D[i][1]:
		acumulador4.append([1])
	if amx==D[i][2]:
		acumulador4.append([2])
	print('----------------------\n')
s=open('Particle/Tensor/Loads/'+webs[0]+'_Autovectores.pickle','rb')
	
A=pickle.load(s)
s.close()
s=open('Particle/Tensor/Loads/'+webs[1]+'_Autovectores.pickle','rb')

B=pickle.load(s)
s.close()
s=open('Particle/Tensor/Loads/'+webs[2]+'_Autovectores.pickle','rb')

C=pickle.load(s)
s.close()
s=open('Particle/Tensor/Loads/'+webs[3]+'_Autovectores.pickle','rb')

D=pickle.load(s)
s.close()


w1AR=[]
w1DEC=[]
w2AR=[]
w2DEC=[]
w3AR=[]
w3DEC=[]
w4AR=[]
w4DEC=[]
ang1=np.zeros((61,3))
ang2=np.zeros((61,3))
ang3=np.zeros((61,3))
ang4=np.zeros((61,3))

print('#\n AUTOVECTORES PARA CADA SAILIDA\n#')
for i,a in zip(Ruta,range(len(Ruta))):
	print(f'Salida de la simulación {i}\n ł1      ł2      ł3')
	print(webs[0])
	print(np.around(A[i],decimals=4))
	e=acumulador1[a]
	x=A[i][0,e]
	y=A[i][1,e]
	z=A[i][2,e]
	ang1[a]=x,y,z
	if z >0:
		DEC=(np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if z<0:
		DEC=(-np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if y>0 and x>0:
		AR=np.arctan(y/x)

	if y<0 and x>0:
		AR=np.arctan(y/x)

	if x<0:
		AR=np.pi+np.arctan(y/x)
		if AR>np.pi:
			AR=AR-2*np.pi
	if DEC<0:
		if AR<0:
			AR=(AR)+np.pi
		if AR>0:
			AR=(AR)-np.pi
		DEC=np.absolute(DEC)

	w1DEC.append(DEC)
	w1AR.append(AR)

	print(webs[1])
	print(np.around(B[i],decimals=4))
	r=acumulador2[a]

	x=B[i][0,r]
	y=B[i][1,r]
	z=B[i][2,r]
	ang2[a]=x,y,z

	if z >0:
		DEC=(np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if z<0:
		DEC=(-np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if y>0 and x>0:
		AR=math.atan2(y,x)
	if y<0 and x>0:
		AR=math.atan2(y,x)
		#AR=2*np.pi-AR
	if x<0:
		AR=math.atan2(y,x)
	if  DEC<0:
		if AR<0:
			AR=(AR)+np.pi
		if AR>0:
			AR=(AR)-np.pi
		DEC=np.absolute(DEC)
	w2DEC.append(DEC)
	w2AR.append(AR)
	print(webs[2])
	print(np.around(C[i],decimals=4))
	t=acumulador3[a]

	x=C[i][0,t]
	y=C[i][1,t]
	z=C[i][2,t]
	ang3[a]=x,y,z

	if z >0:
		DEC=(np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if z<0:
		DEC=(-np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if y>0 and x>0:
		AR=np.arctan(y/x)

	if y<0 and x>0:
		AR=np.arctan(y/x)

	if x<0:
		AR=np.pi+np.arctan(y/x)
		if AR>np.pi:
			AR=AR-2*np.pi
	if DEC<0:
		if AR<0:
			AR=(AR)+np.pi
		if AR>0:
			AR=(AR)-np.pi
		DEC=np.absolute(DEC)
	w3DEC.append(DEC)
	w3AR.append(AR)
	print(webs[3])
	print(np.around(D[i],decimals=4))
	q=acumulador4[a]
	x=D[i][0,q]
	y=D[i][1,q]
	z=D[i][2,q]
	ang4[a]=x,y,z

	if z >0:
		DEC=(np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if z<0:
		DEC=(-np.pi/2-np.arctan(np.sqrt(x**2+y**2)/z) )
	if y>0 and x>0:
		AR=np.arctan(y/x)

	if y<0 and x>0:
		AR=np.arctan(y/x)

	if x<0:
		AR=np.pi+np.arctan(y/x)
		if AR>np.pi:
			AR=AR-2*np.pi
	if DEC<0:
		if AR<0:
			AR=(AR)+np.pi

		DEC=np.absolute(DEC)
		if AR>0:
			AR=(AR)-np.pi
	w4DEC.append(DEC)
	w4AR.append(AR)
	print('----------------------\n')

print(w2AR)
print(np.array(w2AR)*360/2/np.pi)
print(np.array(w2DEC)*360/2/np.pi)
print(np.array(w1AR)*360/2/np.pi)
print(np.array(w1DEC)*360/2/np.pi)
deltaang=[]
deltaang2=[]
for angl in range(len(ang1)):
	deltaang2.append([])
	deltaang.append([])
	deltaang[angl]=math.acos(np.dot(ang1[angl],ang3[angl])/((np.sqrt(ang1[angl,0]**2+ang1[angl,1]**2+ang1[angl,2]**2))*(np.sqrt(ang3[angl,0]**2+ang3[angl,1]**2+ang3[angl,2]**2))))*180/np.pi
	if deltaang[angl]>90:
		deltaang[angl]=180-deltaang[angl]
	deltaang2[angl]=math.acos(np.dot(ang2[angl],ang3[angl])/((np.sqrt(ang2[angl,0]**2+ang2[angl,1]**2+ang2[angl,2]**2))*(np.sqrt(ang3[angl,0]**2+ang3[angl,1]**2+ang3[angl,2]**2))))*180/np.pi
	if deltaang2[angl]>90:
		deltaang2[angl]=180-deltaang2[angl]


R=plt.figure(figsize=(13.0, 10.0))
ax=R.add_subplot(111)
ax.plot(longit,deltaang,linestyle='--',marker='o',color='red',markersize=4,label='Satélites no persistentes')
ax.plot(longit,deltaang2,linestyle='--',marker='o',color='green',markersize=4, label='Satélites persistentes')
plt.yticks(np.linspace(0,90,10))
plt.xticks(longit,Ruta)
plt.xticks(rotation=75)
rects=ax.patches
labels=[ i for i in elip1]
for rect,label in zip(rects,labels):
	height=rect.get_height()
	ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
ax.set_ylabel('Ángulo formado con Web2')
ax.set_title('Variación entre el ángulo Sat. Persistentes/ Sat. No Persistentes & Web2')
plt.legend(loc='lower right')
plt.grid(True)
plt.savefig('Particle/Tensor/Angulos.png')
plt.show()
plt.close()

ceros=np.zeros(61)
a=input('Aitoff')
colors = cm.rainbow(np.linspace(0, 1, 61))

plt.figure(figsize=(13.0, 10.0))

plt.subplot(221, projection="aitoff")
for hui in range(0,61):
	c=[colors[hui]]
	plt.scatter(w1AR[hui],(w1DEC[hui]),color=c,marker='o',s=12)
plt.title('Satétiles no persistentes',fontsize=20)
plt.grid(True)
plt.legend()

plt.subplot(222, projection="aitoff")
for hui in range(0,61):
	c=[colors[hui]]
	plt.scatter(w2AR[hui],w2DEC[hui],color=c,marker='o',s=12)
plt.title('Satélites persistentes',fontsize=20)
plt.grid(True)
plt.legend()


plt.subplot(223, projection="aitoff")
for hui in range(0,61):
	c=[colors[hui]]
	plt.scatter(w3AR[hui],w3DEC[hui],color=c,marker='o',s=12)

plt.title('Web2',fontsize=20)
plt.grid(True)
plt.legend()

plt.subplot(224, projection="aitoff")
for hui in range(0,61):
	c=[colors[hui]]
	plt.scatter(w4AR[hui],w4DEC[hui],color=c,marker='o',s=12)

plt.title('Web3',fontsize=20)

plt.legend()
plt.grid(True)



plt.savefig('Particle/Tensor/Aitoff.png')
plt.show()

plt.close()

plt.figure(figsize=(13.0, 10.0))

plt.subplot(111, projection="aitoff")
plt.scatter(w1AR[15:],(w1DEC[15:]),color='red',marker='o',s=10,label='Sat. No persistentes')
plt.scatter(w2AR[15:],(w2DEC[15:]),color='green',marker='o',s=10,label='Sat. Persistentes')
plt.scatter(w3AR[15:],(w3DEC[15:]),color='blue',marker='o',s=10, label='Web2')
plt.scatter(w4AR[15:],(w4DEC[15:]),color='magenta',marker='o',s=10, label='Web3')
plt.scatter(w1AR[60],w1DEC[60],color='black', marker='*', s=10, label='z=0, Sat. no persistentes')
plt.scatter(w2AR[60],w2DEC[60],color='black', marker='s', s=10, label='z=0, Sat. Persistentes')
plt.scatter(w3AR[60],w3DEC[60],color='black', marker='P', s=10, label='z=0, Web2')

plt.grid(True)
plt.legend(prop={'size': 12},shadow=True,facecolor='grey')


plt.savefig('Particle/Tensor/ALL_plot.png')
plt.show()

plt.close()


plt.close()
plt.figure(figsize=(13.0, 10.0))
plt.subplot(221)

plt.plot(np.array(w1AR)[15:]*180/np.pi,np.array(w1DEC)[15:]*180/np.pi,linestyle='--',marker='o',color='red', label='SatNoPer')
plt.title('Satéites no persistentes',fontsize=20)
plt.xlabel('AR',fontsize=15)
plt.ylabel('DEC',fontsize=20)
plt.plot(np.array(w1AR)[37]*180/np.pi,np.array(w1DEC)[37]*180/np.pi,c='black',marker='*',markersize=6,label='Salida 41132')
plt.plot(np.array(w1AR)[60]*180/np.pi,np.array(w1DEC)[60]*180/np.pi,c='blue',marker='*',markersize=10,label='Salida final')
plt.plot(np.array(w1AR)[15]*180/np.pi,np.array(w1DEC)[15]*180/np.pi,c='magenta',marker='*',markersize=10,label='Salida inicial')
plt.legend()


plt.grid(True)

plt.subplot(222)
plt.plot(np.array(w2AR)[15:]*180/np.pi,np.array(w2DEC)[15:]*180/np.pi,linestyle='--',marker='o',color='green',label='SatPer')
plt.title('Satéites persistentes',fontsize=20)
plt.xlabel('AR',fontsize=15)
plt.ylabel('DEC',fontsize=20)
plt.plot(np.array(w2AR)[37]*180/np.pi,np.array(w2DEC)[37]*180/np.pi,marker='*',c='black',markersize=6,label='Salida 41132')
plt.plot(np.array(w2AR)[60]*180/np.pi,np.array(w2DEC)[60]*180/np.pi,c='blue',marker='*',markersize=10,label='Salida final')
plt.plot(np.array(w2AR)[15]*180/np.pi,np.array(w2DEC)[15]*180/np.pi,c='magenta',marker='*',markersize=10,label='Salida inicial')
plt.legend()
plt.grid(True)

plt.subplot(223)
plt.plot(np.array(w3AR)[15:]*180/np.pi,np.array(w3DEC)[15:]*180/np.pi,linestyle='--',marker='o',color='blue',label='Web2')
plt.plot(np.array(w3AR)[60]*180/np.pi,np.array(w3DEC)[60]*180/np.pi,c='blue',marker='*',markersize=10,label='Salida final')
plt.plot(np.array(w3AR)[15]*180/np.pi,np.array(w3DEC)[15]*180/np.pi,c='magenta',marker='*',markersize=10,label='Salida inicial')

plt.title('Web2',fontsize=20)
plt.xlabel('AR',fontsize=15)
plt.ylabel('DEC',fontsize=20)
plt.grid(True)
plt.legend()

plt.subplot(224)

plt.plot(np.array(w4AR)[15:]*180/np.pi,np.array(w4DEC)[15:]*180/np.pi,linestyle='--',marker='o',color='magenta',label='Web3')
plt.plot(np.array(w4AR)[60]*180/np.pi,np.array(w4DEC)[60]*180/np.pi,c='blue',marker='*',markersize=10,label='Salida final')
plt.plot(np.array(w4AR)[15]*180/np.pi,np.array(w4DEC)[15]*180/np.pi,c='magenta',marker='*',markersize=10,label='Salida inicial')

plt.title('Web3',fontsize=20)
plt.xlabel('AR',fontsize=15)
plt.ylabel('DEC',fontsize=20)
plt.grid(True)
plt.legend()

plt.savefig('Particle/Tensor/NoAitoff.png')
plt.show()
plt.close()

plt.figure(figsize=(13.0, 10.0))
plt.scatter(np.array(w1AR)[15:]*180/np.pi,np.array(w1DEC)[15:]*180/np.pi,marker='o',s=13,color='red', label='SatNoPer')
plt.scatter(np.array(w2AR)[15:]*180/np.pi,np.array(w2DEC)[15:]*180/np.pi,marker='o',s=13,color='green',label='SatPer')
plt.scatter(np.array(w3AR)[15:]*180/np.pi,np.array(w3DEC)[15:]*180/np.pi,marker='o',s=13,color='blue',label='Web2')
plt.scatter(np.array(w4AR)[15:]*180/np.pi,np.array(w4DEC)[15:]*180/np.pi,marker='o',s=13,color='magenta',label='Web3')
plt.scatter(np.array(w1AR)[37]*180/np.pi,np.array(w1DEC)[37]*180/np.pi,c='black',marker='*',s=13,label='Salida 41132')
plt.scatter(np.array(w2AR)[37]*180/np.pi,np.array(w2DEC)[37]*180/np.pi,marker='*',c='black',s=13,label='Salida 41132')
plt.scatter(np.array(w2AR)[60]*180/np.pi,np.array(w2DEC)[60]*180/np.pi,marker='s',c='black',s=10,label='Salida final')
plt.scatter(np.array(w1AR)[60]*180/np.pi,np.array(w1DEC)[60]*180/np.pi,marker='s',c='black',s=10,label='Salida final')

plt.title('Cambios del eje principal asociado al autovalor menor del T.Inercia',fontsize=20)
plt.xlabel('AR',fontsize=20)
plt.ylabel('DEC',fontsize=20)
plt.legend()
plt.grid(True)

plt.savefig('Particle/Tensor/NoAitoff2.png')
plt.show()
plt.close()



