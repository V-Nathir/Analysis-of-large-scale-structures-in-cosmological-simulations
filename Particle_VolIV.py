from mpl_toolkits.mplot3d import Axes3D
from os import listdir
import pylab as pl
from pylab import *
import sys
import csv
import os
import pickle
import cv2
import matplotlib.animation as animation
import matplotlib.cm as cm

csv.field_size_limit(sys.maxsize)
'''
INFORMACIÓN :
               En las funciones mejor no usar w como variable ni i, son variables reservadas a la Ruta y Wdata. 
               Se usa longitud (ws webs -1) por no discriminar en las webs el primer elemento, que corresponde al número total de partículas. 
               WDAta no padece esto.

               REQUISITOS:
                	#  Necesita una carpeta llamada Particle y en ella una carpeta llAmada Webs/satelites con las webs/satélites correspondientes
                			(/CMSatelite)  (/Webs)
                	#  Una carpeta de Datos en Bruto (/DatosBruto) con los correspondientes datos de la simulación


FUNCIONAMIENTO DEL PROGRAMA

Los datos una vez tratados están preparados para utilizarse en este programa.Por ello no se puede ejecutar si no se ha empleado el programa Volcado.

Por el proceso de construcción del script reune código que no es ni necesario ni útil usar. Un ejemplo está el funcionamiento sin cargar webs que puedes 
seleccionar al iniciar el script. Por otra parte, hay métodos obsoletos que están almacenados en funciones a las que no se llaman. Este script no es una versión 
final, pero si está cerrada en cuanto fallos se refiere. 

Los CM pueden ser calculados por dos rutas. La recomendable y más completa es R, método de reducción, que puede consultarse en su correspondiente función.


Por eso explico aquí la opción w con la forma de cálculo de CM R, reducción. Cargar Webs y Sats.

																####	WEBS    ####

	#Se recluta la información de las partículas de los ficheros webs en la carpeta Webs. Se hace una criba de información y se almacena los datos de las 
	partículas implicadas, de interés. 
	Comienza el cálculo del CM con el método de reducción aplicando tres tolerancias. Se reescrie las variables con las nuevas posiciones y se representan.
	Hemos referido el vector r de posición al correspondiente CM, que será distinto en cada salida o redshift.
	Para comprobar que el cálculo es sensato, podemos apoyarnos en el script SeguimientoCM, que representa todos estos CM, mostrando un patrón suave, el 
	sistema se mueve y se redeistribuye como era de esperar.
	
	Finalmente se representa las estructuras afectadas por el método de reducción para comprobar la agresividad del método e ilustrarla. 

	Por tanto este proceso accede a carpetas con información de datos de la simulación, y genera unas nuevas datas para reducir las búsquedas. Crea directorios
	donde almacenar estas informaciones y volcar las representaciones.
	Este proceso es el eje estructurador del programa. En él se incorporan funciones complementarias y útiles como los diagramas de barras de los satélites

																#####	SATÉLITES ####

	# Es una etapa del código anterior al desarrollo principal de las webs, a nivel de organización del script.
	En esencia toma los valores de los centros de masas y un radio de la estructura de 35 satélites para una salida determinada. Esta lista esta dividida en
	satélites persistentes y no persistentes en lo que sus planos orbitales corresponde. NO implica que la estructura vaya a ser/no ser eliminada, aunque será una 
	tendencia en su categoría. 

	El programa crea un .txt para cada galaxia donde almacena las particulas que forman la estructura en base al radio. Con ellas calculará los CM de cada satelite
	en cada salida y su población en base a la lista primera generada. Hace representación en barras. Representa estos satélites, este es otro código integrado en el
	desarrollo de las WEBS. 

	FInalmente realiza el cálculo de CM interpretando a cada categoria como una web propia. 

		                	 
'''



# -------------------------------------------------------------------------------------------------------------------------------------
# ..........................................................#FUNCIONES#..................................................................
#--------------------------------------------------------------------------------------------------------------------------------------

def tensor(ST,bufDicSat,Satdic):
	'''
	#Tiene dos etapas esta función:
			# Cálculo del Cm en cada salida para estas partículas. Ver la Función CMreducc
			# Representacicón de histogramas para ver las poblaciones
		Esta función solo sirve para visualizar el número de galaxias de las webs (no)persistentes
		pero el objetivo es tratar a gran escala, no indivdualmente

	'''
	import os
	for tn in ST:
		cruta=0
		
		try:
			os.mkdir('Particle/'+tn+'/hist')
		except:
			print('#INFO: Carpeta para almacenar hist creada. ')
		try:
			os.mkdir('Particle/CMs')
		except:
			print('#INFO: Carpeta para almacenar los CMs de las cuatro webs creada. ')



	bufcmSat={}
	radios=[]
	for ten in ST:
		#primer cm para cada satélite en cada salida
		cruta=0
		for pa in ST[ten]:
			
			T=str(ST[ten][cruta,0].replace(',','.'))
			radios.append(ST[ten][cruta,3].replace(',','.'))
			bufcms=[]
			cruta=cruta+1
			cruta2=0
			for ruta in Ruta:
				bufcms.append([])
				cmT=0
				cmX=0
				cmY=0
				cmZ=0
				Q=Satdic[ruta]
				for r in (bufDicSat[T]):
					r=int(r)
					cmT=cmT+float(Q[r,0])
					cmX=cmX+(float(Q[r,1])*float(Q[r,0]))
					cmY=cmY+(float(Q[r,2])*float(Q[r,0]))
					cmZ=cmZ+(float(Q[r,3])*float(Q[r,0]))
				cmX=cmX/cmT
				cmY=cmY/cmT
				cmZ=cmZ/cmT
				cmS=[cmX,cmY,cmZ,cmT]
				bufcms[cruta2]=cmS
				cruta2=cruta2+1
			bufcmSat[T]=bufcms
	
	
	b=0
	cataclan=0
	for ten in bufcmSat:
		print(f'cambio al satélite con CMx = {ten}')
		print(f'Porcentaje del actual proceso de representación del diagrama de poblaciones:  {round(int(cataclan)*100/int(len(bufcmSat)))}% \n')
		cataclan=cataclan+1


		F=bufcmSat[ten]

		
		a=0
		elbucle=True
		bufhist=[]
		for ruta in (Ruta):
			#Aquí necesito hacer el CM 
			partbuf=[]
			bufhist.append([])
			tol=float(radios[b])

			while elbucle==True:
				Q=Satdic[ruta]
				partbuf=[]
				for r in bufDicSat[ten]:
					r=int(r)
					cMx=float(F[a][0])
					cMy=float(F[a][1])
					cMz=float(F[a][2])


					newcmsx=float(Q[r,1])-float(F[a][0])
					newcmsy=float(Q[r,2])-float(F[a][1])
					newcmsz=float(Q[r,3])-float(F[a][2])
					#Aquí abre un sat, una ruta y aplica la primera condición

					if np.absolute(newcmsx)<=tol and np.absolute(newcmsy)<= tol and np.absolute(newcmsz)<=tol:
						partbuf.append(int(r))
				newsatX=0
				newsatY=0
				newsatZ=0
				newsatM=0
				
				bufhist[a]=int(len(partbuf))
				if len(partbuf)==0:
					break

				for y in partbuf:
					newsatM=float(Q[y,0])+newsatM
					newsatX=float(Q[y,1])*float(Q[y,0])+newsatX
					newsatY=float(Q[y,2])*float(Q[y,0])+newsatY
					newsatZ=float(Q[y,3])*float(Q[y,0])+newsatZ
				newsatY=newsatY/newsatM
				newsatZ=newsatZ/newsatM
				newsatX=newsatX/newsatM
				if np.absolute(newsatX-cMx)<0.00099 and np.absolute(newsatY-cMy)<0.00099 and np.absolute(newsatZ-cMz)<0.00099:
					elbucle=False
					F[a][0]=newsatX
					F[a][1]=newsatY
					F[a][2]=newsatZ
				else:
					F[a][0]=newsatX
					F[a][1]=newsatY
					F[a][2]=newsatZ
					tol=tol-tol*0.35



			a=a+1
			elbucle=True
		b=b+1
		longit=np.round(np.linspace(0,61,61))
		R=plt.figure(figsize=(13.0, 10.0))
		ax=R.add_subplot(111)
		ax.bar(longit,bufhist)
		plt.xticks(longit,Ruta)
		plt.xticks(rotation=75)
		rects=ax.patches
		labels=[ i for i in bufhist]
		for rect,label in zip(rects,labels):
			height=rect.get_height()
			ax.text(rect.get_x() + rect.get_width() / 2, height + 5, label,ha='center', va='bottom',rotation=90)
		ax.set_ylabel('Número de Partículas')
		ax.set_title('Población del satélite de CMx para la salida 41132: '+ten)
		if b>21:
			plt.savefig('Particle/sat_persistentes_41132.dat/hist/Vida_'+str(ten)+'.png')
		if b<=21:
			plt.savefig('Particle/sat_nopersistentes_41132.dat/hist/Vida_'+str(ten)+'.png')
		plt.close()	
	print('Representación de las poblaciones de las estructuras locales, galaxias : COMPLETADO\n -------------------------------------------\n')
	MegaSat={}
	for tn in ST:
		megabuf=[]
		F=open('Particle/'+str(tn)+'/web.txt','r')
		fsat=F.read()
		fsat=np.array(fsat.split())
		for t in range((len(fsat))):
			megabuf.append([])
			megabuf[t]=fsat[t]
		MegaSat[tn]=np.array(megabuf)
	for tn in MegaSat:
		print('Creando .txt con los CM de la web creada con los'+ tn)
		recopcmsat=[]
		recop=0
		for ruta in Ruta:
			print(f'Proceso:  {round((recop+1)*100/int(len(Ruta)))}%')
			print(ruta)
			cmMsx=0
			cmMsy=0
			cmMsz=0
			cmMsm=0
			F=np.array(Satdic[ruta])
			Ffin=True
			estre=0
			gass=0
			otros=0
			for n in MegaSat[tn]:
				n=np.int(n)

				if int(F[n,7])==1:
					gass=gass+1
				if int(F[n,7])==-1:
					estre=estre+1
				if int(F[n,7])==2:
					otros=otros+1	
				cmMsm=cmMsm+(float(F[n,0]))
				cmMsx=cmMsx+(float(F[n,1])*float(F[n,0]))
				cmMsy=cmMsy+(float(F[n,2])*float(F[n,0]))
				cmMsz=cmMsz+(float(F[n,3])*float(F[n,0]))
			cmMsz=cmMsz/cmMsm
			cmMsy=cmMsy/cmMsm
			cmMsx=cmMsx/cmMsm
			tolsat=0.05
			var=0.0003
			print(f'SALIDA {ruta}  DE {tn} : \n Estrellas -->{estre} Gas --> {gass} Otros-->{otros}')

			while Ffin==True:
				

				newcmMsm=0
				newcmMsx=0
				newcmMsy=0
				newcmMsz=0
				bsat=[]
				for n in MegaSat[tn]:
					n=np.int(n)
					newcmMsx=float(F[n,1])-cmMsx
					newcmMsy=float(F[n,2])-cmMsy
					newcmMsz=float(F[n,3])-cmMsz
					if np.absolute(newcmMsx)<=tolsat and np.absolute(newcmMsy)<=tolsat and np.absolute(newcmMsz)<=tolsat:
						bsat.append(int(n))
				newcmMsm=0
				newcmMsx=0
				newcmMsy=0
				newcmMsz=0
				print(f'tolerancia tolsat {tolsat}')
				for n in bsat:
					n=np.int(n)
					newcmMsm=newcmMsm+(float(F[n,0]))
					newcmMsx=newcmMsx+(float(F[n,1])*float(F[n,0]))
					newcmMsy=newcmMsy+(float(F[n,2])*float(F[n,0]))
					newcmMsz=newcmMsz+(float(F[n,3])*float(F[n,0]))
				newcmMsx=newcmMsx/newcmMsm
				newcmMsy=newcmMsy/newcmMsm
				newcmMsz=newcmMsz/newcmMsm
				if np.absolute(newcmMsx-cmMsx)<=0.00099 and np.absolute(newcmMsy-cmMsy)<=0.00099 and np.absolute(newcmMsz-cmMsz)<=0.00099 :
					cmMsx=newcmMsx
					cmMsy=newcmMsy
					cmMsz=newcmMsz
					cmMsm=newcmMsm
					tolsat2=tolsat
					tolsat2=tolsat-2.1*var*1.1
					print(f'tolerancia 2 {tolsat2}')
					sx=0
					sy=0
					sz=0
					sm=0
					lol=0

					for n in bsat:
						n=np.int(n)
						newcmMsx=float(F[n,1])-cmMsx
						newcmMsy=float(F[n,2])-cmMsy
						newcmMsz=float(F[n,3])-cmMsz
						if np.absolute(newcmMsx)<=tolsat2 and np.absolute(newcmMsy)<=tolsat2 and np.absolute(newcmMsz)<=tolsat2:
							sm=sm+float(F[n,0])
							lol=lol+1
							sx=sx+float(F[n,1])*float(F[n,0])
							sy=sy+float(F[n,2])*float(F[n,0])
							sz=sz+float(F[n,3])*float(F[n,0])
					print(f'Partículas despreciadas para la cota inferior(tolerancia) {int(len(MegaSat[tn]))-lol} por tanto contribuyen {lol}' )
					sx=sx/sm
					sy=sy/sm
					sz=sz/sm
					if np.absolute(cmMsx-sx)<=0.00099 and np.absolute(cmMsy-sy)<=0.00099 and np.absolute(cmMsz-sz)<=0.00099 :
						Ffin=False
						recopcmsat.append([])
						recopcmsat[recop]= cmMsx, cmMsy, cmMsz , cmMsm, float(int(len(bsat))), float(int(len(MegaSat[tn]))) 
						recop=recop+1
						print(f'Partículas despreciadas para la salida {ruta}: {int(len(MegaSat[tn])-int(len(bsat)))} de {len(MegaSat[tn])} partículas \n Por tanto la población para el cálculo es {int(len(bsat))} partículas\n ')
					else:
						var=var*1.1
						tolsat=0.08-var
						
				else:

					tolsat=0.08-var
					var=var*1.1
					cmMsx=newcmMsx
					cmMsy=newcmMsy
					cmMsz=newcmMsz
		recopcmsat=np.array(recopcmsat)
		print(f'Centros de masas para {tn} "COMPLETADO".\n ')
		np.savetxt('Particle/CMs/'+str(tn)+'_web.txt',recopcmsat,fmt='%s')
	cabecera=['|CMx|','|CMy|','|CMz|','|CMm|','|Número de partículas implicadas en CM|','|Total Web|']
	cabe=open('Particle/CMs/Cabecera.txt','w+')
	for cab in cabecera:
		cabe.write('%s  ' %  cab )
	cabe.close()
	print('Cálculo y almacenamiento de los CMs para los satélites: COMPLETADO\n -----------------------------------------------')
	

def PlotSat(ST):
	'''
	Se cargan los valores posicionales, masa, tipo y velocidades en un diccionario. Se incorpora mediante ST los valores de CM y radios para averiguar 
	para la salida 41132 , las partículas que se encuentran en cada región. Se guardan sus valores que los identifican, el número de partícula. A su vez
	se guarda los valores para dicha salida de esas partículas.
	Se usará para PLot

	Solo se usa si no hay datos previos.
	'''
	import os
	import shutil
	global Satdic
	Satdic={}
	
	for rut in range(len(Ruta)):
		F=np.load('d5004/'+ Ruta[rut]+'/Data_mrv.npy')
		print(f'Salida: {Ruta[rut]} || Creando diccionario con los datos en bruto -> {round((rut+1)*100/int(len(Ruta)))}%')
		Satdic[Ruta[rut]]=np.array(F)
	lon=len(Satdic)
	for sat in (ST):
		print('Rastreando partículas para los satélites de: '+sat +' ...'  )
		WebsSat=[]
		websatcon=0
		for s in range(lon):
			if 'd5004.41132'==str(Ruta[s]):
				F=(Satdic[Ruta[s]])
				for a in range(len(ST[sat][:,0])):
					S=ST[sat]
					try:
						os.mkdir('Particle/'+str(sat)+'/CMx_'+(str(S[a,0]).replace(',','.')))
					except:
						shutil.rmtree('Particle/'+str(sat)+'/CMx_'+(str(S[a,0]).replace(',','.')))
						os.mkdir('Particle/'+str(sat)+'/CMx_'+(str(S[a,0]).replace(',','.')))
					satbuf=[]
					o=0
					satcon=0
					particles=[]
					print(f'Analizando...{round((a+1)*100/int(len(ST[sat][:,0])))}%')
					for t in range(len(F[:,0])):

						if int(F[t,7])!=0 and np.absolute(float(F[t,1])-float(str(S[a,0]).replace(',','.')))<=float(str(S[a,3]).replace(',','.')) and np.absolute(float(F[t,2])-float(str(S[a,1]).replace(',','.')))<=float(str(S[a,3]).replace(',','.')) and np.absolute(float(F[t,3])-float(str(S[a,2]).replace(',','.')))<=float(str(S[a,3]).replace(',','.')):
								particles.append([])
								WebsSat.append([])
								satbuf.append([])
								particles[o]=str(t)
								WebsSat[websatcon]=int(t)
								satbuf[o]=F[t,:]
								o=o+1
								websatcon=websatcon+1
					satbuf=np.array(satbuf)
					np.save('Particle/'+str(sat)+'/CMx_'+((str(S[a,0])).replace(',','.'))+'/'+str(Ruta[s])+'_'+str(satcon)+'_NP'+str(o),satbuf, allow_pickle=True,fix_imports=True)
					file=open("Particle/"+str(sat)+'/CMx_'+((str(S[a,0])).replace(',','.'))+"/ListaParticles_CMx_"+((str(S[a,0])).replace(',','.'))+"_Sat"+str(satcon)+"_NP"+str(o)+".txt",'w+')
					for elemento in particles:
						file.write('%s \n' % elemento)
					file.close()
					satcon=satcon+1
		file=open('Particle/'+str(sat)+'/web.txt','w+')
		for elemento in WebsSat:
			file.write('%s \n' % elemento)
		file.close()
		print('Proceso para los satélites /' +str(sat)+ '/ terminado.')
		print('')
	print('RASTREO COMPLETADO')

	print('\n############################################' +'\n$:')
	del o, a , s,  t, satbuf, satcon

 

def CMreducc(WData,Webs,bufpartSat):
	'''
	Consiste en calcular un centro de masas imponiendo unas tolerancias y reduciendo el número de partículas en base a ese criterio para obtener un mejor CM. 
	Una vez conseguido un CM estable, se recolocan las partículas y se representan. Esta función contiene partes sacadas de Centromasas
	Variables prohibidas i y w

	Se plotea: las partículas discriminadas para el cálculo de CM, las webs. 

	EN DETALLE: Hay tres controles o tolerancias. Se calcula un CM primero con el que se recalculan unas nuevas posiciones de las partículas (newxyz)
				respecto al CM.

		# Tolerancia variable: Una función refina el valor de esta tolerancia si las (newxyz), que no han superado esta tolerancia (la cumplen), no generan un 
		nuevo CM que difiere menos del cuarto decimal (tolerancia fija .00099)respecto al primer ( o anterior calculado) CM. Empieza desde el 10% de la caja.
		Es una función suficientemente suave para no provocar un desprecio de partículas muy elevado ni muy corto.Se ha visto su comportamiento y probado con funciones prueba. 

		# Tolerancia fija: Tolerancia que debe de cumplir un centro de masas para coonvertirse en candidato del CM real-estable.Ese CM se calcula con las newxyz que han respetado la tol-var. 
		Es del orden del cuarto decmal (.00099).

		# Tolerancia de cota inferior: Cuando se tiene un CM candidato se debe someter a una tolerancia que se  calcula con dos refinamientos más de la tol-var., es decir, 
		incorpora una cota inferior, dos ciclos más pequeña=más partículas despreciadas, que nuevamente no puede diferir del cuarto decimal(tol-fija) respecto al candidato. Éste da el sello de un CM estable. 
		
		SI NO SE CONSIGUE UN CM VÁLIDO, YA SEA PORQUE NO CONSIGA SER CANDIDATO, O PORQUE EL CANDIDATO NO CUMPLE LA COTA INFERIOR, SE REDUCE LA TOL-VAR. AUMENTA EL 
		NÚMERO DE PART. DESPRECIADAS.
	Recordatorio Ws[Webs[w]] contiene las particulas de cad web , pero también tiene el valor del número de part en la primera fila. De ahí len-1 y Ws[Webs[w]][0]


	'''
	global CM
	global CMzero

	Mtot=0
	mrx=0
	mry=0
	mrz=0
	for j in range(len(Ws[Webs[w]])-1):
		Mtot= float(WData[j,0])+ Mtot
		mrx= (float(WData[j,1])*float(WData[j,0]))+mrx
		mry=(float(WData[j,2])*float(WData[j,0]))+mry
		mrz=(float(WData[j,3])*float(WData[j,0]))+mrz
	CMx=mrx/Mtot
	CMy=mry/Mtot
	CMz=mrz/Mtot
	CM=[CMx, CMy, CMz]
	print('Centro de masas '+ str(CM))

	
	bufglobal=True
	tol=0.1
	tolfix=0.1

	while bufglobal==True:
		'''
		partbuf sirve para almacenar las partículas que no se discriminan. Partbuf2 tiene almecenadas las elimandas. Partbuf3 es igual que partbuf pero con un grado 
		de profundidad en la tolerancia más, es desechable. 
		'''
		partbuf=[]
		partbuf2=[]
		partbuf3=[]
		print(f'Reajustando el CM con una tol de posición del {round(tol,6)}' )

		for y in range(len(Ws[Webs[w]])-1):
			bufdel=False
			newx=(float(WData[y,1])-CMx)
			newy=(float(WData[y,2])-CMy)
			newz=(float(WData[y,3])-CMz)

			if np.absolute(float(newx))>tol:
				bufdel=True
				partbuf2.append(int(y))
			if np.absolute(float(newy))>tol:
				bufdel=True
				partbuf2.append(int(y))
			if np.absolute(float(newz))>tol:
				bufdel=True
				partbuf2.append(int(y))
			if bufdel==False:
				partbuf.append(int(y))	
		print(f'Partículas despreciadas {int((Ws[Webs[w]][0]))-int(len(partbuf))} de {int(Ws[Webs[w]][0])}')
		newCMx=0
		newCMy=0
		newCMz=0
		newMtot=0
		for y in partbuf:
			y=int(y)
			newMtot= float(WData[y,0])+ newMtot
			newCMx= (float(WData[y,1])*float(WData[y,0]))+newCMx
			newCMy=(float(WData[y,2])*float(WData[y,0]))+newCMy
			newCMz=(float(WData[y,3])*float(WData[y,0]))+newCMz
		newCMx=newCMx/newMtot
		newCMy=newCMy/newMtot
		newCMz=newCMz/newMtot

		if np.absolute(newCMx-CMx)<0.00099 and np.absolute(newCMy-CMy)<0.00099 and np.absolute(newCMz-CMz)<0.00099:
			print('Candidato de CM: CMx '+str(newCMx)+' CMy ' +str(newCMy)+' CMz ' + str(newCMz) +' Pendiente de ajuste con una cota de tolerancia inferior.')

			CMx=newCMx
			CMy=newCMy
			CMz=newCMz
			CMM=newMtot
			cmlon=float(int(len(partbuf)))
			cmlong=float(int(Ws[Webs[w]][0]))
			CM=[CMx, CMy, CMz]
			CMzero=[CMx, CMy, CMz,CMM,cmlon,cmlong]
			'''
			Tol2 es un grado más profundo que tol, y así poder ajustar por debajo el CM y comprobar que es estable realmente.
			'''
			tol2=tol-tolfix*0.005
			tol2=tol2-tolfix*0.005

			for y in partbuf:
				'''
				Se hace una nueva criba de particulas para conseguir ajuste por debajo
				'''
				y=int(y)
				newx=(float(WData[y,1])-CMx)
				newy=(float(WData[y,2])-CMy)
				newz=(float(WData[y,3])-CMz)
				if float(newx)<tol2 and float(newy)<tol2 and float(newz)<tol2:
					partbuf3.append(int(y))
			print(f'Partículas despreciadas con la cota inferior de tolerancia {round(tol2,6)}:  {int(Ws[Webs[w]][0])-int(len(partbuf3))}')
			newCMx=0
			newCMy=0
			newCMz=0
			newMtot=0
			for y in partbuf3:
				y=int(y)
				newMtot= float(WData[y,0])+ newMtot
				newCMx= (float(WData[y,1])*float(WData[y,0]))+newCMx
				newCMy=(float(WData[y,2])*float(WData[y,0]))+newCMy
				newCMz=(float(WData[y,3])*float(WData[y,0]))+newCMz
			newCMx=newCMx/newMtot			
			newCMy=newCMy/newMtot
			newCMz=newCMz/newMtot	
			print('CM con dos grado más de profundidad en la tolerencia: CMx '+str(newCMx)+' CMy ' +str(newCMy)+' CMz ' + str(newCMz))
			if np.absolute(newCMx-CMx)<0.00099 and np.absolute(newCMy-CMy)<0.00099 and np.absolute(newCMz-CMz)<0.00099:
				bufglobal=False
				print(f'Reajustado el CM a {CM} con una tolerancia de 0.00099 y tolerancia de posición de {round(tol,6)}')


			else:
				print('Fallo en el ajuste de cota inferior, con dos grados de tolerancia mayor. Reajustando')
				bufglobal=True
				tol=tol-tolfix*0.005
		else:
			CMx=newCMx
			CMy=newCMy
			CMz=newCMz
			CM=[CMx, CMy, CMz]

			bufglobal=True
			#Ver gráficas de la tol
			tol=tol-tolfix*0.005



	'''
	Una vez reajustado el CM se reajustan las posiciones y se corrigen las desplazadas por llegar al límite de la caja
	'''
	
	print()
	print('Reajustando posiciones de las partículas con el CM')
	newx=[0 for x in range(len(Ws[Webs[w]]))]
	newy=[0 for x in range(len(Ws[Webs[w]]))]
	newz=[0 for x in range(len(Ws[Webs[w]]))]
	red=True
	contadoreajuste=0
	while red ==True:
		red=False

		for j in range(len(Ws[Webs[w]])-1):
			newx[j]=(float(WData[j,1])-CMx)
			newy[j]=(float(WData[j,2])-CMy)
			newz[j]=(float(WData[j,3])-CMz)
			if np.absolute(newx[j])>np.absolute((float(WData[j,1])+1)-CMx):
				red=True
				WData[j,1]=float(WData[j,1]+1)
				contadoreajuste=contadoreajuste+1
			if np.absolute(newy[j])>np.absolute((float(WData[j,2])+1)-CMy):
				red=True
				WData[j,2]=float(WData[j,2]+1)
				contadoreajuste=contadoreajuste+1

			if np.absolute(newz[j])>np.absolute((float(WData[j,3])+1)-CMz):
				red=True
				WData[j,3]=float(WData[j,3]+1)
				contadoreajuste=contadoreajuste+1


			#para poder conocer si hay partículas que por la periocidad de la caja se han cambiado en pos, siendo <1. Cm porque esta en la región 
			if np.absolute(newx[j])>np.absolute((float(WData[j,1])-1)-CMx):
				red=True
				WData[j,1]=float(WData[j,1]-1)
				contadoreajuste=contadoreajuste+1

			if np.absolute(newy[j])>np.absolute((float(WData[j,2])-1)-CMy):
				red=True
				WData[j,2]=float(WData[j,2]-1)
				contadoreajuste=contadoreajuste+1

			if np.absolute(newz[j])>np.absolute((float(WData[j,3])-1)-CMz):
				red=True
				WData[j,3]=float(WData[j,3]-1)
				contadoreajuste=contadoreajuste+1
	print(f'Ha habido {contadoreajuste} reajustes de posiciones.')
	del contadoreajuste


	
	#Vamos a representar las partículas discriminadas y las contribuyentes al CM. ALso el CM
	partbufx=[]
	partbufy=[]
	partbufz=[]
	partX=[]
	partY=[]
	partZ=[]
	for p in partbuf:
		k=int(p)
		partbufx.append(float(newx[k]))
		partbufy.append(float(newy[k]))
		partbufz.append(float(newz[k]))
	for p in partbuf2:
		k=int(p)
		partX.append(float(newx[k]))
		partY.append(float(newy[k]))
		partZ.append(float(newz[k]))   



	R=plt.figure(figsize=(13.0, 10.0))
	ax=R.add_subplot(111,projection='3d')

	ax.scatter(partbufx,partbufy,partbufz, marker='o',c='green',s=6, alpha=0.1,label='Partículas que contribuyen')
	ax.scatter(partX,partY,partZ,  marker='o',c='red',s=6, alpha=0.3,label='Partículas discriminadas')
	try:

		CMX=[]
		CMY=[]
		CMZ=[]
		lon=int(len(CMtot))
		A=np.linspace(0,(lon-3),(lon/3))
		B=np.linspace(1,(lon-2),(lon/3))
		C=np.linspace(2,(lon-1),(lon/3))
		for x in range(len(A)):
			a=round(float(A[x]))
			b=round(float(B[x]))
			c=round(float(C[x]))
			CMX.append(float(CMtot[a])-float(CMx))
			CMY.append(float(CMtot[b])-float(CMy))
			CMZ.append(float(CMtot[c])-float(CMz))
		ax.scatter(CMX,CMY,CMZ,marker='*',c='magenta',s=50,label='CM anteriores')
		ax.scatter(0,0,0, marker='*',c='black',s=60,label='CM de la salida '+ str(Ruta[i]))
	except:
		ax.scatter(0,0,0, marker='*',c='black',s=50,label='CM de la salida '+ str(Ruta[i]))
	ax.legend()
	ax.view_init(elev=30, azim=70)

	ax.set_xlabel('X')
	ax.set_title(Ruta[i])
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')

	#plt.show()
	plt.savefig('Particle/CM/CMReduc/'+Webs[w]+'/A_'+Ruta[i]+'.png')

	ax.view_init(elev=30,azim=250)
	plt.savefig('Particle/CM/CMReduc/'+Webs[w]+'/B_'+Ruta[i]+'.png')


	plt.close()
	del partbufx, partbufy, partbufz , partZ, partX, partY
	
	#Con el bufpartsat, que reune las partículas de los satélites, vamos a reposicionarlas con el CM del satélite central. Y luego se las representará
	newsatx=[]
	newsaty=[]
	newsatz=[]
	newsatx2=[]
	newsaty2=[]
	newsatz2=[]
	newconsat2=0
	newconsat=0
	cambiosat=0
	for p in bufpartSat:
		if cambiosat>20:
			for pp in range(len(bufpartSat[p][:,0])):
				newsatx2.append([])
				newsaty2.append([])
				newsatz2.append([])
				newsatx2[newconsat2]=float(bufpartSat[p][pp,1])-float(CMx)
				newsaty2[newconsat2]=float(bufpartSat[p][pp,2])-float(CMy)
				newsatz2[newconsat2]=float(bufpartSat[p][pp,3])-float(CMz)
				newconsat2=newconsat2+1

		if cambiosat<=20:
			for pp in range(len(bufpartSat[p][:,0])):
					newsatx.append([])
					newsaty.append([])
					newsatz.append([])
					newsatx[newconsat]=float(bufpartSat[p][pp,1])-float(CMx)
					newsaty[newconsat]=float(bufpartSat[p][pp,2])-float(CMy)
					newsatz[newconsat]=float(bufpartSat[p][pp,3])-float(CMz)
					newconsat=newconsat+1
		cambiosat=cambiosat+1


	#Vamos a hacer paquetes con el tipo de partículas para el código de color. Representando toda la estructura 

	con1=0
	con_1=0
	con0=0

	q=0
	k=0
	p=0
	gen=False
	for y in range(len(Ws[Webs[w]])-1):
		if gen==False:
			#Para ver la extensión de cada vector en base al tipo de partícula. Tipo 1 , -1 y 2 
			for j in range(len(Ws[Webs[w]])-1):
				if int(round(float(WData[j,7])))==int(1):
					con1=con1+1
				if int(round(float(WData[j,7])))==int(-1):
					con_1=con_1+1
				if int(round(float(WData[j,7])))==int(2):
					con0=con0+1
			print(f'Número de Particulas de gas : {con1}')
			print(f'Número de Estrellas : {con_1}')
			print(f'Número de Transformaciones : {con0}')
			print(f'Total Partículas : {con1+con_1+con0}')
			gen=True
			x1=[0 for x in range((con1))]
			y1=[0 for x in range((con1))]
			z1=[0 for x in range((con1))]
			x_1=[0 for x in range((con_1))]
			y_1=[0 for x in range((con_1))]
			z_1=[0 for x in range((con_1))]
			x0=[0 for x in range(con0)]
			y0=[0 for x in range(con0)]
			z0=[0 for x in range(con0)]


			
		if int(round(float(WData[y,7])))==int(1):
			x1[q]=float(newx[y])
			y1[q]=float(newy[y])
			z1[q]=float(newz[y])
			q=q+1
		if int(round(float(WData[y,7])))==int(-1):
			x_1[p]=(float(newx[y]))
			y_1[p]=(float(newy[y]))
			z_1[p]=(float(newz[y]))
			p=p+1
		if int(round(float(WData[y,7])))==int(2):
			x0[k]=float(newx[y])
			y0[k]=float(newy[y])
			z0[k]=float(newz[y])
			k=k+1

	fig= plt.figure(figsize=(13.0, 10.0))
	ax=fig.add_subplot(111,projection='3d')
	
	ax.scatter(x1,y1,z1,marker='o',c='blue',s=6, alpha=0.1,label='Gas')
	ax.scatter(x_1,y_1,z_1,marker='o',c='magenta',alpha=0.2,s=6,label='Estrellas')
	ax.scatter(x0,y0,z0, marker='o',c='black',s=6,alpha=0.1,label='Transformaciones')
	ax.legend(framealpha=1,edgecolor='black', markerscale=4)
	ax.set_xlabel('X')
	ax.set_title(Ruta[i])
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')

	#plt.show()
	plt.savefig('Particle/Plot/'+Webs[w]+'/'+Ruta[i]+'.png')

	ax.view_init(elev=30,azim=250)
	#ax.dist=3
	plt.savefig('Particle/Plot/'+Webs[w]+'/B_'+Ruta[i]+'.png')
	plt.close()
	fig= plt.figure(figsize=(13.0, 10.0))
	ax=fig.add_subplot(111,projection='3d')
	ax.scatter(0,0,0,marker='o',c='black',s=6,label=CM)
	ax.scatter(newsatx2,newsaty2,newsatz2, marker='^', c='green',alpha=0.05,s=6, label='Satélites persistentes')
	ax.scatter(newsatx,newsaty,newsatz, marker='^', c='red',alpha=0.05,s=6, label='Satélites no persistentes')
	#ax.legend(framealpha=1,edgecolor='black', markerscale=4)
	ax.dist=3
	#if Ruta[i]=='d5004.68502':
		#plt.show()
	plt.savefig('Particle/Plot/'+Webs[w]+'/Sat_'+Ruta[i]+'.png')
	ax.dist=6
	plt.savefig('Particle/Plot/'+Webs[w]+'/Sat2_'+Ruta[i]+'.png')

	plt.close()



def plotCM():
	'''
	Se nutre de las otras funciones, utiliza un fichero guardado, y representa las modificaciones de la posicion del CM. Recopila sus diferentes posiciones. 
	Se ha tenido que ajustar a mano la division de A B C.
	Representa todas las posiciones del CM para ver su comportamiento a lo largo de los redshift. Así verificar que su comportamiento es estable y suave.
	'''
	G=np.load('Particle/CM/Memo_CM'+str(w)+'.npy','r+')
	CM=np.array(G)
	CMX=[]
	CMY=[]
	CMZ=[]
	A=np.linspace(0,180,61)
	B=np.linspace(1,181,61)
	C=np.linspace(2,182,61)

	for x in range(len(A)):
		a=round(float(A[x]))
		b=round(float(B[x]))
		c=round(float(C[x]))
		CMX.append(CM[a])
		CMY.append(CM[b])
		CMZ.append(CM[c])

	
	figCM= plt.figure(figsize=(13.0, 10.0))
	ax=figCM.add_subplot(111,projection='3d')
	colors = cm.rainbow(np.linspace(0, 1, len(CMX)))
	
	ax.view_init(30, 45)

	for x in range(len(CMX)):
		c=[colors[x]]
		ax.scatter(CMX[x],CMY[x],CMZ[x],marker='o',c=c)
	ax.scatter(CMX[0],CMY[0],CMZ[0],marker='*',s=50,c='black',label='Inicio')
	ax.scatter(CMX[60],CMY[60],CMZ[60],marker='*',s=50,c='gray',label='Final')
	ax.legend()
	ax.set_xlabel('X')
	ax.set_title('Variación del Centro de Masas para el método de reajuste.' )
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')

	#plt.show()
	plt.savefig('Particle/CM/CM_Reajuste_web'+str(w)+'.png')
	plt.close()
	


def plots(WData,Webs):
	'''
	Primera función para representar las webs. Carece de control sobre el CM y no reajusta la posicion de partículas
	Para activarse debe de ser desde código, está muteada. NO sirve. Meramente ilustrativa.
	'''
	fig=plt.figure()
	ax=fig.add_subplot(111,projection='3d')
	#for j in range(len(Webs[w])):
	
	
	ax.scatter(WData[:,1],WData[:,2],WData[:,3])
	ax.set_xlabel('X')
	ax.set_title(Ruta[i])
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	plt.savefig('Particle/Plot/'+Webs[w]+'/'+Ruta[i]+'.png')
	#plt.show()
	plt.close()


def Id(j,Ws):
	'''
	Crea un diccionario con las nuevas numeraciones de las partículas en los ficheros webs. Actúa una vez. Utilidad poca.
	'''
	
	asig={}
	asig_con=0

	if j!=int(Ws[Webs[w]][0]):
		if asig_con==0:
			asig={j:int(asig_con)}
			asig_con=+1
		else:
			asig[j]=int(asig_con)
		asig_con=asig_con+1
	t=csv.writer(open('Particle/ID/Particles_ID '+ str(Webs[w]),"w"))
	for key, val in asig.items():
		t.writerow([key,val])


def centromasas(WData,Ws):
	'''
	Calcula, reajusta y reposicina tanto el CM como las partículas de manera paralela. Según calcula un nuevo CM bajo una tolerancia, se reposiciona 
	las partículas, con esos nuevos valores se obtiene un candidato a CM. Si no cumple una tolerancia, se vuelve a calcular el CM y las respectivas posiciones nuevas
	que servirán para el próximo candidato.
	La tolerancia esta basada en la unidad. Es decir, suponemos un comportamiento suficientemente estable, porque sino este método no funciona, no 
	se consigue un CM estable, casos satélites no permanentes.El comportamiento estable solo nos hace preocuparnos por las particulas deslocalizadas por
	superar los valores de la caja, 0 a 1, que por periocidad de box se recolocan.Modo pacman.
	Aquí vemos si su posición actual o la modificada en +- unidad, en alguna dirección, está más cerca del CM. Si es necesario se modifica y recalcula el CM, mejorándolo.
	Los resultados arrojan que no hay muchas partículas. El resultado es válido con el método por reducción. 
	Variables prohibidas i y w
	'''
	global CM

	Mtot=0
	mrx=0
	mry=0
	mrz=0
	for j in range(len(Ws[Webs[w]])-1):
		Mtot= float(WData[j,0])+ Mtot
		mrx= (float(WData[j,1])*float(WData[j,0]))+mrx
		mry=(float(WData[j,2])*float(WData[j,0]))+mry
		mrz=(float(WData[j,3])*float(WData[j,0]))+mrz
	CMx=mrx/Mtot
	CMy=mry/Mtot
	CMz=mrz/Mtot
	CM=[CMx, CMy, CMz]
	print('Centro de masas '+ str(CM))
	newx=[0 for x in range(len(Ws[Webs[w]]))]
	newy=[0 for x in range(len(Ws[Webs[w]]))]
	newz=[0 for x in range(len(Ws[Webs[w]]))]
	red=True
	CMrecupx=[]
	CMrecupy=[]

	CMrecupz=[]

	while red ==True:

		red=False
		print('Reajustando CM')
		for j in range(len(Ws[Webs[w]])-1):
			CMrecupx.append(CMx)
			CMrecupy.append(CMy)
			CMrecupz.append(CMz)
			newx[j]=float(WData[j,1])-CMx
			newy[j]=(float(WData[j,2])-CMy)
			newz[j]=(float(WData[j,3])-CMz)
			#para poder conocer si hay partículas que por la periocidad de la caja se han cambiado en pos, siendo >1. Cm porque esta en la región 
			if np.absolute(newx[j])>np.absolute((float(WData[j,1])+1)-CMx):
				CMx=CMx-(float(WData[j,1])*float(WData[j,0]))/Mtot+((float(WData[j,1])+1)*float(WData[j,0]))/Mtot
				red=True
				WData[j,1]=float(WData[j,1]+1)
			if np.absolute(newy[j])>np.absolute((float(WData[j,2])+1)-CMy):
				CMy=CMy-(float(WData[j,2])*float(WData[j,0]))/Mtot+((float(WData[j,2])+1)*float(WData[j,0]))/Mtot
				red=True
				WData[j,2]=float(WData[j,2]+1)

			if np.absolute(newz[j])>np.absolute((float(WData[j,3])+1)-CMz):
				red=True
				CMx=CMx-(float(WData[j,3])*float(WData[j,0]))/Mtot+((float(WData[j,3])+1)*float(WData[j,0]))/Mtot
				WData[j,3]=float(WData[j,3]+1)

			#para poder conocer si hay partículas que por la periocidad de la caja se han cambiado en pos, siendo <1. Cm porque esta en la región 
			if np.absolute(newx[j])>np.absolute((float(WData[j,1])-1)-CMx):
				CMx=CMx-(float(WData[j,1])*float(WData[j,0]))/Mtot+((float(WData[j,1])-1)*float(WData[j,0]))/Mtot
				red=True
				WData[j,1]=float(WData[j,1]-1)

			if np.absolute(newy[j])>np.absolute((float(WData[j,2])-1)-CMy):
				CMy=CMy-(float(WData[j,2])*float(WData[j,0]))/Mtot+((float(WData[j,2])-1)*float(WData[j,0]))/Mtot
				red=True
				WData[j,2]=float(WData[j,2]-1)


			if np.absolute(newz[j])>np.absolute((float(WData[j,3])-1)-CMz):
				red=True
				CMx=CMx-(float(WData[j,3])*float(WData[j,0]))+((float(WData[j,3])-1)*float(WData[j,0]))/Mtot
				WData[j,3]=float(WData[j,3]-1)





	CM=[CMx, CMy, CMz]
	print('Centro de masas ajustado : '+ str(CM))

	
	#Vamos a hacer paquetes con el tipo de partícula para asignar un código de color

	con1=0
	con_1=0
	con0=0

	q=0
	k=0
	p=0
	gen=False
	for y in range(len(Ws[Webs[w]])-1):
		if gen==False:
			#Para ver la extensión de cada vector en base al tipo de partícula. Tipo 1 , -1 y 2 
			for j in range(len(Ws[Webs[w]])-1):
				if int(round(float(WData[j,7])))==int(1):
					con1=con1+1
				if int(round(float(WData[j,7])))==int(-1):
					con_1=con_1+1
				if int(round(float(WData[j,7])))==int(2):
					con0=con0+1
			print(f'Número de Particulas de gas : {con1}')
			print(f'Número de Estrellas : {con_1}')
			print(f'Número de Transformaciones : {con0}')
			print(f'Total Partículas : {con1+con_1+con0}')
			gen=True
			x1=[0 for x in range((con1))]
			y1=[0 for x in range((con1))]
			z1=[0 for x in range((con1))]
			x_1=[0 for x in range((con_1))]
			y_1=[0 for x in range((con_1))]
			z_1=[0 for x in range((con_1))]
			x0=[0 for x in range(con0)]
			y0=[0 for x in range(con0)]
			z0=[0 for x in range(con0)]


			
		if int(round(float(WData[y,7])))==int(1):
			x1[q]=float(newx[y])
			y1[q]=float(newy[y])
			z1[q]=float(newz[y])
			q=q+1
		if int(round(float(WData[y,7])))==int(-1):
			x_1[p]=(float(newx[y]))
			y_1[p]=(float(newy[y]))
			z_1[p]=(float(newz[y]))
			p=p+1
		if int(round(float(WData[y,7])))==int(2):
			x0=float(newx[y])
			y0=float(newx[y])
			z0=float(newx[y])
			k=k+1
		print(q)
		print(p)
		print(k)
	fig1=plt.figure(figsize=(13.0,10.0))
	ax=fig1.add_subplot(111,projection='3d')
	ax.plot(CMrecupx,CMrecupy,CMrecupz,marker='o')
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	ax.set_title('Desviación del Centro de masas. Salida: '+ Ruta[i])
	plt.savefig('Particle/CM/'+Webs[w]+'/'+Ruta[i]+'.png')
	plt.close()

	fig= plt.figure(figsize=(13.0, 10.0))
	ax=fig.add_subplot(111,projection='3d')
	
	ax.scatter(x1,y1,z1,marker='.',c='blue',s=6, alpha=0.3,label='Gas')
	ax.scatter(x_1,y_1,z_1,marker='o',c='magenta',s=6,label='Estrellas')
	ax.scatter(x0,y0,z0, marker='o',c='green',s=6,label='Transformaciones')
	ax.legend()

	ax.set_xlabel('X')
	ax.set_title(Ruta[i])
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')

	#plt.show()
	plt.savefig('Particle/Plot/'+Webs[w]+'/'+Ruta[i]+'.png')
	plt.close()


#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------











#                   ------------------------------------------      Cabecera      --------------------------------------------
Webs=[]
SAT=[]
Ws={}
ST={}
terminador=False
funciones=[]

bif=input('Este es el programa para hacer el seguimiento de la posición de partículas a diferentes redshift.\n Para cargar una web de partículas introduzca "w". \n Para elegir un conjunto de partículas a mano introduzca "p" : ')
if bif!='w' and bif!='p':
	print('NO has introducido una opción válida. Terminador activado')
	terminador=True
print()
if bif=='w':
	funciones=input('Elija un tipo de representación en base al cálculo del CM.\n Para CM por reducción, se reduce la estructura de partículas hasta obtener un CM estable, "R" (recomendado). \n Para CM reajustado paralelamente con las posiciones de las particulas, "r" :')
	print('\n')
	if funciones!='R' and funciones!='r':
		print('NO has introducido una opción válida. Terminador activado')
		terminador=True
#Busca las webs y muestras algunas particulas.Crea un diccionario con las webs mostrando solo las particulas implicadas.
#WS diccionario almacena el nombre de las webs como clave y su contenido como valor. NO SE HA DESCRIMINADO EL PRIMER ELEMENTO INFORMATIVO 
#También cargamos el nombre de los ficheros satélites y sus CMxyz y el radio para englobar las particulas 


if bif=='w' and funciones=='r' or funciones=='R':
	for web in listdir('./Particle/Webs'):
		Webs.append(web)

	for sat in listdir('./Particle/CMSatelite'):
		SAT.append(sat)
	
	if Webs:
		print('Se han encontrado las siguientes webs: ')
		print(Webs)
		for i in range(len(Webs)):
			F=open('Particle/Webs/'+Webs[i],'r')
			Web1=F.read()
			Web1=np.array((Web1.split()))
			Ws[Webs[i]]=Web1
			print('Algunas partículas seleccionadas de :' + Webs[i] + '\n')
			print(Ws[Webs[i]][1:50])
	if SAT:
		print('Se han encontrado ficheros con CM de satélites:\n // {NOTA}: Estos ficheros tienen el vector [CM_x,CM_y,CM_z] y el radio de la región de selección de particuas.\n')
		for i in range(len(SAT)):
			'''
			Vamos a crear el diccionario ST que almacene las posiciones de los cm de los sotelites y sus radios. 
			'''
			k=open('Particle/CMSatelite/'+SAT[i],'r')
			ksat=k.read()
			ksat=np.array((ksat.split()))
			SATCM=[]
			for t in range(int(len(ksat)/4)):
				t=t*4
				CMx_sat, CMy_sat, CMz_sat, radio = ksat[t:t+4]
				SATCM.append([])
				SATCM[int(t/4)]=CMx_sat, CMy_sat, CMz_sat, radio
			SATCM=np.array(SATCM)
			ST[SAT[i]]=SATCM
			del SATCM
			print('Centromasas de los satélites y el radio para la región de selección pertenecientes al fichero: '+ SAT[i])



	else:	
		bif=input('No se han encontrado la información necesaria. Revisa las carpetas de Webs y CMSatelite. \n Cargar partículas a su elección "y" or "any_key": ')
		if bif!='y':
			terminador=True




Ruta=[]
for out in listdir("./DatosBruto"):
	Ruta.append(out)
Ruta.sort()


rep=True
c=0
ls=[]
CMtot=[]
satcontrol=False
uncontrolador=0
if bif=='w' and terminador==False:

	for s in ST:
		try:
			os.mkdir('Particle/'+s)
		except:
			print('  #INFO: La carpeta para los CM de los satélites //'+ s+'// ya está creada.')
	
		for sar in range(len(ST[s][:,0])):
			if os.path.isdir('Particle/'+s+'/CMx_'+(str(ST[s][sar,0]).replace(',','.'))):
				satcontrol=satcontrol
			else:
				satcontrol=True
	if satcontrol==False:
		Satdic={}
		for rut in range(len(Ruta)):
			F=np.load('d5004/'+ Ruta[rut]+'/Data_mrv.npy')
			print(f'Salida: {Ruta[rut]} || Creando diccionario con los datos en bruto-> {round((rut+1)*100/int(len(Ruta)))}%')
			Satdic[Ruta[rut]]=np.array(F)
	if satcontrol==True:
		PlotSat(ST)

		#Se ha definido Satdicc como variable global

tensorf=True	
if bif=='w' and terminador==False:
	for w in range(len(Webs)):
		try: 
			os.mkdir('Particle/Plot')
		except:
			print('  #INFO: Directorio "Particle/Plot" -> True.')
		try:
			os.mkdir('Particle/Plot/'+Webs[w])
		except:
			print(f'  #INFO: Directorio "Particle/Plot/"{Webs[w]} -> True.')

		try:
			os.mkdir('Particle/CM/')
		except:
			print('  #INFO: Directorio Particle/CM -> True.')
		
		try:
			os.mkdir('Particle/CM/CMReduc')
		except:
			print('  #INFO: Directorio "Particle/CM/CMReduc"-> True')
		try:
			os.mkdir('Particle/CM/CMReduc/'+Webs[w])
		except:
			print(f'  #INFO: Directorio "CM/CMReduc/"{Webs[w]} -> True')
		try:
			os.remove('Particle/CM/Memo_CM'+str(w)+'.npy')

		except:
			print('  #INFO: Se creará un fichero con las posiciones del CM  para cada salida en "Particle/CM".')
		try:
			os.mkdir('Particle/ID')
		except:
			print('  #INFO: Carpeta ID creada. Esta carpeta almacena un diccionario CVS meramente informativo. Asocia el número de párticula con la posición dentro de las varables del programa.\n ------------------')
		try:
			os.mkdir('Particle/'+Webs[w])
		except:
			print('#Particle/'+str(Webs[w])+' ya creado')

		
		print('\n //////////// \n Cargando y Representando partículas.')
		'''
		Este primer bucle crea un diccionario bufDicSat que guarda las partículas obtenidas en la 
		función PLotSat. ENtra en cada carpeta de sat persis. y no persis. a cada directorio y carga solo 
		la extensión .txt, la almaecana y se puede comprobar con un print(bufDicSat['0.278060'][5:10])
		 # El objetivo es aprovechar el bucle en Ruta para ir accedidendo al dicc Satdic y filtrar particulas
		 	crear con ello un nuevo diccionario con el que acabremos representando en la función CMreducc
		 	cada satelite con las webs 
		'''
		bufDicSat={}
		conbuf=0
		bufscm=[]
		for s in ST:
			for satt in range(len(ST[s][:,0])):
				txt_files = [fsa for fsa in os.listdir('./Particle/'+s+'/CMx_'+((ST[s][satt,0]).replace(',','.'))) if fsa.endswith('.txt')]
				Datasat1=open('Particle/'+s+'/CMx_'+((ST[s][satt,0]).replace(',','.'))+'/'+txt_files[0],'r')
				Datasat=Datasat1.read()
				Datasat=np.array(Datasat.split())
				bufDicSat[str(((ST[s][satt,0]).replace(',','.')))]=Datasat
				Datasat1.close()
		if tensorf==True:
			tensor(ST,bufDicSat,Satdic)
			tensorf=False
		print('\n #INFO: INFORMACIÓN DEL LISTADO DE PARTÍCULAS DE LOS SATÉLITES ALMACENADO.\n ')
		for i in range(len(Ruta)):
			#Este primer bucle crea un diccionario con las particulas de los satelites, pero no de manera permanente,es un buffer. Cada resultado del diccionario
			#se pasa a la función CMreducc y se realiza su representación. Sí tenemos el seguimiento de las particulas en cada salida representado. Acomplado al func del programa
			bufpartSat={}
			for s in bufDicSat:
				bufsintetic=[]
				contadorsintetic=0
				for posicion in bufDicSat[s]:
					bufsintetic.append([])
					bufsintetic[contadorsintetic]=Satdic[Ruta[i]][int(posicion),:]
					contadorsintetic=contadorsintetic+1

				bufpartSat[s]=np.array(bufsintetic)
			try:
				C=np.load('Particle/'+Webs[w]+'/'+Ruta[i]+'.npy')
				WData=np.array(C)
				del C

			except:
				#Si no existe archivos np array con la información de las partículas se crea 
				if uncontrolador<4:
					y=input(f'Error.NO existe los datos  de la  web {Webs[w]} de {Ruta[i]}. Pulsa para generarlos.')
				else:
					print('Generando de manera automática los ficheros')
				contador=0
				uncontrolador=uncontrolador+1
				V=np.load('d5004/'+ Ruta[i]+'/Data_mrv.npy')
				memo=[]
				#Toma el valor de la partícula y lo encuentra en los datos formados con los brutos
				for j in Ws[Webs[w]]:
					j=int(j)
					#si es disintinto del primer valor de la lista, que no es una particula
					if j!= int(Ws[Webs[w]][0]):
						memo.append([])
						memo[contador][:]=V[j-1][:]
						contador=contador+1
				del V
				memo=np.array(memo)
				np.save('Particle/'+Webs[w]+'/'+Ruta[i],memo, allow_pickle=True,fix_imports=True)
				onebucle=True
				WData=memo
				del memo

			print('------------------------------------------------')
			print(f'Progreso al {round((i+1)*100/int(len(Ruta)))}%')
			print(Ruta[i])
			print()
			#plots(WData,Webs)
			if funciones=='r':
				centromasas(WData,Ws)
			if funciones=='R':

				CMreducc(WData,Webs,bufpartSat)
			bufscm.append([])
			bufscm[conbuf]=CMzero
			conbuf=conbuf+1

			print()
			
			try:
				CMtot=np.append(CMtot,CM,axis=0)
			except:
				CMtot=CM
			del WData, CM
		bufscm=np.array(bufscm)
		np.savetxt('Particle/CMs/'+str(Webs[w])+'_web.txt',bufscm,fmt='%s')
		np.save('Particle/CM/Memo_CM'+str(w)+'.npy',CMtot,allow_pickle=True,fix_imports=True)
		del CMtot
		plotCM()

		print('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||')
		print('Cargado y representado '+ Webs[w]+' .En la carpeta Plots' )


		#Aquí se intenta cargar o se crea un diccionario para poder saber posición-partícula en WData
		Particle={}
		try:
			for key,val in csv.reader(open('Webs/ID/Particle_ID '+str(Webs[w]))):
				Particle[key]=val
		except:
			for j in Ws[Webs[w]]:
				j=int(j)
				Id(j,Ws)
		print('||||||||||||||||||||||||||||||||||||||||||||||||||||||||||')
		print('Se ha cargado los ID de las partículas de la '+ Webs[w] )





		

			







			

#----------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------Casco Antiguo--------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------
#Programa alternativo. Representación de partículas a elección del usuario. Aún sin modificar. No útil.
b=False
if bif=='p' or bif=='y' and terminador==False:
	while rep==True:

		l=(input('Elige una partícula.Min 1 . Max 625994: '))
		try:
			ls.append(str(int(l)))
		except ValueError:
			print('No es un entero.')
			b=True
			break
		memo=[]
		for i in range(len(Ruta)):
			memo.append([])
			F=np.load('d5004/'+Ruta[i]+'/Data_mrv.npy' )
			memo[i][:]=F[int(l)-1][:]
			
		memo = np.array(memo)
		if c==0:
			Historial={'Partícula '+str(l):memo}
		rep_control=''
		while rep_control!='y' and  rep_control!='n':	
			rep_control=input('¿Desea cargar otra partícula? Introduce "y" or "n": ')
		if rep==True and c!=0:
			Historial['Partícula '+str(l)]=memo
		c=+1
		if rep_control=='n':
			rep=False

	if b==False:
		for i in range(len(ls)):
			print('||||||||||||||||||||||||||||')
			print(f'Partícula {ls[i]}')
			print('||||||||||||||||||||||||||||')
			print(Historial['Partícula '+ls[i]])


		con_graf=input('Imprimir "y" . Para no Imprimir "any_key": ')
		c=0
		if con_graf=='y':
			fig = pl.figure()
			mpl.rcParams['legend.fontsize'] = 10
			ax = Axes3D(fig)
			ax.set_title('Partículas y sus trayectorias')

			ax.set_xlabel('X')
			ax.set_ylabel('Y')
			ax.set_zlabel('Z')
			for i in range(len(ls)):
				RX=np.array(Historial['Partícula '+ls[i]][:,1])
				RY=np.array(Historial['Partícula '+ls[i]][:,2])
				RZ=np.array(Historial['Partícula '+ls[i]][:,3])
				ax.plot(RX,RY,RZ, label='Posiciones para la partícula ' + str(ls[i]))
				ax.legend()
		pl.show()





