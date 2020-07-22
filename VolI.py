import struct
import numpy as np
import pickle
import csv
import os


def volcado(Ruta,SalidaSimulacion):

	valores=open(SalidaSimulacion,'rb')


	dummy=[0 for x in range(100)]
	ibuf=[0 for x in range(100)]
	ibuf1=[0 for x in range(100)]
	ibuf2=[0 for x in range(100)]
	Datos=[]
	Data=[]
	Data1=[]
	Data2=[]
	dummy2=[0 for x in range(100)]

	for i in range(100):
		Data.append([])
		buf=valores.read(20)
		dummy[i], ibuf[i], ibuf1[i], ibuf2[i], dummy2[i]=(struct.unpack('>5i',buf))
		Datos = dummy[i],ibuf[i],ibuf1[i] ,ibuf2[i],dummy2[i]
		Data[i].extend(Datos)
		

	np.save((Ruta+'/Data_ibuf') ,Data, allow_pickle=True,fix_imports=True)

	print('              ')
	print('/////////////////////////////////////////////////////////////////////////////')
	print('-----------------------Volcado y Partición de Datos-----------------------------------')
	print('//////////////////////////////////////////////////////////////////////////////')
	print('Nos interesan las variables: De IBUF1/ iteme - atime- padding - rtot . De IBUF2 / irun - nobj - ngas - ndark - h100- box100 - omega0 - xlambda0 - h0t0 -L ')
	print('           ')

	'''
	Para ibuf 1 se extraen las variables para el diccionario
	'''

	itime=ibuf1[0]
	atime=(struct.unpack('>f', struct.pack('>i',ibuf1[5] )))
	atime=atime[0]
	padding=(struct.unpack('>f', struct.pack('>i',ibuf1[19] )))
	padding=padding[0]
	rmtot=(struct.unpack('>f', struct.pack('>i',ibuf1[32] )))
	rmtot=rmtot[0]
	htime=(struct.unpack('>f', struct.pack('>i',ibuf1[6] )))
	htime=htime[0]


	print(f'Valor de itime {itime} | Valor de atime {atime} |Valor de padding {padding} | Valor de rtot {rmtot} | Valor de htime {htime}')
	'''

	Para ibuf2 se extraen las variables para el diccionario

	'''

	irun=ibuf2[0]
	nobj=ibuf2[1]
	ngas=ibuf2[2]
	ndark=ibuf2[3]
	L=ibuf2[4]
	h100=struct.unpack('>f', struct.pack('>i',ibuf2[12] ))
	h100=h100[0]

	box100=struct.unpack('>f', struct.pack('>i',ibuf2[13] ))
	box100=box100[0]

	omega0=struct.unpack('>f', struct.pack('>i',ibuf2[21] ))
	omega0=omega0[0]

	xlambda0=struct.unpack('>f', struct.pack('>i',ibuf2[22] ))
	xlambda0=xlambda0[0]

	h0t0=struct.unpack('>f', struct.pack('>i',ibuf2[23] ))
	h0t0=h0t0[0]

	print(f'Valor de irun {irun} | Valor de nobj {nobj} | Valor de ngas {ngas} | Valor de ndark {ndark} | Valor de h100 {h100} |'
	f' Valor de box100 {box100} | Valor de omega0 {omega0} | Valor de 	xlambda0 {xlambda0} | Valor de h0t0 {h0t0}')

	'''
	Creamos el diccionario con los valores y  los guardamos para su futuro uso en Conversion. 
	'''
	Values={'itime': itime, 'atime': atime, 'padding':padding , 'rmtot': rmtot, 'irun':irun, 'nobj':nobj, 'ngas':ngas, 'ndark':ndark,'h100':h100,'box100': box100, 'omega0':omega0,'xlambda0':xlambda0, 'h0t0':h0t0, 'L':L, 'htime':htime}
	w = csv.writer(open(Ruta+"/Values.csv", "w"))
	for key, val in Values.items():
	    w.writerow([key, val])


	print('                                                                                                                  ')

	print('--------------------------IBUF1--------------------------')
	print(ibuf1)

	print('--------------------------IBUF2--------------------------')
	print(ibuf2)

	print('--------------------------IBUF--------------------------')
	print(ibuf)

	'''
	print('--------------------------dummy--------------------------')
	print(dummy2)
	print(dummy)
	'''
	print('------------------------------------------------------------------------------------------------------------------------------')
	del Data, Datos, dummy, dummy2, ibuf


	m=[0 for x in range(nobj)]
	rx=[0 for x in range(nobj)]
	ry=[0 for x in range(nobj)]
	rz=[0 for x in range(nobj)]
	vx=[0 for x in range(nobj)]
	vy=[0 for x in range(nobj)]
	vz=[0 for x in range(nobj)]
	itype=[0 for x in range(nobj)]
	dummy=[]
	dummy2=[]
	Datos=[]
	buf=0
	con=0
	#Particle=''
	for i in range(nobj):
		Data1.append([])
		#Particle='P'+str(i+1)
		buf=valores.read(40)
		dummy,m[i], rx[i], ry[i], rz[i],vx[i],vy[i],vz[i],itype[i],dummy2=(struct.unpack('>i7f2i',buf))
		Datos = m[i], rx[i], ry[i], rz[i],vx[i],vy[i],vz[i],itype[i]
		Data1[i].extend(Datos)

	np.save(Ruta+'/Data_mrv',Data1, allow_pickle=True,fix_imports=True)
	del Data1, Datos

	h=[0 for x in range(nobj-ndark)]
	e=[0 for x in range(nobj-ndark)]
	dn=[0 for x in range(nobj-ndark)]
	t_star=[0 for x in range(nobj-ndark)]
	f_g=[0 for x in range(nobj-ndark)]

	for i in range((nobj-ndark)):
		Data2.append([])
		buf=valores.read(28)
		dummy,h[i],e[i],dn[i],t_star[i],f_g[i], dummy2= ((struct.unpack('>i5fi',buf)))
		Datos=h[i],e[i],dn[i],t_star[i],f_g[i]
		Data2[i].extend(Datos)

	np.save(Ruta+'/Data_bar',Data2, allow_pickle=True,fix_imports=True)

	del Data2, Datos
	valores.close()
