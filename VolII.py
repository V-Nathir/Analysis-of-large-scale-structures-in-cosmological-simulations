import numpy as np 
import csv
from os import listdir
import math 

def Conversion():

	Mpc=3.086e24
	G=6.67e-8
	Mpc=3.086e24
	H0=1e7/Mpc
	yr=3.16e7
	kb=1.38e-16
	mum=1.22*1.672649e-24
	Msun=2e33

	Ruta=[]
	for out in listdir("./DatosBruto"):
		Ruta.append(out)
	Ruta.sort()
	for j in range(len(Ruta)):
		print('|||||||||||||||||||||||||Cambio de salida|||||||||||||||||||||||')
		print(Ruta[j])
		Values = {}
		for key, val in csv.reader(open('d5004/'+Ruta[j]+ "/Values.csv")):
		    Values[key] = val
		print(Values)

		rhoc=3*(H0**2)*(float(Values['h100'])**2)/(8.*math.pi*G)
		box=float((Values['box100']))/(float(Values['h100']))
		fl=float(Values['L'])-2.*float(Values['padding'])
		lunit=box*Mpc/fl
		munit=(rhoc*Mpc)*(Mpc/Msun)*Mpc*box**3/(float(Values['rmtot'])*float(Values['omega0']))
		tunit=float(Values['h0t0'])/H0/float(Values['h100'])
		vunit=lunit/tunit
		eunit=vunit**2
		dnuit=lunit*(lunit/Msun)*(lunit/munit)
		Kunit=(eunit*2.*mum)/3./kb

		convL=lunit*fl*float(Values['atime'])/Mpc  #Conversion a Mpc
		convM=munit/1e10          #Conversion a 10^10 M_o
		convV=vunit*fl*float(Values['atime'])/100000. #Conversion a kms^-1
		convT=Kunit*(fl*float(Values['atime']))**2 #Conversion a K
		g43=G*(Msun/Mpc)    #*(1e10/1e10)
		gup=g43/convL*convM/convV**2
		rhoc28=(rhoc*Mpc)*(Mpc/Msun)*(Mpc/1e10)
		rhocup0=rhoc28*convL**3/convM
		rhocup=rhocup0*(fl*float(Values['htime']))**2/float(Values['h0t0'])**2


		SI_Val={'rhoc':rhoc,'box':box,'fl':fl,'lunit':lunit, 'munit':munit, 'tunit':tunit,'vunit':vunit,'eunit':eunit,'dnuit':dnuit,'Kunit':Kunit,
		'convL':convL,'convM':convM,'convV':convV,'g43':g43,'gup':gup, 'rhoc28':rhoc28,'rhocup0':rhocup0,'rhocup':rhocup  }
		w = csv.writer(open('d5004/'+Ruta[j]+"/SI_Val.csv", "w"))
		for key, val in SI_Val.items():
		    w.writerow([key, val])
		print('Parámetros con el cambio de unidades:		')
		print(SI_Val)
Conversion()
