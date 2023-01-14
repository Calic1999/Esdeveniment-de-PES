#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 12:50:19 2022

@author: rogerprat
"""

import matplotlib.pyplot as plt
from importlib import reload
import numpy as np
import math

mp = 1.67262192369E-27
kB = 1.3806504E-23

#Reading files
filename = "/Users/rogerprat/Documents/SynologyDrive/Màster/Space-based astronomy and Space Weather/2-Meteo/Esdeveniment/dades/wind_swe_2m_bhw5zaBDsN.lst.txt"
WIND = open(filename,"r")
WIND_DATA = []
for line in WIND:
    stripped = line.strip()
    DATA = stripped.split(' ')
    while("" in DATA):
        DATA.remove("")
    DATA[0] = int(DATA[0])
    for i in range(1,len(DATA)):
        DATA[i] = float(DATA[i])
    DATA[3] = (mp*(DATA[3]*1000)**2/(2*kB))/100000
    WIND_DATA.append(DATA)
WIND.close()

filename = "/Users/rogerprat/Documents/SynologyDrive/Màster/Space-based astronomy and Space Weather/2-Meteo/Esdeveniment/dades/ACE_EPAM_Data-3.txt"
ACE = open(filename,"r")
ACE_DATA=[]
read = False
for line in ACE:
    if read:
        stripped = line.strip()
        DATA = stripped.split(' ')
        while("" in DATA):
            DATA.remove("")
        for i in range(len(DATA)):
            DATA[i] = float(DATA[i])
        ACE_DATA.append(DATA)
    if line == 'BEGIN DATA\n':
        read = True
ACE.close()

filename = "/Users/rogerprat/Documents/SynologyDrive/Màster/Space-based astronomy and Space Weather/2-Meteo/Esdeveniment/dades/SEPEM_RD_24February2014SEPevent.txt"
SEPEM = open(filename,"r")
SEPEM_DATA=[]
read = False
read2 = False
for line in SEPEM:
    if read and read2:
        stripped = line.strip()
        splitted = stripped.split(',')
        DATA=[]
        for value in splitted:
            DATA.append(value.strip())
        DATE= DATA[0].split(' ')
        Any,Mes,Dia = DATE[0].split('-')
        Hora,Minut,Segon = DATE[1].split(':')
        DOY = 31 + int(Dia) + int(Hora)/24 + int(Minut)/(24*60) + int(Segon)/(24*3600)
        if Mes=='03':
            DOY += 28
        DATA[0] = DOY
        for i in range(1,len(DATA)):
            DATA[i] = float(DATA[i])
        SEPEM_DATA.append(DATA)
    if read:
        read2 = True
    if line == '200.0 -289.2  MeV\n':
        read = True
SEPEM.close()


year,WIND_doy,v,T,d,Bx,By,Bz = zip(*WIND_DATA)
ACE_doy,P1,P2,P3,P4,P5,P6,P7,P8 = zip(*ACE_DATA)
SEPAM_doy,interval,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19 = zip(*SEPEM_DATA)


#Remove bad data from ACE
def BadACE(DOY,PF):
    DOY = list(DOY)
    PF = list(PF)
    N = len(DOY)
    index = []
    for i in range(N):
        if PF[i]==-999.9:
            index.append(i)
    j = 0
    for i in range(len(index)):
        index[i] -= j
        DOY.pop(index[i])
        PF.pop(index[i])
        j += 1
    return tuple(DOY),tuple(PF)

ACE_doy_1,P1 = BadACE(ACE_doy,P1)
ACE_doy_2,P2 = BadACE(ACE_doy,P2)
ACE_doy_3,P3 = BadACE(ACE_doy,P3)
ACE_doy_4,P4 = BadACE(ACE_doy,P4)
ACE_doy_5,P5 = BadACE(ACE_doy,P5)
ACE_doy_6,P6 = BadACE(ACE_doy,P6)
ACE_doy_7,P7 = BadACE(ACE_doy,P7)
ACE_doy_8,P8 = BadACE(ACE_doy,P8)

#Remove bad data from WIND
def BadWIND(DOY,v,T):
    DOY_v = list(DOY)
    DOY_T = list(DOY)
    v = list(v)
    T = list(T)
    N = len(DOY)
    index_v = []
    index_T = []
    for i in range(N):
        if v[i]==99999.9:
            index_v.append(i)
        if T[i]>=3:
            index_T.append(i)
    j = 0
    for i in range(len(index_v)):
        index_v[i] -= j
        DOY_v.pop(index_v[i])
        v.pop(index_v[i])
        j += 1
    j = 0
    for i in range(len(index_T)):
        index_T[i] -= j
        DOY_T.pop(index_T[i])
        T.pop(index_T[i])
        j += 1
    return tuple(DOY_v),tuple(v),tuple(DOY_T),tuple(T)

DOY_v,v,DOY_T,T = BadWIND(WIND_doy,v,T)


#Plotting
fig, axs = plt.subplots(7, 1, sharex=True, height_ratios=[4,1,1,1,1,1,1])
fig.subplots_adjust(hspace=0) # Remove horizontal space between axes

axs[0].plot([],[],' ', label='ACE/EPAM')
axs[0].plot(ACE_doy_1,P1, label='0.047-0.068 MeV')
axs[0].plot(ACE_doy_2,P2, label='0.068-0.115 MeV')
axs[0].plot(ACE_doy_3,P3, label='0.115-0.195 MeV')
axs[0].plot(ACE_doy_4,P4, label='0.195-0.321 MeV')
axs[0].plot(ACE_doy_5,P5, label='0.321-0.580 MeV')
axs[0].plot(ACE_doy_6,P6, label='0.587-1.06 MeV')
axs[0].plot(ACE_doy_7,P7, label='1.06-1.90 MeV')
axs[0].plot(ACE_doy_8,P8, label='1.90-4.80 MeV')
axs[0].plot([],[],' ', label='GOES/EPS')
axs[0].plot(SEPAM_doy,P9, label='5.00-7.23 MeV')
axs[0].plot(SEPAM_doy,P10, label='7.23-10.46 MeV')
axs[0].plot(SEPAM_doy,P11, label='10.46-15.12 MeV')
axs[0].plot(SEPAM_doy,P12, label='15.12- 21.87 MeV')
axs[0].plot(SEPAM_doy,P13, label='21.87-31.62 MeV')
axs[0].plot(SEPAM_doy,P14, label='31.62-45.73 MeV')
axs[0].plot(SEPAM_doy,P15, label='45.73-66.13 MeV')
axs[0].plot(SEPAM_doy,P16, label='66.13-95.64 MeV')
axs[0].plot(SEPAM_doy,P17, label='95.64-138.3  MeV')
axs[0].plot(SEPAM_doy,P18, label='138.3-200.0  MeV')
axs[0].plot(SEPAM_doy,P19, label='200.0-289.2  MeV')
axs[0].set_ylabel('Flux de protons p/(cm$^{2}$ s sr MeV)')
axs[0].set_xlim(55,62)
axs[0].set_yscale('log')
axs[1].plot(DOY_v,v,c='k')
axs[1].set_ylabel('Velocitat [km/s]')
axs[2].plot(DOY_T,T,c='k')
axs[2].set_ylabel('T [$10^5$ K]')
axs[3].plot(WIND_doy,d,c='k')
axs[3].set_ylabel('Densitat [cm$^{-3}$]')
axs[4].plot(WIND_doy,Bx,c='k')
axs[4].set_ylabel("$B_x$ [nT]")
axs[5].plot(WIND_doy,By,c='k')
axs[5].set_ylabel("$B_y$ [nT]")
axs[6].plot(WIND_doy,Bz,c='k')
axs[6].set_ylabel("$B_z$ [nT]")
axs[6].set_xlabel("Dia de l'any 2014")

axs[0].legend(bbox_to_anchor=(1, 0.84), loc='upper left')

fig.set_size_inches(13, 20)
#fig.set_dpi(1000)

pathname='/Users/rogerprat/Documents/SynologyDrive/Màster/Space-based astronomy and Space Weather/2-Meteo/Esdeveniment/plot.pdf'
fig.savefig(pathname, bbox_inches='tight')

#Interplanetary shock passage
IPS_P = 58.658

axs[0].axvline(x=IPS_P,color='gray',linestyle='dashed')
axs[1].axvline(x=IPS_P,color='gray',linestyle='dashed')
axs[2].axvline(x=IPS_P,color='gray',linestyle='dashed')
axs[3].axvline(x=IPS_P,color='gray',linestyle='dashed')
axs[4].axvline(x=IPS_P,color='gray',linestyle='dashed')
axs[5].axvline(x=IPS_P,color='gray',linestyle='dashed')
axs[6].axvline(x=IPS_P,color='gray',linestyle='dashed')
plt.show()
#pathname='/Users/rogerprat/Documents/SynologyDrive/Màster/Space-based astronomy and Space Weather/2-Meteo/Esdeveniment/plot_IPSP.png'
#fig.savefig(pathname, bbox_inches='tight')

hora = int((IPS_P - 58)*24)
minut = int(((IPS_P - 58)*24 - hora)*60)
segon = int((((IPS_P - 58)*24 - hora)*60 - minut)*60)
data = '27-02-2014 ' + str(hora) + ':' + str(minut) + ':' + str(segon)
print('Passatge del xoc interplanetari: ', data)
print("Dia de l'any: ", IPS_P)


def FindShockE(PF,DOY):
    for i in DOY:
        if i>IPS_P:
            get_index = DOY.index(i)
            break
    E = PF[get_index]
    return E

Energies = [0.0575,
            0.0915,
            0.155,
            0.258,
            0.4505,
            0.8235,
            1.48,
            3.35,
            6.115,
            8.845,
            12.79,
            18.495,
            26.745,
            38.675,
            55.93,
            80.885,
            116.97,
            169.15,
            244.6]



Intensitat = []
Intensitat.append(FindShockE(P1,ACE_doy_1))
Intensitat.append(FindShockE(P2,ACE_doy_2))
Intensitat.append(FindShockE(P3,ACE_doy_3))
Intensitat.append(FindShockE(P4,ACE_doy_4))
Intensitat.append(FindShockE(P5,ACE_doy_5))
Intensitat.append(FindShockE(P6,ACE_doy_6))
Intensitat.append(FindShockE(P7,ACE_doy_7))
Intensitat.append(FindShockE(P8,ACE_doy_8))
Intensitat.append(FindShockE(P9,SEPAM_doy))
Intensitat.append(FindShockE(P10,SEPAM_doy))
Intensitat.append(FindShockE(P11,SEPAM_doy))
Intensitat.append(FindShockE(P12,SEPAM_doy))
Intensitat.append(FindShockE(P13,SEPAM_doy))
Intensitat.append(FindShockE(P14,SEPAM_doy))
Intensitat.append(FindShockE(P15,SEPAM_doy))
Intensitat.append(FindShockE(P16,SEPAM_doy))
Intensitat.append(FindShockE(P17,SEPAM_doy))
Intensitat.append(FindShockE(P18,SEPAM_doy))
Intensitat.append(FindShockE(P19,SEPAM_doy))

#Regressió
x = ([])
y = ([])
for i in Energies:
    x = np.append(x,np.log(i))
for i in Intensitat:
    y = np.append(y,np.log(i))

N = len(x)
xquadrat = x**2
sumadexquadrat = xquadrat.sum()
Nsuma = sumadexquadrat*N
sumadex = x.sum()
sumadexalquadrat = sumadex**2

Delta = Nsuma - sumadexalquadrat

xy = x*y
sumadexy = xy.sum()
Nsumadexy = N*sumadexy
sumadey = y.sum()
xsumadey = x*sumadey
sumaxsumay = xsumadey.sum()

a2 = (Nsumadexy - sumaxsumay)/(Delta)

sumaxquadratsumay = sumadexquadrat*sumadey
xsumadexy = x*sumadexy
sumaxsumaxy = xsumadexy.sum()

b2 = (sumaxquadratsumay - sumaxsumaxy)/(Delta)

recta = y - a2*x - b2
rectaquadrat = recta**2
suma = rectaquadrat.sum()
dyreg = (suma/(N - 2))**(1/2)

arrela = (N/Delta)**(1/2)
arrelb = (sumadexquadrat/Delta)**(1/2)

da2 = dyreg*arrela
db2 = dyreg*arrelb

xbarra = sumadex/N
xmenysxbarra = x - xbarra
xmenysxbarraalquadrat = xmenysxbarra**2
sumaxmenysxbarraalquadrat = xmenysxbarraalquadrat.sum()

ybarra = sumadey/N
ymenysybarra = y - ybarra
ymenysybarraalquadrat = ymenysybarra**2
sumaymenysybarraalquadrat = ymenysybarraalquadrat.sum()


factor=0
for i in range(0,N):
    factor = factor + (x[i]-xbarra)*(y[i]-ybarra)

r = factor/(np.sqrt(sumaxmenysxbarraalquadrat*sumaymenysybarraalquadrat))
R = r**2

dif=1-R
def orderOfMagnitude(number):
    return math.floor(math.log(number, 10))
order = orderOfMagnitude(dif)

def ajust(z):
    return a2*z + b2

oa=orderOfMagnitude(da2)
ob=orderOfMagnitude(db2)

ar2 = a2.round(abs(oa))
br2 = b2.round(abs(ob))
dar2 = da2.round(abs(oa))
dbr2 = db2.round(abs(ob))

dR2 = R.round(abs(order)+1)
print('')
print("Ajust de l'índex espectral:")
print('a =', ar2, ' ± ', dar2, ' b =', np.exp(br2), ' ± ', np.exp(dbr2))
#print('dyreg = ',dyreg)
print('R^2 =', dR2)
print('')

plt=reload(plt)
plt.figure()
plt.scatter(Energies,Intensitat,label='Dades',c='k')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('E (MeV)')
plt.ylabel('Intensitat p/(cm$^{2}$ s sr MeV)')
plt.annotate("$I = (221.4\pm1.2) E^{-2.16\pm0.06}$\n $R^2$=0.987",(0.11,0.01))
plt.grid(True)




#DSA theory
nd = 42.3
nu = 17
dnu = 0.9
dnd = 1.1
r = nd/nu
dr = np.sqrt((dnd/nu)**2+(nd/(nu**2))**2*dnu**2)
q=(3*r)/(r-1)
dq = 3*dr*(1-r)/((r-1)**2)
gamma = (1-q)/2
dgamma = -dq/2
print('Acceleració difusiva del xoc:')
print(round(gamma,2), ' ± ', round(dgamma,2))
print('')

def ajust2(z):
    return gamma*z + 5

plt.plot(Energies, [np.exp(ajust2(np.log(i))) for i in Energies],label='ADX', color="blue")
plt.plot(Energies, [np.exp(ajust(np.log(i))) for i in Energies],label='Ajust', color="gray")
plt.legend(loc='upper right')
#fig.set_dpi(1000)
#pathname='/Users/rogerprat/Documents/SynologyDrive/Màster/Space-based astronomy and Space Weather/2-Meteo/Esdeveniment/ajust.pdf'
#plt.savefig(pathname)

