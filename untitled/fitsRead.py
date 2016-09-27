import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from astropy.utils.data import download_file
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from astropy import units as u
from astropy import constants as const
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import math
import pickle
from sortedcontainers import SortedDict, SortedList, SortedSet

def R_Bound(Lambda):
    M200 = 10**((14.499) + (1.33)* math.log10(Lambda/40))
    R200m = ((1.432*10**(-11))*M200)**(1/3)*u.Mpc
    return R200m*3

def isPair(j, k, sortedList, clusterData, coordDict, thetaBoundDict, pairDict, pairIndice):
    if math.fabs(clusterData[j]['Z_LAMBDA']-clusterData[k]['Z_LAMBDA']) <= 0.01:
        #coord1 = SkyCoord(clusterData[j]['RA'], clusterData[j]['DEC'], frame='icrs', unit="deg")
        #coord2 = SkyCoord(clusterData[k]['RA'], clusterData[k]['DEC'], frame='icrs', unit="deg")
        #thetaBound = R_Bound(clusterData[j]['Z_LAMBDA']) / cosmo.angular_diameter_distance(clusterData[j]['Z_LAMBDA'])
        #separationAngle = u.radian.to(coord1.separation(coord2), equivalencies=u.dimensionless_angles())
        if u.radian.to(coordDict[j].separation(coordDict[k]), equivalencies=u.dimensionless_angles()) < thetaBoundDict[j]:
            pairDict[pairIndice].add((j,k, clusterData[j]['ID'], clusterData[k]['ID']))
            #print("cluster bound pair found! ID 1 = %d ID 2 = %d, Thetabound = %f\n" % (clusterData[j]['ID'], clusterData[k]['ID'], thetaBoundDict[j]))

cosmo = LambdaCDM(H0= 70, Om0= 0.3, Ode0= 0.7)
clusters = fits.open('input1.fits', memmap=True)
clusterData = clusters[1].data
print(clusters[1].columns)

galaxyDict = SortedDict()
coordDict = SortedDict()
thetaBoundDict = SortedDict()
bucket = 0.00
for i in range(0,25000):
    print(i)
    coordDict[i] = SkyCoord(clusterData[i]['RA'], clusterData[i]['DEC'], frame='icrs', unit="deg")
    thetaBoundDict[i] = R_Bound(clusterData[i]['Z_LAMBDA']) / cosmo.angular_diameter_distance(clusterData[i]['Z_LAMBDA'])
    bucket = round(clusterData[i]['Z_LAMBDA'], 2)
    if bucket in galaxyDict:
        galaxyDict[bucket].add(i)
    else:
        galaxyDict[bucket] = SortedList()
        galaxyDict[bucket].add(i)

dictKeyList = galaxyDict.keys()
pairDict = SortedDict()
for i in range(0, len(dictKeyList)):
    print(len(galaxyDict[dictKeyList[i]]))
    pairDict[dictKeyList[i]] = SortedList()

#for i in range(0,len(dictKeyList)):
for i in range(27,29):
    print(i)
    bucketList = galaxyDict[dictKeyList[i]]
    for j in range(0, len(bucketList)-1):
        for k in range(j+1, len(bucketList)):
            #is j paired with k from this bucket?
            isPair(bucketList[j], bucketList[k], pairDict, clusterData, coordDict, thetaBoundDict, pairDict, i)
        if i < len(dictKeyList)-1:
            for k in range(0, len(galaxyDict[dictKeyList[i+1]])):
                #is j paired with k from the next bucket?
                isPair(bucketList[j], galaxyDict[dictKeyList[i+1]][k], pairDict, clusterData, coordDict, thetaBoundDict, pairDict, i)

barGraph = plt.bar()
xList = list()
yList = list()
eList = list()
for i in range(0, len(dictKeyList)):
    if len(pairDict[dictKeyList[i]]) > 0:
        xList.append(dictKeyList[i])
        yList.append(len(pairDict[dictKeyList[i]])/len(galaxyDict[dictKeyList[i]]))
        eList.append(math.sqrt(1/len(pairDict[dictKeyList[i]] + 1/len(galaxyDict[dictKeyList[i]]))))
    else:
        xList.append(dictKeyList[i])
        yList.append(0.00)
        eList.append(0.00)

plt.errorbar(xList, yList, eList)
plt.savefig("outputTpair.png")

#all of the i's that apply from the clusterData array

#print(galaxyDict.__len__())


#for i in range(0, 2000):
#    coord1 = SkyCoord(clusterData[i]['RA'], clusterData[i]['DEC'], frame='icrs', unit="deg")
#    thetaBound = R_Bound(clusterData[i]['Z_LAMBDA']) / cosmo.angular_diameter_distance(clusterData[i]['Z_LAMBDA'])
#    print(i)
#    for j in range(i+1, 2000):
#        if math.fabs(clusterData[i]['Z_LAMBDA']-clusterData[j]['Z_LAMBDA']) <= 0.01:
#            coord2 = SkyCoord(clusterData[j]['RA'], clusterData[j]['DEC'], frame='icrs', unit="deg")
#            separationAngle = u.radian.to(coord1.separation(coord2), equivalencies=u.dimensionless_angles())
            #thetaAngle = u.degree.to(thetaBound, u.dimensionless_angles())
            #print(thetaBound)
            #print(separationAngle)
            #print()
#            if separationAngle < thetaBound:
#                print("cluster bound pair found! ID 1 = %d ID 2 = %d, Thetabound = %f\n" % (clusterData[i]['ID'], clusterData[j]['ID'], thetaBound))




