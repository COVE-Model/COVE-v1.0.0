# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 10:01:51 2014

@author: mhurst
"""

#modules
import numpy as np
import matplotlib.pyplot as plt

#X and Y
try:
    X,Y = np.loadtxt("./../driver_files/XY.txt", unpack=True)
    f = open("./../driver_files/XY2.txt","w")
    for i in range(0,len(X)):
        f.write(str(X[i]) + " ")
    f.write("\n")
    for i in range(0,len(Y)):
        f.write(str(Y[i]) + " ")
    f.close()
        
except:
    print "Cannot find XY file!"
    
NoNodes = len(X)

#load Nodes
try:
    f = open("./../driver_files/Nodes.txt",'r')
    Lines = f.readlines()
    NoCells = len(Lines)
except:
    print "Cannot find Nodes file!"

plt.figure(1,figsize=(3,11))
plt.axis('Equal')
plt.plot(X,Y,'k.-',zorder=10)
[plt.text(X[i]-300,Y[i],str(i),fontsize="small") for i in range(0,NoNodes)]

for i in range(0,NoCells):
    Line = Lines[i].strip().split(" ")
    NoNodes = (len(Line)-2)/2
    print i, NoNodes
    X = np.array(Line[2:NoNodes*2+2:2],dtype='float64')
    Y = np.array(Line[3:NoNodes*2+2:2],dtype='float64')
    X = np.append(X,X[0])
    Y = np.append(Y,Y[0])
    if i == 89:
        plt.plot(X,Y,'r.-', lw=2,zorder=10)
        Sum = 0
        for j in range(0,len(X)):
            if (j<(len(X)-1)):
                Sum += (X[j]*Y[j+1] - X[j+1]*Y[j])
            else:
                Sum += X[j]*Y[0] - Y[j]*X[0]
        Area = np.abs(Sum/2.)
        
    else:
        plt.plot(X,Y,'-',color=[0.5,0.5,0.5], lw=2)
        plt.plot(X,Y,'k.')
    #plt.text(X[0],Y[0], str(i))


plt.show()