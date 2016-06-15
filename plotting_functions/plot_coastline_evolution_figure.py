# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 19:04:58 2013

@author: mhurst
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

##########################
# Customise figure style #
##########################
from matplotlib import rc
rc('font',size=8)
rc('ytick.major',pad=5)
rc('xtick.major',pad=5)
padding = 5

#Define wave climate params
Theta = 35
ThetaStd = 25

#open coastline evolution file
filename = "../driver_files/Spiral_6_1_1._0.1_" + str(Theta) + "_" + str(ThetaStd) + ".xy"

#read the header info and get the lines
f = open(filename,"r")
Lines = f.readlines()
NoLines = len(Lines)-1
Header = Lines[0]

#setup the figure space, axis and line to animate
fig = plt.figure(1, figsize=(5,4))
plt.subplots_adjust(0.2,0.15,0.9,0.85)
plt.axis('equal')
X = np.array(Lines[1].strip().split(" ")[1:],dtype='float64')
Y = np.array(Lines[2].strip().split(" ")[1:],dtype='float64')
plt.plot(X,Y,'k.',ms=1,label="Initial Coastline")
plt.xlabel('X-coordinate (m)')
plt.ylabel('Y-coordinate (m)')
plt.xlim(np.min(X)-200,np.max(X)+200)
plt.ylim(np.min(Y),np.max(Y))

TimeIntervals = [0.5,1,2,4,8,16,24.5]

Time = np.zeros(NoLines-1)
for i in range(0,NoLines-1,2):
    Time[i] = float(Lines[i+1].strip().split(" ")[0])

for TimeInt in TimeIntervals:
    Diff = np.abs(Time-TimeInt)
    Line = np.argmin(Diff)

    X = np.array(Lines[Line+1].strip().split(" ")[1:],dtype='float64')
    Y = np.array(Lines[Line+2].strip().split(" ")[1:],dtype='float64')

    plt.plot(X,Y,'k--',lw=0.5)

xmin = np.min(X)-500.
xmax = np.max(X)+500
ymin = np.min(Y)-500
ymax = np.max(Y)+200

SeaX = np.append(X,[xmax,xmax,X[0]])
SeaY = np.append(Y,[ymin,ymax,Y[0]])
BeachX = np.append(X,[xmin,xmin,X[0]])
BeachY = np.append(Y,[ymin,ymax,Y[0]])
    
plt.fill(SeaX, SeaY, color=[0.7,0.9,1.0])
plt.fill(BeachX, BeachY, color=[1.0,0.9,0.6])
plt.plot(X,Y,'k--',lw=0.5,label="1/2, 1, 2, 4, 8, 16 Year Coastlines")
plt.plot(X,Y,'k-',lw=1,label="32 Year Coastline")
plt.plot(X,Y,'ko',ms=3)
plt.plot(X[0:2],Y[0:2],'k-',lw=2)
plt.plot(X[-2:],Y[-2:],'k-',lw=2)

#Display legend
plt.rcParams.update({'legend.labelspacing':0.1}) 
plt.rcParams.update({'legend.numpoints':1}) 
plt.rcParams.update({'legend.frameon':False}) 
plt.rcParams.update({'legend.handlelength':1.5}) 
plt.legend(loc=3)
leg = plt.gca().get_legend()
#sett fontsize to small
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=8) 

#bin all wave data and plot!
Waves = np.random.normal(Theta,ThetaStd,10000)
widths = 15*np.pi/180.0
ax = plt.axes([0.65,0.55,0.25,0.25],polar=True)
hist = np.histogram(Waves*np.pi/180.0,bins=np.arange(0,361*np.pi/180.0,widths))
plt.bar(hist[1][:-1],hist[0],widths,color='white',edgecolor='k')
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_axis_bgcolor('none')
plt.axis('off')
plt.yticks([])
plt.xticks([])
plt.show()
        