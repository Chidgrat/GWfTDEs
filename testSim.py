# -*- coding: utf-8 -*-
"""
Created on Sun Nov  2 20:07:42 2025

@author: Alex
"""

import sim
from scipy.constants import G
from numpy import sqrt, arccos, cos, sin, arange
import matplotlib.pyplot as plt

vals = arange(0.05, 1.51, 0.05)

for i in vals:
    bh = sim.BlackHole(4.297E6)
    #bh = sim.BlackHole(0.1)
    #bh = sim.BlackHole(0.1)
    D = i # Radius of periastron/tidal radius
    print(f"Depth = {D:.2f}")
    
    sun1 = sim.Body(3.2, 2.5, "Primary Star")
    sun2 = sim.Body(0.8, 0.75, "Secondary Star")
    sunSum = sun1.m + sun2.m
    
    a = 7 * 1.496e11 #semi-major 
    theta = 0
    
    xp = (sun2.m/sunSum) * a # Primary star position in binary
    xs = -(sun1.m/sunSum) * a # Secondary star position in binary
    
    omega = sqrt((2 * G * sunSum)/a**3) #angular velocity of binary
    
    vp = xp * omega #velocity of primary in binary
    vs = -xs * omega #velocity of secondary in binary
    
    rTidal = ((2 * bh.m/sunSum)**(1/3)) * a #tidal radius
    r0 = 10 * rTidal #initial radius
    
    rp = rTidal * D #radius periastron
    print(f"{rTidal:e}")
    
    f0 = -arccos(-1 + (D/5)) #initial true anomaly
    
    xcm = r0 * cos(f0) # binary CoM in x
    ycm = r0 * sin(f0) # binary CoM in y
    
    fDot = sqrt(G * bh.m/(rp**3)) * (sqrt(2)/4) * ((1 + cos(f0))**2) #time derivative of true anomaly
    rDot = -2 * rp * (1/((1 + cos(f0))**2)) * (-sin(f0)) * fDot #time derivative of radius
    
    vxcm = rDot * cos(f0) - r0 * fDot * sin(f0) # x 
    vycm = rDot * sin(f0) + r0 * fDot * cos(f0)
    
    sun1.pos = [xcm, xp + ycm]
    sun2.pos = [xcm, xs + ycm]
    
    sun1.vel = [vp + vxcm, vycm]
    sun2.vel = [vs + vxcm, vycm]
    
    timeStep = sqrt((a**3)/(G * sunSum))
    
    program = sim.Sim(bh, [sun1, sun2], timeStep/1000, timeStep*200, rTidal, a = a)
    out = program.propagate(False, False, True, D = D)
    
    val = float(f"{i:.2f}")
    out.savefig(f"C:/Users/doguy/Main/Uni/Liverpool/Y4/Project/Images/270/fig_{str(val).replace('.','')}")
    #out.clear()
    
    #print(out[0][0][1][0])