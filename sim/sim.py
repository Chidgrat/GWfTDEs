# -*- coding: utf-8 -*-
"""
Created on Sun Oct 26 17:37:05 2025

@author: Alex
"""

import numpy as np
from scipy import constants as sp
import matplotlib.pyplot as plt
import time
from .body import Body

from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'qt')

class Sim:
    def __init__(self, bh, bodies, dt, tfinal, rTidal, length = 1, a = 1):
        """
        Constructor function for the Sim class

        Parameters
        ----------
        bh : BlackHole object
            Reference to the "black hole" object fixed at simulation centre
        bodies : List of object Body
            List of bodies in the sim (excluding the black hole)
        length : float
            Length element of the simulation
        dt : float
            timestep
        tfinal : float
            length of time (in seconds) for sim to run for

        Returns
        -------
        None.

        """
        self._dt = dt
        self._tfinal = tfinal
        self._t = 0
        
        self.bh = bh
        self.bodies = bodies

        self.__st = time.time()

        self.d = 0
        
        self._rTidal = rTidal
        
        self._sunSum = 0
        self._suns = []
        for obj in self.bodies:
            if type(obj) == Body:
                self._sunSum += obj.m
                self._suns.append(obj)
        self.a = a
        
        self.trip = False
                
    
    def propagate(self, plot3d = True, plotEn = False, pltout = False, D = 1.00):
        xarr = []
        yarr = []
        zarr = []
        s1en = []
        s2en = []
        tArr = []
        
        while self._t <= self._tfinal:
            y = self._dostep()
            for i in range(len(y)):
                y[i][0].pos = y[i][1][0:3]
                y[i][0].vel = y[i][1][3:6]
                exists = False
                for j in xarr:
                    try:
                        if j[0] == y[i][0]:
                            exists = True
                    except IndexError:
                        pass
                if exists:
                    for j in xarr:
                        if j[0] == y[i][0]:
                            j[1].append(y[i][0].x)
                            yarr[xarr.index(j)][1].append(y[i][0].y)
                            zarr[xarr.index(j)][1].append(y[i][0].z)
                else:
                    xarr.append([y[i][0], [y[i][0].x]])
                    yarr.append([y[i][0], [y[i][0].y]])
                    zarr.append([y[i][0], [y[i][0].z]])
            
            flop = False
            for sun in self._suns:
                for obj in self._suns:
                    if obj != sun:
                        other = obj
                sun.orbEn = 0.5 * sun.m * (sun.velAbs ** 2) \
                    - ((sp.G * (self.bh.m * sun.m))/(sun.rCube(self.bh))) \
                        - ((sp.G * self._sunSum)/(sun.rCube(other)))
                sun.orbEn /= ((sp.G * self._sunSum)/self.a)
                #sun.orbEn -= 0.75e42
                if flop:
                    s2en.append(sun.orbEn)
                else:
                    s1en.append(sun.orbEn)
                flop = not flop
            tArr.append(self._t - 1.10e9)
            #print(f"{s1en}, {s2en}, {tArr}")
                       
            self._t += self._dt
        
        if not plotEn:
            if plot3d:
                ax = plt.figure(figsize=(8,8)).add_subplot(projection='3d')
                ax.view_init(elev=90, azim=-90, roll=0)
                for i in range(len(xarr)):
                    ax.plot(xarr[i][1], yarr[i][1], zarr[i][1], label = xarr[i][0].name)
                
                u = np.linspace(0, 2 * np.pi, 100)
                v = np.linspace(0, np.pi, 50)
                
                u, v = np.meshgrid(u, v)
                
                r = self.bh.r
                
                bhx = r * np.sin(v) * np.cos(u)
                bhy = r * np.sin(v) * np.sin(u)
                bhz = r * np.cos(v)
            
                sol = ax.plot_surface(bhx, bhy, bhz, color = "black", label = "SMBH")
                sol._edgecolors2d = sol._edgecolor3d
                sol._facecolors2d = sol._facecolor3d
                
                ax.legend()
                
                ax.set_xlim(-abs(xarr[0][0].x), abs(xarr[0][0].x))
                ax.set_ylim(abs(-xarr[0][0].x), abs(xarr[0][0].x))
            
            else:
                ax = plt.gca()
                ax.clear()
                for i in range(len(xarr)):
                    ax.plot(xarr[i][1], yarr[i][1], label = xarr[i][0].name)
                
                ax.scatter(xarr[0][1][0], yarr[0][1][0], label = "Start", linewidths = 2, c = "m", marker = 'x')
                
                bhCircle = plt.Circle((0,0), self.bh.r, color = "black")
                tdCircle = plt.Circle((0,0), self._rTidal, color = "black", fill = False, label = "Tidal Radius")
                
                ax.add_patch(bhCircle)
                ax.add_patch(tdCircle)
                
                ax.text(0.75, 0.1, f"Depth = {D:.2f}", transform=ax.transAxes)
                ax.grid()
                
                ax.legend()
                
                #ax.set_xlim(-abs(xarr[0][0].x), abs(xarr[0][0].x))
                #ax.set_ylim(-abs(-xarr[0][0].x), abs(xarr[0][0].x))
                ax.set_xlim(-2.5e15, 2.5e15)
                ax.set_ylim(-2.5e15, 2.5e15)
                
                ax.set_xlabel("Distance from Black Hole (m)")
                ax.set_ylabel("Distance from Black Hole (m)")
        else:
            plt.figure(figsize=(8,8))
            plt.plot(tArr, s1en, label = "Primary", lw = 3)
            plt.plot(tArr, s2en, label = "Secondary", lw = 3)
            
            #plt.hlines(0, -1.1e9, 9302249075 - 1.09e9, color = '0')
            #plt.vlines(0, -0.5e43, 3e43, color = '0')
            
            plt.xlabel("Time (s)")
            plt.ylabel("Relative Orbital Energy (J)")
            
            #plt.xlim(-1.1e9, 9302249075 - 1.1e9)
            #plt.ylim(-0.25e43, 3e43)
            
            plt.legend()
            
        et = time.time()
        ct = et - self.__st
        print(f"Sim execution time = {time.strftime('%H:%M:%S', time.gmtime(ct))}")
        
        if not pltout:
            plt.show()
        else:
            return plt
        
        return [xarr, yarr, zarr]
    
    def _dostep(self):
        y = []
        for body in self.bodies:
            y.append([body, np.concatenate((body.pos, body.vel))])
        y = np.asarray(y, dtype = "object")
        
        k1 = self._differential(y)
        v1 = y[:,1] + k1[:,1] * self._dt/2
        u1 = []
        for i in range(len(k1)):
            u1.append([k1[i][0], v1[i]])
        u1 = np.asarray(u1, dtype = "object")
        
        k2 = self._differential(u1)
        v2 = y[:,1] + k2[:,1] * self._dt/2
        u2 = []
        for i in range(len(k2)):
            u2.append([k2[i][0], v2[i]])
        u2 = np.asarray(u2, dtype = "object")
        
        k3 = self._differential(u2)
        v3 = y[:,1] + k3[:,1] * self._dt
        u3 = []
        for i in range(len(k3)):
            u3.append([k3[i][0], v3[i]])
        u3 = np.asarray(u3, dtype = "object")
        
        k4 = self._differential(u3)
        
        step = (1/6) * (k1[:,1] + 2*k2[:,1] + 2*k3[:,1] + k4[:,1]) * self._dt
        
        out = []
        for i in range(len(step)):
            out.append([k4[i][0], step[i]])
        out = np.asarray(out, dtype = "object")
        return out
    
    def _differential(self, y):
        f = np.empty_like(y)
        for i in range(len(y)):
            f[i][0] = y[i][0]
            f[i][1] = np.concatenate((y[i][1][3:6], np.zeros(3)))

        for h in range(len(y)):
            body = y[h][0]
            if (body.rCube(self.bh)**(1/3)) > self._rTidal:
                for i in range(len(y)):
                    if y[i][0] is body:
                        pass
                    else:
                        grav = -sp.G * y[i][0].m * (body.pos - y[i][0].pos)/(body.rCube(y[i][0]))
                        f[h][1][3:6] += grav
            else:
#                if not self.trip:
#                    print(self._t)
#                    self.trip = True
                pass
            grav = -sp.G * self.bh.m * (body.pos - self.bh.pos)/(body.rCube(self.bh))
            f[h][1][3:6] += grav
        return f
                            
