# -*- coding: utf-8 -*-

import numpy as np

#solar mass and solar radius
sm = 1.988416E30
sr = 6.96340E8

class Body:
    """
    Class for dynamic objects beyond the black hole
    """
    def __init__(self, mass, radius, name = "Unnamed", length = 1):
        self._pos = np.asarray([0,0,0], dtype = np.float64)
        self._vel = np.asarray([0,0,0], dtype = np.float64)
        
        self._mass = mass
        if mass < 10000:
            self._mass *= sm
            
        self._radius = radius
        if radius < 10000:
            self._radius *= sr
        if length != 1:
            self._radius /= length
        
        self.name = name
        self.orbEn = 0
        
        
    def __str__(self):
        return f"{self.name.title()} is a body at {self.pos} with velocity {self.vel}, with mass {self.m} and radius {self.r}"
        
    def rCube(self, obj):
        rx, ry, rz = self.pos-obj.pos
        rCube = (rx**2 + ry**2 + rz**2)**1.5
        
        return rCube
        
    @property
    def x(self):
        return self._pos[0]
    @x.setter
    def x(self, val):
        self._pos[0] += val
    
    @property
    def y(self):
        return self._pos[1]
    @y.setter
    def y(self, val):
        self._pos[1] += val
        
    @property
    def z(self):
        return self._pos[2]
    @z.setter
    def z(self, val):
        self._pos[2] += val
        
    @property
    def pos(self):
        return self._pos
    @pos.setter
    def pos(self, arr):
        if len(arr) == 2:
            self.x = arr[0]
            self.y = arr[1]
        else:
            self._pos += np.asarray(arr)
        
    @property
    def vel(self):
        return self._vel
    @vel.setter
    def vel(self, arr):
        if len(arr) == 2:
            self._vel[0] += arr[0]
            self._vel[1] += arr[1]
        else:
            self._vel += np.asarray(arr)
            
    @property
    def velAbs(self):
        #return (self._vel[0]**2 + self._vel[1]**2 + self._vel[2]**2)**0.5
        return np.sqrt(self._vel[0]**2 + self._vel[1]**2)
        
    @property
    def m(self):
        return self._mass
    @m.setter
    def m(self, val):
        self._mass = val
    
    @property
    def r(self):
        return self._radius
    @r.setter
    def r(self, val):
        self._radius = val
