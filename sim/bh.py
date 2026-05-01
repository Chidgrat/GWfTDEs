# -*- coding: utf-8 -*-
from scipy.constants import G as G
from scipy.constants import c as c
import numpy as np

#solar mass
sm = 1.988416E30

class BlackHole:
    """
    Black hole class (fixed position)
    
    Mass in solar masses
    
    """
    def __init__(self, mass, length = 1):
        """
        Constructor for BlackHole

        Parameters
        ----------
        mass : float
            Mass of SMBH in solar masses.

        Returns
        -------
        None.

        """
        self._mass = mass * sm
        self._srad = (2 * G * self._mass)/c**2
        
        self._pos = np.asarray([0,0,0], dtype = np.float64)
        if length != 1:
            self._srad /= length
    
    @property
    def r(self):
        return self._srad
    @r.setter
    def r(self, val):
        self._srad = val
    
    @property
    def pos(self):
        return self._pos
    @pos.setter
    def pos(self, arr):
        self._pos += np.asarray(arr)
        
    @property
    def m(self):
        return self._mass
    @m.setter
    def m(self, val):
        self._mass = val