"""Definitions for the `ShockCoolingPiro` class."""
from math import isnan

import numpy as np
from astrocats.catalog.source import SOURCE

from mosfit.constants import DAY_CGS, C_CGS
from mosfit.modules.engines.engine import Engine

# Important: Only define one ``Module`` class per file.


class ShockCoolingPiro(Engine):
    """
    Shock cooling model for extended mass.
    """

    _REFERENCES = [{SOURCE.BIBCODE: 'https://doi.org/10.3847/1538-4357/abe2b1'}]

    def process(self, **kwargs):
        """Process module."""
        self._times = kwargs[self.key('dense_times')]
        self._Re = kwargs[self.key('Re')] * 6.957e10 #goes from units of solar radii to cm
        self._vt = kwargs[self.key('vt')] * 1e5 #this takes us into cm/sec from km/sec
        self._Me = kwargs[self.key('Me')] * 2e33 #goes from solar mass units to grams
        self._kappa = kwargs[self.key('kappa')]
        n = 10.0
        delta = 1.1
        K = 0.119
        
        #print("shock cooling, Re = "+str(self._Re))
        #print("shock cooling, vt = "+str(self._vt))
        #print("shock cooling, Me = "+str(self._Me))
        #print("shock cooling, kappa = "+str(self._kappa))
        
        #print("vt: "+str(self._vt))
        
        #v_t = (()/())**0.5*()**0.5

        self._rest_t_explosion = kwargs[self.key('resttexplosion')]


        td = np.sqrt(3. * self._kappa * K * self._Me / ((n - 1.0) * self._vt * C_CGS)) / 86400.0
        
        #print("td = "+str(td)+" days")
        
        Eth_td = np.pi * (n - 1.0) * C_CGS * self._Re * self._vt**2 / (3. * (n - 5.0) * self._kappa)
        
        #print("Eth_td: "+str(Eth_td))
        
        #print("self times: "+str(self._times))

        ts = [
            np.inf if self._rest_t_explosion > x else
            (x - self._rest_t_explosion) for x in self._times
        ]
        
        #print("ts: "+str(ts))

        luminosities = [
            np.pi * (n - 1.0) * C_CGS * self._Re * self._vt**2 / (3.0 * (n - 5.0) * self._kappa) * (td / t)**(4. / (n - 2.0)) if t - td < 0 else Eth_td * np.exp(-0.5 * (t**2/td**2 - 1.0)) for t in ts
        ]
        luminosities = [0.0 if isnan(x) else x for x in luminosities]
        
        #print("luminosities: "+str(luminosities))

        return {self.dense_key('luminosities'): luminosities}
