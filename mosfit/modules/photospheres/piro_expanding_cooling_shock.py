"""Definitions for the `PiroExpandingCoolingShock` class."""
import numpy as np
from astrocats.catalog.source import SOURCE
from astropy import constants as c

from mosfit.constants import DAY_CGS, FOUR_PI, KM_CGS
from mosfit.modules.photospheres.photosphere import Photosphere


# Important: Only define one ``Module`` class per file.


class PiroExpandingCoolingShock(Photosphere):
    """Photosphere of Piro's expanding thermal material
    From equations 7-10 of https://doi.org/10.3847/1538-4357/abe2b1
    """

    _REFERENCES = [
        {SOURCE.BIBCODE: 'https://doi.org/10.3847/1538-4357/abe2b1'}
    ]

    STEF_CONST = (FOUR_PI * c.sigma_sb).cgs.value
    RAD_CONST = KM_CGS * DAY_CGS

    def process(self, **kwargs):
        """Process module."""
        kwargs = self.prepare_input(self.key('luminosities'), **kwargs)
        self._rest_t_explosion = kwargs[self.key('resttexplosion')]
        self._times = kwargs[self.key('rest_times')]
        
        self._Re = kwargs[self.key('Re')] * 6.957e10 #goes from units of solar radii to cm
        self._vt = kwargs[self.key('vt')] * 1.e5 #I'm not really sure, but I have to think this takes us into cm/sec from km.sec
        self._Me = kwargs[self.key('Me')] * 2.e33 #goes from solar mass units to grams
        self._kappa = kwargs[self.key('kappa')]
        self._luminosities = kwargs[self.key('luminosities')]
        #print("Luminosities: "+str(self._luminosities))
        n = 10.0
        delta = 1.1
        K = 0.119
        
        t_ph = ((3.*self._kappa * K * self._Me)/(2*(n-1.)*self._vt**2.))**0.5
        ts_day = [x - self._rest_t_explosion  for x in self._times]
        
        rphot = []
        for i in range(len(ts_day)):
            t_day = ts_day[i]
            t = t_day * 86400.
            if t_day < 0:
                r = 0.0
            elif t < t_ph:
                r = (t_ph/t)**(2./(n-1))*self._vt*t
            else:
                #r = (((delta-1)/(n-1))*((t/t_ph)**2.-1)+1)**(-1./(delta-1))*self._vt*t
                r = ((((delta-1)/(n-1))*(t**2/t_ph**2-1))+1)**(-1./(delta-1))*self._vt*t
            rphot.append(r)
        
        Tphot = []
        
        
        for li, lum in enumerate(self._luminosities):

            if lum == 0.0:
                temperature = 0.0
            else:
                temperature = (lum / (self.STEF_CONST * rphot[li]**2. * 4 * np.pi)) ** 0.25

            

            Tphot.append(temperature)
        #print("Luminosity: "+str(self._luminosities))
        #print("Temperatures: "+str(Tphot))
        #print("radius: "+str(rphot))

        return {self.key('radiusphot'): rphot,
                self.key('temperaturephot'): Tphot}
