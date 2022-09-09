"""Definitions for the `Simplefallback` class."""
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

    _REFERENCES = [{SOURCE.BIBCODE: 'TBD'}]

    def process(self, **kwargs):
        """Process module."""
        self._times = kwargs[self.key('dense_times')]
        self._Re = kwargs[self.key('Re')] * 6.957e10
        self._ve = kwargs[self.key('ve')] * 1e9 
        self._Me = kwargs[self.key('Me')] * 2e33
        self._kappa = kwargs[self.key('kappa')]
        n = 10.0
        delta = 1.1
        K = 0.119

        self._rest_t_explosion = kwargs[self.key('resttexplosion')]


        td = np.sqrt(3. * self._kappa * K * self._Me / ((n - 1.0) * self._ve * C_CGS)) / 86400.0
        Eth = np.pi * (n - 1.0) * C_CGS * self._Re * self._ve**2 / (3. * (n - 5.0) * self._kappa)

        ts = [
            np.inf if self._rest_t_explosion > x else
            (x - self._rest_t_explosion) for x in self._times
        ]

        luminosities = [
            np.pi * (n - 1.0) * C_CGS * self._Re * self._ve**2 / (3.0 * (n - 5.0) * self._kappa) * (td / t)**(4. / (n - 2.0)) if t - td < 0 else
            Eth / td * np.exp(-0.5 * (t**2/td**2 - 1.0)) for t in ts
        ]
        luminosities = [0.0 if isnan(x) else x for x in luminosities]

        return {self.dense_key('luminosities'): luminosities}
