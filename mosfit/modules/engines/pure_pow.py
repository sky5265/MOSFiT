"""Definitions for the `ExpPow` class."""
from math import isnan

import numpy as np

from mosfit.modules.engines.engine import Engine


# Important: Only define one ``Module`` class per file.


class pure_exp(Engine):
    """A simple analytical engine."""

    def process(self, **kwargs):
        """Process module."""
        self._times = kwargs[self.key('dense_times')]
        self._t_peak = kwargs[self.key('tpeak')]
        self._lum_scale = kwargs[self.key('lumscale')]
        self._rest_t_explosion = kwargs[self.key('resttexplosion')]
        self._alpha = kwargs[self.key('alpha')]

        ts = [
            np.inf
            if self._rest_t_explosion > x else (x - self._rest_t_explosion)
            for x in self._times
        ]

        luminosities = [
            self._lum_scale *  (t / self._t_peak)**(-1.0*self._alpha) for t in ts
        ]
        luminosities = [0.0 if isnan(x) else x for x in luminosities]

        return {self.dense_key('luminosities'): luminosities}
