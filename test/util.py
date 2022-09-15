#!/usr/bin/env python3

import numpy as np

import pyspeckit
from astropy import units as u
from pyspeckit.spectrum.models import ammonia, ammonia_constants


# Update pyspeckit constants to match implementation
ammonia.aval_dict["oneone"] = 1.67524303e-7
ammonia.aval_dict["twotwo"] = 2.24162441e-7
ammonia.TCMB = 2.72548
ammonia_constants.Brot = 298192.92e6
ammonia_constants.Crot = 186695.86e6


class PskAmmoniaSpec:
    def __init__(self, vchan=0.158):
        self.vchan = vchan  # km/s
        freqs = ammonia_constants.freq_dict
        Axis = pyspeckit.spectrum.units.SpectroscopicAxis
        vaxis = np.arange(30, -30, -vchan) * u.km / u.s
        xa11 = Axis(vaxis, velocity_convention="radio",
                refX=freqs["oneone"]).as_unit("Hz")
        xa22 = Axis(vaxis, velocity_convention="radio",
                refX=freqs["twotwo"]).as_unit("Hz")
        self.xa11 = xa11
        self.xa22 = xa22
        self.varr = vaxis.value

    def get_spec11(self, xoff_v=0.0, trot=15.0, tex=6.0, ntot=15.0, width=0.6,
            fortho=0.0):
        return ammonia.ammonia(
                self.xa11, xoff_v=xoff_v, trot=trot, tex=tex, ntot=ntot, width=width,
                fortho=fortho, line_names=["oneone"],
        )

    def get_spec22(self, xoff_v=0.0, trot=15.0, tex=6.0, ntot=15.0, width=0.6,
            fortho=0.0):
        return ammonia.ammonia(
                self.xa22, xoff_v=xoff_v, trot=trot, tex=tex, ntot=ntot, width=width,
                fortho=fortho, line_names=["twotwo"],
        )

    def get_spectra(self, xoff_v=0.0, trot=15.0, tex=6.0, ntot=15.0, width=0.6,
            fortho=0.0):
        return self.get_spec11(), self.get_spec22()


