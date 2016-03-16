#!/usr/bin/env python
"""
AWAP weather generator functions
- testing before writing them in C for GDAY
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos, sin, exp, sqrt


__author__  = "Martin De Kauwe"
__version__ = "1.0 (16.03.2016)"
__email__   = "mdekauwe@gmail.com"


def main():

    tmin = 2.0
    tmax = 24.0
    day_length = 8.0

    tday = estimate_diurnal_temp(tmin, tmax, day_length)
    tday2 = maestra_diurnal_func(tmin, tmax, day_length)
    #par = estimate_diurnal_par(1500.0, day_length)
    hours = np.arange(48) / 2.0
    plt.plot(hours, tday, "r-")
    plt.plot(hours, tday2, "k-")
    plt.show()

def maestra_diurnal_func(tmin, tmax, day_length):
    """ Not sure where this function original comes from... """
    tav = (tmax + tmin) / 2.0
    tampl = (tmax - tmin) / 2.0

    tday = np.zeros(48)
    for i in xrange(48):
        hrtime = i - 0.5
        time = i + day_length * 0.5 - 48.0 / 2.0
        if time < 0.0 or time > day_length:
            if time < 0.0:
                hrtime += 48
            tday[i] = tav - (tav - tmin) * (hrtime - day_length * 0.5 - \
                                            (48.0 / 2.0)) / (48.0 - day_length)
        else:
            tday[i] = tav - tampl * cos(1.5 * pi * time / day_length)

    return (tday)

def estimate_diurnal_temp(tmin, tmax, day_length):
    """
    Calculate diurnal temperature following Parton and Logan
    the day is divided into two segments and using a truncated sine wave
    in the daylight and an exponential decrease in temperature
    at night.

    References:
    ----------
    * Parton and Logan (1981) A model for dirunal variation in soil and
       air temperature. Agricultural Meteorology, 23, 205--216.
    * Kimball and Bellamy (1986) Energy in Agriculture, 5, 185-197.
    """
    # 1.5 m air temperature values from Parton and Logan, table 1
    a = 1.86
    b = 2.2     # nighttime coeffcient
    c = -0.17   # lag of the min temp from the time of runrise


    night_length = 24 - day_length

    sunrise = 12.0 - day_length / 2.0 + c
    sunset = 12.0 + day_length / 2.0

    # temperature at sunset
    m = sunset - sunrise + c
    tset = (tmax - tmin) * sin(pi * m / (day_length + 2.0 * a)) + tmin

    tday = np.zeros(48)
    for i in xrange(48):
        hour = i / 2.0

        # hour - time of the minimum temperature (accounting for lag time)
        m = hour - sunrise + c
        if hour >= sunrise and hour <= sunset:
            tday[i] = tmin + (tmax - tmin) * \
                        sin((pi * m) / (day_length + 2.0 * a))
        else:
            if hour > sunset:
                n = hour - sunset
            elif hour < sunrise:
                n = (24.0 + hour) - sunset


            d = (tset - tmin) / (exp(b) - 1.0)

            # includes missing displacement to allow T to reach Tmin, this
            # removes a discontinuity in the original Parton and Logan eqn.
            # See Kimball and Bellamy (1986) Energy in Agriculture, 5, 185-197
            tday[i] = (tmin -d) + (tset - tmin - d) * \
                        exp(-b * n / (night_length + c))


    return (tday)

if __name__ == "__main__":

    main()
