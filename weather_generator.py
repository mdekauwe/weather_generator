#!/usr/bin/env python
"""
AWAP weather generator functions
- testing before writing them in C for GDAY
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos, sin, exp, sqrt
import random

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
    plt.plot(hours, tday, "r-", label="Parton & Logan")
    plt.plot(hours, tday2, "k-", label="MAESPA")
    plt.legend(numpoints=1, loc="best")
    plt.ylabel("Air Temperature (deg C)")
    plt.xlabel("Hour of day")
    plt.show()

    rain = 10.0
    ppt = disaggregate_rainfall(rain)
    plt.plot(hours, ppt, "ro")
    plt.ylabel("PPT (mm)")
    plt.xlabel("Hour of day")
    plt.ylim(0, 3)
    plt.show()

    vpd09 = 1.4
    vpd09_next = 1.8
    vpd15 = 2.3
    vpd15_prev = 3.4
    vpd = estimate_diurnal_vpd(vpd09, vpd15, vpd09_next, vpd15_prev)

    plt.plot(hours, vpd, "ro")
    plt.ylabel("VPD (kPa)")
    plt.xlabel("Hour of day")
    plt.ylim(0, 3)
    plt.show()

def estimate_diurnal_vpd(vpd09, vpd15, vpd09_next, vpd15_prev):
    """
    Interpolate VPD between 9am and 3pm values to generate diurnal VPD
    following the method of Haverd et al.

    Reference:
    ---------
    * Haverd et al. (2013) Multiple observation types reduce uncertainty in
      Australia's terrestrial carbon and water cycles. Biogeosciences, 10,
      2011-2040.
    """
    # number of hours gap, i.e. 3pm to 9am the next day
    gap = 18.0
    
    vpd = np.zeros(48)
    for i in xrange(48):
        hour = i / 2.0
        if hour <= 9.0:
           vpd[i] = vpd15_prev + (vpd09 - vpd15_prev) * (9.0 + hour) / gap
        elif hour > 9.0 and hour <= 15.0:
           vpd[i] = vpd09 + (vpd15 - vpd09) * (hour - 9.0) / (15.0 - 9.0)
        elif hour > 15.0:
            vpd[i] =  vpd15 + (vpd09_next - vpd15) * (hour - 15.0) / gap

    return vpd

def disaggregate_rainfall(rain):
    """
    Assign daily PPT total to hours of the day, following MAESTRA, which follows
    algorithm from GRAECO (model of D. Loustau).

    Reference:
    * Loustau, D., F. Pluviaud, A. Bosc, A. Porte, P. Berbigier, M. Deque
      and V. Perarnaud. 2001. Impact of a regional 2 x CO2 climate scenario
      on the water balance, carbon balance and primary production
      of maritime pine in southwestern France. In Models for the Sustainable
      Management of Plantation Forests. Ed. M. Tome. European
      Cultivated Forest Inst., EFI Proc. No. 41D, Bordeaux, pp 45-58.
    """
    ppt = np.zeros(48)

    # All rain falls in one hour for light storms (<2 mm)
    if rain <= 2.0:
        hour_index = randint(0,48)
        ppt[hour_index] = rain

    # All rain falls in 24 hours for storms >46 mm
    elif rain > 46.0:
        for i in xrange(48):
            ppt[i] = rain / 48.0
    # All rain falls at 2mm/hour at a random time of the day
    else:
        #num_hrs_with_rain = min(int((rain / 2.0) * 48. / 24.), 48)
        num_hrs_with_rain = int(rain / 2.0)
        rate = rain / float(num_hrs_with_rain)
        # sample without replacement
        random_hours = random.sample(range(0, 48), num_hrs_with_rain)
        for i in xrange(num_hrs_with_rain):
            ppt[random_hours[i]] = rate

    return ppt

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

    TO DO:
    - Hours between 00:00 and sunrise should be modelled using the previous
      days information.

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
