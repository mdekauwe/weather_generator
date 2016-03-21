#!/usr/bin/env python
"""
AWAP weather generator functions
- testing before writing them in C for GDAY
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from math import pi, cos, sin, exp, sqrt, acos, asin
import random

__author__  = "Martin De Kauwe"
__version__ = "1.0 (16.03.2016)"
__email__   = "mdekauwe@gmail.com"


def main():

    SW_2_PAR = 2.3
    PAR_2_SW = 1.0 / SW_2_PAR
    J_TO_MJ = 1E-6
    MJ_TO_J = 1E6
    SEC_2_DAY = 86400.0
    DAY_2_SEC = 1.0 / SEC_2_DAY
    UMOLPERJ = 4.57     # Conversion from J to umol quanta
    lat = 51.5
    lon = -0.13
    hours = np.arange(48) / 2.0
    tmin = 2.0
    tmax = 24.0
    doy = 180.0
    sw_rad_day = 500.0 # mj m-2 d-1

    # MJ m-2 d-1 -> J m-2 s-1 = W m-2 -> MJ m-2 d-1
    # all the other units conv cancel.
    par_day = sw_rad_day * SW_2_PAR / UMOLPERJ

    day_length = calc_day_length(doy, 365, lat)
    print day_length
    cos_zenith = calculate_solar_geometry(doy, lat, lon)
    diffuse_frac = spitters(doy, par_day, cos_zenith)
    par = estimate_dirunal_par(par_day, cos_zenith, diffuse_frac)

    # check disaggregation
    print sw_rad_day, np.sum(par * PAR_2_SW * J_TO_MJ * SEC_2_DAY), "MJ m-2 d-1"
    print sw_rad_day * MJ_TO_J * DAY_2_SEC * SW_2_PAR, np.sum(par), "umol m-2 s-1"


    plt.plot(hours, par, "b-")
    plt.ylabel("par (umol m$^{-2}$ s$^{-1}$)")
    plt.xlabel("Hour of day")
    #plt.legend(numpoints=1, loc="best")
    plt.show()

    tday = estimate_diurnal_temp(tmin, tmax, day_length)
    tday2 = maestra_diurnal_func(tmin, tmax, day_length)

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

def estimate_dirunal_par(par_day, cos_zenith, diffuse_frac):
    """
    Calculate daily course of incident PAR from daily totals using routinue
    from MAESTRA

    Arguments:
    ----------
    par_day : double
        daily total photosynthetically active radiation (MJ m-2 d-1)
    cos_zenith : float
        cosine of zenith angle (radians)
    diffuse_frac : float
        fraction of par which is diffuse [0-1]
    """

    # Transmissivity of atmosphere
    tau = 0.76
    direct_frac = 1.0 - diffuse_frac
    MJ_TO_J = 1E6
    UMOLPERJ = 4.57     # Conversion from J to umol quanta
    SEC_2_DAY = 86400.0
    DAY_2_SEC = 1.0 / SEC_2_DAY

    cos_bm = np.zeros(48)
    cos_df = np.zeros(48)

    sum_bm = 0.0
    sum_df = 0.0
    for i in xrange(1,48+1):

        hrtime = float(i) - 0.5

        if cos_zenith[i-1] > 0.0:
            zenith = acos(cos_zenith[i-1])

            # set FBM = 0.0 for ZEN > 80 degrees
            if zenith < 80.0 * pi / 180.0:
                cos_bm[i-1] = cos_zenith[i-1] * tau**(1.0 / cos_zenith[i-1])
            else:
                cos_bm[i-1] = 0.0

            cos_df[i-1] = cos_zenith[i-1]
            sum_bm += cos_bm[i-1]
            sum_df += cos_df[i-1]

    par = np.zeros(48)
    for i in xrange(1,48+1):

        if sum_bm > 0.0:
            rdbm = par_day * direct_frac * cos_bm[i-1] / sum_bm
        else:
            rdbm = 0.0

        if sum_df > 0.0:
            rddf = par_day * diffuse_frac * cos_df[i-1] / sum_df
        else:
            rddf = 0.0

        # MJ m-2 d-1 -> J m-2 s-1 -> umol m-2 s-1
        par[i-1] = (rddf + rdbm) * MJ_TO_J * DAY_2_SEC * UMOLPERJ

    return par

def estimate_diurnal_vpd(vpd09, vpd15, vpd09_next, vpd15_prev):
    """
    Interpolate VPD between 9am and 3pm values to generate diurnal VPD
    following the method of Haverd et al. This seems reasonable, vapour pressure
    plotted aginst time of day often does not reveal consistent patterns, with
    small fluctuations (see Kimball and Bellamy, 1986).


    Reference:
    ---------
    * Haverd et al. (2013) Multiple observation types reduce uncertainty in
      Australia's terrestrial carbon and water cycles. Biogeosciences, 10,
      2011-2040.
    """
    # number of hours gap, i.e. 3pm to 9am the next day
    gap = 18.0

    vpd = np.zeros(48)
    for i in xrange(1, 48+1):
        hour = i / 2.0
        if hour <= 9.0:
           vpd[i-1] = vpd15_prev + (vpd09 - vpd15_prev) * (9.0 + hour) / gap
        elif hour > 9.0 and hour <= 15.0:
           vpd[i-1] = vpd09 + (vpd15 - vpd09) * (hour - 9.0) / (15.0 - 9.0)
        elif hour > 15.0:
            vpd[i-1] =  vpd15 + (vpd09_next - vpd15) * (hour - 15.0) / gap

    return vpd

def disaggregate_rainfall(rain_day):
    """
    Assign daily PPT total to hours of the day, following MAESTRA, which follows
    algorithm from GRAECO (model of D. Loustau).

    Arguments:
    ----------
    rain_day : double
        daily rainfall total (mm)

    Returns:
    ----------
    rain : array
        dirunal time course of rainfall (mm)

    Reference:
    ----------
    * Loustau, D., F. Pluviaud, A. Bosc, A. Porte, P. Berbigier, M. Deque
      and V. Perarnaud. 2001. Impact of a regional 2 x CO2 climate scenario
      on the water balance, carbon balance and primary production
      of maritime pine in southwestern France. In Models for the Sustainable
      Management of Plantation Forests. Ed. M. Tome. European
      Cultivated Forest Inst., EFI Proc. No. 41D, Bordeaux, pp 45-58.
    """
    rain = np.zeros(48)

    # All rain falls in one hour for light storms (<2 mm)
    if rain_day <= 2.0:
        hour_index = np.random.randint(low=0, high=47)
        rain[hour_index] = rain_day

    # All rain falls in 24 hours for storms >46 mm
    elif rain_day > 46.0:
        for i in xrange(48):
            rain[i] = rain_day / 48.0

    # Aim if for all rain to fall at 2mm/hour at a random time of the day.
    # If we generate the same random number, then we increase rainfall
    # for this hour
    else:
        #num_hrs_with_rain = min(int((rain / 2.0) * 48. / 24.), 48)
        num_hrs_with_rain = int(rain_day / 2.0)
        rate = rain_day / float(num_hrs_with_rain)
        # sample without replacement
        #random_hours = random.sample(range(0, 48), num_hrs_with_rain)
        random_hours = np.random.randint(low=0, high=47, size=num_hrs_with_rain)
        print random_hours
        for i in xrange(num_hrs_with_rain):
            rain[random_hours[i]] += rate

    return rain

def maestra_diurnal_func(tmin, tmax, day_length):
    """ Not sure where this function original comes from... """
    tav = (tmax + tmin) / 2.0
    tampl = (tmax - tmin) / 2.0

    dayl = day_length * 2.

    tday = np.zeros(48)
    for i in xrange(1, 48+1):
        hrtime = i - 0.5
        time = i + dayl * 0.5 - 48.0 / 2.0
        if time < 0.0 or time > dayl:
            if time < 0.0:
                hrtime += 48

            arg1 = tav
            arg2 = (tav - tmin) * (hrtime - dayl * 0.5 - (48.0 / 2.0))
            arg3 = 48.0 - dayl

            tday[i-1] = arg1 - arg2 / arg3
        else:
            tday[i-1] = tav - tampl * cos(1.5 * pi * time / dayl)

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

    Arguments:
    ----------
    tmin : double
        day minimum temp (deg C)
    tmax : float
        day maximum temp (deg C)
    day_length : float
        length of day (hours)


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
    for i in xrange(1, 48+1):
        hour = i / 2.0

        # hour - time of the minimum temperature (accounting for lag time)
        m = hour - sunrise + c
        if hour >= sunrise and hour <= sunset:
            tday[i-1] = tmin + (tmax - tmin) * \
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
            tday[i-1] = (tmin -d) + (tset - tmin - d) * \
                        exp(-b * n / (night_length + c))


    return (tday)

def calculate_solar_geometry(doy, latitude, longitude):

    """
    The solar zenith angle is the angle between the zenith and the centre
    of the sun's disc. The solar elevation angle is the altitude of the
    sun, the angle between the horizon and the centre of the sun's disc.
    Since these two angles are complementary, the cosine of either one of
    them equals the sine of the other, i.e. cos theta = sin beta. I will
    use cos_zen throughout code for simplicity.

    Arguments:
    ----------
    doy : double
        day of year
    latitude : float
        latitude (degrees)
    longitude : float
        longitude (degrees)

    References:
    -----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    """

    cos_zenith = np.zeros(48)
    for i in xrange(1, 48+1):
        # need to convert 30 min data, 0-47 to 0-23.5
        hod = i / 2.0

        gamma = day_angle(doy)
        dec = calculate_solar_declination(doy, gamma)
        et = calculate_eqn_of_time(gamma)
        t0 = calculate_solar_noon(et, longitude)
        h = calculate_hour_angle(hod, t0)
        rlat = latitude * pi / 180.0

        # A13 - De Pury & Farquhar
        sin_beta = sin(rlat) * sin(dec) + cos(rlat) * cos(dec) * cos(h)
        cos_zenith[i-1] = sin_beta; # The same thing, going to use throughout
        if cos_zenith[i-1] > 1.0:
            cos_zenith[i-1] = 1.0
        elif cos_zenith[i-1] < 0.0:
            cos_zenith[i-1] = 0.0

        zenith = 180.0 / pi * acos(cos_zenith[i-1])
        elevation = 90.0 - zenith

    return cos_zenith

def calculate_solar_noon(et, longitude):
    """
    Calculation solar noon - De Pury & Farquhar, '97: eqn A16

    Arguments:
    ----------
    et : float
        equation of time (radians)
    longitude : float
        longitude (degrees)

    Returns:
    ---------
    t0 - solar noon (hours).

    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    """
    # all international standard meridians are multiples of 15deg east/west of
    # greenwich
    Ls = round_to_value(longitude, 15.)
    t0 = 12.0 + (4.0 * (Ls - longitude) - et) / 60.0

    return (t0)

def round_to_value(number, roundto):
    return round(number / roundto) * roundto

def calculate_solar_declination(doy, gamma):
    """
    Solar Declination Angle is a function of day of year and is indepenent
    of location, varying between 23deg45' to -23deg45'

    Arguments:
    ----------
    doy : int
        day of year, 1=jan 1
    gamma : double
        fractional year (radians)

    Returns:
    --------
    dec: float
        Solar Declination Angle [radians]

    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    * Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.
    * J. W. Spencer (1971). Fourier series representation of the position of
      the sun.
    """

    # declination (radians) */
    #decl = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - \
    #       0.006758 * cos(2.0 * gamma) + 0.000907 * sin(2.0 * gamma) -\
    #       0.002697 * cos(3.0 * gamma) + 0.00148 * sin(3.0 * gamma);*/

    # (radians) A14 - De Pury & Farquhar
    decl = -23.4 * (pi / 180.) * cos(2.0 * pi * (doy + 10) / 365)

    return (decl)

def calculate_hour_angle(t, t0):
    """
    Calculation solar noon - De Pury & Farquhar, '97: eqn A15
    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    Returns:
    ---------
    h - hour angle (radians).
    """
    return (pi * (t - t0) / 12.0)


def day_angle(doy):
    """
    Calculation of day angle - De Pury & Farquhar, '97: eqn A18

    Arguments:
    ----------
    doy : int
        day of year

    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    * J. W. Spencer (1971). Fourier series representation of the position of
      the sun.

    Returns:
    ---------
    gamma - day angle in radians.
    """

    return (2.0 * pi * (float(doy) - 1.0) / 365.0)

def calculate_eqn_of_time(gamma):
    """
    Equation of time - correction for the difference btw solar time
    and the clock time.

    Arguments:
    ----------
    gamma : double
        fractional year (radians)

    References:
    -----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    * Campbell, G. S. and Norman, J. M. (1998) Introduction to environmental
      biophysics. Pg 169.
    * J. W. Spencer (1971). Fourier series representation of the position of
      the sun.
    * Hughes, David W.; Yallop, B. D.; Hohenkerk, C. Y. (1989),
      "The Equation of Time", Monthly Notices of the Royal Astronomical
      Society 238: 1529-1535
    """

    # minutes - de Pury and Farquhar, 1997 - A17
    et = (0.017 + 0.4281 * cos(gamma) - 7.351 * sin(gamma) - 3.349 *
          cos(2.0 * gamma) - 9.731  * sin(gamma))

    return (et)

def calc_extra_terrestrial_rad(doy, cos_zen):
    """
    Solar radiation incident outside the earth's atmosphere, e.g.
    extra-terrestrial radiation. The value varies a little with the earths
    orbit.
    Using formula from Spitters not Leuning!

    Arguments:
    ----------
    doy : double
        day of year
    cos_zen : double
        cosine of zenith angle (radians)

    Returns:
    --------
    So : float
        solar radiation normal to the sun's bean outside the Earth's atmosphere
        (J m-2 s-1)

    Reference:
    ----------
    * Spitters et al. (1986) AFM, 38, 217-229, equation 1.
    """

    # Solar constant (J m-2 s-1)
    Sc = 1370.0

    if cos_zen > 0.0:
        # remember sin_beta = cos_zenith; trig funcs are cofuncs of each other
        # sin(x) = cos(90-x) and cos(x) = sin(90-x).
        So = Sc * (1.0 + 0.033 * cos(float(doy) / 365.0 * 2.0 * pi)) * cos_zen
    else:
        So = 0.0

    return So

def spitters(doy, par, cos_zenith):
    """
    Spitters algorithm to estimate the diffuse component from the total daily
    incident radiation.

    NB. Eqns. 2a-d, not 20a-d

    Parameters:
    ----------
    doy : int
        day of year
    par : double
        daily total photosynthetically active radiation (MJ m-2 d-1)
    cos_zenith : float
        cosine of zenith angle (radians)

    Returns:
    -------
    diffuse : double
        diffuse component of incoming radiation

    References:
    ----------
    * Spitters, C. J. T., Toussaint, H. A. J. M. and Goudriaan, J. (1986)
      Separating the diffuse and direct component of global radiation and its
      implications for modeling canopy photosynthesis. Part I. Components of
      incoming radiation. Agricultural Forest Meteorol., 38:217-229.
    """

    # Fraction of global radiation that is PAR
    fpar = 0.5
    SEC_2_HFHR = 1800.0
    J_TO_MJ = 1E-6
    CONV = SEC_2_HFHR * J_TO_MJ

    # Calculate extra-terrestrial radiation
    S0 = 0.0
    for i in xrange(48):
        cos_zen = cos(acos(cos_zenith[i]) * 180.0 / pi) # degrees
        S0 += calc_extra_terrestrial_rad(doy, cos_zen) * CONV

    # Spitter's formula (Eqns. 2a-d)
    trans = (par / fpar) / S0
    if trans < 0.07:
        diffuse_frac = 1.0
    elif trans < 0.35:
        diffuse_frac = 1.0 - 2.3 * (trans - 0.07)**2
    elif trans < 0.75:
        diffuse_frac = 1.33 - 1.46 * trans
    else:
        diffuse_frac = 0.23

    return (diffuse_frac)

def calc_day_length(doy, yr_days, latitude):

    """
    Daylength in hours
    Eqns come from Leuning A4, A5 and A6, pg. 1196

    Reference:
    ----------
    Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.

    Parameters:
    -----------
    doy : int
        day of year, 1=jan 1
    yr_days : int
        number of days in a year, 365 or 366
    latitude : float
        latitude [degrees]

    Returns:
    --------
    dayl : float
        daylength [hrs]
    """

    deg2rad = pi / 180.0;
    rlat = latitude * deg2rad;
    sindec = -sin(23.5 * deg2rad) * cos(2.0 * pi * (doy + 10.0) / yr_days);
    a = sin(rlat) * sindec;
    b = cos(rlat) * cos(asin(sindec));

    return 12.0 * (1.0 + (2.0 / pi) * asin(a / b))


if __name__ == "__main__":

    main()
