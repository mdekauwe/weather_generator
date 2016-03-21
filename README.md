# Weather generator

- Calculate the dirurnal air temperature from Tmin & Tmax following Parton and Logan (1981). I have added the missing displacement factor in the original Parton paper, which removes the discontinuity and forces the temp to reach Tmin.
The match is very close to the routinue in MAESTRA.

- Disaggregate daily rainfall following function in MAESTRA, which comes from
 Loustau et al. (2001).

- Interpolate between VPD 9am and 3pm values from AWAP data follow Haverd et al. (2013)

- Estimate diurnal radiation time course using spitters and routinue from
  MAESTRA
