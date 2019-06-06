import math
import mpmath

mpmath.dps = 50
mpmath.pretty = True
# Scintillometer constants
alpha11 = 4.9491e-2  # Equivalent to 4.48 * D^(7/3), BLS manual (App A.5)
lamda = 880e-9  # BLS wavelength 880 nm
m1_opt = 7.806109e-7  # Needed for AT and Aq, from Owens (1967)
m2_opt = 6.561286e-7  # Needed for AT and Aq, from Owens (1967)
At_opt = -270e-6  # AT coefficient for 880 nm and typical atmos conditions

Aq_opt = -0.685e-6  # Aq coefficient for 880 nm and typical atmos conditions

# Physical constants
Rd = 287.04  # Specific gas constant for dry air [J K^-1 kg^-1]
Rv = 461.5  # Specific gas contstant for water vapour [J K^-1 kg^-1]
RatioRMM = 0.622  # Ratio of molecular masses of water vapour and dry air
Cp = 1004.67  # Specific heat capacity of air at constant pressure [J K^-1
# kg^-1]
Lv = 2.45e6  # Latent heat of vapourisation (at 20 deg C) [J kg^-1]
Rho = 1.225  # Density of air at STP [kg m^-3]
Kv = 0.4  # von K?rm?n's constant
g = 9.81  # Acceleration due to gravity [m s^-2]


# Useful functions
# -------------------------------------------------------------
# ------------------------------------------------------------------------------
# FUNCTION: Calculates specific humidity (q) in kg kg^-3
#   from relative humidity (RH) in %
#        absolute temperature (Tabs) in K
#        atmospheric pressure (P) in Pa
def qfromRH(RH, Tabs, P):
    esat = 611.2 * mpmath.exp(
        17.67 * (Tabs - 273.16) / (Tabs - 29.66))  # Stull
    Qabs = (RatioRMM / (Rd * Tabs)) * RH / 100 * esat
    R = Rd / (1 - Qabs * Tabs / P * (Rv - Rd))
    qfromRH = Qabs * R * Tabs / P
    return qfromRH


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# FUNCTION: Integrated stability function for momentum for unstable conditions
#   apply when Lob < 0
def Psi_m_unstable(Lob, z):
    # Lob   Obukhov Length / m
    # z     Height / m
    x = (1 - 16 * (z / Lob)) ^ (1 / 4)
    pmu = (2 * math.log((1 + x) / 2) + math.log((1 + x ^ 2) / 2) -
                      2 * mpmath.atan(x) + math.pi / 2)
    return pmu


# -------------------------------------------------------------------------------
# FUNCTION: Integrated stability function for momentum for stable conditions
#   apply when Lob > 0
def Psi_m_stable(Lob, z):
    # Lob   Obukhov Length / m
    # z     Height / m
    pms = (-5) * z / Lob
    return pms


# -------------------------------------------------------------------------------
# FUNCTION: Integrated stability function for either stability
#   Uses stable function for Lob > 0; unstable otherwise (i.e. Lob < 0)
def Psi_m(Lob, z):
    # Lob   Obukhov Length / m
    # z     Height / m (e.g. (Zu - d), z0)
    if Lob > 0:
        pm = Psi_m_stable(Lob, z)
    else:
        pm = Psi_m_unstable(Lob, z)
    return pm


# -------------------------------------------------------------------------------


# FUNCTION: Calculates friction velocity in m s-1
def calc_u_star(u, Zu, z0, Lob):
    # u     Wind speed / m s^-1
    # Zu    Height of wind speed measurement including displacement (z-d)
    # zO    Roughness length for momentum / m
    # Lob   Obukhov length / m

    cus = Kv * u / (
                math.log(Zu / z0) - Psi_m(Lob, z=Zu) + Psi_m(Lob, z=z0))
    return cus
# ------------------------------------------------------------------------------
