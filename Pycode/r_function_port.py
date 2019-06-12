#!/usr/bin/env python3
"""Python port of R functions used for processing scintillometer data.
Original R functions by Helen Ward, port and minor editing by Nicolas
Gampierakis (2019)."""

import time

# Scintillometer constants
alpha11 = 4.9491e-2  # Equivalent to 4.48 * D^(7/3), BLS manual (App A.5)
lamda = 880e-9  # BLS wavelength 880 nm
m1_opt = 7.806109e-7  # Needed for AT and Aq, from Owens (1967)
m2_opt = 6.561286e-7  # Needed for AT and Aq, from Owens (1967)
At_opt = -270e-6  # AT coefficient for 880 nm and typical atmos conditions

Aq_opt = -0.685e-6  # Aq coefficient for 880 nm and typical atmos conditions

# Physical constants
R_dry = 287.04  # Specific gas constant for dry air [J K^-1 kg^-1]
R_vapour = 461.5  # Specific gas contstant for water vapour [J K^-1 kg^-1]
ratio_rmm = 0.622  # Ratio of molecular masses of water vapour and dry air
cp = 1004.67  # Specific heat capacity of air at constant pressure [J K^-1
# kg^-1]
latent_vapour = 2.45e6  # Latent heat of vapourisation (at 20 deg C) [J kg^-1]
rho = 1.225  # Density of air at STP [kg m^-3]
k = 0.4  # von K?rm?n's constant
g = 9.81  # Acceleration due to gravity [m s^-2]


def convert_rh_to_q(rh, temp_abs, pressure):
    """Calculates specific humidity (q) in kg kg^-3

    Args:
        rh (float): relative humidity, %
        temp_abs (float): absolute temperature, K
        pressure (float): atmospheric pressure, P

    Returns:
        specific_humidity (float): specific humidity, %

    """

    esat = 611.2 * math.exp(
        17.67 * (temp_abs - 273.16) / (temp_abs - 29.66))  # Stull
    q_abs = (ratio_rmm / (R_dry * temp_abs)) * rh / 100 * esat
    r = R_dry / (1 - q_abs * temp_abs / pressure * (R_vapour - R_dry))
    specific_humidity = q_abs * r * temp_abs / pressure
    return specific_humidity


def psi_m_unstable(obukhov, z):
    """Integrated stability function for momentum for unstable
        conditions, for when obukhov <0

    Args:
        obukhov (float): Obukhov length, m
        z (float): height, m

    Returns:
        pmu (float): value for instability function psi in unstable
            conditions

    """

    x = (1 - 16 * (z / obukhov)) ** (1 / 4)
    pmu = (2 * math.log((1 + x) / 2) + math.log((1 + x ** 2) / 2) -
           2 * math.atan(x) + math.pi / 2)
    return pmu


def psi_m_stable(obukhov, z):
    """Integrated stability function for momentum for stable conditions,

    Args:
        obukhov (float): Obukhov length, m
        z (float): height, m

    Returns:
        pmu (float): value for instability function psi in stable
            conditions

    """

    pms = (-5) * z / obukhov
    return pms


def psi_m(obukhov, z):
    """Calculates stability function psi according to stability
    conditions

    Args:
        obukhov (float): Obukhov length, m
        z (float): height, m

    Returns:
        pm (float): psi value for instability function

    """

    if obukhov > 0:
        pm = psi_m_stable(obukhov, z)
    else:
        pm = psi_m_unstable(obukhov, z)
    return pm


def calc_u_star(u, z_u, z0, obukhov):
    """Calculates friction velocity u* in ms^-1

    Args:
        u (float): wind speed, ms^-1
        z_u (float): height of wind speed measurement including
            displacement (z-d), m
        z0 (float): roughness length for momentum, m
        obukhov (float): Obukhov length, m

    Returns:
        cus (float): friction velocity u*, ms^-1
    """
    cus = k * u / (math.log(z_u / z0) - psi_m(obukhov, z=z_u)
                   + psi_m(obukhov, z=z0))
    return cus


def ward_iteration(dataframe, index, zm_bls):
    """Calculates Obukhov length, friction velocity, and temperature
    scale according to the method given in Helen Ward's scintillometry
    lecture.

    Args:
        dataframe (pandas.DataFrame): dataframe containing variables
            necessary for calculating returned variables
        index (int): index of row being iterated over, i.e. time step
            in time series
        zm_bls (float): effective path height

    Returns:
        obukhov_new (float): computed obukhov length for the given time
            and conditions

    """

    con_T1 = 4.9  # 4.9
    con_T2 = 6.1  # 9.0
    # con_T3 = 2.2  # 0
    iteration_step = 0
    z0 = zm_bls * 0.1

    while iteration_step <= 5:
        f_mo = con_T1 * (1 - con_T2 * (zm_bls / dataframe["obukhov"].loc[
            index])) ** (-2 / 3)
        dataframe["u_star"].loc[index] = calc_u_star(u=dataframe[
            "windspeed"].loc[index], z_u=zm_bls, z0=z0, obukhov=dataframe[
            "obukhov"].loc[index])
        dataframe["theta_star"].loc[index] = (-1) * math.sqrt(
            dataframe["CT2"].loc[index] * (zm_bls ** (2 / 3)) / f_mo)
        # Calculate new Obukhov length
        obukhov_new = (dataframe["temperature"].loc[index] * dataframe[
            "u_star"].loc[index] ** 2) / (
                              g * k * dataframe["theta_star"].loc[index])
        # Store difference between new and current Obukhov length
        dataframe["obukhov_diff"].loc[index] = (obukhov_new - dataframe[
            "obukhov"].loc[index])
        dataframe.loc[index, "l_ob"] = obukhov_new
        # Calculate new heat flux
        h_new = (-1) * dataframe["rho_air"].loc[index] * cp * dataframe[
            "u_star"].loc[index] * dataframe["theta_star"].loc[index]
        # Store difference between new and current heat flux
        #         cdf["H_iter_diff"].loc[index] = (H_new - cdf[
        #         "H_iter"].loc[index])
        dataframe["shf"].loc[index] = h_new
        iteration_step += 1
        return obukhov_new


def ward_method(dataframe, eff_h, switch_time):
    """Overarching function for calculating Obukhov lengths and sensible
    heat fluxes as a time series.

    Args:
        dataframe (pandas.DataFrame): dataframe containing data
            necessary for deriving fluxes
        eff_h (dict): contaings effective path heights for stable and
            unstable conditions
        switch_time (string): the time at which the atmosphere switches
            from stable to unstable conditions. This must be estimated
            manually.

    Returns:
        dataframe (pandas.DataFrame): dataframe with additional data for
            Obukhov lengths, sensible heat fluxes, frictional velocity,
            and temperature scale

    """
    start = time.time()
    z_eff = eff_h
    print("\nIteration started")
    # Initialise dataframe
    dataframe["obukhov"] = -100  # unstable conditions
    dataframe["shf"] = float
    dataframe["u_star"] = float
    dataframe["theta_star"] = float
    dataframe["obukhov_diff"] = float
    dataframe["shf_diff"] = float

    dataframe_stable = dataframe.iloc[
        dataframe.index.indexer_between_time("00:00", switch_time)]
    for index, row in dataframe_stable.iterrows():
        obukhov = ward_iteration(dataframe_stable, index, z_eff["stable"])
        dataframe_stable.loc[index, "obukhov"] = obukhov
    dataframe_unstable = dataframe.iloc[
        dataframe.index.indexer_between_time(switch_time, "23:59")]
    for index, row in dataframe_unstable.iterrows():
        obukhov = ward_iteration(dataframe_stable, index, z_eff["unstable"])
        dataframe_unstable.loc[index, "obukhov"] = obukhov
    dataframe = dataframe_stable.append(dataframe_unstable)
    end = time.time()
    print(str(end - start) + "s")
    print("\nIteration completed!")
    return dataframe
