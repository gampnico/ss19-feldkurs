#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Derives CT^2.
"""

import re
import data_parser as dp
import path_weighting as pw


def derive_ct2(data_file):
    """Derives the CT^2 structure parameter from measured Cn^2,
    and external weather data.

    Args:
        data_file (str): the name of the data file, without the .mnd
            extension, located in ../data/SRun/ folder
    Returns:
        derived (pandas.dataFrame): data frame containing derived CT^2
            alongside weather conditions
    """
    scint_data = dp.scintillometer_parse(data_file)

    # scrub extra characters from file name
    day = re.sub("[^0-9\-]", "", data_file)
    # placeholder until we can access acinn data
    acinn_data = dp.weather_download(day, "off")

    derived = scint_data.filter(["Cn2"], axis=1)
    kelvin = 273.15
    # merge weather and scintillometer data according to datetimes
    derived = derived.join(acinn_data[["t", "ldred", "wg"]]).rename(columns={
        "t": "temperature", "ldred": "pressure", "wg": "windspeed"})
    # adjust values
    derived["temperature"] = derived["temperature"] + kelvin
    derived["windspeed"] = derived["windspeed"] / 3.6  # convert to ms^-1

    transmit_lambda = 880 * (10 ** -9)  # m
    lambda_2 = 7.53 * (10 ** -3)  # micron^2
    lambda_2 = lambda_2 / (10 ** 6) ** 2
    alpha_factor_2 = 77.6 * (10 ** -6)  # K hPa^-1
    alpha_factor_1 = alpha_factor_2 * (1 + (lambda_2 / transmit_lambda))

    # This equation is accurate to 10^-5, but doesn't account for humidity
    # fluctuations - errors should be within 3% of inverse Bowen ratio
    # (Moene 2003)

    derived["CT2"] = derived["Cn2"] * (derived["temperature"] ** 4) * (
            (alpha_factor_1 * derived["pressure"]) ** -2)

    derived = derived.tz_convert("CET")

    return derived


def kinematic_shf(dataframe, z_eff):
    k = 0.4  # von Karman constant
    g = 9.81
    dataframe["Q_0"] = 1.165 * k * z_eff * (
            dataframe["CT2"] ** (3 / 4)) * (
                               g / dataframe["temperature"]) ** (1 / 2)
    return dataframe


def compute_fluxes(file_name, z_eff):
    r_d = 287.05  # J kg^-1 K^-1, specific gas constant for dry air
    cp = 1004  # J kg^-1 K^-1, heat capacity of air
    dataframe = derive_ct2(file_name)
    # Calculate kinematic surface heat flux
    dataframe = kinematic_shf(dataframe, z_eff)
    # Air density
    dataframe["rho_air"] = 100 * dataframe["pressure"] / (
            r_d * dataframe["temperature"])
    # Surface sensible heat flux under free convection (little wind shear,
    # high instability)
    dataframe["H_free"] = dataframe["Q_0"] * cp * dataframe[
        "rho_air"]
    return dataframe


def data_processor(filename, station):
    # Get effective path height
    z_eff = pw.return_z_effective(station)
    # Derive CT2, merge weather data, calculate free convection fluxes
    dataframe = compute_fluxes(filename, z_eff)
    processed_data = {"computed": dataframe, "effective_height": z_eff}
    return processed_data
