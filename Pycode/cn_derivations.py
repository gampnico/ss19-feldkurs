#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Derives CT^2 from parsed scintillometer
data and calculates fluxes.
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
    """Calculate kinematic surface heat flux.

    Args:
        dataframe (pandas.DataFrame): contains weather conditions and
            scintillometer data

    Returns:
        dataframe (pandas.DataFrame): includes additional column for
            kinematic surface heat flux
    """
    k = 0.4  # von Karman constant
    g = 9.81
    dataframe["Q_0"] = 1.165 * k * z_eff * (
            dataframe["CT2"] ** (3 / 4)) * (
                               g / dataframe["temperature"]) ** (1 / 2)
    return dataframe


def compute_fluxes(file_name, z_eff):
    """Merges weather data, derives CT2, computes various fluxes incl.
    free convection.

    Args:
        file_name (str): the name of the file, without extension. Should
            be .mnd format with a date within the file name
        z_eff (dict): dictionary containing effective path heights for
            stable and unstable conditions

    Returns:
        dataframe (pandas.DataFrame): modified dataframe with columns
            containing the fluxes calculated for the time series

    """

    r_d = 287.05  # J kg^-1 K^-1, specific gas constant for dry air
    cp = 1004  # J kg^-1 K^-1, heat capacity of air

    switch_time = input("\nPlease enter the time at which stable"
                        " switches to unstable: ")

    dataframe = derive_ct2(file_name)

    # Calculate kinematic surface heat flux for stable conditions

    dataframe_stable = dataframe.iloc[
        dataframe.index.indexer_between_time("00:00", switch_time)]
    dataframe_stable = kinematic_shf(dataframe_stable, z_eff["stable"])
    dataframe_unstable = dataframe.iloc[
        dataframe.index.indexer_between_time(switch_time, "23:59")]
    dataframe_unstable = kinematic_shf(dataframe_unstable, z_eff["unstable"])

    dataframe = dataframe_stable.append(dataframe_unstable)

    # Air density
    dataframe["rho_air"] = 100 * dataframe["pressure"] / (
            r_d * dataframe["temperature"])
    # Surface sensible heat flux under free convection (little wind shear,
    # high instability)
    dataframe["H_free"] = dataframe["Q_0"] * cp * dataframe[
        "rho_air"]
    return dataframe


def data_processor(filename, station):
    """Overarching function handling data post-processing.

    Args:
        filename (str): the name of the file, without extension. Should
            be .mnd format with a date within the file name
        station (str): the name of the transmitter station

    Returns:
        processed_data (dict): contains dataframe with calculated fluxes
        and effective path heights under various conditions
    """
    # Get effective path height
    z_eff = pw.return_z_effective(station)
    # Derive CT2, merge weather data, calculate free convection fluxes

    # input hour in format str(2204)

    dataframe = compute_fluxes(filename, z_eff)
    processed_data = {"computed": dataframe, "z_eff_stable": z_eff[
        "stable"], "z_eff_unstable": z_eff["unstable"]}
    return processed_data
