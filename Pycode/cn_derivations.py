#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Derives CT^2.
"""

import re
import data_parser as dp


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
    scint_data = dp.data_parsing(data_file)

    # scrub extra characters from file name
    day = re.sub("[^0-9\-]", "", data_file)
    # placeholder until we can access acinn data
    acinn_data = dp.weather_parsing(day, "off")

    derived = scint_data.filter(["Cn2", "temp", "pressure"], axis=1)
    kelvin = 273.15
    derived["temp"] = acinn_data["t"][8] + kelvin
    derived["pressure"] = acinn_data["ldred"][8]

    transmit_lambda = 880 * (10 ** -9)  # m
    lambda_2 = 7.53 * (10 ** -3)  # micron^2
    lambda_2 = lambda_2 / (10 ** 6) ** 2
    alpha_factor_2 = 77.6 * (10 ** -6)  # K hPa^-1
    alpha_factor_1 = alpha_factor_2 * (1 + (lambda_2 / transmit_lambda))

    # This equation is accurate to 10^-5, but doesn't account for humidity
    # fluctuations - errors should be within 3% of inverse Bowen ratio
    # (Moene 2003)

    derived["CT2"] = derived["Cn2"] * (derived["temp"] ** 4) * (
            (alpha_factor_1 * derived["pressure"]) ** -2)

    return derived
