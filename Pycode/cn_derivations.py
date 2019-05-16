#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Parses SRun data files.
"""

import data_parser as dp

scint_data = dp.data_parsing()
day = "2019-04-19"
# placeholder until we actually get acinn data
acinn_data = dp.weather_parsing(day, "on")

derived = scint_data.filter(["Cn2", "temp", "pressure"], axis=1)
T_0 = 273.15

derived["temp"] = acinn_data["t"][8] + T_0
derived["pressure"] = acinn_data["ldred"][8]

transmit_lambda = 880 * (10 ** -9)  # m
lambda_2 = 7.53 * (10 ** -3)  # micron^2
lambda_2 = lambda_2 / (10 ** 6) ** 2
alpha_factor_2 = 77.6 * (10 ** -6)  # K hPa^-1
alpha_factor_1 = alpha_factor_2 * (1 + (lambda_2 / transmit_lambda))

# This equation is accurate to 10^-5, but doesn't account for humidity
# fluctuations - errors should be within 3% of inverse Bowen ratio (Moene 2003)
derived["CT2"] = derived["Cn2"] * (derived["temp"] ** 4) * (
        (alpha_factor_1 * derived["pressure"]) ** -2)
print(derived["pressure"])
print(derived["CT2"] - scint_data["CT2"])
