#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Derives sensible heat flux, H.
"""

import math

import cn_derivations as cn
import numpy as np
import path_weighting as pw
import mpmath


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
    dataframe = cn.derive_ct2(file_name)
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
