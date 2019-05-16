#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Scintillometer path weighting and
effective height calculator from path position and normalised height.
"""

import math

import numpy as np
import pandas as pd
import scipy as sp
from scipy import integrate, special


def bessel_second(x):
    """Calculates the bessel function for the specific path position.

    Args:
        x (float): path position, m

    Returns:
        y (float): bessel function of path point
    """
    bessel_variable = 2.283 * math.pi * (x - 0.5)
    if bessel_variable == 0:
        y = 1
    else:
        y = 2 * (sp.special.jv(1, bessel_variable)) / bessel_variable
    return y


def pwf(path_position):
    """Path weighting function for effective path height

    Args:
        path_position (Series): input data array
    Returns:
        weight (float): weight
    """
    weight = []
    for i in path_position:
        weight.append(2.163 * bessel_second(i))
    return weight


def effective_z(path_height, path_position):
    """Calculates the effective path height across the entire
    scintillometer path based on actual path height and position.

    Args:
        path_height (Series): actual path height, m
        path_position (Series): the position along the
            scintillometer path.
    Returns:
        z_eff (float): effective path height, m

    """
    # b = 1  # No height dependency
    b = - 2 / 3  # Stable
    # b = - 4 / 3  # Unstable
    ph_list = (np.multiply(path_height ** b, pwf(path_position)))

    # For path-building intersections
    ph_list[np.isnan(ph_list)] = 0
    z_eff = ((sp.integrate.trapz(ph_list)) / (sp.integrate.trapz(pwf(
        path_position)))) ** (1 / b)
    test = (sp.integrate.trapz(ph_list)) / (sp.integrate.trapz(pwf(
        path_position)))
    print(test)

    return z_eff


path_height_data = pd.read_csv(
    "../MATLAB/path_height_schiessstand.csv", header=None,
    names=["path_height", "norm_position"])

effective_path_height = effective_z(
    path_height_data["path_height"], path_height_data["norm_position"])

mean_path_height = np.mean(path_height_data["path_height"])
path_weight = pwf(path_height_data["norm_position"])
path_height_data["path_weight"] = path_weight

print("Mean path height: " + str(mean_path_height) + "m")
print("Effective path height: " + str(effective_path_height) + "m")
# For Matlab export
# path_height_data.to_csv(path_or_buf="../MATLAB/hungerburg_sim.csv",
#                         index=None)
