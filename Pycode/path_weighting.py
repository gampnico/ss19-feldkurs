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


def effective_z(path_height, path_position, b):
    """Calculates the effective path height across the entire
    scintillometer path based on actual path height and position.

    Args:
        path_height (Series): actual path height, m
        path_position (Series): the position along the
            scintillometer path
        b (str) = stability conditions, either stable, unstable, or no
            height dependency
    Returns:
        z_eff (float): effective path height, m

    """

    if b.lower() in ["stable", "s"]:
        b = -2 / 3
    elif b.lower() in ["unstable", "u"]:
        b = -4 / 3
    else:
        b = 1
        print("No height dependency selected")

    ph_list = (np.multiply(path_height ** b, pwf(path_position)))

    # For path-building intersections
    ph_list[np.isnan(ph_list)] = 0
    z_eff = ((sp.integrate.trapz(ph_list)) / (sp.integrate.trapz(pwf(
        path_position)))) ** (1 / b)

    return z_eff


def return_z_effective(location_name):
    """Returns effective path height.

    Args:
        location_name (str): the name of the transmitter station

    Returns:
        effective_path_heights (dict): effective path heights for stable
            and unstable conditions.
    """
    if location_name.lower() == "h":
        location_name = "hungerburg"
    elif location_name.lower() in ["s", "schießstand"]:
        location_name = "schiessstand"

    path_string = "../MATLAB/path_height_" + location_name.lower() + ".csv"
    path_height_data = pd.read_csv(
        path_string, header=None,
        names=["path_height", "norm_position"])

    effective_path_height_stable = effective_z(
        path_height_data["path_height"], path_height_data["norm_position"],
        b="stable")
    effective_path_height_unstable = effective_z(
        path_height_data["path_height"], path_height_data["norm_position"],
        b="unstable")

    effective_path_heights = {"stable": effective_path_height_stable,
                              "unstable": effective_path_height_unstable}

    mean_path_height = np.mean(path_height_data["path_height"])
    path_weight = pwf(path_height_data["norm_position"])
    path_height_data["path_weight"] = path_weight

    print("Mean path height: " + str(mean_path_height) + "m")
    print("Effective path height, stable: " + str(effective_path_heights[
                                                      "stable"]) + "m")
    print("Effective path height, unstable: " + str(effective_path_heights[
                                                        "unstable"]) + "m")
    return effective_path_heights
