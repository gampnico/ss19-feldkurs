#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Derives sensible heat flux, H.
"""

import path_weighting as pw
import cn_derivations as cn

# Preliminary derivation of H under free convection conditions i.e
# little wind shear, strong instability

k = 0.4  # von Karman constant
cp = 1004  # J kg^-1 K^-1, heat capacity of air
r_d = 287.05  # J kg^-1 K^-1, specific gas constant for dry air
location = input("Please enter the location name (Hungerburg or "
                 "Schießstand):\n")

# file_name = input("Please enter the name of the file (without extension), "
#                   "eg 2019-04-09a")

file_name = "2018-05-06"
computed_data = cn.derive_ct2(file_name)

# Calculate kinematic surface heat flux

computed_data["Q_0"] = 1.165 * k * pw.return_z_effective(location) * (
        computed_data["CT2"] ** (3 / 4)) * (
                               9.81 / computed_data["temp"]) ** (1 / 2)

# Air density
computed_data["rho_air"] = 100 * computed_data["pressure"] / (
        r_d * computed_data["temp"])
# Surface sensible heat flux
computed_data["H"] = computed_data["Q_0"] * cp * computed_data["rho_air"]
print(computed_data)
