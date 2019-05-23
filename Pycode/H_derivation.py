#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Derives sensible heat flux, H.
"""

import path_weighting as pw
import cn_derivations as cn
import numpy as np
import math

# Preliminary derivation of H under free convection conditions i.e
# little wind shear, strong instability
k = 0.4  # von Karman constant
cp = 1004  # J kg^-1 K^-1, heat capacity of air
r_d = 287.05  # J kg^-1 K^-1, specific gas constant for dry air
beta_1 = 0.86

location = input("Please enter the location name (Hungerburg or "
                 "Schießstand):\n")

file_name = input("Please enter the name of the file (without extension), "
                  "eg 2019-04-09a")

file_name = "2018-05-06"
computed_data = cn.derive_ct2(file_name)

# Calculate kinematic surface heat flux
z_eff = pw.return_z_effective(location)
computed_data["Q_0"] = 1.165 * k * z_eff * (
        computed_data["CT2"] ** (3 / 4)) * (
                               9.81 / computed_data["temp"]) ** (1 / 2)

# Air density
computed_data["rho_air"] = 100 * computed_data["pressure"] / (
        r_d * computed_data["temp"])
# Surface sensible heat flux
computed_data["H"] = computed_data["Q_0"] * cp * computed_data["rho_air"]
print(computed_data)

# Iterative scheme

# Start calculation with initial Obukhov length
L_Ob = -1000
# Calculate corresponding value for f_H for unstable conditions according
zeta = z_eff / L_Ob
if zeta < 0:
    f_H = 4 * (k ** (-2 / 3)) * beta_1 * (1 - 7 * zeta + 75 * zeta ** 2) ** (
            -1 / 3)
else:
    f_H = 4 * (k ** (-2 / 3)) * beta_1 * (1 + 7 * zeta + 20 * zeta ** 2) ** (
            -1 / 3)

# Calculate temperature scale θ* with measured CT2
computed_data["theta_star"] = np.sqrt((computed_data["CT2"]
                                       * (z_eff ** (2 / 3))) / f_H)

# Calculate corresponding value for Ψ_M for unstable conditions using L_Ob and
# the height of the wind sensor z_u
z_u = 579  # from ZAMG
x = (1 - 15 * zeta) ** (1 / 4)

if zeta < 0:
    psi_M = -2 * math.log((1 + x) / 2) - math.log(
        (1 + x ** 2) / 2) + 2 * math.atan(x) - \
            math.pi / 2
else:
    psi_M = 4.7 * zeta

# Calculate friction velocity u*. Horizontal wind speed U(z_u) roughness
# length z_0 are required.
z_0 = 0.1 * z_u


