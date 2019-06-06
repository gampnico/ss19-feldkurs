#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Derives sensible heat flux, H.
"""

import math

import cn_derivations as cn
import numpy as np
import path_weighting as pw
import mpmath

# Preliminary derivation of H under free convection conditions i.e
# little wind shear, strong instability
k = 0.4  # von Karman constant
cp = 1004  # J kg^-1 K^-1, heat capacity of air
r_d = 287.05  # J kg^-1 K^-1, specific gas constant for dry air
beta_1 = 0.86
g = 9.81


def free_flux_compute(file_name):
    computed_df = cn.derive_ct2(file_name)
    # Calculate kinematic surface heat flux
    computed_df["Q_0"] = 1.165 * k * z_eff * (
            computed_df["CT2"] ** (3 / 4)) * (
                                 g / computed_df["temperature"]) ** (1 / 2)

    # Air density
    computed_df["rho_air"] = 100 * computed_df["pressure"] / (
            r_d * computed_df["temperature"])
    # Surface sensible heat flux under free convection
    computed_df["H_free"] = computed_df["Q_0"] * cp * computed_df[
        "rho_air"]
    return computed_df


def obukhov_iteration(row_index, inital_length=-1000.0):
    """Calculates Obukhov length using procedure in SRun Theory Manual
    Args:
        inital_length (float): inital guess for Obukhov length, default
            initial_length = -1000.0
    Returns:
        l_ob (pandas.core.series.Series): dataframe of computed Obukhov
            lengths for each time interval.
    """

    # Calculate corresponding value for f_H for unstable conditions
    zeta = z_eff / inital_length
    if zeta < 0:
        f_H = 4 * (k ** (-2 / 3)) * beta_1 * (
                1 - 7 * zeta + 75 * zeta ** 2) ** (
                      -1 / 3)
    else:
        f_H = 4 * (k ** (-2 / 3)) * beta_1 * (
                1 + 7 * zeta + 20 * zeta ** 2) ** (
                      -1 / 3)

    # Calculate corresponding value for Ψ_M for unstable conditions
    # using L_Ob and the height of the wind sensor z_u.
    z_u = 579  # from ZAMG
    mpmath.dps = 30
    mpmath.pretty = True
    x = mpmath.power((1 - 15 * zeta), (1 / 4))
    # x = (1 - 15 * zeta) ** (1 / 4)

    if zeta < 0:
        psi_M = -2 * math.log((1 + x) / 2) - math.log(
            (1 + x ** 2) / 2) + 2 * math.atan(x) - \
                math.pi / 2
    else:
        psi_M = 4.7 * zeta

    # roughness length z_0 are required.
    z_0 = 0.1 * z_u

    # Calculate temperature scale θ* with measured CT2
    computed_data["theta_star"].loc[row_index] = np.sqrt(
        (computed_data["CT2"].loc[row_index] * (z_eff ** (2 / 3))) / f_H)

    # Calculate friction velocity u*. Horizontal wind speed U(z_u)
    computed_data["u_star"].loc[row_index] = k * computed_data[
        "windspeed"].loc[row_index] / (math.log(z_u / z_0) - psi_M)

    # Calculate Obukhov length
    l_obu = (computed_data["temperature"].loc[row_index] *
             computed_data["u_star"].loc[row_index] ** 2) / (
                    k * g * computed_data["theta_star"].loc[row_index])
    return l_obu


station = "schiessstand"
filename = "2019-05-24"

z_eff = pw.return_z_effective(station)
computed_data = free_flux_compute(filename)
computed_data["theta_star"] = 0
computed_data["u_star"] = 0
print("Iteration started...\n")
for index, row in computed_data.iterrows():
    # Initialise first iteration step
    ob_ini = -1000
    l_ob = obukhov_iteration(index, ob_ini)
    computed_data.loc[index, "l_ob"] = l_ob
    step = 0
    while abs(l_ob - ob_ini) > 0.1:
        step += 1
        ob_prev = ob_ini
        ob_ini = l_ob
        l_ob = obukhov_iteration(index, ob_ini)
        # avoid infinite loops
        if math.isclose(l_ob, ob_prev, rel_tol=0.1) or step > 100:
            break
    print(l_ob)

print("Iteration completed!")
# Calculate sensible heat flux
computed_data["H"] = -computed_data["rho_air"] * cp \
                     * computed_data["theta_star"] * computed_data["u_star"]
# Calculate momentum flux
computed_data["momentum_flux"] = computed_data["rho_air"] * computed_data[
    "u_star"] ** 2

print(computed_data["H"])
