import r_function_port as rfp
import pandas as pd
from shf_derivation import computed_data as cdf

cdf["LOb_iter"] = 100  # stable conditions
cdf["H_iter"] = 0
NRuns = 0
# Do the iteration
while NRuns <= 5:
    fMO = ConT1 * (1 - ConT2 * (zm_BLS / cdf["LOb_iter"])) ** (-2 / 3)
    cdf["ustar_iter"] = rfp.calc_u_star(u=cdf["windspeed"], Zu=zm_BLS, z0=z0,
                                        Lob=cdf["LOb_iter"])
    cdf["Tstar_iter"] = (-1) * mpmath.sqrt(
        cdf["CT2"] * (zm_BLS ** (2 / 3)) / fMO)
    # Calculate new Obukhov length
    Lob_new = (cdf["ustar_iter"] ** 2 * cdf["T_K"]) / (
                g * Kv * Data["Tstar_iter"])
    # Store difference between new and current Obukhov length
    cdf["LOb_iter_diff"] = (Lob_new - cdf["LOb_iter"])
    cdf["LOb_iter"] < - Lob_new
    # Calculate new heat flux
    H_new < - (-1) * Rho * Cp * cdf["ustar_iter"] * cdf["Tstar_iter"]
    # Store difference between new and current heat flux
    cdf["H_iter_diff"] < - (H_new - cdf["H_iter"])
    cdf["H_iter"] < - H_new
    NRuns = NRuns + 1
