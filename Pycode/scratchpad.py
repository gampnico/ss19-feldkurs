"""Ward's iteration"""

import cn_derivations as cn
import data_parser as dp
import path_weighting as pw
import r_function_port as rfp

filename = "2019-05-24"
k = 0.4  # von Karman constant
cp = 1004  # J kg^-1 K^-1, heat capacity of air
r_d = 287.05  # J kg^-1 K^-1, specific gas constant for dry air
beta_1 = 0.86
g = 9.81
scint_data = dp.scintillometer_parse(filename)
station = "schiessstand"
filename = "2019-05-24"
z_eff = pw.return_z_effective(station)
cdf = cn.derive_ct2(filename)
cdf = cdf.tz_convert("CET")

# Calculate kinematic surface heat flux
cdf["Q_0"] = 1.165 * k * z_eff * (cdf["CT2"] ** (3 / 4)) * (
        g / cdf["temperature"]) ** (1 / 2)

# Air density
cdf["rho_air"] = 100 * cdf["pressure"] / (
        r_d * cdf["temperature"])
# Surface sensible heat flux under free convection
cdf["H_free"] = cdf["Q_0"] * cp * cdf["rho_air"]

cdf = rfp.ward_method(cdf, z_eff)
