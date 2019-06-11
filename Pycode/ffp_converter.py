"""Converts data into FFP readable format

CSV file format:
yyyy = Year
mm = month
day = day
HH = hour, UTC
MM = minutes
zm = measurement height above ground
d = displacement height
z0 = roughness length
u_mean = mean wind speed at zm
L = Obukhov length
sigma_v = Standard deviation of lateral velociity fluctuations (urms)
u_star = friction_velocity
wind_dir = wind direction in degrees
"""
import data_parser as dp
import cn_derivations as cn
import r_function_port as rfp

filename = "2019-05-24"
station = "schiessstand"
sdf = dp.scintillometer_parse(filename)
pdf = cn.data_processor(filename, station)
cdf = pdf["computed"]
z_eff = pdf["effective_height"]
cdf = rfp.ward_method(cdf, z_eff)
cdf = cdf.tz_convert("UTC")
cdf.index = cdf.reset_index
print(cdf)

# path_height_data.to_csv(path_or_buf="../MATLAB/hungerburg_sim.csv",
#                         index=None)
