#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Parses SRun data files.
"""

import numpy as np
import pandas as pd

header_list = ["Time", "Cn2", "CT2", "Heat Flux (Free)", "Crosswind",
               "SD Crosswind", "Pressure", "Temperature", "Rel H.",
               "Path Length", "Path Height", "Cn2 Extinction",
               "Cn2 Saturation", "Cn2 Covariance Correction",
               "Data Counter", "Av XA", "Av YA", "Norm SD XA", "Norm SD YA",
               "Correl XA/YA", "VDSP", "VDSP Crosswind", "Total VDSP",
               "CCF", "ADR Failure", "Error Code"]
scintill_data = pd.read_csv("../../SRun/2019-04-09a.mnd", header=None,
                            skiprows=35, names=header_list, sep="\t")

# Remove timestamp fluff
scintill_data["Time"] = scintill_data["Time"].str.replace("PT00H00M59S/", "")
scintill_data["Time"] = scintill_data["Time"].str.replace("PT00H00M29S/", "")

scintill_data['Time'] = pd.to_datetime(scintill_data['Time'])

scintill_data = scintill_data.set_index("Time")
print(scintill_data)
