#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Parses SRun data files.
"""

import numpy as np
import pandas as pd

header_list = ["Time", "Cn2", "CT2", "H_convection", "crosswind",
               "sigCrosswind", "pressure", "temp", "humidity", "pathLength",
               "pathHeight", "correctCn2EO", "correctCn2Sat", "correctCn2Cov",
               "mndCounter", "<XA>(c)", "<YA>(c)", "nSigXA(c)", "nSigYA(c)",
               "corXAYA(c)", "numDgnValid", "numDgnValidCrosswind",
               "numDgnTotal", "channelFlagsCombined", "error"]

scint_data = pd.read_csv("../../SRun/2019-04-09a.mnd", header=None,
                         skiprows=35, names=header_list, sep="\t")

# Remove timestamp fluff
scint_data["Time"] = scint_data["Time"].str.replace("PT00H00M59S/", "")
scint_data["Time"] = scint_data["Time"].str.replace("PT00H00M29S/", "")

scint_data['Time'] = pd.to_datetime(scint_data['Time'])

scint_data = scint_data.set_index("Time")
print(scint_data)
