#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Parses SRun data files.
"""

import pandas as pd
import get_weather as gw
import pathlib


def data_parsing():
    header_list = ["Time", "Cn2", "CT2", "H_convection", "crosswind",
                   "sigCrosswind", "pressure", "temp", "humidity",
                   "pathLength", "pathHeight", "correctCn2EO",
                   "correctCn2Sat", "correctCn2Cov", "mndCounter",
                   "<XA>(c)", "<YA>(c)", "nSigXA(c)", "nSigYA(c)",
                   "corXAYA(c)", "numDgnValid", "numDgnValidCrosswind",
                   "numDgnTotal", "channelFlagsCombined", "error"]

    scint_data = pd.read_csv("../../SRun/2019-04-09a.mnd", header=None,
                             skiprows=35, names=header_list, sep="\t")

    # Remove timestamp fluff
    scint_data["Time"] = scint_data["Time"].str.replace("PT00H00M59S/", "")
    scint_data["Time"] = scint_data["Time"].str.replace("PT00H00M29S/", "")

    scint_data['Time'] = pd.to_datetime(scint_data['Time'])

    scint_data = scint_data.set_index("Time")
    return scint_data


def weather_parsing(day, switch):
    variable_list = ["t", "rf", "wr", "regen", "ldred", "ldstat", "sonne"]
    header_list = ["station_no", "station_name", "number", "date", "time",
                   "variable", "unit", "datetime"]
    weather_data = pd.DataFrame()
    for variable in variable_list:
        client_path = pathlib.Path.cwd().joinpath("weather_data/" + variable
                                                  + ".csv")
        if switch == "on":
            server_path = "http://at-wetter.tk/api/v1/station/11121/" \
                      + variable + "/" + day
            gw.download_database(server_path, client_path)
        else:
            pass
        temp_data = pd.read_csv(client_path, header=None,
                                names=header_list, sep="';'")

        temp_data["datetime"] = pd.to_datetime(temp_data["datetime"])
        temp_data = temp_data.set_index("datetime")
        weather_data[variable] = temp_data["variable"]
    return weather_data
