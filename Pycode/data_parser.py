#!/usr/bin/env python3
"""Nicolas Gampierakis (2019). Parses SRun data files.
"""

import pathlib

import download_weather as dw
import pandas as pd


def scintillometer_parse(filename):
    """Parses .mnd files that follow the format specified in header_list
    Args:
        filename (str): the name of the .mnd data file

    Returns:
        scint_data (pandas.DataFrame): correctly indexed scintillometry
            data
    """
    header_list = ["Time", "Cn2", "CT2", "H_convection", "crosswind",
                   "sigCrosswind", "pressure", "temp", "humidity",
                   "pathLength", "pathHeight", "correctCn2EO",
                   "correctCn2Sat", "correctCn2Cov", "mndCounter",
                   "<XA>(c)", "<YA>(c)", "nSigXA(c)", "nSigYA(c)",
                   "corXAYA(c)", "numDgnValid", "numDgnValidCrosswind",
                   "numDgnTotal", "channelFlagsCombined", "error"]

    path_string = "../../data/SRun/" + filename + ".mnd"
    scint_data = pd.read_csv(path_string, header=None,
                             skiprows=35, names=header_list, sep="\t")

    # Remove timestamp fluff
    scint_data["Time"] = scint_data["Time"].str.replace("PT00H00M21S/", "")
    scint_data["Time"] = scint_data["Time"].str.replace("PT00H00M23S/", "")
    scint_data["Time"] = scint_data["Time"].str.replace("PT00H00M29S/", "")
    scint_data["Time"] = scint_data["Time"].str.replace("PT00H00M59S/", "")
    scint_data["Time"] = pd.to_datetime(scint_data["Time"])
    scint_data = scint_data.set_index("Time")
    # scint_data = scint_data[scint_data["error"] == 0]
    scint_data = scint_data.tz_convert("CET")
    # scint_data.index = scint_data.index + pd.DateOffset(hours=2)
    return scint_data


def weather_download(day, connect="off"):
    """Downloads and parses external weather data for the university
    rooftop.

    Args:
        day (str): the date we wish to get observations for, YYY-MM-DD
        connect (str): switch for downloading new datafiles. Default on

    Returns:
        weather_data (pandas.DataFrame): data frame containing hourly
        weather observations for university rooftop.
    """

    # Matches format used in the weather service files
    variable_list = ["t", "rf", "wr", "wg", "regen", "ldred", "ldstat",
                     "sonne"]
    header_list = ["station_no", "station_name", "number", "date", "time",
                   "variable", "unit", "datetime"]
    weather_data = pd.DataFrame()

    # Download observations for each variable
    for variable in variable_list:
        client_path = pathlib.Path.cwd().joinpath("weather_data/" + variable
                                                  + ".csv")
        if connect == "on":
            server_path = "http://at-wetter.tk/api/v1/station/11121/" \
                          + variable + "/" + day
            dw.download_database(server_path, client_path)
        else:
            pass

        # Parse data
        temp_data = pd.read_csv(client_path, header=None, names=header_list,
                                sep="';'", engine="python")
        temp_data["datetime"] = pd.to_datetime(temp_data["datetime"])
        temp_data = temp_data.set_index("datetime")
        # append to returned data
        weather_data[variable] = \
            temp_data["variable"].resample("T").pad().fillna(
                method="bfill").tz_localize("UTC")

    return weather_data


def weather_parsing(day):
    """Downloads and parses university-provided rooftop weather data.

        Args:
            day (str): the date we wish to get observations for, YYY-MM-DD

        Returns:
            weather_data (pandas.DataFrame): data frame containing
            minute-interval weather observations for university rooftop.
        """

    header_list = ["Datetime", "airtemp", "airtemp_hygro", "soiltemp+0.5m",
                   "soiltemp-.1m", "soiltemp-0.2m", "soiltemp-0.5m",
                   "dewpoint", "relhum", "relhum_hygro", "prcp_weight",
                   "prcp", "pressure", "sunshine_duration", "global_rad",
                   "diffuse_sky_rad", "windspeed_scalar", "windspeed_vector",
                   "winddir", "windspeed_max", "gustdir", "shm30"]

    weather_data = pd.read_csv("../../data/weather_data/data_TAWES_UIBK.csv",
                               header=None, skiprows=2, index_col=0,
                               parse_dates=True, names=header_list, sep=";")
    weather_data = weather_data.loc[day]

    return weather_data
