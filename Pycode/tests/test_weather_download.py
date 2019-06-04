from unittest import TestCase
import pathlib
import pandas as pd


def weather_download(day):
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
        client_path = pathlib.Path.cwd().joinpath("test_data/"
                                                  + variable + ".csv")
        # Parse data
        temp_data = pd.read_csv(client_path, header=None, names=header_list,
                                sep="';'", engine="python")
        temp_data["datetime"] = pd.to_datetime(temp_data["datetime"])
        temp_data = temp_data.set_index("datetime")
        # append to returned data
        weather_data[variable] =\
            temp_data["variable"].resample("T").pad().fillna(method="bfill")
    return weather_data


def test_weather_download():
    wd = weather_download("2019-05-24")
    assert wd["t"][0]


wd = weather_download("2019-05-24")
print(wd)
