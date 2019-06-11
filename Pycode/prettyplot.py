import matplotlib.pyplot as plt
import matplotlib.dates as dates


def label_selector(dependent):
    if "shf" or "H" in str(dependent):
        plt.title(r"Sensible Heat Flux $Q_{H}$", fontweight="bold")
        plt.ylabel("Sensible Heat Flux, [W$\cdot$m$^{-2}$]")
    elif "theta_star" in str(dependent):
        plt.title(r"Temperature Scale $\theta^{*}$", fontweight="bold")
        plt.ylabel("Temperature Scale, [K]")
    elif "u_star" in str(dependent):
        plt.title(r"Friction Velocity $u^{*}$", fontweight="bold")
        plt.ylabel("Friction Velocity, [m$\cdot$s$^{-2}$]")
    elif "obukhov" in str(dependent):
        plt.title(r"Obukhov Length $L_{Ob}$", fontweight="bold")
        plt.ylabel("Obukhov Length, [m]")
    elif "temperature" in str(dependent):
        plt.title(r"Temperature $T$", fontweight="bold")
        plt.ylabel("Temperature, [K]")
    elif "pressure" in str(dependent):
        plt.title(r"Pressure $P$", fontweight="bold")
        plt.ylabel("Pressure, [mbar]")
    elif "Cn2" in str(dependent):
        plt.title(r"Structure Parameter of Refractive Index $C_{n}^{2}$",
                  fontweight="bold")
        plt.ylabel("$C_{n}^{2}$")
    elif "CT2" in str(dependent):
        plt.title(r"Structure Parameter of Temperature $C_{T}^{2}$",
                  fontweight="bold")
        plt.ylabel("$C_{T}^{2}$")


def plot_free_convection(scint, processed):
    hourly_sdf_mean = scint.resample("H").mean()
    hourly_sdf_mean["H_std"] = scint["H_convection"].resample("H").std()

    hourly_cdf_mean = processed.resample("H").mean()
    hourly_cdf_mean["H_free_std"] = processed["H_free"].resample("H").std()

    fig = plt.figure(figsize=(26, 6))
    scint["H_convection"].plot(color="black", label="Scintillometer")
    processed["H_free"].plot(color="red", label="Free Convection")

    hourly_sdf_mean["H_convection"].plot(label="Scintillometer hourly mean",
                                         color="black",
                                         linestyle="dashed")
    hourly_cdf_mean["H_free"].plot(label="Free convection hourly mean",
                                   color="red",
                                   linestyle="dashed")

    plt.legend(loc="upper left")
    plt.title(
        r"Sensible Heat Fluxes $Q_{H}$ from Scintillometer and for Free "
        r"Convection",
        fontweight="bold")
    plt.xlabel("Time")
    plt.ylabel("Sensible Heat Flux, [W$\cdot$m$^{-2}$]")
    ax = plt.gca()
    ax.xaxis.set_major_formatter(
        dates.DateFormatter('%H:%M'))  # hours and minutes
    return fig


def plot_generic(dataframe, name, colour="grey"):
    name_std = name + "_std"
    name_avg = name + "_avg"
    dataframe[name_std] = dataframe[name].astype("float64").resample("H").std()
    dataframe[name_avg] = dataframe[name].astype("float64").resample(
        "H").mean()

    fig = plt.figure(figsize=(26, 6))
    dataframe[name].astype("float64").plot(color=colour, label="Time series")
    dataframe[name_avg].astype("float64").plot(
        color="black", linestyle="dashed", label="Time series")

    plt.legend(loc="upper left")
    label_selector(name)
    plt.xlabel("Time, CET")
    ax = plt.gca()
    ax.xaxis.set_major_formatter(
        dates.DateFormatter('%H:%M'))  # hours and minutes
    plt.show()
    return fig


def plot_comparison(dataframe, name_1, name_2):
    fig = plt.figure(figsize=(26, 6))
    dataframe[name_1].plot(color="black", label="Scintillometer")
    dataframe[name_2].plot(color="red", label="Iterated Flux")

    plt.legend(loc="upper left")
    plt.title(
        r"Sensible Heat Fluxes $Q_{H}$ from Scintillometer and from Iteration",
        fontweight="bold")
    plt.xlabel("Time")
    plt.ylabel("Sensible Heat Flux, [W$\cdot$m$^{-2}$]")
    ax = plt.gca()
    ax.xaxis.set_major_formatter(
        dates.DateFormatter('%H:%M'))  # hours and minutes
    return fig
