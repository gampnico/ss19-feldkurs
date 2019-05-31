# -*- coding: utf-8 -*-
"""
Created on Thu May 30 22:21:19 2019

@author: josef
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import Constants_and_functions as C

# Temperature
T = pd.read_csv('HATPRO_temp.csv', sep=';', comment='#', header=0, names=C.Z,
                index_col=0, date_parser=pd.to_datetime)

# Humidity
a = pd.read_csv('HATPRO_humidity.csv', sep=';', comment='#', header=0,
                names=C.Z, index_col=0, date_parser=pd.to_datetime)
a = a*10**(-3)  # to get kg/m^-3

# water vapour pressure
e = C.e(T, a)

# pressure of dry air
p_d = pd.DataFrame(index=T.index, columns=T.columns)
p_d[0] = 1000*10**2  # initial condition, have to replaced by weather data
for i in np.arange(1, len(p_d.columns)):
    p_d.iloc[:, i] = C.p_z2(p_d.iloc[:, i-1], p_d.columns[i-1], p_d.columns[i],
            T.iloc[:, i-1])

# mixing ratio
w = C.w(e, p_d)

# virtual Temperature
Tv = C.Tv(T, w)


# pressure
p = pd.DataFrame(index=T.index, columns=T.columns)
p[0] = 1000*10**2  # initial condition
for i in np.arange(1, len(p.columns)):
    p.iloc[:, i] = C.p_z2(p.iloc[:, i-1], p.columns[i-1], p.columns[i],
          Tv.iloc[:, i-1])

# potential Temperature
theta = C.Theta(T, p)

# Gradient of potential Temperature
Gradient = pd.DataFrame(index=T.index, columns=T.columns)
for i in np.arange(0, len(theta.columns)-1):
    Gradient.iloc[:, i] = (theta.iloc[:, i+1]- theta.iloc[:, i])/(theta.columns[i+1]-theta.columns[i])


def Plot_theta(index):
    plt.plot(theta.loc[index], theta.columns)
    plt.xlabel('Potential Temperature in K')
    plt.ylabel('Height in m')

def Plot_Gradient(index):
    plt.plot(Gradient.loc[index], Gradient.columns)
    plt.xlabel('Gradient of Potential Temperature in K')
    plt.ylabel('Height in m')

Plot_Gradient('2019-05-23 22:00:00')