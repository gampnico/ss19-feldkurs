# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:15:19 2019

@author: josef
"""

import numpy as np

# Vertical scanning levels of Hatpro
Z = np.array([0, 10, 30, 50, 75, 100, 125, 150, 200, 250, 325, 400, 475, 550,
              625, 700, 800, 900, 1000, 1150, 1300, 1450, 1600, 1800, 2000,
              2200, 2500, 2800, 3100, 3500, 3900, 4400, 5000, 5600, 6200, 7000,
              8000, 9000, 10000])

# constants
R_d = 287  # J K^-1 kg^-1, gas constant of dry air
R_v = 461  # J K^-1 kg^-1, gas constant of water vapor
g = 9.81  # m s^-2 gravitational constant
p_0 = 1000*10**2  # Pa reference pressure
c_p = 1005  # J kg^-1 K^-1, specific heat capacity of air


def e(a, T):
    ''' water vapor pressure.

    Parameters:
    ---------
    a: absolute humidity
    T: Temperature
    '''
    return a*R_v*T


def w(e, p_d):
    ''' mixing ratio.

    Parameters:
    ---------
    e: water vapor pressure
    p_d: pressure of dry air
    '''
    return (R_d*e)/(R_v*p_d)


def Tv(T, w):
    ''' virtual Temperature.
    Parameters:
    ---------
    T: Temperature
    w: mixing ratio '''

    return T*(1+0.61*w)


def p_z2(p_z1, z1, z2, Tv):
    ''' Pressure at height z2.
    Parameters:
    ---------
    p_z1: pressure at height z1
    z1: height one
    z2: height two
    Tv: virtual Temperature '''

    return p_z1*np.exp(-g*(z2-z1)/(R_d*Tv))


def Theta(T, p):
    ''' potential Temperature.
    Parameters:
    ---------
    T: Temperature
    p: pressure '''

    return T*(p_0/p)**(R_d/c_p)
