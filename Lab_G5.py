import numpy as np, matplotlib.pyplot as plt, pandas as pd, scipy.interpolate as interp
import scipy.optimize as opt, scipy.integrate as integrate, os
from gekko import GEKKO
from scipy.optimize import curve_fit

D_column = 80 / 1000  # Column diameter (Source: https://www.gunt.de/en/products/adsorptive-air-drying/083.54000/ce540/glct-1:pa-148:pr-53)
Z_column = 800 / 1000 # Column height (Source: https://www.gunt.de/en/products/adsorptive-air-drying/083.54000/ce540/glct-1:pa-148:pr-53)

dporekin_4x8 = 0.9E-09       # Kinetic pore diameter in m
dpore_4x8 = 1.0E-09          # Pore diameter in m
dbead_4x8 = 4.0 / 1000       # Bead diameter in m
rho_bulk_4x8 = 650           # Bulk density in kg/m3
As_4x8 = 950                 # Specific surface area in m2/g

dporekin_8x12 = 0.9E-09      # Kinetic pore diameter in m
dpore_8x12 = 1.0E-09         # Pore diameter in m
dbead_8x12 = 2.0 / 1000      # Bead diameter in m
rho_bulk_8x12 = 650          # Bulk density in kg/m3
As_8x12 = 950                # Specific surface area in m2/g

# Note Pequil = P0 - deltaP cause P drops during adsorption

def experiment():
    df_4x8 = pd.read_excel('G10_Lab_G5.xlsx', sheet_name='4x8')
    df_8x12 = pd.read_excel('G10_Lab_G5.xlsx', sheet_name='8x12')

    # t_4x8 = df_4x8['Time'].values
    # CbyCo_4x8 = df_4x8['C/Co'].values
    # Co_4x8 = 
    # m_4x8 = 

    # t_8x12 = df_8x12['Time'].values
    # CbyCo_8x12 = df_8x12['C/Co'].values
    # Co_8x12 = 
    # m_8x12 = 

    area_4x8 = integrate.simpson(1 - CbyCo_4x8, t_4x8)
    area_8x12 = integrate.simpson(1 - CbyCo_8x12, t_8x12)
    qstar_4x8 = area_4x8 * Co_4x8 / m_4x8
    qstar_8x12 = area_8x12 * Co_8x12 / m_8x12


def isotherms(input, P, T):
    if input == "Langmuir":
        qm = interp.interp1d([298, 308, 318, 328], [2.7, 1.94, 1.801, 0.789], kind='cubic')
        b = interp.interp1d([298, 308, 318, 328], [6.9E-5, 20.8E-5, 6.35E-5, 6.02E-5], kind='cubic')
        return (qm * b * P) / (1 + b * P)
    elif input == "freundlich":
        k = interp.interp1d([298, 308, 318, 328], [1.4, 1.73, 0.6, 0.55], kind='cubic')
        n = interp.interp1d([298, 308, 318, 328], [1.84, 1.1, 1.09, 1.01], kind='cubic')
    elif input == "langmuirfreundlich":
        qm = interp.interp1d([298, 308, 318, 328], [7.05, 6.85, 6.15, 5.9], kind='cubic')
        b = interp.interp1d([298, 308, 318, 328], [4.6E-5, 2.3E-5, 1.01E-5, 0.95E-5], kind='cubic')
        n = interp.interp1d([298, 308, 318, 328], [1.84, 1.1, 1.09, 1.01], kind='cubic')
        return (qm * (b * P)**(1/n)) / (1 + (b * P)**(1/n))
    elif input == "toth":
        qm = interp.interp1d([298, 308, 318, 328], [5.25, 4.75, 3.89, 3.35], kind='cubic')
        b = interp.interp1d([298, 308, 318, 328], [5.58E-5, 3.91E-5, 2.01E-5, 1.85E-5], kind='cubic')
        n = interp.interp1d([298, 308, 318, 328], [1.84, 1.1, 1.09, 1.01], kind='cubic')        # Double check this??
        return (qm * b * P) / (1 + (b * P)**(1/n))**(n)
    elif input == "BET":
        qm = interp.interp1d([298, 308, 318, 328], [4.25, 4.05, 3.75, 3.35], kind='cubic')
        b = interp.interp1d([298, 308, 318, 328], [4.58E-5, 3.01E-5, 2.57E-5, 1.35E-5], kind='cubic')
        Psat = 10**(6.81228 - (1301.679 / (T - 3.494)))     # Constants A, B, and C from: https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI&Mask=4#Thermo-Phase
        return (qm * b * P) / ((Psat - P) * (1 + (b - 1) * P / Psat))
    

##### Pressure Drop Correlations #####

U = G / (np.pi * (D_column / 2)**2)
deltaP = Z_column * (((150 * (1 - eta)**2 * mu) / (eta**3 * (Dp * phi_s)**2)) * U + ((1.75 * (1 - eta) * rhog) / (eta**3 * Dp * phi_s)) * U**2)

##### Mass Based Models #####

def bohartadam(t, y):
    kprime, qs = y
    return (np.exp(kprime * Co * (t - z_bed / U)) / (np.exp(kprime * Co * (t - z_bed / U)) + np.exp(kprime * (qs * z_bed / U) * ((1 - eta) / eta)) - 1))

def matbal_eqns():
    r = (k / Co) / (1 + k / Co)
    h = (C / Co) / (r + (1 + r) * (C / Co))
    beta = (qo_star / C_To) * ((1 + Yo) / Yo)
    Ni = (Kia * Li * beta) / Ui_star
    Mi = (Kia * tc) / rhoc
    return [Ni, Mi]


def solvingpdes():
    m = GEKKO()



    # first segment
    m.Equation(rho*A*L_seg*cp*T[0].dt() == \
                keff*A*(Th-T[0])/L_seg \
                - keff*A*(T[0]-T[1])/L_seg \
                - heff*As*(T[0]-Ts))
    # middle segments
    m.Equations([rho*A*L_seg*cp*T[i].dt() == \
                keff*A*(T[i-1]-T[i])/L_seg \
                - keff*A*(T[i]-T[i+1])/L_seg \
                - heff*As*(T[i]-Ts) for i in range(1,seg-1)])
    # last segment
    m.Equation(rho*A*L_seg*cp*T[seg-1].dt() == \
                keff*A*(T[seg-2]-T[seg-1])/L_seg \
            - heff*(As+A)*(T[seg-1]-Ts))


