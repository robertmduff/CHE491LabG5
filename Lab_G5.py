import numpy as np, matplotlib.pyplot as plt, pandas as pd, scipy.interpolate as interp, scipy.optimize as opt

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
        n = interp.interp1d([298, 308, 318, 328], [1.84, 1.1, 1.09, 1.01], kind='cubic')
        return (qm * b * P) / (1 + (b * P)**(1/n))**(n)
    elif input == "BET":
        qm = interp.interp1d([298, 308, 318, 328], [4.25, 4.05, 3.75, 3.35], kind='cubic')
        b = interp.interp1d([298, 308, 318, 328], [4.58E-5, 3.01E-5, 2.57E-5, 1.35E-5], kind='cubic')
        Psat = 10**(6.81228 - (1301.679 / (T - 3.494)))     # Constants A, B, and C from: https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI&Mask=4#Thermo-Phase
        return (qm * b * P) / ((Psat - P) * (1 + (b - 1) * P / Psat))