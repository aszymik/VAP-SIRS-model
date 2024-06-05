import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from model import modified_vap_sirs_model
from app import *

def generate_synthetic_data(initial_conditions, days, params, seasonal_amplitude, noise_std):
    t = np.linspace(0, days, days+1)
    y0 = initial_conditions
    result = odeint(modified_vap_sirs_model, y0, t, args=params)

    df = pd.DataFrame(result, columns=['Sd', 'Sn', 'Smn', 'Smd', 'S1', 'S2', 'Sm1', 'Sm2', 'V', 'Vm', 'Id', 'In', 'Imn', 'Imd', 'I1', 'I2', 'Im1', 'Im2', 'Rd', 'Rn', 'Rmn', 'Rmd', 'Rv', 'Rmv'])
    df['time'] = t
    df.to_csv('data/synthetic_data_without_noise.tsv', sep='\t', index=False)

    # Fluktuacje sezonowe
    # seasonal_variation = seasonal_amplitude * np.sin(2 * np.pi * t / 365)
    # result += seasonal_variation[:, np.newaxis]

    # Szum Gaussowski
    noise = np.random.normal(scale=noise_std, size=result.shape)
    result += noise

    # Zapisujemy do df
    df = pd.DataFrame(result, columns=['Sd', 'Sn', 'Smn', 'Smd', 'S1', 'S2', 'Sm1', 'Sm2', 'V', 'Vm', 'Id', 'In', 'Imn', 'Imd', 'I1', 'I2', 'Im1', 'Im2', 'Rd', 'Rn', 'Rmn', 'Rmd', 'Rv', 'Rmv'])
    df['time'] = t
    df.to_csv('data/synthetic_data.tsv', sep='\t', index=False)

    return result


# Wartości początkowe
initial_conditions = [(1-I_percn)*Sd, (1-I_percn)*Sn, (1-I_percn)*Smn, (1-I_percn)*Smd, S1, S2, Sm1, Sm2, V, Vm, Id, In, Imn, Imd, I1, I2, Im1, Im2, Rd, Rn, Rmn, Rmd, Rv, Rmv]
days = 730
# beta_0, beta_m0, f, f_v, kappa, upsilon, upsilon_r, upsilon_m, upsilon_mr, omega, omega_m, a, a_m, gamma, und_inf
params=(0.1, 0.2, 0.77, 0.55, 1/400, 0.005, 0.005, 0.01, 0.01, 1/365, 1/365, 0.95, 0.8, 7, 1.0)

if __name__ == '__main__':
    generate_synthetic_data(initial_conditions, days, params, seasonal_amplitude=0.2, noise_std=1)
