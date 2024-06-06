from scipy.integrate import odeint
import numpy as np
from scipy.optimize import minimize


def modified_vap_sirs_model(y, _, beta_0, beta_m0, f, f_v, kappa, upsilon, upsilon_r, upsilon_m, upsilon_mr, omega, omega_m, a, a_m, gamma, und_inf):
    # S1 myślą ze są odporni
    # S2 wiedzą ze stracili odporność

    # Sn planują się zaszczepić
    # Sd nie planują się zaszczepić
    
    Sd, Sn, Smn, Smd, S1, S2, Sm1, Sm2, V, Vm, Id, In, Imn, Imd, I1, I2, Im1, Im2, Rd, Rn, Rmn, Rmd, Rv, Rmv = y

    I = Id + In 
    Im = Imn + Imd
    Iv = I1 + I2
    Imv = Im1 + Im2

    beta = beta_0 * (1-f)
    beta_v = beta_0 * (1-f_v)
    beta_m = beta_m0 * (1-f)
    beta_mv = beta_m0 * (1-f_v)

    # Susceptible
    dSd = - beta * und_inf * (I + Im + Iv + Imv) * Sd + kappa * Rd
    dSn = - beta * und_inf * (I + Im + Iv + Imv) * Sn - upsilon * Sn + kappa * Rn
    dS1 = upsilon_r * (1 - a) * S2 + upsilon * (1 - a) * Sn + upsilon *(1-a)*Rn + upsilon_r*(1-a)*Rv - omega * S1 - (beta * (I + Im) + beta_v * (Iv + Imv)) * und_inf * S1
    dS2 = - upsilon_r * S2 + omega * V + omega * S1 - (beta * (I + Im) + beta_v * (Iv + Imv))* und_inf * S2 + kappa * Rv


    # More susceptible
    dSmd = - beta_m * und_inf *(I + Im + Iv + Imv) * Smd + kappa * Rmd
    dSmn = - beta_m * und_inf *(I + Im + Iv + Imv) * Smn - upsilon_m * Smn + kappa * Rmn

    dSm1 = upsilon_mr * (1 - a_m) * Sm2 + upsilon_m * (1 - a_m) * Smn + upsilon_m *(1-a_m)*Rmn + upsilon_mr*(1-a_m)*Rmv - omega_m * Sm1 - (beta_m * (I + Im) + beta_mv * (Iv+Imv)) * und_inf * Sm1
    dSm2 = - upsilon_mr * Sm2 + omega_m * Vm + omega_m * Sm1 - (beta_m * (I + Im) + beta_mv * (Iv+Imv)) * und_inf * Sm2 + kappa * Rmv

    
    # Vaccinated
    dV = upsilon * a * Sn + upsilon_r * a * S2  - omega * V + upsilon_r * a * Rv  + upsilon * a * Rn 
    dVm = upsilon_m * a_m * Smn + upsilon_mr * a_m * Sm2 - omega_m * Vm + upsilon_mr* a_m * Rmv + upsilon_m* a_m * Rmn

    # Infected
    dId = beta * (I + Im + Iv + Imv) * und_inf * Sd - gamma * Id
    dIn = beta * (I + Im + Iv + Imv) * und_inf * Sn - gamma * In
    dImd = beta_m *(I + Im + Iv + Imv) * und_inf * Smd - gamma * Imd
    dImn = beta_m *(I + Im + Iv + Imv) * und_inf * Smn - gamma * Imn

    dI1 = (beta * (I + Im) + beta_v * (Iv + Imv)) * und_inf * S1 - gamma * I1
    dI2 = (beta * (I + Im) + beta_v * (Iv + Imv)) * und_inf * S2 - gamma * I2
    dIm1 = (beta_m * (I + Im) + beta_mv * (Iv+Imv)) * und_inf * Sm1 - gamma * Im1
    dIm2 = (beta_m * (I + Im) + beta_mv * (Iv+Imv)) * und_inf * Sm2 - gamma * Im2

    # Recovered
    dRv = gamma * (I1 + I2) - kappa * Rv - upsilon_r * Rv
    dRmv = gamma * (Im1 + Im2) - kappa * Rmv - upsilon_mr * Rmv

    dRd = gamma * Id - kappa * Rd
    dRn = gamma * In - kappa * Rn - upsilon * Rn
    dRmd = gamma * Imd - kappa * Rmd
    dRmn = gamma * Imn - kappa * Rmn - upsilon_m * Rmn

    
    return [dSd, dSn, dSmn, dSmd, dS1, dS2, dSm1, dSm2, dV, dVm, dId, dIn, dImn, dImd, dI1, dI2, dIm1, dIm2, dRd, dRn, dRmn, dRmd, dRv, dRmv]


def simulate_vap_sirs_model(initial_conditions, beta_0, beta_m0, f, f_v, kappa, upsilon, upsilon_r, upsilon_m, upsilon_mr, omega, omega_m, a, a_m, gamma, und_inf, days):
    t = np.linspace(0, days, days)  # time points
    y0 = initial_conditions
    return odeint(modified_vap_sirs_model, y0, t, args=(beta_0, beta_m0, f, f_v, kappa, upsilon, upsilon_r, upsilon_m, upsilon_mr, omega, omega_m, a, a_m, gamma, und_inf))


def run_model_with_seasonal_variations(initial_conditions, beta_0, beta_m0, f, f_v, kappa, upsilon, upsilon_r, upsilon_m, upsilon_mr, omega, omega_m, a, a_m, gamma, und_inf, days, beta_values):
   
    results = np.array([initial_conditions])  # macierz do przechowywania wyników
    interval_length = days // len(beta_values) * 2 # długość sezonu

    for i in range(len(beta_values)):
        new_result = simulate_vap_sirs_model(results[-1], beta_0, beta_m0, f, f_v, kappa, upsilon, upsilon_r, upsilon_m, upsilon_mr, omega, omega_m, a, a_m, gamma, und_inf, interval_length)
        beta_0, beta_m0 = beta_values[i]  # zmieniamy bety
        results = np.concatenate((results, new_result), axis=0)  # dodajemy wynik symulacji do listy

    return results


def generate_synthetic_data(initial_conditions, beta_0, beta_m0, f, f_v, kappa, upsilon, upsilon_r, upsilon_m, upsilon_mr, omega, omega_m, a, a_m, gamma, und_inf, days, seasonal_amplitude, noise_std):
    t = np.linspace(0, days, days)
    y0 = initial_conditions
    result = odeint(modified_vap_sirs_model, y0, t, args=(beta_0, beta_m0, f, f_v, kappa, upsilon, upsilon_r, upsilon_m, upsilon_mr, omega, omega_m, a, a_m, gamma, und_inf))

    # Fluktuacje sezonowe
    seasonal_variation = seasonal_amplitude * np.sin(2 * np.pi * t / 365)
    result += seasonal_variation[:, np.newaxis]

    # Szum Gaussowski
    noise = np.random.normal(scale=noise_std, size=result.shape)
    result += noise

    return result

def fit_to_real_data(initial_condition, noised_data, init_params):

    days = len(noised_data)
    t = np.linspace(0, days, days) 
    bounds = [(0,1),(0,1),(0,1),(0,1), (1/700, 1/30), (0,1),(0,1),(0,1),(0,1), (1/700, 1/30), (1/700, 1/30),(0,1),(0,1), (1/21, 1), (0,1)]
    def mse(params):
        pred = odeint(modified_vap_sirs_model, initial_condition, t, args=tuple(params))
        return np.sum(np.square(np.subtract(noised_data, pred)), axis=1).mean()
    
    # Minimalizacja funkcji straty
    fitted = minimize(mse, init_params, bounds=bounds)
    # Wyniki
    optimal_params = tuple(fitted.x)

    pred_opt = odeint(modified_vap_sirs_model, initial_condition, t, args=optimal_params)

    return pred_opt, optimal_params


def add_noise(data, seasonal_amplitude, noise_std):
    days = len(data)
    t = np.linspace(0, days, days)
     # Fluktuacje sezonowe
    seasonal_variation = seasonal_amplitude * np.sin(2 * np.pi * t / 365)
    data += seasonal_variation[:, np.newaxis]

    # Szum Gaussowski
    noise = np.random.normal(scale=noise_std, size=data.shape)
    data += noise

    return data
