from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np


def vap_sirs_model(y, _, beta, beta_v, kappa, upsilon, upsilon_r, omega, a, gamma):
    S_D, S_N, S_1, S_2, V, I_D, I_N, I_1, I_2, R_V, R_D, R_N = y
    I = I_D + I_N
    I_V = I_1 + I_2

    # Susceptible
    dS_D = - (beta * I + beta * I_V) * S_D + kappa * R_D
    dS_N = - (beta * I + beta * I_V) * S_N - upsilon * S_N + kappa * R_N
    dS_1 = upsilon_r * (1 - a) * S_2 + upsilon * (1 - a) * S_N - omega * S_1 - (beta * I + beta_v * I_V) * S_1
    dS_2 = - upsilon_r * S_2 + omega * V + omega * S_1 - (beta * I + beta_v * I_V) * S_2 + kappa * R_V

    # Vaccinated
    dV = upsilon * a * S_N + upsilon_r * a * S_2 - omega * V + upsilon_r * R_V + upsilon * R_N

    # Infected
    dI_D = (beta * I + beta * I_V) * S_D - gamma * I_D
    dI_N = (beta * I + beta * I_V) * S_N - gamma * I_N
    dI_1 = (beta * I + beta_v * I_V) * S_1 - gamma * I_1
    dI_2 = (beta * I + beta_v * I_V) * S_2 - gamma * I_2

    # Recovered
    dR_V = gamma * I_V - kappa * R_V - upsilon_r * R_V
    dR_D = gamma * I_D - kappa * R_D
    dR_N = gamma * I_N - kappa * R_N - upsilon * R_N

    return [dS_D, dS_N, dS_1, dS_2, dV, dI_D, dI_N, dI_1, dI_2, dR_V, dR_D, dR_N]


def simulate_vap_sirs_model(initial_conditions, beta, beta_v, kappa, upsilon, upsilon_r, omega, a, gamma, days):
    t = np.linspace(0, days, days)  # time points
    y0 = initial_conditions
    return odeint(vap_sirs_model, y0, t, args=(beta, beta_v, kappa, upsilon, upsilon_r, omega, a, gamma))


def plot_absolute_values(result):
    days = len(result)
    t = np.linspace(0, days, days)
    S_D, S_N, S_1, S_2, V, I_D, I_N, I_1, I_2, R_V, R_D, R_N = [result[:, i] for i in range(len(result[0]))]
    
    fig, axes = plt.subplots()
    plt.plot(t, S_D, label='S_D')
    plt.plot(t, S_N, label='S_N')
    plt.plot(t, S_1, label='S_1')
    plt.plot(t, S_2, label='S_2')
    plt.plot(t, V, label='V')
    plt.plot(t, I_D, label='I_D')
    plt.plot(t, I_N, label='I_N')
    plt.plot(t, I_1, label='I_1')
    plt.plot(t, I_2, label='I_2')
    plt.plot(t, R_V, label='R_V')
    plt.plot(t, R_D, label='R_D')
    plt.plot(t, R_N, label='R_N')
    plt.xlabel('Time (days)')
    plt.ylabel('Population')
    plt.title('VAP-SIRS Model Simulation')
    plt.legend()
    plt.grid(True)
    return fig

def plot_changes_in_infected(result):
    days = len(result)
    t = np.linspace(0, days, days)
    diff_result = np.diff(result, axis=0)  # calculate the differences

    I_D, I_N, I_1, I_2 = [diff_result[:, i] for i in range(5, 9)]
    I_sigma = I_N + I_D + I_1 + I_2
    I = I_N + I_D

    fig, axes = plt.subplots()
    plt.plot(t[:-1], I, label='I')
    plt.plot(t[:-1], I_sigma, label='I_sigma')
    plt.plot(t[:-1], I_1, label='I_1')
    plt.plot(t[:-1], I_2, label='I_2')
    plt.xlabel('Time (days)')
    plt.ylabel('Change in population')
    plt.title('VAP-SIRS Model Simulation (Differences)')
    plt.legend()
    plt.grid(True)
    return fig
