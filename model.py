from scipy.integrate import odeint
import plotly.graph_objects as go
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from itertools import cycle


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
    interval_length = days // len(beta_values) # długość sezonu

    for i in range(len(beta_values)):
        new_result = simulate_vap_sirs_model(results[-1], beta_0, beta_m0, f, f_v, kappa, upsilon, upsilon_r, upsilon_m, upsilon_mr, omega, omega_m, a, a_m, gamma, und_inf, interval_length)
        beta_0, beta_m0 = beta_values[i]  # zmieniamy bety
        results = np.concatenate((results, new_result), axis=0)  # dodajemy wynik symulacji do listy

    return results


def plot_absolute_values(result, N, m):
    # Assuming 'result' is a numpy array
    days = len(result)
    t = np.linspace(0, days, days)
    
    Sd, Sn, Smn, Smd, S1, S2, Sm1, Sm2, V, Vm, Id, In, Imn, Imd, I1, I2, Im1, Im2, Rd, Rn, Rmn, Rmd, Rv, Rmv = [result[:, i] for i in range(len(result[0]))]

    n_scale = 100 / (N * (1-m))
    m_scale = 100 / (N * m)
    # Create a plotly graph object
    palette = cycle(px.colors.qualitative.Pastel)
    fig1 = go.Figure()

    # Add traces
    fig1.add_trace(go.Scatter(x=t, y=Sd*n_scale, mode='lines', name='Sd', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=Sn*n_scale, mode='lines', name='Sn', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=S1*n_scale, mode='lines', name='S1', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=S2*n_scale, mode='lines', name='S2', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=V*n_scale, mode='lines', name='V', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=Id*n_scale, mode='lines', name='Id', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=In*n_scale, mode='lines', name='In', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=I1*n_scale, mode='lines', name='I1', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=I2*n_scale, mode='lines', name='I2', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=Rd*n_scale, mode='lines', name='Rd', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=Rn*n_scale, mode='lines', name='Rn', line_color=next(palette)))

    # Set layout properties
    fig1.update_layout(
        title='VAP-SIRS Model Simulation for normal people',
        xaxis_title='Time (days)',
        yaxis_title='Population',
        legend_title='Variables',
        hovermode="x",
        autosize=False, width=1000, height=600,
    )


    # Variables with 'm'
    palette = cycle(px.colors.qualitative.Pastel)
    fig2 = go.Figure()

    # Add traces
    fig2.add_trace(go.Scatter(x=t, y=Smn*m_scale, mode='lines', name='Smn', line_color=next(palette)))
    fig2.add_trace(go.Scatter(x=t, y=Smd*m_scale, mode='lines', name='Smd', line_color=next(palette)))
    fig2.add_trace(go.Scatter(x=t, y=Sm1*m_scale, mode='lines', name='Sm1', line_color=next(palette)))
    fig2.add_trace(go.Scatter(x=t, y=Sm2*m_scale, mode='lines', name='Sm2', line_color=next(palette)))
    fig2.add_trace(go.Scatter(x=t, y=Vm*m_scale, mode='lines', name='Vm', line_color=next(palette)))
    fig2.add_trace(go.Scatter(x=t, y=Imn*m_scale, mode='lines', name='Imn', line_color=next(palette)))
    fig2.add_trace(go.Scatter(x=t, y=Imd*m_scale, mode='lines', name='Imd', line_color=next(palette)))
    fig2.add_trace(go.Scatter(x=t, y=Im1*m_scale, mode='lines', name='Im1', line_color=next(palette)))
    fig2.add_trace(go.Scatter(x=t, y=Im2*m_scale, mode='lines', name='Im2', line_color=next(palette)))
    fig2.add_trace(go.Scatter(x=t, y=Rmn*m_scale, mode='lines', name='Rmn', line_color=next(palette)))
    fig2.add_trace(go.Scatter(x=t, y=Rmd*m_scale, mode='lines', name='Rmd', line_color=next(palette)))

    # Set layout properties
    fig2.update_layout(
        title='VAP-SIRS Model Simulation for more susceptible people',
        xaxis_title='Time (days)',
        yaxis_title='Population',
        legend_title='Variables',
        hovermode="x",
        autosize=False, width=1000, height=600,
    )
    return fig1, fig2


def plot_changes_in_infected(result):
    days = len(result)
    t = np.linspace(0, days, days)
    diff_result = np.diff(result, axis=0)  # calculate the differences

    Id, In, Imn, Imd, I1, I2, Im1, Im2 = [diff_result[:, i] for i in range(10, 18)]
    susceptible_not_vaccinated = Imn + Imd
    susceptible_vaccinated = Im1 + Im2
    normal_not_vaccinated = In + Id
    normal_vaccinated = I1 + I2

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=t[:-1], y=susceptible_not_vaccinated, mode='lines', name='More susceptible not vaccinated'))
    fig.add_trace(go.Scatter(x=t[:-1], y=susceptible_vaccinated, mode='lines', name='More susceptible vaccinated'))
    fig.add_trace(go.Scatter(x=t[:-1], y=normal_not_vaccinated, mode='lines', name='Normal not vaccinated'))
    fig.add_trace(go.Scatter(x=t[:-1], y=normal_vaccinated, mode='lines', name='Normal vaccinated'))
    fig.update_layout(title='VAP-SIRS Model Simulation (Differences)',
                      xaxis_title='Time (days)',
                      yaxis_title='Change in infected population',
                      hovermode="x",
                      autosize=False, width=1000, height=600,
                      )

    return fig


def plot_infected_seasons(result, beta_values):
    days = len(result)
    interval_length = days // len(beta_values)
    t = np.linspace(0, days, days)
    
    Sd, Sn, Smn, Smd, S1, S2, Sm1, Sm2, V, Vm, Id, In, Imn, Imd, I1, I2, Im1, Im2, Rd, Rn, Rmn, Rmd, Rv, Rmv = [result[:, i] for i in range(len(result[0]))]

    Id, In, Imn, Imd, I1, I2, Im1, Im2 = [result[:, i] for i in range(10, 18)]
    susceptible_not_vaccinated = Imn + Imd
    susceptible_vaccinated = Im1 + Im2
    normal_not_vaccinated = In + Id
    normal_vaccinated = I1 + I2
    
    
    palette = cycle(px.colors.qualitative.Pastel)
    fig1 = go.Figure()

    fig1.add_trace(go.Scatter(x=t, y=susceptible_not_vaccinated, mode='lines', name='More susceptible not vaccinated', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=susceptible_vaccinated, mode='lines', name='More susceptible vaccinated', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=normal_not_vaccinated, mode='lines', name='Normal not vaccinated', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=normal_vaccinated, mode='lines', name='Normal vaccinated', line_color=next(palette)))

    # Kolory pór roku
    season_colors = ['green', 'yellow', 'orange', 'blue']

    # Pionowe linie i zmiana koloru tła
    for i in range(8):
        start = i * interval_length
        end = start + interval_length
        color = season_colors[i % 4]
        fig1.add_shape(type="rect", xref="x", yref="paper", x0=start, y0=0, x1=end, y1=1, fillcolor=color, opacity=0.2, layer="below", line_width=0)
        fig1.add_shape(type="line", xref="x", yref="paper", x0=start, y0=0, x1=start, y1=1, line=dict(color="Black", width=1))

    fig1.update_layout(title='Number of infected at different seasons of the year',
                      xaxis_title='Time (days)',
                      yaxis_title='Total',
                      hovermode="x",
                      autosize=False, width=1000, height=600,)    

    return fig1



def plot_compartments(result):
    days = len(result)
    t = np.linspace(0, days, days)
    Sd, Sn, Smn, Smd, S1, S2, Sm1, Sm2, V, Vm, Id, In, Imn, Imd, I1, I2, Im1, Im2, Rd, Rn, Rmn, Rmd, Rv, Rmv = [result[:, i] for i in range(len(result[0]))]

    # Summing all compartments
    S_total = np.sum([Sd, Sn, Smn, Smd, S1, S2, Sm1, Sm2], axis=0)
    I_total = np.sum([Id, In, Imn, Imd, I1, I2, Im1, Im2], axis=0)
    R_total = np.sum([Rd, Rn, Rmn, Rmd, Rv, Rmv], axis=0)
    V_total = np.sum([V, Vm], axis=0)

    # Creating the plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=t, y=S_total, mode='lines', name='S'))
    fig.add_trace(go.Scatter(x=t, y=I_total, mode='lines', name='I'))
    fig.add_trace(go.Scatter(x=t, y=R_total, mode='lines', name='R'))
    fig.add_trace(go.Scatter(x=t, y=V_total, mode='lines', name='V'))

    fig.update_layout(title='Sum of Compartments Over Time',
                      xaxis_title='Time (days)',
                      yaxis_title='Total',
                      hovermode="x",
                      autosize=False, width=1000, height=600,)
    return fig
