import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from itertools import cycle


def plot_absolute_values(result, N, m):
    days = len(result)
    t = np.linspace(0, days, days+1)
    
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
    fig1.add_trace(go.Scatter(x=t, y=Rv*n_scale, mode='lines', name='Rn', line_color=next(palette)))

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
    fig2.add_trace(go.Scatter(x=t, y=Rmv*m_scale, mode='lines', name='Rmd', line_color=next(palette)))

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
    t = np.linspace(0, days, days+1)
    diff_result = np.diff(result, axis=0)  # calculate the differences

    Id, In, Imn, Imd, I1, I2, Im1, Im2 = [diff_result[:, i] for i in range(10, 18)]
    susceptible_not_vaccinated = Imn + Imd
    susceptible_vaccinated = Im1 + Im2
    normal_not_vaccinated = In + Id
    normal_vaccinated = I1 + I2

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=t[:-1], y=susceptible_not_vaccinated, mode='lines', name='More susceptible not vaccinated'))
    fig.add_trace(go.Scatter(x=t[:-1], y=susceptible_vaccinated, mode='lines', name='More susceptible vaccinated'))
    fig.add_trace(go.Scatter(x=t[:-1], y=normal_not_vaccinated, mode='lines', name='Standard not vaccinated'))
    fig.add_trace(go.Scatter(x=t[:-1], y=normal_vaccinated, mode='lines', name='Standard vaccinated'))
    fig.update_layout(title='VAP-SIRS Model Simulation (Differences)',
                      xaxis_title='Time (days)',
                      yaxis_title='Change in infected population',
                      hovermode="x",
                      autosize=False, width=1000, height=600,
                      )

    return fig


def plot_infected_seasons(result, beta_values):
    days = len(result)
    interval_length = days // len(beta_values) * 2
    t = np.linspace(0, days, days+1)
    
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
    fig1.add_trace(go.Scatter(x=t, y=normal_not_vaccinated, mode='lines', name='Standard not vaccinated', line_color=next(palette)))
    fig1.add_trace(go.Scatter(x=t, y=normal_vaccinated, mode='lines', name='Standard vaccinated', line_color=next(palette)))

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
    t = np.linspace(0, days, days+1)
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


def plot_from_file(filepath):
    df = pd.read_csv(filepath, sep='\t')
    days = len(df)
    t = np.linspace(0, days, days+1)
    
    # Summing all compartments
    S_total = df[['Sd', 'Sn', 'Smn', 'Smd', 'S1', 'S2', 'Sm1', 'Sm2']].sum(axis=1)
    I_total = df[['Id', 'In', 'Imn', 'Imd', 'I1', 'I2', 'Im1', 'Im2']].sum(axis=1)
    R_total = df[['Rd', 'Rn', 'Rmn', 'Rmd', 'Rv', 'Rmv']].sum(axis=1)
    V_total = df[['V', 'Vm']].sum(axis=1)


    # Creating the plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=t, y=S_total, mode='lines', name='S'))
    fig.add_trace(go.Scatter(x=t, y=I_total, mode='lines', name='I'))
    fig.add_trace(go.Scatter(x=t, y=R_total, mode='lines', name='R'))
    fig.add_trace(go.Scatter(x=t, y=V_total, mode='lines', name='V'))

    fig.update_layout(title='Synthetic data',
                      xaxis_title='Time (days)',
                      yaxis_title='Total',
                      hovermode="x",
                      autosize=False, width=1000, height=600,)
    return fig


def plot_all(result, label):

    days = len(result)
    t = np.linspace(0, days, days) 

    palette = cycle(px.colors.qualitative.Pastel)
    fig = go.Figure()

    Sd, Sn, Smn, Smd, S1, S2, Sm1, Sm2, V, Vm, Id, In, Imn, Imd, I1, I2, Im1, Im2, Rd, Rn, Rmn, Rmd, Rv, Rmv = [result[:, i] for i in range(len(result[0]))]

    # Add traces
    fig.add_trace(go.Scatter(x=t, y=Sd, mode='lines', name='Sd', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Sn, mode='lines', name='Sn', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=S1, mode='lines', name='S1', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=S2, mode='lines', name='S2', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=V, mode='lines', name='V', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Id, mode='lines', name='Id', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=In, mode='lines', name='In', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=I1, mode='lines', name='I1', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=I2, mode='lines', name='I2', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Rd, mode='lines', name='Rd', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Rn, mode='lines', name='Rn', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Rv, mode='lines', name='Rn', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Smn, mode='lines', name='Smn', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Smd, mode='lines', name='Smd', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Sm1, mode='lines', name='Sm1', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Sm2, mode='lines', name='Sm2', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Vm, mode='lines', name='Vm', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Imn, mode='lines', name='Imn', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Imd, mode='lines', name='Imd', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Im1, mode='lines', name='Im1', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Im2, mode='lines', name='Im2', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Rmn, mode='lines', name='Rmn', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Rmd, mode='lines', name='Rmd', line_color=next(palette)))
    fig.add_trace(go.Scatter(x=t, y=Rmv, mode='lines', name='Rn', line_color=next(palette)))

    # Set layout properties
    fig.update_layout(
        title=label,
        xaxis_title='Time (days)',
        yaxis_title='Population',
        legend_title='Variables',
        hovermode="x",
        autosize=False, width=1000, height=600,
    )
    return fig