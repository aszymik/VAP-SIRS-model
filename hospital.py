from model import *
import plotly.graph_objects as go
from json import load
import numpy as np

def hospital_stat():
    t = np.linspace(0, 150, 150)

    # b=? bm=0.2, f=0.95, fv=0.8, v=vr=0.1, vm=vmr=0.2, hid_case=0.55
    fig_beta = go.Figure() # 1,2,3
    fig_beta.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s1.json'), mode='lines', name='0.10'))
    fig_beta.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s2.json'), mode='lines', name='0.12'))
    fig_beta.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s3.json'), mode='lines', name='0.14'))

    fig_beta.update_layout(title='Impact of the Beta Coefficient on the Total Case Count',
                        xaxis_title='Time (days)',
                        yaxis_title='% of population',
                        legend_title_text='Beta',
                        hovermode="x",
                        autosize=False, width=1000, height=600,)

    # b=0.1 bm=0.2, f=?, fv=0.8, v=vr=0.1, vm=vmr=0.2, hid_case=0.55
    fig_f = go.Figure() # 1, 6, 7
    fig_f.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s6.json'), mode='lines', name='0.94'))
    fig_f.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s1.json'), mode='lines', name='0.95'))
    fig_f.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s7.json'), mode='lines', name='0.96'))

    fig_f.update_layout(title='Impact of Restrictions Enforced on Unvaccinated Individuals on the Total Case Count',
                        xaxis_title='Time (days)',
                        yaxis_title='% of population',
                        legend_title_text='f',
                        hovermode="x",
                        autosize=False, width=1000, height=600,)

    # b=0.1 bm=0.2, f=0.95, fv=?, v=vr=0.1, vm=vmr=0.2, hid_case=0.55
    fig_fv = go.Figure() # 1, 4, 5
    fig_fv.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s5.json'), mode='lines', name='0.50'))
    fig_fv.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s1.json'), mode='lines', name='0.8'))
    fig_fv.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s4.json'), mode='lines', name='0.9'))

    fig_fv.update_layout(title='Impact of Restrictions Imposed on Vaccinated Individuals on the Total Case Count',
                        xaxis_title='Time (days)',
                        yaxis_title='% of population',
                        legend_title_text='fv',
                        hovermode="x",
                        autosize=False, width=1000, height=600,)

    # b=0.1 bm=0.2, f=0.95, fv=0.8, v=vr=?, vm=vmr=?, hid_case=0.55
    fig_v = go.Figure() # 1, 8, 9
    fig_v.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s9.json'), mode='lines', name='0.005')) # \t0.01
    fig_v.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s1.json'), mode='lines', name='0.01')) # \t0.02
    fig_v.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s8.json'), mode='lines', name='0.02')) # \t0.03' 

    fig_v.update_layout(title='Impact of Vaccination Rate on Infections',
                        xaxis_title='Time (days)',
                        yaxis_title='% of population',
                        legend_title_text='v',
                        hovermode="x",
                        autosize=False, width=1000, height=600,)
    

    # b=0.1 bm=0.2, f=0.95, fv=0.8, v=vr=0.1, vm=vmr=0.2, hid_case=?
    fig_hid = go.Figure() # 1, 10, 11
    fig_hid.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s9.json'), mode='lines', name='45%'))
    fig_hid.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s1.json'), mode='lines', name='55%'))
    fig_hid.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s8.json'), mode='lines', name='65%'))

    fig_hid.update_layout(title='Impact of the Size of the Fraction of Undetected Cases on the Total Case Count',
                        xaxis_title='Time (days)',
                        yaxis_title='% of population',
                        legend_title_text='% of hidden cases',
                        hovermode="x",
                        autosize=False, width=1000, height=600,)
    return fig_beta, fig_f, fig_fv, fig_v, fig_hid


def beta_m_stat():
    t = np.linspace(0, 150, 150)

    # b=? bm=0.2, f=0.95, fv=0.8, v=vr=0.1, vm=vmr=0.2, hid_case=0.55
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s12.json'), mode='lines', name='0.2'))
    fig.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s13.json'), mode='lines', name='0.4'))
    fig.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s14.json'), mode='lines', name='0.6'))
    fig.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s15.json'), mode='lines', name='0.8'))
    fig.add_trace(go.Scatter(x=t, y=run_from_file('data/params_s16.json'), mode='lines', name='1.0'))

    fig.update_layout(title='Impact of the Beta Coefficient for More Susceptible on the Total Case Count',
                        xaxis_title='Time (days)',
                        yaxis_title='% of population',
                        legend_title_text='Beta',
                        hovermode="x",
                        autosize=False, width=1000, height=600,)

    return fig

N = 100
V, Vm = 0, 0

I_percn = 0.01

Sd = 0.12 * (1-0.1) * N
Sn = (1-0.12)* (1-0.1) * N
Smn = (1-0.12) * 0.1 * N
Smd = 0.12 * 0.1 * N
S1, S2, Sm1, Sm2 = 0, 0, 0, 0

Id = I_percn * Sd
In = I_percn * Sn
Imd = I_percn *Smd
Imn = I_percn * Smn

I1, I2, Im1, Im2 = 0, 0, 0, 0
Rd, Rn, Rmn, Rmd, Rv, Rmv = 0, 0, 0, 0, 0, 0

initial_conditions = [(1-I_percn)*Sd, (1-I_percn)*Sn, (1-I_percn)*Smn, (1-I_percn)*Smd, S1, S2, Sm1, Sm2, V, Vm, Id, In, Imn, Imd, I1, I2, Im1, Im2, Rd, Rn, Rmn, Rmd, Rv, Rmv]

def run_from_file(filename):
    with open(filename, 'r') as file:
        data = load(file)
    result = simulate_vap_sirs_model(initial_conditions, **data)
    I_total = np.sum([result[:, i] for i in range(10, 18)], axis=0)
    return I_total
