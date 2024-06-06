import streamlit as st
import json
from model import *
from plots import *
from hospital import *

from fit import *


st.title('VAP-SIRS Model Simulation')

with st.sidebar:

    params = {
        'beta_0': st.slider('Beta', min_value=0.0, max_value=1.0, value=0.1),
        'beta_m0': st.slider('Beta m', min_value=0.0, max_value=1.0, value=0.2),
        'f': st.slider('f', min_value=0.0, max_value=1.0, value=0.77),
        'f_v': st.slider('fv', min_value=0.0, max_value=1.0, value=0.55),
        'kappa': 1/st.slider('Number of days after which recovered people lose their immunity', min_value=30, max_value=700, value=400),
        'upsilon': st.slider('Upsilon', min_value=0.0, max_value=1.0, value=0.005),
        'upsilon_r': st.slider('Upsilon R', min_value=0.0, max_value=1.0, value=0.005),
        'upsilon_m': st.slider('Upsilon m', min_value=0.0, max_value=1.0, value=0.01),
        'upsilon_mr': st.slider('Upsilon m R', min_value=0.0, max_value=1.0, value=0.01),
        'omega': 1/st.slider('Number of days after which vaccinated normal people lose their immunity', min_value=30, max_value=700, value=365),
        'omega_m': 1/st.slider('Number of days after which vaccinated more susceptible people lose their immunity', min_value=30, max_value=700, value=365),
        'a': st.slider('Vaccine effectiveness for normal people', min_value=0.0, max_value=1.0, value=0.95),
        'a_m': st.slider('Vaccine effectiveness for more susceptible people', min_value=0.0, max_value=1.0, value=0.80),
        'gamma': 1/st.slider('Duration of the disease', min_value=1, max_value=21, value=7),
        'und_inf': st.slider('% of undiagnosed infected people', min_value=0.0, max_value=1.0, value=0.55),
        'days': st.number_input('Days', min_value=1, max_value=2000, value=150, step=1)
    }

    d = st.slider("% of people who won't vaccinate", min_value=0.0, max_value=1.0, value=0.12)
    m = st.slider('% of more susceptible people', min_value=0.0, max_value=1.0, value=0.05)
    # Text input for filename
    filename = st.text_input('Enter filename', 'params.json')

    # When the "Save" button is pressed
    if st.button('Save'):
        # Write parameters to a JSON file
        with open(f'data/{filename}', 'w') as f:
            json.dump(params, f)
        st.success(f'Parameters saved to {filename}')


# Wartości początkowe
N = 100
V = 0
Vm = 0

I_percn = 0.01

Sd = d * (1-m) * N
Sn = (1-d)* (1-m) * N
Smn = (1-d) * m * N
Smd = d * m * N
S1 = 0
S2 = 0
Sm1 = 0
Sm2 = 0


Id = I_percn * Sd
In = I_percn * Sn
Imd = I_percn *Smd
Imn = I_percn * Smn
I1 = 0
I2 = 0
Im1 = 0
Im2 = 0

Rd = 0
Rn = 0
Rmn = 0
Rmd = 0
Rv = 0
Rmv = 0


# Paramtery do wariacji sezonowych (wiosna, lato, jesień, zima przez dwa lata)
years = 2
beta_values = [[0.02, 0.06], [0.02, 0.06],
               [0.01, 0.05], [0.02, 0.06], 
               [0.04, 0.2], [0.04, 0.2],
               [0.4, 0.6], [0.4, 0.6]] * years


# Przycisk do uruchomienia symulacji
if st.button('Run simulation'):
    initial_conditions = [(1-I_percn)*Sd, (1-I_percn)*Sn, (1-I_percn)*Smn, (1-I_percn)*Smd, S1, S2, Sm1, Sm2, V, Vm, Id, In, Imn, Imd, I1, I2, Im1, Im2, Rd, Rn, Rmn, Rmd, Rv, Rmv]
    # print(f'{sum(initial_conditions)=}')
    result = simulate_vap_sirs_model(initial_conditions, 
                                    params['beta_0'], 
                                    params['beta_m0'], 
                                    params['f'], 
                                    params['f_v'], 
                                    params['kappa'], 
                                    params['upsilon'], 
                                    params['upsilon_r'], 
                                    params['upsilon_m'], 
                                    params['upsilon_mr'], 
                                    params['omega'], 
                                    params['omega_m'], 
                                    params['a'], 
                                    params['a_m'], 
                                    params['gamma'], 
                                    params['und_inf'],
                                    params['days']
                                    )

    total_plot = plot_compartments(result)
    st.markdown('## Sum of Compartments Over Time')
    st.plotly_chart(total_plot)

    abs_plot1, abs_plot2 = plot_absolute_values(result, N, m)
    st.markdown('## Group populations over time')
    st.plotly_chart(abs_plot1)
    st.plotly_chart(abs_plot2)

    inf_change_plot = plot_changes_in_infected(result)
    st.markdown('## Changes in infected populations over time')
    st.plotly_chart(inf_change_plot)

    result_s = run_model_with_seasonal_variations(initial_conditions,
                                                   params['beta_0'], 
                                                   params['beta_m0'], 
                                                   params['f'], 
                                                   params['f_v'], 
                                                   params['kappa'], 
                                                   params['upsilon'], 
                                                   params['upsilon_r'], 
                                                   params['upsilon_m'], 
                                                   params['upsilon_mr'], 
                                                   params['omega'], 
                                                   params['omega_m'], 
                                                   params['a'], 
                                                   params['a_m'], 
                                                   params['gamma'], 
                                                   params['und_inf'],
                                                   730,
                                                   beta_values
                                                   )
    
    inf_seasons_plot = plot_infected_seasons(result_s, beta_values)
    st.markdown('## Seasonal variation')
    st.plotly_chart(inf_seasons_plot)

    st.markdown('## Hospital capacity')
    fig_beta, fig_f, fig_fv, fig_v, fig_hid = hospital_stat()
    st.markdown('#### Transmition rate')
    st.plotly_chart(fig_beta)
    st.markdown('#### Restrictions')
    st.plotly_chart(fig_f)
    st.plotly_chart(fig_fv)
    st.markdown('#### Vaccination rate')
    st.plotly_chart(fig_v)
    st.markdown('#### Fraction of hidden infection cases')
    st.plotly_chart(fig_hid)

    st.markdown('## Transmition rate for More Susceptible')

    plot_ims, plot_im, plot_i = v_m_stat()
    st.plotly_chart(plot_im)
    st.plotly_chart(plot_ims)
    st.plotly_chart(plot_i)

    st.markdown('### Influence on Infected number')
    fig_beta_m = beta_m_stat()
    st.plotly_chart(fig_beta_m)

    st.markdown('## Real-world data simulation')
    st.markdown('### Data based on Hospital scenario')
    plot_syn, plot_fit, params_plot = fitted_scenario('data/params_s1.json', initial_conditions)
    st.plotly_chart(plot_syn)
    st.plotly_chart(plot_fit)
    st.plotly_chart(params_plot)

    st.markdown('### Data based on seasonalvariations scenario')
    seasonal_noise = add_noise(result_s, seasonal_amplitude=0.2, noise_std=1)
    plot_syn, plot_fit, params_plot = fit_seasonal(seasonal_noise, params, initial_conditions)
    st.plotly_chart(plot_syn)
    st.plotly_chart(plot_fit)
    st.plotly_chart(params_plot)




    
