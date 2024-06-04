import streamlit as st
from model import *
import json


st.title('VAP-SIRS Model Simulation')

with st.sidebar:

    params = {
        'beta': st.slider('Beta', min_value=0.0, max_value=1.0, value=0.2),
        'beta_m': st.slider('Beta m', min_value=0.0, max_value=1.0, value=0.4),
        'f': st.slider('f', min_value=0.0, max_value=1.0, value=0.77),
        'fv': st.slider('fv', min_value=0.0, max_value=1.0, value=0.55),
        'upsilon': st.slider('Upsilon', min_value=0.0, max_value=1.0, value=0.004),
        'upsilon_r': st.slider('Upsilon R', min_value=0.0, max_value=1.0, value=0.003),
        'upsilon_m': st.slider('Upsilon m', min_value=0.0, max_value=1.0, value=0.006),
        'upsilon_mr': st.slider('Upsilon m R', min_value=0.0, max_value=1.0, value=0.006),
        'omega': 1/st.slider('Number of days after which vaccinated normal people lose their immunity', min_value=30, max_value=700, value=365),
        'omega_m': 1/st.slider('Number of days after which vaccinated more susceptible people lose their immunity', min_value=30, max_value=700, value=365),
        'kappa': 1/st.slider('Number of days after which recovered people lose their immunity', min_value=30, max_value=700, value=400),
        'a': st.slider('Vaccine effectiveness for normal people', min_value=0.0, max_value=1.0, value=0.79),
        'a_m': st.slider('Vaccine effectiveness for more susceptible people', min_value=0.0, max_value=1.0, value=0.60),
        'gamma': 1/st.slider('Duration of the disease', min_value=1, max_value=21, value=7),
        'und_inf': st.slider('% of undiagnosed infected people', min_value=0.0, max_value=1.0, value=1.0),
        'd': st.slider("% of people who won't vaccinate", min_value=0.0, max_value=1.0, value=0.12),
        'm': st.slider('% of more susceptible people', min_value=0.0, max_value=1.0, value=0.05),
        'days': st.number_input('Days', min_value=1, max_value=2000, value=150, step=1)
    }

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

Sd = params['d'] * (1-params['m']) * N
Sn = (1-params['d'])* (1-params['m']) * N
Smn = (1-params['d']) * params['m'] * N
Smd = params['d'] * params['m'] * N
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


# Przycisk do uruchomienia symulacji
if st.button('Run simulation'):
    initial_conditions = [(1-I_percn)*Sd, (1-I_percn)*Sn, (1-I_percn)*Smn, (1-I_percn)*Smd, S1, S2, Sm1, Sm2, V, Vm, Id, In, Imn, Imd, I1, I2, Im1, Im2, Rd, Rn, Rmn, Rmd, Rv, Rmv]
    # print(f'{sum(initial_conditions)=}')
    result = simulate_vap_sirs_model(initial_conditions, 
                                    params['beta'], 
                                    params['beta_m'], 
                                    params['f'], 
                                    params['fv'], 
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
                                    params['days'],
                                    params['und_inf']
                                    )

    total_plot = plot_compartments(result)
    st.markdown('### Sum of Compartments Over Time')
    st.plotly_chart(total_plot)

    abs_plot1, abs_plot2 = plot_absolute_values(result, N, params['m'])
    st.markdown('### Group populations over time')
    st.plotly_chart(abs_plot1)
    st.plotly_chart(abs_plot2)

    inf_change_plot = plot_changes_in_infected(result)
    st.markdown('### Changes in infected populations over time')
    st.plotly_chart(inf_change_plot)
    
