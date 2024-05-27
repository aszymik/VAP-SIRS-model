import streamlit as st
from model import *


st.title('VAP-SIRS Model Simulation')

with st.sidebar:
    # Suwaki do wyboru parametrów
    beta = st.slider('Beta', min_value=0.0, max_value=1.0, value=0.23)
    beta_v = st.slider('Beta V', min_value=0.0, max_value=1.0, value=0.45)
    kappa = st.slider('Kappa', min_value=0.0, max_value=1.0, value=1/500)
    upsilon = st.slider('Upsilon', min_value=0.0, max_value=1.0, value=1/250)
    upsilon_r = st.slider('Upsilon R', min_value=0.0, max_value=1.0, value=1/250)
    omega = st.slider('Omega', min_value=0.0, max_value=1.0, value=1/500)
    a = st.slider('A', min_value=0.0, max_value=1.0, value=0.79)
    gamma = st.slider('Gamma', min_value=0.0, max_value=1.0, value=1/6)
    d = st.slider('D', min_value=0.0, max_value=1.0, value=0.12)
    days = st.number_input('Days', min_value=1, max_value=1000, value=750, step=1)

# Wartości początkowe
V = 0
I_1 = 0
I_2 = 0
I_D = d * 10e-6
I_N = (1-d) * 10e-6
S = 1 - I_D - I_N
S_D = d * S
S_N = (1-d) * S
S_1 = 0
S_2 = 0
R_V = 0
R_D = 0
R_N = 0

# Przycisk do uruchomienia symulacji
if st.button('Run simulation'):
    initial_conditions = [S_D, S_N, S_1, S_2, V, I_D, I_N, I_1, I_2, R_V, R_D, R_N]
    result = simulate_vap_sirs_model(initial_conditions, beta, beta_v, kappa, upsilon, upsilon_r, omega, a, gamma, days)

    absolute_values_plot = plot_absolute_values(result)
    st.markdown('### Group populations over time')
    st.pyplot(absolute_values_plot)

    changes_in_infected_plot = plot_changes_in_infected(result)
    st.markdown('### Changes in infected populations over time')
    st.pyplot(changes_in_infected_plot)
