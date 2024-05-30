import streamlit as st
from model import *


st.title('VAP-SIRS Model Simulation')

with st.sidebar:
    # Suwaki do wyboru parametrów
    beta_0 = st.slider('Beta', min_value=0.0, max_value=5.0, value=1.0)
    beta_m0 = st.slider('Beta m', min_value=0.0, max_value=5.0, value=2.0)

    f = st.slider('f', min_value=0.0, max_value=1.0, value=0.77)
    f_v = st.slider('fv', min_value=0.0, max_value=1.0, value=0.55)

    kappa = st.slider('Kappa', min_value=0.0, max_value=1.0, value=1/500)

    upsilon = st.slider('Upsilon', min_value=0.0, max_value=1.0, value=1/250)
    upsilon_r = st.slider('Upsilon R', min_value=0.0, max_value=1.0, value=1/250)
    upsilon_m = st.slider('Upsilon m', min_value=0.0, max_value=1.0, value=2/250)
    upsilon_mr = st.slider('Upsilon m R', min_value=0.0, max_value=1.0, value=2/250)
    
    omega = st.slider('Omega', min_value=0.0, max_value=1.0, value=1/500)
    omega_m = st.slider('Omega m', min_value=0.0, max_value=1.0, value=1/400)

    a = st.slider('A', min_value=0.0, max_value=1.0, value=0.79)
    a_m = st.slider('Am', min_value=0.0, max_value=1.0, value=0.60)

    gamma = st.slider('Gamma', min_value=0.0, max_value=1.0, value=1/6)
    d = st.slider("% of people who won't vaccinate", min_value=0.0, max_value=1.0, value=0.12)
    m = st.slider('% of more susceptible people', min_value=0.0, max_value=1.0, value=0.12)
    days = st.number_input('Days', min_value=1, max_value=1000, value=750, step=1)

    

# Wartości początkowe
V = 0

Id = d * 10e-6 * (1-m)
In = (1-d) * 10e-6 * (1-m)
Imd = d * 10e-6 * m
Imn = (1-d) * 10e-6 * m
I1 = 0
I2 = 0
Im1 = 0
Im2 = 0

S = 1 - Id - In - Imn - Imd
Sd = d * S * (1-m)
Sn = (1-d) * S * (1-m)
Smn = (1-d) * S * m
Smd = d * S * m
S1 = 0
S2 = 0
Sm1 = 0
Sm2 = 0


Rd = 0
Rn = 0
Rmn = 0
Rmd = 0
Rv = 0
Rmv = 0


# Przycisk do uruchomienia symulacji
if st.button('Run simulation'):
    initial_conditions = [Sd, Sn, Smn, Smd, S1, S2, Sm1, Sm2, V, Id, In, Imn, Imd, I1, I2, Im1, Im2, Rd, Rn, Rmn, Rmd, Rv, Rmv]
    result = simulate_vap_sirs_model(initial_conditions, beta_0, beta_m0, f, f_v, kappa, upsilon, upsilon_r, upsilon_m, upsilon_mr, omega, omega_m, a, a_m, gamma, days)

    abs_plot1, abs_plot2 = plot_absolute_values(result)
    st.markdown('### Group populations over time')
    st.pyplot(abs_plot1)
    st.pyplot(abs_plot2)

    # changes_in_infected_plot = plot_changes_in_infected(result)
    # st.markdown('### Changes in infected populations over time')
    # st.pyplot(changes_in_infected_plot)
