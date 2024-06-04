import streamlit as st
from model import *


st.title('VAP-SIRS Model Simulation')

with st.sidebar:
    # Suwaki do wyboru parametrów
    beta_0 = st.slider('Beta', min_value=0.0, max_value=1.0, value=0.2)
    beta_m0 = st.slider('Beta m', min_value=0.0, max_value=1.0, value=0.4)

    f = st.slider('f', min_value=0.0, max_value=1.0, value=0.77)
    f_v = st.slider('fv', min_value=0.0, max_value=1.0, value=0.55)

    upsilon = st.slider('Upsilon', min_value=0.0, max_value=1.0, value=0.004)
    upsilon_r = st.slider('Upsilon R', min_value=0.0, max_value=1.0, value=0.003)
    upsilon_m = st.slider('Upsilon m', min_value=0.0, max_value=1.0, value=0.006)
    upsilon_mr = st.slider('Upsilon m R', min_value=0.0, max_value=1.0, value=0.006)
    
    omega = 1/st.slider('Number of days after which vacinated normal people lose their immunity', min_value=30, max_value=700, value=365)
    omega_m = 1/st.slider('Number of days after which vacinated more susceptible people lose their immunity', min_value=30, max_value=700, value=365)
    kappa = 1/st.slider('Number of days after which recovered people lose their immunity', min_value=30, max_value=700, value=400)

    a = st.slider('Vaccine effectiveness for normal people', min_value=0.0, max_value=1.0, value=0.79)
    a_m = st.slider('vaccine effectiveness for more susceptible people', min_value=0.0, max_value=1.0, value=0.60)

    gamma = 1/st.slider('Duration of the disease', min_value=1, max_value=21, value=7)
    d = st.slider("% of people who won't vaccinate", min_value=0.0, max_value=1.0, value=0.12)
    m = st.slider('% of more susceptible people', min_value=0.0, max_value=1.0, value=0.10)
    days = st.number_input('Days', min_value=1, max_value=2000, value=100, step=1)

    

# Wartości początkowe
N = 100
V = 0
Vm = 0

I_percn = 0.1

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


# Przycisk do uruchomienia symulacji
if st.button('Run simulation'):
    initial_conditions = [(1-I_percn)*Sd, (1-I_percn)*Sn, (1-I_percn)*Smn, (1-I_percn)*Smd, S1, S2, Sm1, Sm2, V, Vm, Id, In, Imn, Imd, I1, I2, Im1, Im2, Rd, Rn, Rmn, Rmd, Rv, Rmv]
    # print(f'{sum(initial_conditions)=}')
    result = simulate_vap_sirs_model(initial_conditions, beta_0, beta_m0, f, f_v, kappa, upsilon, upsilon_r, upsilon_m, upsilon_mr, omega, omega_m, a, a_m, gamma, days)

    abs_plot1, abs_plot2 = plot_absolute_values(result)
    st.markdown('### Group populations over time')
    st.plotly_chart(abs_plot1)
    st.plotly_chart(abs_plot2)

    # changes_in_infected_plot = plot_changes_in_infected(result)
    # st.markdown('### Changes in infected populations over time')
    # st.pyplot(changes_in_infected_plot)
