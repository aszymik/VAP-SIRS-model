from json import load
from model import *
from plots import *
import pandas as pd


def fitted_scenario(filename, initial_conditions):


    with open(filename, 'r') as file:
        params = load(file)

    del params['days']

    result_syn = generate_synthetic_data(initial_conditions, 
                                            days = 730, 
                                            seasonal_amplitude=0.2, 
                                            noise_std=1,
                                            **params)


    original_params = list(params.values())
    init_params = [0.15, 0.4, 0.6, 0.45, 0.0050, 0.01, 0.01, 0.02, 0.02, 0.005, 0.0035, 0.87, 0.7, 0.2, 0.45]

    result_fit, fit_params = fit_to_real_data(initial_conditions, result_syn, init_params)

    plot_syn = plot_all(result_syn, label='Real-world data simulation')
    plot_fit = plot_all(result_fit, label='Fitted model data simulation')

    params_labels = ["beta_0", "beta_m0", "f", "fv", "kappa", "v", "v_r", "v_m", "v_mr", "omega", "omega_m", "alpha", "alpha_m", "gamma", "% hidden cases"]
    params_df = pd.DataFrame({'params': params_labels, 'Real-world data params': original_params, 'Fitted model params': fit_params, 'Initial params for optimalization': init_params},)
    params_df = pd.melt(params_df, id_vars='params', var_name='Data_type', value_name='value')

    params_plot = px.scatter(params_df, x='params', y='value', color='Data_type')
    params_plot.update_layout(
        title='Parameters comparison',
        xaxis_title='Value',
        yaxis_title='Parameters',
        legend_title='Data Type',
        hovermode="y",
        autosize=False, width=1000, height=600,
    )

    return plot_syn, plot_fit, params_plot
    
def fit_seasonal(data_noise, params, initial_conditions):

    del params['days']
    original_params = list(params.values())
    init_params = [0.15, 0.4, 0.6, 0.45, 0.0050, 0.01, 0.01, 0.02, 0.02, 0.005, 0.0035, 0.87, 0.7, 0.2, 0.45]

    result_fit, fit_params = fit_to_real_data(initial_conditions, data_noise, init_params)

    plot_syn = plot_all(data_noise, label='Real-world data simulation')
    plot_fit = plot_all(result_fit, label='Fitted model data simulation')

    params_labels = ["beta_0", "beta_m0", "f", "fv", "kappa", "v", "v_r", "v_m", "v_mr", "omega", "omega_m", "alpha", "alpha_m", "gamma", "% hidden cases"]
    params_df = pd.DataFrame({'params': params_labels, 'Real-world data params': original_params, 'Fitted model params': fit_params, 'Initial params for optimalization': init_params},)
    params_df = pd.melt(params_df, id_vars='params', var_name='Data_type', value_name='value')

    params_plot = px.scatter(params_df, x='params', y='value', color='Data_type')
    params_plot.update_layout(
        title='Parameters comparison',
        xaxis_title='Value',
        yaxis_title='Parameters',
        legend_title='Data Type',
        hovermode="y",
        autosize=False, width=1000, height=600,
    )

    return plot_syn, plot_fit, params_plot