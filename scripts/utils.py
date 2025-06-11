import itertools

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import pearsonr

# Models
def dose_time_response_model(params, dt):
    k_alpha, k_beta, k_gamma, k_delta = params
    d, t = dt
    return (2**(k_alpha*t) - 2**(k_delta*t))/(1 + 10**(k_beta*np.abs(2**(k_alpha*t) - 2**(k_delta*t))*(d - k_gamma))) + 2**(k_delta*t)

def dose_time_response_model_vus(params, dt):
    k_alpha, k_beta, k_gamma, k_delta = params
    d, t = dt
    
    # To treat drug effects based on relative change in cell proliferation rates and not absolute change such that cell lines with high proliferation rate do not get better VUS values than those with low proliferation rate when the drug effect is the same
    if k_alpha >= k_delta:
        k_alpha, k_delta = 1, k_delta/k_alpha
    else: # avoids numerical issues if both k_alpha and k_delta are between 0 and 1, and the relative change from k_alpha to k_delta is still the same
        k_alpha, k_delta = k_alpha/k_delta, 1
    
    return (2**(k_alpha*t) - 2**(k_delta*t))/(1 + 10**(k_beta*np.abs(2**(k_alpha*t) - 2**(k_delta*t))*(d - k_gamma))) + 2**(k_delta*t)

def dose_response_model(params, d):
    alpha, beta, gamma, delta = params
    return (alpha - delta)/(1 + 10**(beta*(d - gamma))) + delta

def residuals_dose_time_response_model(params, dt, z):
    return z - dose_time_response_model(params, dt)

def residuals_dose_response_model(params, d, z):
    return z - dose_response_model(params, d)

def dose_response_model_auc(params, d_min, d_max):
    alpha, beta, gamma, delta = params
    norm_factor = dose_response_model(params, d_min)*(d_max - d_min)
    auc = alpha*(d_max - d_min) + (alpha - delta)/beta*np.log10((10**(beta*d_min) + 10**(beta*gamma))/(10**(beta*d_max) + 10**(beta*gamma)))
    auc_norm = auc/norm_factor
    return auc, auc_norm

# Metrics
def scoring_function(z, z_pred, metric):
    z = np.array(z)
    z_pred = np.array(z_pred)
    if metric == "rmse":
        score = np.sqrt(np.sum(((z_pred - z)**2))/len(z))
    if metric == "mae":
        score = np.sum(np.abs(z_pred - z))/len(z)
    if metric == "mape":
        score = np.sum(np.abs((z_pred - z)/z))/len(z)
    return score

def concordance_correlation_coefficient(z, z_pred):
    mean_true = np.mean(z)
    mean_pred = np.mean(z_pred)
    std_true = np.std(z)
    std_pred = np.std(z_pred)
    pearson_correlation_coefficient = pearsonr(z, z_pred)[0]
    ccc = (2*pearson_correlation_coefficient*std_true*std_pred)/(std_true**2 + std_pred**2 + (mean_true - mean_pred)**2)
    return ccc

# Fitting
def fit_model(model_function, residual_function, param_guesses, function_input, function_output, bounds, metric):
    best_score = np.inf
    best_params = [np.nan]*len(param_guesses)
    for guess in itertools.product(*param_guesses):
        result = least_squares(residual_function, guess, args=(function_input, function_output), bounds=bounds)
        params = result["x"]
        score = scoring_function(model_function(params, function_input), function_output, metric)
        if best_score > score:
            best_score = score
            best_params = params
    return best_params, best_score