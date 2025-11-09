import itertools

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import pearsonr

# Models
def dose_time_response_model(params, dt):
    k_alpha, a_alpha, k_beta, k_gamma, k_delta, a_delta = params
    d, t = dt
    
    alpha_t = (a_alpha*2**(k_alpha*t))/((a_alpha - 1) + 2**(k_alpha*t))
    delta_t = (a_delta*2**(k_delta*t))/((a_delta - 1) + 2**(k_delta*t))
    beta_t = k_beta*np.abs(alpha_t - delta_t)
    gamma_t = k_gamma
    
    return (alpha_t - delta_t)/(1 + 10**(beta_t*(d - gamma_t))) + delta_t

def dose_time_response_model_grivus(params, dt):
    k_alpha, a_alpha, k_beta, k_gamma, k_delta, a_delta = params
    d, t = dt
    
    if k_alpha > k_delta:
        k_alpha, k_delta = 1, k_delta/k_alpha
    else:
        k_alpha, k_delta = k_alpha/k_delta, 1
    
    alpha_t = (a_alpha*2**(k_alpha*t))/((a_alpha - 1) + 2**(k_alpha*t))
    delta_t = (a_delta*2**(k_delta*t))/((a_delta - 1) + 2**(k_delta*t))
    beta_t = k_beta*np.abs(alpha_t - delta_t)
    gamma_t = k_gamma
    return (alpha_t - delta_t)/(1 + 10**(beta_t*(d - gamma_t))) + delta_t

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

# Bootstrap confidence interval for CCC
def bootstrap_ccc(y_true, y_pred, n_boot=100, alpha=0.05, random_state=None):
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    
    rng = np.random.default_rng(random_state)
    n = len(y_true)
    boot_stats = []
    
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        boot_y_true = y_true[idx]
        boot_y_pred = y_pred[idx]
        boot_stats.append(concordance_correlation_coefficient(boot_y_true, boot_y_pred))

    lower = np.percentile(boot_stats, 100 * alpha/2)
    upper = np.percentile(boot_stats, 100 * (1 - alpha/2))
    return np.mean(boot_stats), (lower, upper)

# For scipy's least_squares: if trf algorithm throws an error, stick with last found solution
last_x = None
def callback(x, *args):
    global last_x
    last_x = x.copy()

# Fitting
def fit_model(model_function, residual_function, param_guesses, function_input, function_output, bounds, metric):
    best_score = np.inf
    best_params = [np.nan]*len(param_guesses)
    for guess in itertools.product(*param_guesses):
        try:
            result = least_squares(residual_function, guess, args=(function_input, function_output), bounds=bounds, callback=callback)
            params = result["x"]
        except ValueError:
            # Use last_x as best-guess fallback
            params = guess if last_x is None else last_x
        score = scoring_function(model_function(params, function_input), function_output, metric)
        if best_score > score:
            best_score = score
            best_params = params
    return best_params, best_score