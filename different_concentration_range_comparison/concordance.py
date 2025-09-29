import numpy as np
import pandas as pd

from scipy.stats import pearsonr

def concordance_correlation_coefficient(y_true, y_pred):
    mean_true = np.mean(y_true)
    mean_pred = np.mean(y_pred)
    std_true = np.std(y_true)
    std_pred = np.std(y_pred)
    pearson_correlation_coefficient, _ = pearsonr(y_true, y_pred)
    ccc = (2*pearson_correlation_coefficient*std_true*std_pred)/(std_true**2 + std_pred**2 + (mean_true - mean_pred)**2)
    return ccc

# Bootstrap confidence interval for CCC
def bootstrap_ccc(y_true, y_pred, n_boot=100, alpha=0.05, random_state=0):
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

relative_vus_decrease = pd.read_csv("fit_surface_results_sorted_by_relative_vus_decrease_after_72_daily_False.csv", index_col=0)
relative_vus_decrease_no_min_conc = pd.read_csv("fit_surface_results_sorted_by_relative_vus_decrease_after_72_daily_False_no_min_conc.csv", index_col=0)
relative_vus_decrease_no_max_conc = pd.read_csv("fit_surface_results_sorted_by_relative_vus_decrease_after_72_daily_False_no_max_conc.csv", index_col=0)
relative_vus_decrease_no_minmax_conc = pd.read_csv("fit_surface_results_sorted_by_relative_vus_decrease_after_72_daily_False_no_minmax_conc.csv", index_col=0)

relative_vus_decrease = relative_vus_decrease.sort_values(by=["cell_line", "drug"])["relative_vus_decrease"]
relative_vus_decrease_no_min_conc = relative_vus_decrease_no_min_conc.sort_values(by=["cell_line", "drug"])["relative_vus_decrease"]
relative_vus_decrease_no_max_conc = relative_vus_decrease_no_max_conc.sort_values(by=["cell_line", "drug"])["relative_vus_decrease"]
relative_vus_decrease_no_minmax_conc = relative_vus_decrease_no_minmax_conc.sort_values(by=["cell_line", "drug"])["relative_vus_decrease"]

pearsonr_no_min, pval_no_min = pearsonr(relative_vus_decrease, relative_vus_decrease_no_min_conc)
ccc_no_min = concordance_correlation_coefficient(relative_vus_decrease, relative_vus_decrease_no_min_conc)
bootstrap_ccc_no_min = bootstrap_ccc(relative_vus_decrease, relative_vus_decrease_no_min_conc)
#print(pearsonr_no_min, pval_no_min)
print(ccc_no_min, bootstrap_ccc_no_min)

pearsonr_no_max, pval_no_max = pearsonr(relative_vus_decrease, relative_vus_decrease_no_max_conc)
ccc_no_max = concordance_correlation_coefficient(relative_vus_decrease, relative_vus_decrease_no_max_conc)
bootstrap_ccc_no_max = bootstrap_ccc(relative_vus_decrease, relative_vus_decrease_no_max_conc)
#print(pearsonr_no_max, pval_no_max)
print(ccc_no_max, bootstrap_ccc_no_max)

pearsonr_no_minmax, pval_no_minmax = pearsonr(relative_vus_decrease, relative_vus_decrease_no_minmax_conc)
ccc_no_minmax = concordance_correlation_coefficient(relative_vus_decrease, relative_vus_decrease_no_minmax_conc)
bootstrap_ccc_no_minmax = bootstrap_ccc(relative_vus_decrease, relative_vus_decrease_no_minmax_conc)
#print(pearsonr_no_minmax, pval_no_minmax)
print(ccc_no_minmax, bootstrap_ccc_no_minmax)

pearsonr_no_min_no_max, pval_no_min_no_max = pearsonr(relative_vus_decrease_no_max_conc, relative_vus_decrease_no_min_conc)
ccc_no_min_no_max = concordance_correlation_coefficient(relative_vus_decrease_no_max_conc, relative_vus_decrease_no_min_conc)
bootstrap_ccc_no_min_no_max = bootstrap_ccc(relative_vus_decrease_no_max_conc, relative_vus_decrease_no_min_conc)
#print(pearsonr_no_min_no_max, pval_no_min_no_max)
print(ccc_no_min_no_max, bootstrap_ccc_no_min_no_max)