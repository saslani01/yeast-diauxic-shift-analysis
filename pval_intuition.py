import random 
import matplotlib.pyplot as plt
from scipy.stats import binom, hypergeom
import os

"""
N : population size
n : sample size
K : positive counts in the population
k : positive counts in the sample
t : number of times we sample
"""

def population(N, K):
    return [1] * K + [0] * (N - K)

def sample(population, n):
    sample = random.sample(population, n)
    k = sum(sample)
    return k

def repeated_sample(population, t, n):
    results = []
    for _ in range(t):
        k = sample(population, n)
        results.append(k)
    return results

def empirical_pval(N, K, n, t, observed_val):
    results = repeated_sample(population(N, K), t, n)
    return sum(1 for k in results if k >= observed_val) / len(results)

def create_p_val_trial_count_convergence_plot(N, K, n, observed_val, t_vals):
    pvals = []
    for t in t_vals:
        pvals.append(empirical_pval(N, K, n, t, observed_val))
    
    pval_changes = [abs(pvals[i] - pvals[i - 1]) for i in range(1, len(pvals))]

    plt.plot(t_vals[1:], pval_changes, marker='o') # skipping first trial to match the dimension of pval_changes (also useless because of small n)
    plt.xscale('log')
    plt.xlabel("Number of Trials (Upper Bound)")
    plt.xticks(t_vals[1:])
    plt.minorticks_off()
    plt.ylabel("P-value Difference")
    plt.title("Convergence of P-Value Differences Across Trials")
    os.makedirs("outputs", exist_ok=True)
    plt.savefig("outputs/pval_delta_convergence.png")
    plt.close()
    return pvals, pval_changes

def p_vals_using_different_models(N, K, n, observed_val, t):
    emp_pval = empirical_pval(N, K, n, t, observed_val)
    bin_pval = 1 - binom.cdf(observed_val - 1, n, K / N)
    hyp_pval = 1 - hypergeom.cdf(observed_val - 1, N, K, n)

    return emp_pval, bin_pval, hyp_pval

def is_significant(p_val, thresh):
    return p_val < thresh

def find_sucess_sample_and_population_params_for_pval_range(N, n, p_min, p_max):
    for K in range(0, N + 1, 20):
        for k in range(0, n + 1):
            hyp_p = 1 - hypergeom.cdf(k - 1, N, K, n)
            if p_min <= hyp_p <= p_max:
                return (K, k, hyp_p)
    return None
            