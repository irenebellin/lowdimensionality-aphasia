"""
Sensitivity analysis: the smallest effect the study could detect.

No a priori power analysis was possible (retrospective Pure Aphasia Cohort,
General Stroke Cohort from a previous study). This script reports, for a
two-sided alpha and 80% power, the smallest standardized difference (Cohen's d)
detectable in an independent-samples comparison between a group of n patients
and the General Stroke Cohort, and the smallest detectable correlation.

    python sensitivity_analysis.py
"""
import numpy as np
from scipy import stats
from statsmodels.stats.power import TTestIndPower

N_GENERAL = 155        # General Stroke Cohort
POWER = 0.80
ALPHA = 0.05
N_COMPARISONS = 9      # nine syndrome groups (Bonferroni)
GROUPS = {"Pure Aphasia Cohort (all)": 21,
          "largest syndrome groups (n = 4)": 4,
          "n = 3": 3,
          "n = 2": 2}

power = TTestIndPower()

def min_d(n_small, alpha):
    return power.solve_power(nobs1=n_small, ratio=N_GENERAL / n_small,
                             alpha=alpha, power=POWER, alternative="two-sided")

def min_r(n, alpha):
    z = stats.norm.ppf(1 - alpha / 2) + stats.norm.ppf(POWER)
    return np.tanh(z / np.sqrt(n - 3))

print(f"Smallest detectable effect (power {POWER:.0%}, two-sided)\n")
print(f"{'comparison with General Stroke Cohort':40s} {'d (alpha=0.05)':>15s} {'d (Bonferroni)':>15s}")
for label, n in GROUPS.items():
    print(f"{label:40s} {min_d(n, ALPHA):15.2f} {min_d(n, ALPHA / N_COMPARISONS):15.2f}")
print(f"\nSmallest detectable correlation, n = {N_GENERAL}: r = {min_r(N_GENERAL, ALPHA):.2f}")
print("Groups with n = 1 cannot be analysed inferentially.")
