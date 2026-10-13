"""
Comparison of the two cohorts from the summary statistics of Tables 1 and 2.

Raw demographic data were not available, so the tests are computed from group
means, SDs and percentages (percentages are converted to counts, which makes
the p-values approximate). Welch t-tests for continuous variables, Fisher's
exact tests for proportions.

    python cohort_comparison_from_summary.py
"""
from scipy import stats
import numpy as np

N_G, N_P = 155, 21          # General Stroke Cohort, Pure Aphasia Cohort

def welch(m1, s1, m2, s2, label):
    r = stats.ttest_ind_from_stats(m1, s1, N_G, m2, s2, N_P, equal_var=False)
    d = (m1 - m2) / np.sqrt(((N_G - 1) * s1**2 + (N_P - 1) * s2**2) / (N_G + N_P - 2))
    print(f"{label:34s} t = {r.statistic:6.2f}  p = {r.pvalue:.4f}  d = {d:5.2f}")

def fisher(a, b, c, d, label):
    odds, p = stats.fisher_exact([[a, b], [c, d]])
    print(f"{label:34s} Fisher p = {p:.4f}")

print("Continuous variables (General Stroke Cohort vs Pure Aphasia Cohort)")
welch(68.4, 8.0, 68.8, 11.6, "Age (years)")
welch(9.6, 4.0, 10.6, 5.0, "Education (years)")
welch(7.0, 5.2, 4.0, 3.0, "NIHSS")

print("\nProportions (counts derived from the reported percentages)")
fisher(round(0.54 * N_G), N_G - round(0.54 * N_G), round(0.71 * N_P), N_P - round(0.71 * N_P), "Sex (males)")
fisher(round(0.90 * N_G), N_G - round(0.90 * N_G), N_P, 0, "Stroke subtype (ischemic)")
fisher(round(0.48 * N_G), round(0.44 * N_G), round(0.86 * N_P), round(0.14 * N_P),
       "Lesion side (left vs right)")
