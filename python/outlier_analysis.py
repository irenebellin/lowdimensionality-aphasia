"""
Multivariate outlier analysis (Mahalanobis distance) of the Pure Aphasia Cohort
relative to the General Stroke Cohort.

Input: two plain-text files, whitespace-separated, one row per patient,
columns = PC1 PC2 PC3 (the behavioral PC scores used for Figure 2).
    python outlier_analysis.py general_stroke_scores.txt pure_aphasia_scores.txt
Runs the analysis in the PC1-PC3 plane (as in Figure 2) and in the full PC1-PC3 space.
"""
import sys
import numpy as np
from scipy import stats

G = np.loadtxt(sys.argv[1]); P = np.loadtxt(sys.argv[2])
print(f"General Stroke Cohort n = {len(G)}, Pure Aphasia Cohort n = {len(P)}")

def outlier_analysis(cols, label, reps=20000, seed=1):
    g, p = G[:, cols], P[:, cols]
    mu, cov_inv = g.mean(0), np.linalg.inv(np.cov(g.T))
    d2 = lambda X: np.einsum('ij,jk,ik->i', X - mu, cov_inv, X - mu)
    dg, dp = d2(g), d2(p)
    cut = np.percentile(dg, 95)             # empirical cutoff from the General Stroke Cohort itself
    n_out = int((dp > cut).sum())
    rng = np.random.default_rng(seed)
    obs, n = dp.mean(), len(dp)
    count = sum(dg[rng.choice(len(dg), n, replace=False)].mean() >= obs for _ in range(reps))
    print(f"\n[{label}] cutoff D2 = {cut:.1f}; beyond cutoff: {n_out} of {n} (expected {0.05*n:.1f})")
    print(f"  mean D2 {obs:.2f} vs {dg.mean():.2f}; permutation p = {(count+1)/(reps+1):.3f}")
    print(f"  median D2 {np.median(dp):.2f} vs {np.median(dg):.2f}; one-sided Mann-Whitney p = "
          f"{stats.mannwhitneyu(dp, dg, alternative='greater').pvalue:.3f}")

outlier_analysis([0, 2], "PC1-PC3 plane (Figure 2)")
outlier_analysis([0, 1, 2], "PC1-PC3 space (sensitivity analysis)")

print("\nComponent-wise comparisons (two-sided Mann-Whitney)")
for k in range(3):
    u = stats.mannwhitneyu(P[:, k], G[:, k])
    print(f"  PC{k+1}: median {np.median(P[:,k]):.2f} vs {np.median(G[:,k]):.2f}; p = {u.pvalue:.4f}")
imp = G[:, 0] < np.median(G[:, 0])   # post hoc: more impaired half (PC1 below median; check sign convention)
u = stats.mannwhitneyu(P[:, 0], G[imp, 0])
print(f"  PC1 vs more impaired half of General Stroke Cohort: median {np.median(P[:,0]):.2f} vs "
      f"{np.median(G[imp,0]):.2f}; p = {u.pvalue:.3f}")
