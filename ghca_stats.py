"""Statistics harness for Track 3a (see docs/stats_sweeps_plan.md).

Small, dependency-light, and seeded like everything else in the repo (never the
global RNG — the perturb_tau lesson). Provides:

- `bootstrap_ci`  : percentile bootstrap CI for any statistic of a sample.
- `wilson_ci`     : CI for a proportion (for rates like E3 joint-success).
- `effect_size`   : standardised gap between two samples (Cohen's d, pooled SD).
- `bimodality`    : Sarle's bimodality coefficient + a verdict flag (the E3
                    discipline, automated so no headline hides a two-mode mean).
- `summarise`     : one-line row (mean, CI, n, bimodality) for the master table.
- `strip_ci`      : matplotlib strip-of-seeds + mean±CI plot helper.

Publication-grade default is n=50 seeds and 10 000 bootstrap resamples.
"""

import numpy as np

BOOT_DEFAULT = 10_000
# Sarle's bimodality coefficient exceeds this for a uniform / two-mode split;
# 5/9 is the value for a uniform distribution (the standard cut).
BC_UNIFORM = 5.0 / 9.0


def bootstrap_ci(data, statistic=np.mean, n_boot=BOOT_DEFAULT, ci=95, seed=0):
    """Percentile bootstrap CI. Returns (point, lo, hi).

    `data` is a 1-D array of per-seed values; `statistic` is applied to each
    resample. Seeded explicitly.
    """
    data = np.asarray(data, float)
    data = data[np.isfinite(data)]
    n = data.size
    point = float(statistic(data)) if n else float("nan")
    if n < 2:
        return point, point, point
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    boot = statistic(data[idx], axis=1)
    lo, hi = np.percentile(boot, [(100 - ci) / 2, 100 - (100 - ci) / 2])
    return point, float(lo), float(hi)


def wilson_ci(k, n, ci=95):
    """Wilson score interval for a proportion k/n. Returns (p, lo, hi).

    Preferred over bootstrap for a rate (e.g. E3 joint-success across seeds),
    which is what a small-n binomial needs.
    """
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    from scipy.stats import norm
    z = norm.ppf(1 - (1 - ci / 100) / 2)
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = (z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))) / denom
    return float(p), float(centre - half), float(centre + half)


def effect_size(a, b):
    """Cohen's d between samples a and b (pooled SD). Sign: mean(a) - mean(b)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    na, nb = a.size, b.size
    if na < 2 or nb < 2:
        return float("nan")
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    if sp == 0:
        return float("inf") if a.mean() != b.mean() else 0.0
    return float((a.mean() - b.mean()) / sp)


def bimodality(data):
    """Sarle's bimodality coefficient and a verdict.

    BC = (skew^2 + 1) / (kurtosis + 3*(n-1)^2/((n-2)(n-3))).
    BC > 5/9 flags a distribution at least as split as uniform (possibly
    bimodal). Returns dict(bc, bimodal, n). Robust to tiny/degenerate samples.
    """
    data = np.asarray(data, float)
    data = data[np.isfinite(data)]
    n = data.size
    if n < 4 or data.std() == 0:
        return {"bc": float("nan"), "bimodal": False, "n": int(n)}
    m = data - data.mean()
    s = data.std()  # population moments
    skew = np.mean(m**3) / s**3
    kurt = np.mean(m**4) / s**4 - 3.0  # excess kurtosis
    corr = 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    bc = (skew**2 + 1.0) / (kurt + corr)
    return {"bc": float(bc), "bimodal": bool(bc > BC_UNIFORM), "n": int(n)}


def summarise(name, data, seed=0, n_boot=BOOT_DEFAULT):
    """One row for the master table: mean, 95% bootstrap CI, SD, n, bimodality."""
    data = np.asarray(data, float)
    finite = data[np.isfinite(data)]
    point, lo, hi = bootstrap_ci(data, seed=seed, n_boot=n_boot)
    bim = bimodality(data)
    return {
        "name": name,
        "n": int(finite.size),
        "mean": point,
        "ci_lo": lo,
        "ci_hi": hi,
        "sd": float(finite.std(ddof=1)) if finite.size > 1 else float("nan"),
        "bc": bim["bc"],
        "bimodal": bim["bimodal"],
    }


def fmt_row(row):
    """Render a summarise() row as a markdown-table line fragment."""
    flag = " ⚠bimodal" if row["bimodal"] else ""
    return (f"{row['name']}: {row['mean']:.3f} "
            f"[{row['ci_lo']:.3f}, {row['ci_hi']:.3f}] (n={row['n']})"
            f"{flag}")


def strip_ci(ax, samples, labels, seed=0):
    """Strip plot of per-seed points + mean and 95% bootstrap CI per group.

    `samples` is a list of 1-D arrays; `labels` names them. Deterministic jitter
    (seeded), so figures reproduce.
    """
    rng = np.random.default_rng(seed)
    for i, (s, lab) in enumerate(zip(samples, labels)):
        s = np.asarray(s, float)
        s = s[np.isfinite(s)]
        jit = rng.uniform(-0.08, 0.08, size=s.size)
        ax.scatter(np.full(s.size, i) + jit, s, s=18, alpha=0.5, color="0.5", zorder=1)
        pt, lo, hi = bootstrap_ci(s, seed=seed + i)
        ax.errorbar(i, pt, yerr=[[pt - lo], [hi - pt]], fmt="o", color="crimson",
                    capsize=5, zorder=3, ms=7)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    return ax


if __name__ == "__main__":
    # Self-test: unimodal vs bimodal detection, CI sanity, effect size.
    rng = np.random.default_rng(0)
    uni = rng.normal(0.9, 0.02, 50)
    bim = np.concatenate([rng.normal(0.0, 0.02, 25), rng.normal(1.0, 0.02, 25)])
    print("unimodal :", fmt_row(summarise("uni", uni)))
    print("bimodal  :", fmt_row(summarise("bim", bim)))
    assert not bimodality(uni)["bimodal"], "false positive on unimodal"
    assert bimodality(bim)["bimodal"], "missed a clear bimodal split"
    p, lo, hi = wilson_ci(1, 5)  # E3's historical 1/5 joint-success
    print(f"wilson 1/5: {p:.2f} [{lo:.2f}, {hi:.2f}]")
    print(f"effect size uni-vs-bim: {effect_size(uni, bim):.2f}")
    print("ghca_stats self-test OK")
