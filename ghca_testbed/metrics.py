"""Method-agnostic causal metrics used across the C-series.

All take samples/arrays, not a live substrate. Formulas match the C-series
(C2 fat-hand band, C4 macro-sufficiency + outcome matrix, C7 mediation/screening).
"""

from dataclasses import dataclass
import numpy as np


def binarize_median(x):
    """Median split to a balanced binary label (C-series behaviour B)."""
    x = np.asarray(x, float)
    return (x >= np.median(x)).astype(int)


def fat_hand_band(effect_by_policy, baseline_std=None):
    """Fat-handedness band in sigma units (C2/C5).

    `effect_by_policy`: {policy: outcome-samples-or-mean}. band = across-policy
    range of the per-policy means, divided by a spread: `baseline_std` if given
    (C2/C3 observational-sigma units), else the pooled within-policy std (C5).
    Returns exactly 0.0 when all policy means are equal (the C6 short-circuit).
    """
    vals = list(effect_by_policy.values())
    means = np.array([np.mean(np.atleast_1d(v)) for v in vals])
    rng = float(means.max() - means.min())
    if rng == 0.0:
        return 0.0
    if baseline_std is not None:
        std = float(baseline_std)
    else:
        within = [np.std(np.atleast_1d(v)) for v in vals if np.size(np.atleast_1d(v)) > 1]
        std = float(np.mean(within)) if within else 0.0
    if std <= 0.0:
        return float("inf")
    return rng / std


def _default_decoder(X, y, cv=5):
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
    X = np.asarray(X)
    if X.ndim == 1:
        X = X[:, None]
    return float(cross_val_score(clf, X, y, cv=cv).mean())


def macro_sufficiency(B, W, S, decoder=None, chance=0.5, cv=5):
    """I(B;W)/I(B;S) proxied by decodability above chance (C4).

    suff = max(acc(W->B)-chance, 0) / max(acc(S->B)-chance, eps). May exceed 1
    (the macro can be as sufficient as the micro); NOT clipped.
    """
    dec = decoder or (lambda X, y: _default_decoder(X, y, cv=cv))
    aW = dec(W, B) - chance
    aS = dec(S, B) - chance
    return float(max(aW, 0.0) / max(aS, 1e-6))


@dataclass
class OutcomeResult:
    raw: np.ndarray
    normalized: np.ndarray
    handles: list
    outcomes: list

    def is_diagonal_dominant(self, tol=0.2):
        M = self.normalized
        k = min(M.shape)
        return all(M[i, i] >= M[i, j] - tol for i in range(k) for j in range(M.shape[1]))


def outcome_matrix(sweeps, normalize="column"):
    """(handle, outcome) causal matrix (C4/C7).

    `sweeps`: {handle: {outcome: array-of-outcome-values-over-the-handle-sweep}}.
    effect[i,j] = range (ptp) of outcome j as handle i is swept. `column`
    normalisation divides each outcome column by its max range.
    """
    handles = list(sweeps)
    outcomes = list(next(iter(sweeps.values())))
    raw = np.array([[np.ptp(np.asarray(sweeps[h][o])) for o in outcomes]
                    for h in handles], dtype=float)
    if normalize == "column":
        norm = raw / (raw.max(0, keepdims=True) + 1e-9)
    elif normalize == "row":
        norm = raw / (raw.max(1, keepdims=True) + 1e-9)
    else:
        norm = raw.copy()
    return OutcomeResult(raw=raw, normalized=norm, handles=handles, outcomes=outcomes)


@dataclass
class MediationResult:
    screened: bool
    residual: float          # max generator effect on B within a mediator stratum
    by_mediator: dict


def mediation(generator, mediator, outcome, min_stratum=1):
    """Screening test on a generator -> mediator -> outcome chain (C7-B).

    Does the generator affect the outcome only THROUGH the mediator? Within each
    mediator value, measure how much E[outcome] still varies with the generator;
    if ~0 the mediator screens the generator off. `residual` is that max variation.
    """
    g = np.asarray(generator)
    m = np.asarray(mediator)
    o = np.asarray(outcome, float)
    by = {}
    residual = 0.0
    for mv in np.unique(m):
        sel = m == mv
        if sel.sum() < min_stratum:
            continue
        gvals = np.unique(g[sel])
        means = []
        for gv in gvals:
            s2 = sel & (g == gv)
            if s2.sum() >= min_stratum:
                means.append(o[s2].mean())
        if len(means) >= 2:
            spread = float(max(means) - min(means))
            by[float(mv)] = spread
            residual = max(residual, spread)
    return MediationResult(screened=residual < 0.05, residual=residual, by_mediator=by)
