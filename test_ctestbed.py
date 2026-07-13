"""Regression suite for ghca_testbed.

Two tiers:
  * end-to-end validation of the executable C1 scenarios (verdicts + effects);
  * unit behaviour of the method-agnostic certificate / metric library;
  * a dedicated guard for the C1 interventional-vs-observational subtlety (the
    single highest-risk correctness property).

Substrate scenarios (C2–C7) are not executed here (they need the C-scripts and
minutes of compute); their ground truth is carried in the registry and checked
for well-formedness. Run: `pytest test_ctestbed.py` or `python test_ctestbed.py`.
"""

import numpy as np
import ghca_testbed as tb
from ghca_testbed.scm import build_scm, EXPECTED_VERDICT
from ghca_testbed.certificates import (epiphenomenality_test, theorem1_certificate,
                                       strata_means)


# --- C1 end-to-end --------------------------------------------------------

def test_c1_all_reference_checks_pass():
    rep = tb.score(tasks=("T1",))
    assert rep.total == 18 and rep.ok, str(rep)


def test_c1_theorem1_matches_expected_for_all_six_graphs():
    for name, expected in EXPECTED_VERDICT.items():
        _, _, edges, S = build_scm(name, n=2000, seed=0)
        assert theorem1_certificate(edges, S) == expected, name


def test_c1_confounded_high_association_zero_effect():
    # correlation != causation: strong obs association, epiphenomenal under do
    obs, do, edges, S = build_scm("confounded", n=200_000, seed=0)
    epi = epiphenomenality_test(obs, do, cond=tuple(S))
    assoc = tb.observational_association(obs, tuple(S))
    assert epi.verdict == "epiphenomenal" and epi.effect < 0.03
    assert assoc > 0.5


# --- the subtlety guard (highest value) -----------------------------------

def test_frontdoor_needs_interventional_vs_observational():
    """The correct Definition-1 contrast scores front-door CAUSAL; a naive
    do(1)-vs-do(0) contrast misses it (scores it ~epiphenomenal). This guards
    the C1 catch against a future regression."""
    obs, do, edges, S = build_scm("frontdoor", n=200_000, seed=0)

    correct = epiphenomenality_test(obs, do, cond=tuple(S))
    assert correct.verdict == "causal" and correct.effect > 0.2

    # naive do1-vs-do0: max over strata | E[B|c;do1] - E[B|c;do0] |
    d0 = strata_means(do[0], tuple(S))
    d1 = strata_means(do[1], tuple(S))
    naive = max(abs(d1[c] - d0[c]) for c in set(d0) & set(d1))
    assert naive < 0.03, f"do1-vs-do0 should miss front-door, got {naive}"


# --- metric library units -------------------------------------------------

def test_fat_hand_band_short_circuit_and_ratio():
    assert tb.fat_hand_band({"a": [1, 1], "b": [1, 1]}) == 0.0
    b = tb.fat_hand_band({"a": np.zeros(50), "b": np.ones(50)}, baseline_std=0.1)
    assert abs(b - 10.0) < 1e-6


def test_outcome_matrix_diagonal():
    sweeps = {"do_tau": {"timing": [8, 12, 18, 26], "identity": [0.5] * 4},
              "do_route": {"timing": [14] * 4, "identity": [0.1, 0.4, 0.6, 0.9]}}
    om = tb.outcome_matrix(sweeps)
    assert np.allclose(om.normalized, np.eye(2), atol=1e-6)
    assert om.is_diagonal_dominant()


def test_mediation_screens_when_generator_acts_through_mediator():
    rng = np.random.default_rng(0)
    n = 4000
    gen = rng.integers(0, 2, n)
    med = gen.copy()
    out = (med == 1).astype(float)
    res = tb.mediation(gen, med, out)
    assert res.screened and res.residual < 0.05


def test_macro_sufficiency_high_when_wave_informative():
    rng = np.random.default_rng(0)
    n = 2000
    B = rng.integers(0, 2, n)
    W = (B + 0.01 * rng.standard_normal(n))[:, None]
    S = rng.standard_normal((n, 5))
    assert tb.macro_sufficiency(B, W, S) > 3.0


# --- registry / export ----------------------------------------------------

def test_registry_shape():
    assert len(tb.REGISTRY) == 12
    c1 = [k for k in tb.REGISTRY if tb.REGISTRY[k].track == "C1"]
    assert len(c1) == 6 and all(tb.REGISTRY[k].executable for k in c1)


def test_export_idempotent(tmp_path):
    p1 = tb.export_ground_truth(str(tmp_path / "a"))
    p2 = tb.export_ground_truth(str(tmp_path / "b"))
    import filecmp
    for fn in ("scenarios.json", "ground_truth.json"):
        assert filecmp.cmp(f"{p1}/{fn}", f"{p2}/{fn}", shallow=False), fn


if __name__ == "__main__":
    import sys
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            import inspect
            if "tmp_path" in inspect.signature(fn).parameters:
                import tempfile
                fn(__import__("pathlib").Path(tempfile.mkdtemp()))
            else:
                fn()
            print(f"PASS {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL {fn.__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
