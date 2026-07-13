"""The ground-truth scenario registry + language-agnostic export.

C1 scenarios are executable in-process (`run`); C2–C7 substrate scenarios carry
their ground truth as verified data (from the C-series results docs / plan
table), with `script` naming the experiment that reproduces them in full.
"""

import functools
import json
import os

from .scm import build_scm, EXPECTED_VERDICT, EDGES, COND
from .certificates import (epiphenomenality_test, observational_association,
                           theorem1_certificate)
from .scenario import Scenario, GraphSpec, GroundTruth

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_OUT = os.path.join(_ROOT, "result", "testbed")

# reference effect / association per C1 graph (deterministic MC, N=200k seed 0)
_C1_REF = {
    "fork": (0.001, 0.004), "confounded": (0.001, 0.642), "direct": (0.400, 0.799),
    "mediator_unobs": (0.322, 0.640), "mediator_obs": (0.008, 0.001),
    "frontdoor": (0.258, 0.634),
}


def _run_scm(name, n=200_000, seed=0):
    """Measure a C1 scenario end-to-end: Definition-1 effect + Theorem-1 cert."""
    obs, do, edges, S = build_scm(name, n=n, seed=seed)
    epi = epiphenomenality_test(obs, do, cond=tuple(S))
    return {
        "effect": epi.effect,
        "verdict": epi.verdict,
        "obs_assoc": observational_association(obs, tuple(S)),
        "certificate": theorem1_certificate(edges, S),
    }


def _c1_scenarios():
    scen = {}
    for name in EXPECTED_VERDICT:
        eff, assoc = _C1_REF[name]
        scen[f"c1_{name}"] = Scenario(
            name=f"c1_{name}", track="C1", task="T1", kind="scm",
            graph=GraphSpec(edges=tuple(EDGES[name]), S=tuple(COND[name]),
                            expected_verdict=EXPECTED_VERDICT[name]),
            ground_truth=GroundTruth(
                task="T1", verdict=EXPECTED_VERDICT[name],
                scalars={"effect": eff, "obs_assoc": assoc},
                exact=True,
                notes=("Definition 1 scored interventional-vs-OBSERVATIONAL; "
                       "do(1)-vs-do(0) mis-scores front-door as epiphenomenal.")),
            run=functools.partial(_run_scm, name), script="experiments/c1_graph_certificates.py")
    return scen


def _substrate_scenarios():
    return {
        "c2_cnet": Scenario(
            name="c2_cnet", track="C2", task="T2", kind="cnet",
            ground_truth=GroundTruth(
                task="T2",
                scalars={"band_collective": 0.24, "band_labeled": 33.1},
                ordering=("band_labeled >> band_collective (~140x): do(W) pins a "
                          "collective readout but leaves a labeled-line reader free",),
                notes="fat-handedness of a constituted do(W)."),
            script="experiments/c2_fat_handed.py"),
        "c3_cnet_e3": Scenario(
            name="c3_cnet_e3", track="C3", task="T3", kind="cnet",
            ground_truth=GroundTruth(
                task="T3",
                scalars={"doW_labeled": 33.1, "doW_collective": 0.24,
                         "doTheta": 0.014, "e3_latency_slope": 1.00},
                ordering=("do(theta) << do(W): the generating parameter is the "
                          "well-posed handle",),
                notes="do(theta) realization band ~0 vs do(W) fat-handed."),
            script="experiments/c3_do_theta.py"),
        "c4_e3_cnet": Scenario(
            name="c4_e3_cnet", track="C4", task="T4", kind="e3",
            ground_truth=GroundTruth(
                task="T4",
                scalars={"do_tau_timing": 1.00, "do_tau_identity": 0.00,
                         "do_route_timing": 0.06, "do_route_identity": 1.00,
                         "suff_collective": 1.03, "suff_labeled": 0.11},
                ordering=("outcome matrix is diagonal: do(tau)->timing, do(route)->identity",
                          "macro-sufficiency collective >> labeled (causal emergence)"),
                notes="outcome-relativity; the causal role is (handle,outcome)-relative."),
            script="experiments/c4_outcome_relativity.py"),
        "c5_spiral": Scenario(
            name="c5_spiral", track="C5", task="T5", kind="spiral",
            ground_truth=GroundTruth(
                task="T5",
                scalars={"band_center": 6.2, "band_tracked": 1.0, "band_global": 2.6,
                         "routing_center": 0.55, "routing_tracked": 0.78},
                ordering=("fat-hand band center > global > tracked: chirality is "
                          "well-posed only read topologically",),
                notes="do(chi) fat-handed at a fixed locus."),
            script="experiments/c5_do_chirality.py"),
        "c6_spiral": Scenario(
            name="c6_spiral", track="C6", task="T5", kind="spiral",
            ground_truth=GroundTruth(
                task="T5",
                scalars={"band_theta_chi": 0.0, "switch_intact": 0.85,
                         "switch_ablated": 0.52, "single_intact": 0.90,
                         "single_ablated": 0.89},
                ordering=("do(theta_chi) well-posed (0 sigma) for every reader",
                          "ablating persistence collapses switching, spares single-rule"),
                notes="generative handle well-posed; persistent core necessary."),
            script="experiments/c6_do_theta_chi.py"),
        "c7_spiral": Scenario(
            name="c7_spiral", track="C7", task="T5", kind="spiral",
            ground_truth=GroundTruth(
                task="T5",
                scalars={"chi_rule": 1.00, "chi_content": 0.11,
                         "route_rule": 0.57, "route_content": 1.00,
                         "screen_seedplus_injplus": 0.86, "screen_seedplus_injminus": 0.19},
                ordering=("outcome matrix diagonal-DOMINANT (off-diagonals honest)",
                          "O_rule tracks injected chi, independent of nucleation seed"),
                notes="outcome-relativity + theta->chi->B mediation/screening."),
            script="experiments/c7_outcome_relativity.py"),
    }


REGISTRY = {**_c1_scenarios(), **_substrate_scenarios()}


def _scenario_to_dict(sc):
    d = {
        "name": sc.name, "track": sc.track, "task": sc.task, "kind": sc.kind,
        "seed": sc.seed, "executable": sc.executable, "script": sc.script,
        "ground_truth": {
            "task": sc.ground_truth.task, "verdict": sc.ground_truth.verdict,
            "scalars": sc.ground_truth.scalars,
            "ordering": list(sc.ground_truth.ordering),
            "exact": sc.ground_truth.exact, "notes": sc.ground_truth.notes,
        },
    }
    if sc.graph is not None:
        d["graph"] = {"edges": [list(e) for e in sc.graph.edges],
                      "S": list(sc.graph.S),
                      "expected_verdict": sc.graph.expected_verdict}
    return d


def export_ground_truth(path=_OUT):
    """Write scenarios.json + ground_truth.json. Deterministic / idempotent."""
    os.makedirs(path, exist_ok=True)
    scenarios = [_scenario_to_dict(REGISTRY[k]) for k in sorted(REGISTRY)]
    with open(os.path.join(path, "scenarios.json"), "w") as f:
        json.dump(scenarios, f, indent=1, sort_keys=True)
    rows = []
    for k in sorted(REGISTRY):
        sc = REGISTRY[k]
        gt = sc.ground_truth
        if gt.verdict is not None:
            rows.append({"task": gt.task, "scenario": sc.name,
                         "quantity": "verdict", "value": gt.verdict,
                         "exact": gt.exact})
        for q, v in gt.scalars.items():
            rows.append({"task": gt.task, "scenario": sc.name,
                         "quantity": q, "value": v, "exact": gt.exact})
    with open(os.path.join(path, "ground_truth.json"), "w") as f:
        json.dump(rows, f, indent=1, sort_keys=True)
    return path
