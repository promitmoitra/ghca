"""ghca_testbed — a synthetic-SCM benchmark for spike–wave causality.

A reusable packaging of the C-series (C0–C7): a substrate where the constitution
W = f(S) is explicit, the three do-operators are available, the causal graph is
KNOWN, and every benchmark task ships a ground-truth answer. Others can score
their epiphenomenality / causal-discovery / causal-emergence methods against it —
impossible in vivo, where the graph is unknown (Causal Hierarchy Theorem).

Quick start:
    import ghca_testbed as tb
    tb.score().  # run the reference impls over the registry
    tb.reproduce()

See docs/causal_testbed.md and docs/causal_testbed_spec.md.
"""

from .operators import (  # noqa: F401
    do_S, do_W, do_theta,
    wave_coherence, wave_active_fraction, WAVE_FNS,
    make_observation_mask, read_spikes,
    Network, ring, lattice2d, rgg, smallworld, count_phase_singularities,
)
from .scm import build_scm, SCM_GRAPHS, EXPECTED_VERDICT  # noqa: F401
from .certificates import (  # noqa: F401
    epiphenomenality_test, theorem1_certificate,
    strata_means, observational_association, EpiResult,
)
from .metrics import (  # noqa: F401
    binarize_median, fat_hand_band, macro_sufficiency, outcome_matrix,
    mediation, OutcomeResult, MediationResult,
)
from .scenario import Scenario, GraphSpec, GroundTruth  # noqa: F401
from .registry import REGISTRY, export_ground_truth  # noqa: F401
from .harness import score, reproduce, ScoreReport, ScoreRow  # noqa: F401

__all__ = [
    "do_S", "do_W", "do_theta", "wave_coherence", "wave_active_fraction",
    "WAVE_FNS", "make_observation_mask", "read_spikes", "Network", "ring",
    "lattice2d", "rgg", "smallworld", "count_phase_singularities",
    "build_scm", "SCM_GRAPHS", "EXPECTED_VERDICT",
    "epiphenomenality_test", "theorem1_certificate", "strata_means",
    "observational_association", "EpiResult",
    "binarize_median", "fat_hand_band", "macro_sufficiency", "outcome_matrix",
    "mediation", "OutcomeResult", "MediationResult",
    "Scenario", "GraphSpec", "GroundTruth", "REGISTRY", "export_ground_truth",
    "score", "reproduce", "ScoreReport", "ScoreRow",
]
