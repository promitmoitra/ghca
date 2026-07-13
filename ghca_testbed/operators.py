"""Single-source intervention / observation / substrate layer.

Re-exports (does NOT reimplement) the operators and substrate primitives the
C-series already defines, so the testbed has one canonical home for them. See
`ghca_causal.py` and `ghca_net.py`.
"""

import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# intervention operators + wave variables + observation model
from ghca_causal import (  # noqa: E402,F401
    do_S, do_W, do_theta,
    wave_coherence, wave_active_fraction, WAVE_FNS,
    make_observation_mask, read_spikes,
)

# substrate + topology builders + topological-charge readouts
from ghca_net import (  # noqa: E402,F401
    Network, ring, lattice2d, rgg, smallworld, count_phase_singularities,
)

__all__ = [
    "do_S", "do_W", "do_theta",
    "wave_coherence", "wave_active_fraction", "WAVE_FNS",
    "make_observation_mask", "read_spikes",
    "Network", "ring", "lattice2d", "rgg", "smallworld",
    "count_phase_singularities",
]
