"""Track 3b follow-on — does the E1 learning dissociation hold off the small-world
hidden medium?

E1's recurrent hidden reservoir is a small-world ring (`layered_graph`, hh_topo=
'smallworld'). 3b showed the *raw* excitable dynamics generalise; this asks whether
the *learned* result does: swap the hidden medium to an ordered ring and a
random-geometric graph (matched mean degree ~ hh_k=6), and re-run E1 conditioning
at n=50. Headline: does reward-driven routing (Line A) still learn the
stimulus→action map, and does the A-vs-B dissociation survive, on each medium?

Reuses E1's config, trial, and trial count verbatim (only the H->H topology
changes). The 'smallworld' arm reproduces E1's committed behaviour as a check.
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from ghca_learn import layered_graph, GHLearner  # noqa: E402
import e1_conditioning as e1  # noqa: E402
import ghca_stats as st  # noqa: E402

N = int(os.environ.get("STATS_N", "50"))
TOPOS = ["smallworld", "ring", "rgg"]
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def hidden_degree(W, roles):
    h = roles["hidden"]
    sub = W[np.ix_(h, h)] > 0
    return float(sub.sum(1).mean())


def run_topo(line, seed, hh_topo):
    W, plastic, roles = layered_graph(seed=seed, w_hm=e1.CFG["w_hm"],
                                      w_hh=e1.CFG["w_hh"], hh_topo=hh_topo)
    net = GHLearner(W, plastic, roles, line=line, act=e1.CFG["act"], pas=e1.CFG["pas"],
                    theta=e1.CFG["theta"], p_s=e1.CFG["p_s"], eta_w=e1.CFG["eta_w"],
                    eta_tau=e1.CFG["eta_tau"], seed=seed + 100)
    rng = np.random.default_rng(seed + 7)
    R = np.zeros(e1.N_TRIALS)
    for t in range(e1.N_TRIALS):
        R[t], _ = e1.trial(net, int(rng.integers(net.roles["K"])), rng)
    return float(R[-120:].mean()), (hidden_degree(W, roles) if seed == 0 else None)


def main():
    groups = {}
    degs = {}
    for topo in TOPOS:
        for line in ("A", "B"):
            vals = np.zeros(N)
            for s in range(N):
                acc, deg = run_topo(line, s, topo)
                vals[s] = acc
                if deg is not None and line == "A":
                    degs[topo] = deg
            groups[f"{topo}_{line}"] = vals
            print("  " + st.fmt_row(st.summarise(f"{topo}:{line}", vals)), flush=True)
        eff = st.effect_size(groups[f"{topo}_A"], groups[f"{topo}_B"])
        print(f"  [{topo}] hidden degree={degs.get(topo, float('nan')):.1f}  "
              f"A-vs-B Cohen d={eff:.2f}", flush=True)

    np.savez(os.path.join(OUT, "e1_topology_port.npz"), n=N,
             degrees=np.array([degs[t] for t in TOPOS]),
             **{f"seed_{k}": v for k, v in groups.items()})
    rows = {k: st.summarise(k, v) for k, v in groups.items()}
    with open(os.path.join(OUT, "e1_topology_port.json"), "w") as f:
        json.dump({"n": N, "degrees": degs, "rows": rows}, f, indent=2, default=float)
    print("wrote e1_topology_port.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
