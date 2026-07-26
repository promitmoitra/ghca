"""Topology — cycle-space bound on reentrant-loop (memory) capacity.

E2 builds working memory as a circulating pulse on a directed ring and shows the
duration is tau-controlled: a loop sustains while `tau < L` (loop transit
length) and dies in ~L steps once `tau >= L`. That fixes the duration of ONE
loop but says nothing about HOW MANY independent loops a given topology can hold
at once. This experiment supplies that bound from graph theory.

The method (standard algebraic graph theory; see e.g. Gallier & Quaintance,
"Algebra, Topology, Differential Calculus, and Optimization Theory for CS and
ML", Ch 20-21 on graph Laplacians):

  * For the graph Laplacian `L = D - A` with oriented incidence matrix `B`
    (`L = B B^T`), the CYCLE SPACE -- the space of independent cycles -- has
    dimension the circuit rank (cyclomatic number)

        beta_1 = m - N + c        (m edges, N nodes, c components)

    equivalently `m - rank(L)` or `m - rank(B)`. This is a purely TOPOLOGICAL
    ceiling on the number of independent reentrant loops, computed with no
    dynamics.

  * That ceiling massively overcounts what the DYNAMICS can use, for two
    reasons drawn from E2: (i) a loop only carries a circulating pulse if its
    length exceeds tau, and (ii) two simultaneous circulating waves cannot share
    an edge (they would collide and annihilate in the refractory tail). So the
    usable count is the maximum set of EDGE-DISJOINT cycles of length > tau.
    We report a greedy packing, which is a LOWER bound on that maximum.

Substrate / analysis boundary: this is a pure ANALYSIS of the topology (the
`W` matrices from `ghca_net`), not a dynamical run. It bounds how many loops the
graph could hold; it does not simulate loops coexisting, and K_dyn is therefore
a capacity *bound*, not a measured count. The one dynamical check here is the
E2 gate itself (`tau < L`), re-verified on `ghca_net.Network` so the length
criterion used by the packing is grounded in the actual dynamics.

Note this is a different notion of "capacity" from `continual_capacity.py`
(3c/P4), which measures *representational* capacity of a shared readout. This
one is topological: how many independent loops the medium admits.

Outputs
-------
docs/figures/topology_cycle_capacity.png
result/topology/cycle_capacity.npz
"""

import os
import sys

import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ghca_net import Network, lattice2d, ring, smallworld, rgg  # noqa: E402

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "topology")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

TAUS = np.arange(3, 22)
SEED = 0


def _undirected(W):
    """Symmetrised unweighted support of a (possibly directed) weight matrix."""
    return ((W + W.T) > 0).astype(int)


def cycle_space_dim(W):
    """beta_1 = m - N + c, cross-checked against m - rank(L) and m - rank(B)."""
    A = _undirected(W)
    N = A.shape[0]
    iu = np.triu_indices(N, 1)
    edges = [(i, j) for i, j in zip(*iu) if A[i, j]]
    m = len(edges)
    B = np.zeros((N, m))
    for e, (i, j) in enumerate(edges):
        B[i, e] = 1.0
        B[j, e] = -1.0
    L = np.diag(A.sum(1)) - A.astype(float)
    rank_L = np.linalg.matrix_rank(L)
    c = N - rank_L
    return dict(N=N, m=m, c=c,
                beta1_formula=m - N + c,
                beta1_laplacian=m - rank_L,
                beta1_incidence=m - np.linalg.matrix_rank(B))


def pack_long_cycles(W, tau, max_iter=500, seed=SEED):
    """Greedy edge-disjoint packing of cycles longer than tau.

    Lower bound on the max number of simultaneously-sustainable loops: each
    accepted cycle is longer than tau (E2's sustain gate) and shares no edge
    with any other (no refractory collision). RNG seeded explicitly.
    """
    rng = np.random.default_rng(seed)
    H = nx.from_numpy_array(_undirected(W))
    packed = 0
    for _ in range(max_iter):
        found = None
        nodes = list(H.nodes())
        rng.shuffle(nodes)
        for s in nodes:
            sub = H.subgraph(nx.node_connected_component(H, s))
            longer = [cyc for cyc in nx.cycle_basis(sub) if len(cyc) > tau]
            if longer:
                found = max(longer, key=len)
                break
        if found is None:
            break
        el = [(found[i], found[(i + 1) % len(found)]) for i in range(len(found))]
        H.remove_edges_from([(u, v) for u, v in el if H.has_edge(u, v)])
        packed += 1
    return packed


def verify_e2_gate(L=24, taus=(6, 10, 14, 18, 22, 24, 26, 30), seed=SEED):
    """Re-verify E2's sustain gate on the repo's own Network, directed ring.

    A bare pulse on a BIDIRECTIONAL ring splits and the two fronts annihilate
    head-on, so it never circulates; E2's `two_ring` uses directed rings, and so
    must this check. Returns {tau: persisted?}.
    """
    W = np.zeros((L, L))
    for i in range(L):
        W[i, (i - 1) % L] = 1.0            # directed ring, as in e2 two_ring()
    out = {}
    for tau in taus:
        net = Network(W, act=2, pas=int(tau) - 2, theta=1.0, p_s=0.0, seed=seed)
        net.phi[:] = 0
        net.phi[0] = 1
        A = net.run(400, record=False)["A"]
        out[int(tau)] = bool(A[-30:].mean() > 0)
    return out


def topologies():
    """The repo's own constructors, at comparable size."""
    return {
        "ring (N=60, k=3)":            ring(60, k=3),
        "lattice2d (12x12, r=2)":      lattice2d(12, r=2, periodic=True),
        "smallworld (N=120, k=6)":     smallworld(120, k=6, beta=0.1, seed=SEED),
        "rgg (N=150, radius=0.14)":    rgg(150, radius=0.14, seed=SEED),
    }


def main():
    tops = topologies()

    print("E2 sustain gate on the repo Network (directed ring, L=24):")
    gate = verify_e2_gate()
    for tau, alive in gate.items():
        print(f"   tau={tau:>3}  tau<L={str(tau < 24):>5}  "
              f"{'PERSISTS' if alive else 'dies'}")

    print("\nCycle-space dimension (three equivalent computations must agree):")
    dims = {}
    for name, W in tops.items():
        d = cycle_space_dim(W)
        dims[name] = d
        assert d["beta1_formula"] == d["beta1_laplacian"] == d["beta1_incidence"]
        print(f"   {name:<26} N={d['N']:<5} m={d['m']:<6} c={d['c']}  "
              f"beta_1={d['beta1_formula']}")

    print("\nUsable capacity K_dyn(tau) = greedy edge-disjoint cycles longer than tau:")
    K = {}
    for name, W in tops.items():
        K[name] = np.array([pack_long_cycles(W, int(t)) for t in TAUS])
        print(f"   {name:<26} " + " ".join(
            f"t={t}:{k}" for t, k in zip(TAUS[::4], K[name][::4])))

    # Print the results-doc table verbatim from the data, and the beta_1/K_dyn
    # ratio range, so docs/topology_cycle_capacity.md is transcribed rather than
    # hand-typed (an earlier draft mis-indexed the tau=8 column by hand).
    i3, i8, i21 = [TAUS.tolist().index(t) for t in (3, 8, 21)]
    print("\n-- results-doc table (copy verbatim) --")
    print("| topology | N | m | c | beta_1 | K tau=3 | tau=8 | tau=21 |")
    for name in K:
        d, row = dims[name], K[name]
        print(f"| `{name}` | {d['N']} | {d['m']} | {d['c']} | "
              f"{d['beta1_formula']} | {row[i3]} | {row[i8]} | {row[i21]} |")
    ratios = [dims[n]["beta1_formula"] / K[n][i8] for n in K]
    allr = [dims[n]["beta1_formula"] / k for n in K for k in K[n] if k > 0]
    print(f"beta_1/K_dyn at tau=8: {min(ratios):.1f}x - {max(ratios):.1f}x "
          f"(min over all swept tau: {min(allr):.1f}x)")

    plot(dims, K, gate)

    np.savez(os.path.join(DATADIR, "cycle_capacity.npz"),
             taus=TAUS,
             names=np.array(list(tops.keys())),
             beta1=np.array([dims[n]["beta1_formula"] for n in tops]),
             K=np.array([K[n] for n in tops]),
             gate_taus=np.array(sorted(gate)),
             gate_persist=np.array([gate[t] for t in sorted(gate)]))
    print("\nwrote", os.path.join(DATADIR, "cycle_capacity.npz"))


def plot(dims, K, gate):
    names = list(K)
    cols = ["#6b7280", "#2563eb", "#d97706", "#16a34a"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.7))

    for name, col in zip(names, cols):
        ax1.plot(TAUS, K[name], "-o", color=col, ms=3, lw=1.6,
                 label=name.split(" (")[0])
    ax1.axvspan(4, 18, color="#2563eb", alpha=0.05)
    ax1.set_xlabel(r"local timescale $\tau$")
    ax1.set_ylabel(r"edge-disjoint loops longer than $\tau$  (lower bound)")
    ax1.set_title("Usable loop capacity falls with $\\tau$ on 3 of 4 topologies\n"
                  "(ring is flat: its packed cycles all exceed the $\\tau$ range)",
                  loc="left", fontsize=8)
    ax1.legend(frameon=False, fontsize=6.5, ncol=2, loc="upper right")
    ax1.margins(y=0.12)

    x = np.arange(len(names))
    b1 = [dims[n]["beta1_formula"] for n in names]
    i8 = TAUS.tolist().index(8)
    k8 = [K[n][i8] for n in names]
    ax2.bar(x - 0.2, b1, 0.4, color="#cbd5e1",
            label=r"topological ceiling $\beta_1 = m-N+c$")
    ax2.bar(x + 0.2, k8, 0.4, color=cols, label=r"usable at $\tau=8$")
    for xi, (b, k) in enumerate(zip(b1, k8)):
        ax2.text(xi - 0.2, b + max(b1) * 0.015, str(b), ha="center", fontsize=6)
        ax2.text(xi + 0.2, k + max(b1) * 0.015, str(k), ha="center", fontsize=6,
                 fontweight="bold")
    ax2.set_xticks(x)
    ax2.set_xticklabels([n.split(" (")[0] for n in names], rotation=20,
                        ha="right", fontsize=6.5)
    ax2.set_ylabel("number of independent loops")
    ax2.set_ylim(0, max(b1) * 1.2)
    ax2.set_title("Cycle-space ceiling overcounts usable memory", loc="left",
                  fontsize=8)
    ax2.legend(frameon=False, fontsize=6.5, loc="upper right")

    fig.tight_layout()
    path = os.path.join(FIGDIR, "topology_cycle_capacity.png")
    fig.savefig(path, dpi=110)
    print("wrote", path)


if __name__ == "__main__":
    main()
