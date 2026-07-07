"""
C4 - Outcome-relativity and degeneracy (synthesis cap).

Two closing points for the C-series:

A. Outcome-relativity. The causal role of a variable is a property of the
   (handle, outcome) PAIR, not of the variable alone. On the E3 timed-response
   net we run two well-posed interventions -- do(tau_gate) and do(g_route) (a
   readout-routing parameter) -- and measure two outcomes -- response TIMING
   (latency) and IDENTITY (which channel fires). The causal matrix is diagonal:
   do(tau) drives timing but not identity; do(g_route) drives identity but not
   timing. So "is theta causal?" has no answer without naming the outcome --
   exactly the paper's structure-dependence, and the causal grounding of the E3
   double dissociation.

B. Degeneracy / causal emergence. W = f(S) is a massive many-to-one compression.
   We quantify how much of the behaviourally-relevant information survives at the
   macro level: for a collective code the macro W is (near-)sufficient despite the
   huge S->W degeneracy -- so the macro is the "natural" causal variable (Hoel-
   style causal emergence); for a labeled-line code it is not.

Outputs
-------
docs/figures/c4_outcome_relativity.png
result/c4/c4_data.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from ghca_net import Network, smallworld
from ghca_causal import do_theta, wave_active_fraction
import e3_timed_response as e3

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "c4")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)


# ---------------------------------------------------------------------------
# A. Outcome-relativity on the E3 net: 2x2 causal matrix.
# ---------------------------------------------------------------------------
def e3_net():
    W, plastic, roles, theta, _ = e3.build(w_hm=0.15, w_gate=4.0, theta_m=5.0, seed=0)
    net = Network(W, act=e3.ACT, pas=e3.TAU0 - e3.ACT, theta=theta, p_s=0.0, seed=0)
    relay = np.concatenate([roles["hidden"]] + roles["motor"])
    net.tau[relay] = e3.ACT
    return net, roles, W.copy()


def run_e3(net, roles, W0, tau_gate, g_route, n=40):
    """Return (mean latency, P(channel 0 wins)) under do(tau_gate) and do(g_route)."""
    # do(g_route): scale the two motor channels' readout weights
    net.W[:] = W0
    m0 = roles["motor"][0]; m1 = roles["motor"][1]
    net.W[m0, :] *= (1.0 + g_route)
    net.W[m1, :] *= (1.0 - g_route)
    do_theta(net, tau=tau_gate, nodes=np.array([roles["gate"]]))   # do(tau_gate)
    lat, ch0 = [], 0
    for tr in range(n):
        net.phi[:] = 0
        net.phi[roles["gate"]] = e3.ACT + 1
        x = tr % 2
        ft = -1
        tot = np.zeros(2)                        # integrated channel activity
        for tt in range(80):
            net.phi[roles["src"]] = 1
            net.phi[roles["sensory"][x]] = 1
            net.step(None)
            am = net.active_mask()
            sc = np.array([am[roles["motor"][c]].sum() for c in range(2)])
            tot += sc
            if ft < 0 and sc.sum() > 0:
                ft = tt                          # latency = first-fire time (timing)
        if ft >= 0:
            lat.append(ft)
            ch0 += (int(np.argmax(tot)) == 0)    # identity = integrated winner (routing)
    return (np.mean(lat) if lat else np.nan), ch0 / n


def outcome_matrix():
    net, roles, W0 = e3_net()
    tau_grid = [8, 14, 20, 26]
    g_grid = [-0.4, -0.15, 0.15, 0.4]
    # do(tau) sweep at neutral routing
    lat_tau, id_tau = [], []
    for t in tau_grid:
        la, p0 = run_e3(net, roles, W0, t, 0.0)
        lat_tau.append(la); id_tau.append(p0)
    # do(g_route) sweep at default tau
    lat_g, id_g = [], []
    for g in g_grid:
        la, p0 = run_e3(net, roles, W0, e3.TAU0, g)
        lat_g.append(la); id_g.append(p0)
    # causal effect = range of outcome over the intervention sweep
    eff = np.array([[np.ptp(lat_tau), np.ptp(id_tau)],       # do(tau): timing, identity
                    [np.ptp(lat_g),   np.ptp(id_g)]])        # do(g_route): timing, identity
    return eff, (tau_grid, lat_tau, id_tau), (g_grid, lat_g, id_g)


# ---------------------------------------------------------------------------
# B. Degeneracy / macro-sufficiency on the c-net.
# ---------------------------------------------------------------------------
def decode(X, y):
    clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
    return cross_val_score(clf, X, y, cv=5).mean()


def macro_sufficiency(seed=0, n=1200, N_H=120):
    W = smallworld(N_H, k=6, beta=0.15, seed=seed)
    net = Network(W, act=2, pas=8, theta=1.0, p_s=0.02, seed=seed)
    rng = np.random.default_rng(seed + 1)
    w_coll = 1.0 + 0.1 * rng.standard_normal(N_H)
    w_lab = rng.standard_normal(N_H); w_lab -= w_lab.mean()
    S, Wv, rc, rl = [], [], [], []
    for _ in range(n):
        net.phi[:] = 0
        idx = rng.permutation(N_H)
        net.phi[idx[:N_H // 10]] = rng.integers(1, 3, size=N_H // 10)
        for _ in range(40):
            net.step(None)
        a = net.active_mask().astype(float)
        S.append(a); Wv.append(wave_active_fraction(net))
        rc.append(w_coll @ a); rl.append(w_lab @ a)
    S = np.array(S); Wv = np.array(Wv)[:, None]
    B_coll = (np.array(rc) >= np.median(rc)).astype(int)
    B_lab = (np.array(rl) >= np.median(rl)).astype(int)
    # macro-sufficiency = info(B;W) / info(B;S_full), info ~ (acc - chance)
    def suff(B):
        aW = decode(Wv, B) - 0.5
        aS = decode(S, B) - 0.5
        return max(aW, 0) / max(aS, 1e-6)
    return suff(B_coll), suff(B_lab), N_H


def main():
    eff, (tg, lat_tau, id_tau), (gg, lat_g, id_g) = outcome_matrix()
    # normalize each column (outcome) by its max range for a comparable matrix
    effn = eff / (eff.max(0, keepdims=True) + 1e-9)
    print("outcome-relativity causal matrix (normalized):")
    print(f"            timing   identity")
    print(f"  do(tau)   {effn[0,0]:.2f}     {effn[0,1]:.2f}")
    print(f"  do(route) {effn[1,0]:.2f}     {effn[1,1]:.2f}")
    print(f"  (raw: latency range do(tau)={eff[0,0]:.1f} steps, "
          f"identity range do(route)={eff[1,1]:.2f})")

    suff_c, suff_l, N_H = macro_sufficiency()
    print(f"\nmacro-sufficiency I(B;W)/I(B;S): collective={suff_c:.2f}  labeled={suff_l:.2f}")
    print(f"S->W degeneracy: ~2^{N_H} micro-states compressed onto one scalar W")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    im = axes[0].imshow(effn, cmap="magma", vmin=0, vmax=1)
    axes[0].set_xticks([0, 1]); axes[0].set_xticklabels(["timing\n(latency)", "identity\n(channel)"])
    axes[0].set_yticks([0, 1]); axes[0].set_yticklabels(["do(tau_gate)", "do(g_route)"])
    for i in range(2):
        for j in range(2):
            axes[0].text(j, i, f"{effn[i,j]:.2f}", ha="center", va="center",
                         color="white" if effn[i, j] < 0.5 else "black", fontsize=13)
    axes[0].set_title("A. Outcome-relativity: causal effect of\nhandle on outcome (normalized) -- diagonal")
    fig.colorbar(im, ax=axes[0], fraction=0.046)

    axes[1].bar(["collective", "labeled-line"], [suff_c, suff_l],
                color=["crimson", "steelblue"])
    axes[1].axhline(1.0, ls="--", color="k", alpha=0.5, label="macro fully sufficient")
    axes[1].set_ylabel("macro-sufficiency  I(B;W) / I(B;S)")
    axes[1].set_ylim(0, 1.2)
    axes[1].set_title(f"B. Degeneracy: W compresses ~2^{N_H} states to a scalar;\n"
                      "macro is sufficient only for a collective code")
    axes[1].legend()
    fig.suptitle("C4: the causal role is (handle, outcome)-relative; the wave is "
                 "the natural variable only where behaviour is collective")
    fig.tight_layout()
    p = os.path.join(FIGDIR, "c4_outcome_relativity.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "c4_data.npz"),
             eff=eff, effn=effn, tau_grid=tg, lat_tau=lat_tau, id_tau=id_tau,
             g_grid=gg, lat_g=lat_g, id_g=id_g,
             suff_coll=suff_c, suff_lab=suff_l, N_H=N_H)
    print("wrote", os.path.join(DATADIR, "c4_data.npz"))


if __name__ == "__main__":
    main()
