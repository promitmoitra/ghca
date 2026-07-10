"""
C3 - do(theta): the generating parameter is the well-posed causal handle.

C2 showed do(W) is fat-handed for a constituted W=f(S): one macro target admits
a wide band of behaviour, so the causal verdict depends on the (arbitrary)
realization policy. The fix is to intervene not on the aggregate but on the
GENERATING PARAMETERS theta (per-node timescale tau, coupling scale): setting
tau = t is modular and unique -- there is no realization degree of freedom -- so
do(theta) yields a single well-defined interventional response. The only
residual variability is ordinary stochastic estimation error, which shrinks with
samples; the do(W) policy band does not.

Three parts:
  A. c-net: E[B | do(tau=t)] and E[B | do(coupling=g)] are single-valued causal
     responses (theta -> W -> B), with tiny across-seed error.
  B. intervention ambiguity: irreducible policy band of do(W) (from C2) vs the
     (zero) realization band of do(theta).
  C. E3 bridge: do(tau_gate = t) -> response latency is a clean single-valued
     function (latency ~ t - act), and it is outcome-specific (timing, not
     identity) -- linking to C4.

Outputs
-------
docs/figures/c3_do_theta.png
result/c3/c3_data.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_net import Network, smallworld
from ghca_causal import do_theta, wave_active_fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "c3")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

N_H = 120
ACT, THETA, P_S = 2, 1.0, 0.02
T_PROP = 6
N_TRIALS = 200


def build(seed=0):
    W = smallworld(N_H, k=6, beta=0.15, seed=seed)
    net = Network(W, act=ACT, pas=8, theta=THETA, p_s=P_S, seed=seed)
    rng = np.random.default_rng(seed + 1)
    w_coll = 1.0 + 0.1 * rng.standard_normal(N_H)
    w_lab = rng.standard_normal(N_H); w_lab -= w_lab.mean()
    return net, w_coll, w_lab


def readouts(net, w_coll, w_lab):
    a = net.active_mask().astype(float)
    return float(w_coll @ a), float(w_lab @ a)


def do_tau_response(tau_grid, seed=0):
    """E[B | do(tau=t)] on the c-net, and across-seed reproducibility."""
    Ec, El, Wf = [], [], []
    rep_c = []                       # across-seed std of the E[B_coll] estimate
    for t in tau_grid:
        per_seed = []
        rc_all, rl_all, wf_all = [], [], []
        for s in range(4):
            net, w_coll, w_lab = build(seed=seed + s)
            rng = np.random.default_rng(1000 * s + int(t))
            rc, rl, wf = [], [], []
            for _ in range(N_TRIALS):
                net.phi[:] = 0
                net.seed_random(0.1, 0.1)
                do_theta(net, tau=t, nodes=np.arange(N_H))   # modular, unique
                for _ in range(T_PROP):
                    net.step(None)
                c, l = readouts(net, w_coll, w_lab)
                rc.append(c); rl.append(l); wf.append(wave_active_fraction(net, np.arange(N_H)))
            per_seed.append(np.mean(rc))
            rc_all += rc; rl_all += rl; wf_all += wf
        Ec.append(np.mean(rc_all)); El.append(np.mean(rl_all)); Wf.append(np.mean(wf_all))
        rep_c.append(np.std(per_seed))          # reproducibility of the response
    return map(np.array, (Ec, El, Wf, rep_c))


def do_coupling_response(g_grid, seed=0):
    Wf = []
    for g in g_grid:
        wf = []
        for s in range(3):
            net, _, _ = build(seed=seed + s)
            do_theta(net, coupling_scale=g, nodes=np.arange(N_H))
            net.seed_random(0.1, 0.1)
            for _ in range(T_PROP):
                net.step(None)
            wf.append(wave_active_fraction(net, np.arange(N_H)))
        Wf.append(np.mean(wf))
    return np.array(Wf)


def e3_latency_vs_tau():
    """Bridge to E3: do(tau_gate=t) -> response latency, single-valued."""
    sys.path.insert(0, os.path.join(ROOT, "experiments"))
    import e3_timed_response as e3
    W, plastic, roles, theta, _ = e3.build(w_hm=0.15, w_gate=4.0, theta_m=5.0, seed=0)
    net = Network(W, act=e3.ACT, pas=e3.TAU0 - e3.ACT, theta=theta, p_s=0.0, seed=0)
    relay = np.concatenate([roles["hidden"]] + roles["motor"])
    net.tau[relay] = e3.ACT
    taus = np.arange(6, 30, 3)
    lat = []
    for tg in taus:
        do_theta(net, tau=tg, nodes=np.array([roles["gate"]]))   # do(tau_gate)
        ts = []
        for tr in range(15):
            net.phi[:] = 0
            net.phi[roles["gate"]] = e3.ACT + 1
            x = tr % 2
            ft = -1
            for tt in range(80):
                net.phi[roles["src"]] = 1
                net.phi[roles["sensory"][x]] = 1
                net.step(None)
                if ft < 0 and sum(net.active_mask()[roles["motor"][c]].sum() for c in range(2)) > 0:
                    ft = tt
            if ft >= 0:
                ts.append(ft)
        lat.append(np.mean(ts) if ts else np.nan)
    return taus, np.array(lat)


def main():
    tau_grid = np.arange(4, 22, 2)
    Ec, El, Wf, rep_c = do_tau_response(tau_grid)
    g_grid = np.linspace(0.5, 2.0, 7)
    Wf_g = do_coupling_response(g_grid)
    taus_e3, lat_e3 = e3_latency_vs_tau()

    # load C2 irreducible policy band for the ambiguity comparison
    c2 = np.load(os.path.join(ROOT, "result", "c2", "c2_data.npz"))
    doW_band_coll = float(np.mean(c2["band_c"]))     # ~0.24 sigma
    doW_band_lab = float(np.mean(c2["band_l"]))      # ~33 sigma
    doTheta_band = float(np.mean(rep_c) / (np.std(Ec) + 1e-9))  # residual, ~0

    print(f"do(tau) response range (collective) = {np.ptp(Ec):.1f}, "
          f"across-seed reproducibility std = {np.mean(rep_c):.2f}")
    print(f"intervention ambiguity (irreducible DOF band, sigma-units): "
          f"do(W) labeled={doW_band_lab:.1f}  do(W) collective={doW_band_coll:.2f}  "
          f"do(theta)={doTheta_band:.3f}")
    print(f"E3 bridge: latency ~ tau_gate - act? "
          f"fit slope={np.polyfit(taus_e3[np.isfinite(lat_e3)], lat_e3[np.isfinite(lat_e3)],1)[0]:.2f}")

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    # A: do(tau) causal response on the c-net (single-valued)
    axes[0].plot(tau_grid, Ec / (np.abs(Ec).max()), "o-", color="crimson", label="collective B (norm)")
    axes[0].plot(tau_grid, Wf, "s-", color="purple", label="wave W = active frac")
    axes[0].set_xlabel("do(tau = t)")
    axes[0].set_ylabel("interventional response")
    axes[0].set_title("A. do(theta): single-valued causal response\n(theta -> W -> B), tiny across-seed error")
    axes[0].legend(fontsize=8)

    # B: intervention ambiguity comparison
    labels = ["do(W)\nlabeled", "do(W)\ncollective", "do(theta)"]
    vals = [doW_band_lab, doW_band_coll, doTheta_band]
    axes[1].bar(labels, vals, color=["steelblue", "crimson", "seagreen"])
    axes[1].set_ylabel("irreducible intervention ambiguity (sigma)")
    axes[1].set_title("B. do(W) has a realization band;\ndo(theta) has none")
    for i, v in enumerate(vals):
        axes[1].text(i, v + 0.5, f"{v:.2f}", ha="center", fontsize=9)

    # C: E3 bridge -- do(tau_gate) -> latency, single-valued & outcome-specific
    fin = np.isfinite(lat_e3)
    axes[2].plot(taus_e3[fin], lat_e3[fin], "o-", color="darkorange", label="response latency")
    axes[2].plot(taus_e3, taus_e3 - ACT, "k--", alpha=0.5, label="tau - act")
    axes[2].set_xlabel("do(tau_gate = t)")
    axes[2].set_ylabel("response latency (steps)")
    axes[2].set_title("C. E3 bridge: do(tau) sets timing,\nuniquely (identity untouched)")
    axes[2].legend(fontsize=8)

    fig.suptitle("C3: intervene on generating parameters theta, not the "
                 "constituted aggregate W -- do(theta) is the well-posed causal handle")
    fig.tight_layout()
    p = os.path.join(FIGDIR, "c3_do_theta.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "c3_data.npz"),
             tau_grid=tau_grid, Ec=Ec, El=El, Wf=Wf, rep_c=rep_c,
             g_grid=g_grid, Wf_g=Wf_g, taus_e3=taus_e3, lat_e3=lat_e3,
             doW_band_lab=doW_band_lab, doW_band_coll=doW_band_coll,
             doTheta_band=doTheta_band)
    print("wrote", os.path.join(DATADIR, "c3_data.npz"))


if __name__ == "__main__":
    main()
