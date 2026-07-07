"""
E3 (composition study) - decomposing the A+B interference.

E3 found that jointly training Line A (routing -> identity) and Line B
(timescale -> timing) under a SINGLE shared scalar reward interferes: A+B is
worse than either line alone. C4 showed identity and timing are causally
orthogonal, and the discrete-diffusion analysis (Casado Noguerales et al., 2026)
shows the optimum can be coordinate-free while the parameterization governs
learnability. Both point to FACTORED CREDIT. This experiment decomposes the
interference into its sources and tests each fix:

  1. Reward-conflation. A single scalar r = 0.5*id + 0.5*time credits both lines
     for a trial that may be right on one axis and wrong on the other.
     FIX: factored credit -- delta_A = id-error, delta_B = timing-error, each
     line credited only by the component it controls (each still an
     environment-emitted scalar; no per-node teaching).
  2. Dynamical non-stationarity. Even with factored credit, Line B's ongoing
     timescale exploration keeps shifting WHEN the response fires, so Line A's
     routing target never settles.
     FIX: a slow-first curriculum -- learn the slow variable tau first (E0: tau
     is a slow, near-static variable), then freeze it and learn routing.
  3. Substrate resonance (a caveat, not a credit-assignment issue). In the
     reduced gate-metronome substrate, identity learnability is a jagged
     (aliasing) function of the gate tau (panel C). We therefore evaluate at a
     resonance-compatible operating point (target latency chosen so B's tau
     lands in an identity-learnable zone).

Conditions (channel/identity accuracy and latency error, several seeds):
  A only | B only | AB shared | AB factored | AB factored + curriculum

Outputs
-------
docs/figures/e3_factored.png
result/e3/e3_factored.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))
import e3_timed_response as e3

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e3")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

# resonance-compatible operating point: latency 15-17 -> gate tau in [17,19],
# an identity-learnable zone of the substrate (see panel C).
e3.TARGET_LAT = 16
e3.LAT_TOL = 1
N_SEEDS = 5


def trial(net, roles, x, rng, mode, base, learn=True):
    net.reset_traces()
    net.phi[:] = 0
    if learn:
        net.perturb_tau()
    net.phi[roles["gate"]] = e3.ACT + 1
    first_t, first_ch = -1, -1
    for t in range(e3.TWIN):
        net.phi[roles["src"]] = 1
        net.phi[roles["sensory"][x]] = 1
        net.step_learn(None)
        am = net.active_mask()
        mact = [am[roles["motor"][c]].sum() for c in range(roles["A"])]
        if first_t < 0 and sum(mact) > 0:
            first_t = t
            first_ch = int(np.argmax(np.array(mact) + 1e-6 * rng.standard_normal(roles["A"])))
    ch_ok = 1.0 if first_ch == x else 0.0
    lat_ok = 1.0 if (first_t >= 0 and abs(first_t - e3.TARGET_LAT) <= e3.LAT_TOL) else 0.0
    if learn:
        # order-parameter critics (state-dependent baselines, as in E3), one per
        # factor. Identity and timing carry equal 0.5 weight, so factoring changes
        # only WHICH line sees WHICH component, not the magnitudes.
        f = net.features()
        r_id, r_time = 0.5 * ch_ok, 0.5 * lat_ok
        if mode == "shared":
            d = (r_id + r_time) - base["joint"] @ f
            base["joint"] += net.alpha_v * d * f
            net.learn(d, d)
        else:  # factored
            dA = r_id - base["id"] @ f
            dB = r_time - base["time"] @ f
            base["id"] += net.alpha_v * dA * f
            base["time"] += net.alpha_v * dB * f
            net.learn(dA, dB)
    return ch_ok, first_t


def _fresh_base():
    return {"id": np.zeros(3), "time": np.zeros(3), "joint": np.zeros(3)}


def _evaluate(net, roles, rng):
    net.tau[roles["gate"]] = net.tau_scalar
    net._eps_s = 0.0
    ch, lat = [], []
    b = _fresh_base()
    for _ in range(120):
        c, t = trial(net, roles, int(rng.integers(2)), rng, "factored", b, learn=False)
        ch.append(c)
        if t >= 0:
            lat.append(t)
    return np.mean(ch), (abs(np.mean(lat) - e3.TARGET_LAT) if lat else np.nan)


def run_single(line, mode, seed, ntr=1800):
    net, roles = e3.make(line, seed)
    rng = np.random.default_rng(seed + 3)
    base = _fresh_base()
    for _ in range(ntr):
        trial(net, roles, int(rng.integers(2)), rng, mode, base, learn=True)
    return _evaluate(net, roles, rng)


def run_curriculum(seed, n1=1500, n2=1000):
    net, roles = e3.make("AB", seed)
    rng = np.random.default_rng(seed + 3)
    base = _fresh_base()
    net.line = "B"                                    # phase 1: learn slow tau
    for _ in range(n1):
        trial(net, roles, int(rng.integers(2)), rng, "factored", base, learn=True)
    net.line = "A"; net.tau_sigma = 0.0               # phase 2: freeze tau, learn routing
    net.tau[roles["gate"]] = net.tau_scalar
    for _ in range(n2):
        trial(net, roles, int(rng.integers(2)), rng, "factored", base, learn=True)
    return _evaluate(net, roles, rng)


def resonance_map(taus, seeds=2, ntr=800):
    """Identity accuracy vs fixed gate tau (the substrate aliasing artifact)."""
    acc = np.zeros(len(taus))
    for i, tau in enumerate(taus):
        vals = []
        for s in range(seeds):
            net, roles = e3.make("A", s)
            net.tau[roles["gate"]] = tau; net.tau_scalar = tau
            rng = np.random.default_rng(s + 3)
            base = _fresh_base()
            for _ in range(ntr):
                trial(net, roles, int(rng.integers(2)), rng, "factored", base, learn=True)
            c, _ = _evaluate(net, roles, rng)
            vals.append(c)
        acc[i] = np.mean(vals)
    return acc


def main():
    conds = {
        "A only": lambda s: run_single("A", "factored", s),
        "B only": lambda s: run_single("B", "factored", s),
        "AB shared": lambda s: run_single("AB", "shared", s),
        "AB factored": lambda s: run_single("AB", "factored", s),
        "AB factored\n+curriculum": lambda s: run_curriculum(s),
    }
    ch = {k: [] for k in conds}
    le = {k: [] for k in conds}
    for k, fn in conds.items():
        for s in range(N_SEEDS):
            c, l = fn(s)
            ch[k].append(c); le[k].append(l)
        print(f"{k:24s}: channel={np.mean(ch[k]):.2f}  latency_err={np.nanmean(le[k]):.1f}")

    taus = np.arange(10, 25)
    res = resonance_map(taus)
    print("resonance map (identity acc vs gate tau):",
          {int(t): round(a, 2) for t, a in zip(taus, res)})

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    labels = list(conds.keys())
    colors = ["crimson", "steelblue", "gray", "goldenrod", "seagreen"]
    axes[0].bar(labels, [np.mean(ch[k]) for k in labels],
                yerr=[np.std(ch[k]) / np.sqrt(N_SEEDS) for k in labels],
                color=colors, capsize=4)
    axes[0].axhline(0.5, ls="--", color="k", alpha=0.5)
    axes[0].set_ylim(0, 1); axes[0].set_ylabel("channel (identity) accuracy")
    axes[0].set_title("Identity: factoring + curriculum recover it")
    axes[0].tick_params(axis="x", rotation=20, labelsize=8)
    axes[1].bar(labels, [np.nanmean(le[k]) for k in labels],
                yerr=[np.nanstd(le[k]) / np.sqrt(N_SEEDS) for k in labels],
                color=colors, capsize=4)
    axes[1].axhline(e3.LAT_TOL, ls="--", color="k", alpha=0.5, label=f"tol={e3.LAT_TOL}")
    axes[1].set_ylabel("latency error"); axes[1].set_title("Timing")
    axes[1].tick_params(axis="x", rotation=20, labelsize=8); axes[1].legend()
    axes[2].plot(taus, res, "o-", color="purple")
    axes[2].axvspan(17, 19, color="green", alpha=0.15, label="identity-learnable zone")
    axes[2].set_xlabel("fixed gate tau"); axes[2].set_ylabel("identity accuracy")
    axes[2].set_title("Substrate resonance (caveat):\nidentity learnability is jagged in tau")
    axes[2].legend(fontsize=8)
    fig.suptitle("E3 composition study: factored credit + slow-first curriculum "
                 "beat the shared-reward interference")
    fig.tight_layout()
    p = os.path.join(FIGDIR, "e3_factored.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "e3_factored.npz"),
             labels=np.array(labels),
             ch=np.array([ch[k] for k in labels], dtype=object),
             le=np.array([le[k] for k in labels], dtype=object),
             taus=taus, resonance=res)
    print("wrote", os.path.join(DATADIR, "e3_factored.npz"))


if __name__ == "__main__":
    main()
