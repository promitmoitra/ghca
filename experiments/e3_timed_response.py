"""
E3 - Timed response (reduced sequence task: identity + latency).

The full E3 in the design doc is multi-element sequence reproduction with order
AND inter-element timing. This is a reduced single-event version that tests the
same spatial-vs-temporal double dissociation: the network must fire the CORRECT
motor channel (identity/order -> spatial routing -> Line A) at a TARGET LATENCY
after go (timing -> a tau-controlled metronome gate -> Line B).

Mechanism: a stimulus is routed S->H->M (plastic, Line A). A gate node is reset
to just-refractory at go and driven tonically, so it first fires ~tau-act steps
later and then every ~tau steps (a refractory-limited metronome; period tracks
tau). Motor fires only in coincidence with a gate beat, so the response latency
= first gate beat ~ tau - act. Line B tunes the gate's (shared) tau to hit the
target latency; Line A tunes which channel fires.

Prediction (double dissociation):
  Line A only  -> channel correct, latency wrong (gate tau fixed at default)
  Line B only  -> latency correct, channel wrong (routing not learned)
  Line A+B     -> both correct

Reward is a single scalar, graded: 0.5*[channel correct] + 0.5*[latency within
tolerance]. (Graded so each line can earn its half; still one scalar, no
per-node teaching.)

Outputs
-------
docs/figures/e3_mechanism.png  : response latency vs gate tau
docs/figures/e3_dissociation.png : channel accuracy and latency error by line
result/e3/e3_data.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_net import Network
from ghca_learn import GHLearner

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e3")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

ACT = 2
TAU0 = 10                    # default gate tau -> default latency ~ TAU0-ACT = 8
TARGET_LAT = 18              # target latency (needs tau ~ 20; != default)
LAT_TOL = 2
TWIN = 70
N_S, N_H, N_M = 10, 80, 10


def build(w_hm=0.15, w_gate=4.0, theta_h=2.0, theta_m=5.0,
          fanin_sh=6, fanin_hm=12, channel_bias=0.85, jitter=0.12, seed=0):
    rng = np.random.default_rng(seed)
    K = A = 2
    src = 0
    s0 = 1
    h0 = s0 + K * N_S
    m0 = h0 + N_H
    gate = m0 + A * N_M
    N = gate + 1
    sensory = [np.arange(s0 + c * N_S, s0 + (c + 1) * N_S) for c in range(K)]
    hidden = np.arange(h0, h0 + N_H)
    motor = [np.arange(m0 + c * N_M, m0 + (c + 1) * N_M) for c in range(A)]
    allsen = np.concatenate(sensory)
    pref = rng.integers(K, size=N_H)
    W = np.zeros((N, N))
    plastic = np.zeros((N, N), dtype=bool)
    theta = np.zeros(N)
    theta[src] = 99.0
    for i, h in enumerate(hidden):
        p = sensory[pref[i]]
        npf = int(round(channel_bias * fanin_sh))
        s = np.unique(np.concatenate([
            rng.choice(p, min(npf, len(p)), replace=False),
            rng.choice(allsen, max(fanin_sh - npf, 0), replace=False)]))
        W[h, s] = 1.0 * (1 + jitter * rng.standard_normal(len(s)))
        plastic[h, s] = True
        theta[h] = theta_h
    for c in range(A):
        for m in motor[c]:
            s = rng.choice(hidden, fanin_hm, replace=False)
            W[m, s] = w_hm * (1 + jitter * rng.standard_normal(fanin_hm))
            plastic[m, s] = True
            theta[m] = theta_m
    W[gate, src] = 5.0
    theta[gate] = 1.0
    for c in range(A):
        W[motor[c], gate] = w_gate
    hidden_all = hidden
    roles = {"sensory": sensory, "motor": motor, "hidden": hidden_all,
             "src": src, "gate": gate, "K": K, "A": A, "N": N, "n_m": N_M}
    tau_mask = np.zeros(N, dtype=bool)
    tau_mask[gate] = True
    return np.clip(W, 0, None), plastic, roles, theta, tau_mask


def make(line, seed, eta_w=0.05, eta_tau=0.5, tau_sigma=2.5):
    W, plastic, roles, theta, tmask = build(seed=seed)
    net = GHLearner(W, plastic, roles, line=line, act=ACT, pas=TAU0 - ACT,
                    theta=theta, p_s=0.0, eta_w=eta_w, eta_tau=eta_tau,
                    tau_min=4, tau_max=40, tau_sigma=tau_sigma, tau_mask=tmask,
                    tau_shared=True, seed=seed + 50)
    # gapless relay: hidden and motor are re-excitable every step, so the
    # response reliably fires at the FIRST gate beat (latency ~ gate_tau - act).
    # Only the gate's (shared) timescale is plastic under Line B.
    relay = np.concatenate([roles["hidden"]] + roles["motor"])
    net.tau[relay] = ACT
    net.tau_base[relay] = ACT
    return net, roles


def run_trial(net, roles, x, rng, learn=True):
    net.reset_traces()
    net.phi[:] = 0
    if learn:
        net.perturb_tau()
    net.phi[roles["gate"]] = ACT + 1        # reset gate to just-refractory at go
    first_t, first_ch = -1, -1
    for t in range(TWIN):
        net.phi[roles["src"]] = 1           # tonic drive
        net.phi[roles["sensory"][x]] = 1    # stimulus clamped through trial
        net.step_learn(None)
        am = net.active_mask()
        mact = [am[roles["motor"][c]].sum() for c in range(roles["A"])]
        if first_t < 0 and sum(mact) > 0:
            first_t = t
            first_ch = int(np.argmax(np.array(mact) + 1e-6 * rng.standard_normal(roles["A"])))
    ch_ok = 1.0 if first_ch == x else 0.0
    lat_ok = 1.0 if (first_t >= 0 and abs(first_t - TARGET_LAT) <= LAT_TOL) else 0.0
    r = 0.5 * ch_ok + 0.5 * lat_ok
    if learn:
        V = net.value()
        net.learn(r - V)
        net.update_critic(r - V, net.features())
    return ch_ok, first_t


def run_condition(line, seed, ntr=1500):
    net, roles = make(line, seed)
    rng = np.random.default_rng(seed + 3)
    for _ in range(ntr):
        run_trial(net, roles, int(rng.integers(2)), rng, learn=True)
    net.tau[roles["gate"]] = net.tau_scalar
    net._eps_s = 0.0
    ch, lat = [], []
    for _ in range(120):
        c, t = run_trial(net, roles, int(rng.integers(2)), rng, learn=False)
        ch.append(c)
        if t >= 0:
            lat.append(t)
    return np.mean(ch), (np.mean(lat) if lat else np.nan), net.tau_scalar


def mechanism(seed=0):
    """Latency vs gate tau with routing strengthened so motor fires at 1st beat."""
    W, plastic, roles, theta, _ = build(seed=seed)
    # strengthen a correct routing so the response reliably fires at the 1st beat
    net = Network(W, act=ACT, pas=TAU0 - ACT, theta=theta, p_s=0.0, seed=seed)
    relay = np.concatenate([roles["hidden"]] + roles["motor"])
    net.tau[relay] = ACT                     # gapless relay (see make())
    taus = np.arange(6, 34, 4)
    lat = []
    rng = np.random.default_rng(0)
    for tau in taus:
        net.tau[roles["gate"]] = tau
        ts = []
        for tr in range(20):
            net.phi[:] = 0
            net.phi[roles["gate"]] = ACT + 1
            x = tr % 2
            ft = -1
            for t in range(TWIN):
                net.phi[roles["src"]] = 1
                net.phi[roles["sensory"][x]] = 1
                net.step(None)
                am = net.active_mask()
                if ft < 0 and sum(am[roles["motor"][c]].sum() for c in range(2)) > 0:
                    ft = t
            if ft >= 0:
                ts.append(ft)
        lat.append(np.mean(ts) if ts else np.nan)
    return taus, np.array(lat)


def main():
    print("E3: mechanism (latency vs gate tau) ...")
    taus, lat = mechanism()

    print("E3: learning (lines A / B / AB) ...")
    lines = ["A", "B", "AB"]
    n_seeds = 5
    ch_acc = {L_: np.zeros(n_seeds) for L_ in lines}
    lat_err = {L_: np.zeros(n_seeds) for L_ in lines}
    tau_l = {L_: np.zeros(n_seeds) for L_ in lines}
    for L_ in lines:
        for s in range(n_seeds):
            c, la, ts = run_condition(L_, s)
            ch_acc[L_][s] = c
            lat_err[L_][s] = abs(la - TARGET_LAT) if not np.isnan(la) else np.nan
            tau_l[L_][s] = ts
        print(f"  line {L_}: channel_acc={ch_acc[L_].mean():.2f} "
              f"latency_err={np.nanmean(lat_err[L_]):.1f} tau*={tau_l[L_].mean():.1f}")

    # mechanism figure
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.plot(taus, lat, "o-", color="purple")
    ax.plot(taus, taus - ACT, "k--", alpha=0.5, label=r"latency $=\tau-a$")
    ax.axhline(TARGET_LAT, ls=":", color="green", label=f"target latency={TARGET_LAT}")
    ax.set_xlabel(r"gate timescale $\tau$")
    ax.set_ylabel("response latency (steps)")
    ax.set_title("E3 mechanism: response latency tracks gate tau")
    ax.legend()
    fig.tight_layout()
    p1 = os.path.join(FIGDIR, "e3_mechanism.png")
    fig.savefig(p1, dpi=110)
    print("wrote", p1)

    # dissociation figure
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    colors = {"A": "crimson", "B": "steelblue", "AB": "seagreen"}
    axes[0].bar(lines, [ch_acc[L_].mean() for L_ in lines],
                yerr=[ch_acc[L_].std() / np.sqrt(n_seeds) for L_ in lines],
                color=[colors[L_] for L_ in lines], capsize=4)
    axes[0].axhline(0.5, ls="--", color="k", alpha=0.5)
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel("channel (identity) accuracy")
    axes[0].set_title("Order/identity - learned by Line A")
    axes[1].bar(lines, [np.nanmean(lat_err[L_]) for L_ in lines],
                yerr=[np.nanstd(lat_err[L_]) / np.sqrt(n_seeds) for L_ in lines],
                color=[colors[L_] for L_ in lines], capsize=4)
    axes[1].axhline(LAT_TOL, ls="--", color="k", alpha=0.5, label=f"tolerance={LAT_TOL}")
    axes[1].set_ylabel("latency error |lat - target|")
    axes[1].set_title("Timing - learned by Line B")
    axes[1].legend()
    fig.suptitle("E3 double dissociation: identity (A) vs timing (B)")
    fig.tight_layout()
    p2 = os.path.join(FIGDIR, "e3_dissociation.png")
    fig.savefig(p2, dpi=110)
    print("wrote", p2)

    np.savez(os.path.join(DATADIR, "e3_data.npz"),
             mech_taus=taus, mech_lat=lat,
             **{f"ch_{L_}": ch_acc[L_] for L_ in lines},
             **{f"laterr_{L_}": lat_err[L_] for L_ in lines},
             **{f"tau_{L_}": tau_l[L_] for L_ in lines})
    print("wrote", os.path.join(DATADIR, "e3_data.npz"))


if __name__ == "__main__":
    main()
