"""
E6 - Emergent categories (Horde / GVF readout).

Operationalises the claim that "memory / attention / executive are EMERGENT
CATEGORIES, not modules." See docs/learning_experiments.md section 5, E6.

We build ONE frozen Greenberg-Hastings substrate whose single node-phase stream
simultaneously carries three dynamical motifs (three non-interacting regions of one
weight matrix, sharing nothing but the one phi vector):

  (M) MEMORY   : two stimulus-specific directed rings (the E2 mechanism). Each
                 epoch a random stimulus is cued; with some probability the ring is
                 ignited in a DYING mode (tau >= L) so it dies ~L steps later ->
                 the held-stimulus lifetime varies epoch to epoch.
  (A) ATTENTION: an excitable chain (the E4 arena). Each episode two waves ignite
                 the two ends with a random bias/jitter, travel inward and
                 annihilate; whichever captures the centre wins -> the winner varies
                 and must be FORECAST before the collision resolves.
  (X) EXECUTIVE: a slow rotating ring (a persistent loop) whose phase advances one
                 node/step with period = the context-switch interval, so an IMMINENT
                 switch is encoded in the loop's phase.

Onto this ONE frozen substrate we attach three linear GVF demons, each a tuple
<cumulant c, continuation gamma> learned by linear TD from the SAME feature vector
(the whole-network active mask + bias). Only the cumulant/gamma differ per demon:

  memory demon    : c = 1 while the current stimulus ring is alive, gamma=0.9
                    -> value = discounted future "still-remembering" time.
  attention demon : predict the eventual centre-capture winner (terminal target),
                    gamma=1 within an episode -> value = P(stream 0 wins | state).
  executive demon : c = 1 at a context-switch step, gamma=0.8
                    -> value = discounted proximity to the next switch.

Claims tested:
  1. All three demons predict their target well above a mean-predictor baseline,
     from ONE unchanged substrate read by identical linear probes.
  2. The three learned weight vectors are nearly ORTHOGONAL (low cosine similarity)
     -> genuinely distinct questions, not one signal relabelled.
  3. DISCRIMINATOR: a generic full-substrate demon does about as well as an ORACLE
     demon restricted to only its own region, and a demon restricted to the OTHER
     regions fails -> no function-specific sub-network is required; the "function"
     lives in the probe, not the machine.

Outputs
-------
docs/figures/e6_horde.png   : prediction quality (full/own/other), traces, weight overlap
result/e6/e6_data.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_net import Network

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e6")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

ACT = 2
# memory rings
L_MEM = 16
TAU_MEM_LIVE = 12         # < L -> sustains (stimulus held)
TAU_MEM_DIE = 18          # >= L -> dies ~L steps after cue
T_MEM = 44                # memory epoch length
P_DISRUPT = 0.45          # fraction of epochs ignited in the dying mode
# attention chain
L_ATT = 31
CENTER = L_ATT // 2
TAU_ATT = 6
T_ATT = 40                # attention episode length
# executive slow ring
L_EXE = 20
TAU_EXE = 16              # < L -> sustains; single pulse period ~ L_EXE
T_SWITCH = L_EXE          # switch interval aligned to the loop period


def build_substrate(seed=0):
    """One weight matrix with three non-interacting regions; returns W, roles, tau."""
    m0 = 0                                   # memory ring 0
    m1 = m0 + L_MEM                          # memory ring 1
    a0 = m1 + L_MEM                          # attention chain
    x0 = a0 + L_ATT                          # executive ring
    N = x0 + L_EXE
    W = np.zeros((N, N))
    tau = np.full(N, ACT + 1, dtype=np.int64)

    mem = [np.arange(m0, m0 + L_MEM), np.arange(m1, m1 + L_MEM)]
    for rg in mem:                           # two directed memory rings
        for a in range(L_MEM):
            W[rg[(a + 1) % L_MEM], rg[a]] = 1.0
        tau[rg] = TAU_MEM_LIVE

    att = np.arange(a0, a0 + L_ATT)          # attention chain (bidirectional)
    for i in range(L_ATT - 1):
        W[att[i], att[i + 1]] = 1.0
        W[att[i + 1], att[i]] = 1.0
    tau[att] = TAU_ATT

    exe = np.arange(x0, x0 + L_EXE)          # executive slow ring (directed)
    for a in range(L_EXE):
        W[exe[(a + 1) % L_EXE], exe[a]] = 1.0
    tau[exe] = TAU_EXE

    roles = {"mem": mem, "att": att, "exe": exe, "N": N,
             "center": a0 + CENTER, "att0": a0}
    return W, roles, tau


def run_world(seed=0, T=6000):
    """Drive the frozen substrate and record the phase stream + ground-truth signals.

    Returns
    -------
    mask   : (T, N) float   whole-network active mask each step (the demon features)
    signals: dict of ground-truth arrays for the three demons.
    """
    W, roles, tau = build_substrate(seed)
    theta = np.ones(W.shape[0])
    net = Network(W, act=ACT, pas=0, theta=theta, p_s=0.0, seed=seed)
    net.tau = tau.copy()
    net.phi[:] = 0
    rng = np.random.default_rng(seed + 3)
    N = roles["N"]
    mem, att, exe = roles["mem"], roles["att"], roles["exe"]

    # ignite the executive slow ring once (persistent rotation)
    net.phi[exe[0]] = 1

    mask = np.zeros((T, N), np.float32)
    mem_alive = np.zeros(T)          # is the current epoch's stimulus ring alive
    att_winner = np.full(T, -1)      # eventual winner of the current episode (0/1)
    att_forecast = np.zeros(T, bool) # True in the pre-collision forecast window
    exe_switch = np.zeros(T)         # 1 at a context-switch step

    cur_stim = 0
    ep_att_start = -T_ATT
    winner = -1
    collide_t = 0
    s0 = s1 = -1
    for t in range(T):
        # --- MEMORY epoch management ---
        if t % T_MEM == 0:
            for rg in mem:
                net.phi[rg] = 0
            cur_stim = int(rng.integers(2))
            dying = rng.random() < P_DISRUPT
            net.tau[mem[cur_stim]] = TAU_MEM_DIE if dying else TAU_MEM_LIVE
            net.phi[mem[cur_stim][0]] = 1        # single seed -> rotating pulse
        # --- ATTENTION episode management ---
        if t % T_ATT == 0:
            net.phi[att] = 0
            bias = int(rng.integers(-3, 4))
            j0, j1 = rng.integers(0, 2, size=2)
            s0 = t + max(0, -bias) + int(j0)
            s1 = t + max(0, bias) + int(j1)
            winner = 0 if s0 < s1 else (1 if s1 < s0 else int(rng.integers(2)))
            collide_t = max(s0, s1) + CENTER      # approx centre-capture time
            ep_att_start = t
        # --- EXECUTIVE switch schedule ---
        if t % T_SWITCH == 0 and t > 0:
            exe_switch[t] = 1.0

        # drive the attention ends at their ignition times
        drive = np.zeros(N, bool)
        if t == s0:
            drive[att[0]] = True
        if t == s1:
            drive[att[-1]] = True
        net.step(drive if drive.any() else None)

        mask[t] = net.active_mask().astype(np.float32)
        mem_alive[t] = 1.0 if net.active_mask()[mem[cur_stim]].any() else 0.0
        att_winner[t] = winner
        # forecast window: after both ignitions, before centre capture
        att_forecast[t] = (t >= max(s0, s1)) and (t < collide_t)

    signals = {"mem_alive": mem_alive, "att_winner": att_winner,
               "att_forecast": att_forecast, "exe_switch": exe_switch,
               "roles": roles}
    return mask, signals


# ----------------------------------------------------------------------------
# Linear GVF demons (TD) reading a shared feature vector.
# ----------------------------------------------------------------------------

def features(mask, feat_mask=None):
    """Feature matrix = active mask (optionally region-masked) + bias column."""
    X = mask.copy()
    if feat_mask is not None:
        X = X * feat_mask[None, :]
    return np.concatenate([X, np.ones((len(X), 1), np.float32)], axis=1)


def td_continuing(X, cumulant, gamma, alpha=0.02, n_pass=6):
    """Linear TD(0) for a continuing prediction demon. Returns weight vector."""
    w = np.zeros(X.shape[1])
    for _ in range(n_pass):
        for t in range(len(X) - 1):
            v = w @ X[t]
            delta = cumulant[t] + gamma * (w @ X[t + 1]) - v
            w += alpha * delta * X[t]
    return w


def mc_forecast(X, live, winner, lam=1.0):
    """Episodic forecasting GVF (attention winner), learned by Monte-Carlo.

    The attention demon is an EPISODIC terminal-outcome predictor: at each step of
    the pre-collision forecast window (the states where the question is live) it
    predicts the eventual centre-capture winner, whose realised value is known at
    the episode's resolution. The natural GVF estimator for a terminal outcome is
    Monte-Carlo -- regress the value on the realised return (here the winner label)
    over the live states -- as opposed to the online TD(0) used for the CONTINUING
    memory/executive demons. Both are standard Horde/GVF learning rules; the choice
    follows the GVF type (episodic vs continuing), not the substrate. Ridge
    regularisation (lam) keeps the region-masked control scopes well-posed."""
    v = live & (winner >= 0)
    Xt, yt = X[v], winner[v].astype(float)
    return np.linalg.solve(Xt.T @ Xt + lam * np.eye(X.shape[1]), Xt.T @ yt)


def mc_return(cumulant, gamma, horizon=60):
    """Truncated discounted Monte-Carlo return for evaluation."""
    T = len(cumulant)
    G = np.zeros(T)
    disc = gamma ** np.arange(horizon)
    for t in range(T):
        seg = cumulant[t:t + horizon]
        G[t] = np.sum(disc[:len(seg)] * seg)
    return G


def r2(pred, truth):
    var = np.var(truth)
    return 1.0 - np.mean((truth - pred) ** 2) / var if var > 1e-9 else 0.0


def evaluate(seed=0):
    mask, sig = run_world(seed=seed)
    roles = sig["roles"]
    N = roles["N"]
    T = len(mask)
    split = T // 2
    # region masks (over the N node features; bias handled inside features())
    reg = {"mem": np.zeros(N, bool), "att": np.zeros(N, bool), "exe": np.zeros(N, bool)}
    for rg in roles["mem"]:
        reg["mem"][rg] = True
    reg["att"][roles["att"]] = True
    reg["exe"][roles["exe"]] = True

    out = {}
    for demon in ["memory", "attention", "executive"]:
        res = {}
        for scope in ["full", "own", "other"]:
            reg_of = {"memory": "mem", "attention": "att", "executive": "exe"}[demon]
            if scope == "full":
                fmask = None
            elif scope == "own":
                fmask = reg[reg_of]
            else:  # other
                fmask = ~reg[reg_of]
            Xtr = features(mask[:split], fmask)
            Xte = features(mask[split:], fmask)
            if demon == "memory":
                w = td_continuing(Xtr, sig["mem_alive"][:split], gamma=0.9)
                G = mc_return(sig["mem_alive"], 0.9)[split:]
                V = Xte @ w
                res[scope] = {"metric": r2(V, G), "kind": "R2", "V": V, "G": G}
            elif demon == "executive":
                w = td_continuing(Xtr, sig["exe_switch"][:split], gamma=0.8)
                G = mc_return(sig["exe_switch"], 0.8)[split:]
                V = Xte @ w
                res[scope] = {"metric": r2(V, G), "kind": "R2", "V": V, "G": G}
            else:  # attention -> classification accuracy in forecast window
                w = mc_forecast(Xtr, sig["att_forecast"][:split], sig["att_winner"][:split])
                V = Xte @ w
                fc = sig["att_forecast"][split:]
                winner = sig["att_winner"][split:]
                valid = fc & (winner >= 0)
                acc = np.mean((V[valid] > 0.5).astype(int) == winner[valid]) if valid.any() else np.nan
                res[scope] = {"metric": acc, "kind": "acc", "V": V,
                              "G": winner.astype(float), "valid": valid}
            res[scope]["w"] = w
        out[demon] = res
    return out, mask, sig


def main():
    n_seeds = 5
    demons = ["memory", "attention", "executive"]
    scopes = ["full", "own", "other"]
    metric = {d: {s: [] for s in scopes} for d in demons}
    cos_all = []
    example = None
    for seed in range(n_seeds):
        out, mask, sig = evaluate(seed=seed)
        for d in demons:
            for s in scopes:
                metric[d][s].append(out[d][s]["metric"])
        # inter-demon weight cosine similarity (full-substrate weights, drop bias)
        ws = np.array([out[d]["full"]["w"][:-1] for d in demons])
        ws = ws / (np.linalg.norm(ws, axis=1, keepdims=True) + 1e-9)
        cos_all.append(ws @ ws.T)
        if seed == 0:
            example = out
    cos = np.mean(cos_all, axis=0)

    print("E6 - Horde / GVF readout on one frozen substrate")
    for d in demons:
        kind = example[d]["full"]["kind"]
        line = f"  {d:10s} [{kind}]:"
        for s in scopes:
            line += f"  {s}={np.nanmean(metric[d][s]):.2f}"
        print(line)
    chance = {"memory": "baseline R2=0", "attention": "chance acc=0.50",
              "executive": "baseline R2=0"}
    print("  (baselines:", ", ".join(f"{d} {chance[d]}" for d in demons), ")")
    print("  inter-demon weight cosine similarity (off-diagonal, want ~0):")
    print(f"    mem-att={cos[0,1]:+.2f}  mem-exe={cos[0,2]:+.2f}  att-exe={cos[1,2]:+.2f}")

    # ---------------- figure ----------------
    fig = plt.figure(figsize=(15, 9))
    colors = {"memory": "steelblue", "attention": "crimson", "executive": "seagreen"}

    # (1) prediction quality: full vs own vs other, per demon
    ax = fig.add_subplot(2, 3, 1)
    xg = np.arange(len(demons))
    wbar = 0.26
    for j, s in enumerate(scopes):
        vals = [np.nanmean(metric[d][s]) for d in demons]
        errs = [np.nanstd(metric[d][s]) / np.sqrt(n_seeds) for d in demons]
        ax.bar(xg + (j - 1) * wbar, vals, wbar, yerr=errs, capsize=3,
               label=s, color=["#4c72b0", "#55a868", "#c44e52"][j])
    ax.axhline(0.5, ls=":", color="gray", alpha=0.6)
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xticks(xg)
    ax.set_xticklabels([d + "\n(acc)" if d == "attention" else d + "\n(R2)" for d in demons],
                       fontsize=8)
    ax.set_ylabel("prediction quality")
    ax.set_title("Full-substrate probe ~ own-region oracle >> other-region\n"
                 "(no function-specific wiring needed)")
    ax.legend(fontsize=8)

    # (2) weight cosine-similarity matrix
    ax = fig.add_subplot(2, 3, 2)
    im = ax.imshow(cos, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(3)); ax.set_yticks(range(3))
    ax.set_xticklabels(["mem", "att", "exe"]); ax.set_yticklabels(["mem", "att", "exe"])
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{cos[i,j]:.2f}", ha="center", va="center",
                    color="k" if abs(cos[i, j]) < 0.5 else "w", fontsize=9)
    ax.set_title("Inter-demon weight similarity\n(off-diagonal ~0 = distinct questions)")
    fig.colorbar(im, ax=ax, fraction=0.046)

    # (3) memory demon trace
    ax = fig.add_subplot(2, 3, 4)
    V, G = example["memory"]["full"]["V"], example["memory"]["full"]["G"]
    sl = slice(0, 300)
    ax.plot(G[sl], color="gray", label="true return")
    ax.plot(V[sl], color=colors["memory"], label="demon V")
    ax.set_title(f"Memory demon: discounted 'still-remembering'\nR2={np.nanmean(metric['memory']['full']):.2f}")
    ax.set_xlabel("test step"); ax.legend(fontsize=8)

    # (4) executive demon trace
    ax = fig.add_subplot(2, 3, 5)
    V, G = example["executive"]["full"]["V"], example["executive"]["full"]["G"]
    ax.plot(G[sl], color="gray", label="true return")
    ax.plot(V[sl], color=colors["executive"], label="demon V")
    ax.set_title(f"Executive demon: proximity to switch\nR2={np.nanmean(metric['executive']['full']):.2f}")
    ax.set_xlabel("test step"); ax.legend(fontsize=8)

    # (5) attention demon: predicted P(stream0) vs winner in forecast windows
    ax = fig.add_subplot(2, 3, 6)
    V = example["attention"]["full"]["V"]
    winner = example["attention"]["full"]["G"]
    valid = example["attention"]["full"]["valid"]
    idx = np.where(valid)[0][:200]
    ax.scatter(idx, V[idx], c=winner[idx], cmap="bwr", s=10, alpha=0.7)
    ax.axhline(0.5, ls=":", color="k")
    ax.set_title(f"Attention demon: forecast P(stream0 wins)\n(colour=true winner) acc={np.nanmean(metric['attention']['full']):.2f}")
    ax.set_xlabel("test step (forecast windows)"); ax.set_ylabel("demon V")

    fig.suptitle("E6: three GVF demons, one frozen substrate, identical linear probes "
                 "-- memory / attention / executive are questions, not modules",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(FIGDIR, "e6_horde.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "e6_data.npz"),
             cos=cos,
             **{f"{d}_{s}": np.array(metric[d][s]) for d in demons for s in scopes})
    print("wrote", os.path.join(DATADIR, "e6_data.npz"))


if __name__ == "__main__":
    main()
