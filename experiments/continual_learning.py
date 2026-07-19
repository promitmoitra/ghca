"""Track 3c / P1 — continual learning on one shared substrate (harness + baselines).

See docs/continual_learning_plan.md. v1 uses a remapping task triple on ONE shared
E1-style substrate (three conflicting stimulus->action mappings, K=A=3) rather than
the heterogeneous E1/E2/E5 — this isolates catastrophic interference from
architecture-swapping and reuses E1's validated conditioning machinery. Tasks are
cyclic relabellings that maximally conflict on the shared readout:

    T0: action = stim            T1: action = stim+1 (mod 3)   T2: action = stim+2

P1 establishes the harness, the continual-learning metrics, and the two "what
adapts" baselines:
  - frozen : the S->H representation is fixed (reservoir); only the shared H->M
             head adapts. Reservoir theory predicts little interference.
  - plastic: S->H and H->M both adapt (the shared representation is rewritten by
             each task) — the hard case where interference is expected.

Both use the current correlational eligibility-trace credit (Line A). The causal
do(theta) credit rule is P2.
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

N = int(os.environ.get("STATS_N", "30"))
N_TRAIN = int(os.environ.get("STATS_NTRAIN", "500"))   # trials per task
K = A = int(os.environ.get("STATS_K", "3"))            # P1 uses K=3; P2 uses K=2
CHANCE = 1.0 / A


def _shift(s):
    return lambda x: (x + s) % K                        # cyclic relabelling task


TASKS = [_shift(s) for s in range(K)]                   # K conflicting remappings
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def make(seed, regime, eta_w=None):
    W, plastic, roles = layered_graph(K=K, A=A, seed=seed,
                                      w_hm=e1.CFG["w_hm"], w_hh=e1.CFG["w_hh"])
    if regime == "frozen":
        # freeze the S->H representation; keep only the H->M head plastic
        h0, h1 = roles["hidden"][0], roles["hidden"][-1] + 1
        plastic = plastic.copy()
        plastic[h0:h1, :] = False          # no S->H updates (rows = hidden nodes)
    net = GHLearner(W, plastic, roles, line="A", act=e1.CFG["act"], pas=e1.CFG["pas"],
                    theta=e1.CFG["theta"], p_s=e1.CFG["p_s"],
                    eta_w=e1.CFG["eta_w"] if eta_w is None else eta_w,
                    eta_tau=e1.CFG["eta_tau"], seed=seed + 100)
    return net


def trial(net, x, target, rng, learn=True):
    """One conditioning trial (mirrors e1.trial) with reward vs an explicit target."""
    net.reset_traces()
    for _ in range(e1.SETTLE):
        net.step_learn(None)
    feats = net.features()
    V = net.value()
    for _ in range(e1.CUE):
        net.step_learn(net.sensory_drive(x))
    sc = np.zeros(net.roles["A"])
    for _ in range(e1.WWIN):
        net.step_learn(None)
        sc += net.motor_scores()
    action = int(np.argmax(sc + 1e-6 * rng.standard_normal(len(sc)))) if sc.sum() > 0 else -1
    r = 1.0 if action == target else 0.0
    if learn:
        delta = r - V
        net.learn(delta)
        net.update_critic(delta, feats)
    return r


def evaluate(net, task_idx, rng, n=120):
    m = TASKS[task_idx]
    hits = 0
    for _ in range(n):
        x = int(rng.integers(K))
        hits += trial(net, x, m(x), rng, learn=False)
    return hits / n


# P2 causal-credit hyperparameters (weight-perturbation on the plastic couplings)
CAUSAL_LR = float(os.environ.get("STATS_CAUSAL_LR", "3.0"))
CAUSAL_SIGMA = float(os.environ.get("STATS_CAUSAL_SIGMA", "0.05"))


def causal_train_step(net, x, target, rng, bstate, pidx, lr=CAUSAL_LR, sigma=CAUSAL_SIGMA):
    """One causal-credit update (P2): a do(θ)-style intervention on the plastic
    couplings. Perturb the weights by ε, measure the reward change, and credit each
    weight by its *measured causal effect* Δw ∝ (r − baseline)·ε — a weight-
    perturbation REINFORCE estimator (an interventional do(w+ε)), in contrast to the
    correlational eligibility trace. `bstate` is a 1-element running-baseline holder.
    """
    w = net.W
    eps = sigma * rng.standard_normal(pidx[0].size)
    w[pidx] += eps
    r = trial(net, x, target, rng, learn=False)     # forward only, perturbed couplings
    w[pidx] -= eps                                   # remove the perturbation
    b = bstate[0]
    w[pidx] += lr * (r - b) * eps                    # credit by measured causal effect
    np.clip(w, 0.0, None, out=w)
    bstate[0] = 0.98 * b + 0.02 * r                  # running reward baseline
    return r


CAUSAL_M = int(os.environ.get("STATS_CAUSAL_M", "4"))   # antithetic samples for low-var


def causal_lowvar_step(net, x, target, rng, pidx, lr=CAUSAL_LR, sigma=CAUSAL_SIGMA, M=CAUSAL_M):
    """P3 — a *low-variance* causal-credit update: antithetic central-difference
    weight-perturbation, averaged over M perturbation pairs. The (r+ − r−)·ε
    estimator is baseline-free (the difference cancels the baseline) and much lower
    variance than the P2 single-sided rule — the tractable stand-in for the Mesnard
    low-variance (hindsight) estimator. Tests whether *variance*, not credit quality,
    was what pinned P2 to the correlational frontier.
    """
    w = net.W
    grad = np.zeros(pidx[0].size)
    for _ in range(M):
        eps = sigma * rng.standard_normal(pidx[0].size)
        w[pidx] += eps
        r_plus = trial(net, x, target, rng, learn=False)
        w[pidx] -= 2.0 * eps
        r_minus = trial(net, x, target, rng, learn=False)
        w[pidx] += eps                               # restore
        grad += (r_plus - r_minus) * eps             # central difference (antithetic)
    w[pidx] += lr * grad / M
    np.clip(w, 0.0, None, out=w)
    # report the +eps arm's reward on the last sample as a rough progress signal
    return r_plus


def run_sequence(seed, regime, credit="correlational", eta_w=None, causal_lr=CAUSAL_LR):
    """Train T0->T1->T2 on one net; return accuracy matrix R[after_task, eval_task].

    credit='correlational' uses Line A's eligibility trace (P1); credit='causal'
    uses weight-perturbation do(θ) credit on the plastic couplings (P2). eta_w and
    causal_lr expose the two rules' plasticity knobs for the stability-plasticity
    frontier control.
    """
    net = make(seed, regime, eta_w=eta_w)
    rng = np.random.default_rng(seed + 7)
    T = len(TASKS)
    R = np.zeros((T, T))
    just_trained = np.zeros(T)              # acc on task k right after training it
    pidx = np.where(net.plastic)
    bstate = [CHANCE]
    for k in range(T):
        m = TASKS[k]
        for _ in range(N_TRAIN):
            x = int(rng.integers(K))
            if credit == "causal":
                causal_train_step(net, x, m(x), rng, bstate, pidx, lr=causal_lr)
            elif credit == "causal_lowvar":
                causal_lowvar_step(net, x, m(x), rng, pidx, lr=causal_lr)
            else:
                trial(net, x, m(x), rng, learn=True)
        just_trained[k] = evaluate(net, k, rng)
        for j in range(T):
            R[k, j] = evaluate(net, j, rng)
    return R, just_trained


def cl_metrics(R, just_trained):
    T = R.shape[0]
    avg_acc = float(R[-1].mean())                                   # after full sequence
    # backward transfer: change on task k (k<T-1) from its just-trained level to the end
    bwt = float(np.mean([R[-1, k] - just_trained[k] for k in range(T - 1)]))
    # forward transfer: acc on task k before training it (after training k-1), vs chance
    fwt = float(np.mean([R[k - 1, k] - CHANCE for k in range(1, T)]))
    return avg_acc, bwt, fwt


def main():
    regimes = ["frozen", "plastic"]
    out = {}
    for regime in regimes:
        avg = np.zeros(N); bwt = np.zeros(N); fwt = np.zeros(N)
        Rsum = np.zeros((len(TASKS), len(TASKS)))
        for s in range(N):
            R, jt = run_sequence(s, regime)
            a, b, f = cl_metrics(R, jt)
            avg[s], bwt[s], fwt[s] = a, b, f
            Rsum += R
        out[regime] = {"avg_acc": avg, "bwt": bwt, "fwt": fwt, "Rmean": Rsum / N}
        print(f"\n=== {regime} (n={N}, chance={CHANCE:.3f}) ===", flush=True)
        for lab, arr in [("avg_acc", avg), ("backward_transfer", bwt), ("forward_transfer", fwt)]:
            print("  " + st.fmt_row(st.summarise(lab, arr)), flush=True)
        print("  Rmean[after_task, eval_task]:\n" +
              "\n".join("    " + " ".join(f"{v:.2f}" for v in row) for row in out[regime]["Rmean"]),
              flush=True)

    np.savez(os.path.join(OUT, "continual_p1.npz"), n=N, chance=CHANCE,
             **{f"{r}_{k}": out[r][k] for r in regimes for k in ("avg_acc", "bwt", "fwt")},
             **{f"{r}_Rmean": out[r]["Rmean"] for r in regimes})
    with open(os.path.join(OUT, "continual_p1.json"), "w") as f:
        json.dump({"n": N, "chance": CHANCE,
                   "rows": {f"{r}:{k}": st.summarise(k, out[r][k])
                            for r in regimes for k in ("avg_acc", "bwt", "fwt")}},
                  f, indent=2, default=float)
    print("\nwrote continual_p1.{npz,json}", flush=True)


def main_p2():
    """P2 — causal vs correlational credit on the 2-task reversal (K=2), in the
    frozen-representation regime so only the shared H->M head adapts (both rules act
    on the *same* parameters; isolates the credit rule). Requires STATS_K=2."""
    assert K == 2, "run P2 with STATS_K=2 (2-task reversal)"
    conds = [("frozen", "correlational"), ("frozen", "causal")]
    out = {}
    for regime, credit in conds:
        key = f"{regime}_{credit}"
        avg = np.zeros(N); bwt = np.zeros(N); fwt = np.zeros(N)
        Rsum = np.zeros((len(TASKS), len(TASKS)))
        for s in range(N):
            R, jt = run_sequence(s, regime, credit=credit)
            a, b, f = cl_metrics(R, jt)
            avg[s], bwt[s], fwt[s] = a, b, f
            Rsum += R
        out[key] = {"avg_acc": avg, "bwt": bwt, "fwt": fwt, "Rmean": Rsum / N}
        print(f"\n=== {key} (n={N}, chance={CHANCE:.3f}, lr={CAUSAL_LR}) ===", flush=True)
        for lab, arr in [("avg_acc", avg), ("backward_transfer", bwt), ("forward_transfer", fwt)]:
            print("  " + st.fmt_row(st.summarise(lab, arr)), flush=True)
        print("  Rmean:\n" + "\n".join("    " + " ".join(f"{v:.2f}" for v in row)
                                       for row in out[key]["Rmean"]), flush=True)
    d_bwt = st.effect_size(out["frozen_causal"]["bwt"], out["frozen_correlational"]["bwt"])
    print(f"\n  causal − correlational backward-transfer effect size (Cohen d) = {d_bwt:.2f}",
          flush=True)
    np.savez(os.path.join(OUT, "continual_p2.npz"), n=N, chance=CHANCE, causal_lr=CAUSAL_LR,
             **{f"{k}_{m}": out[k][m] for k in out for m in ("avg_acc", "bwt", "fwt")},
             **{f"{k}_Rmean": out[k]["Rmean"] for k in out})
    with open(os.path.join(OUT, "continual_p2.json"), "w") as f:
        json.dump({"n": N, "chance": CHANCE, "causal_lr": CAUSAL_LR,
                   "rows": {f"{k}:{m}": st.summarise(m, out[k][m])
                            for k in out for m in ("avg_acc", "bwt", "fwt")}},
                  f, indent=2, default=float)
    print("\nwrote continual_p2.{npz,json}", flush=True)


def main_frontier():
    """Stability-plasticity frontier control (K=2, frozen regime): sweep each rule's
    plasticity knob and record (new-task acquisition R[1,1], retention of task 0
    R[1,0]). If causal sits ON the correlational frontier, the P2 effect is just an
    effective learning-rate difference; if causal's frontier dominates (more
    retention at matched acquisition), the causal credit genuinely interferes less.
    """
    assert K == 2, "run frontier with STATS_K=2"
    nf = int(os.environ.get("STATS_NF", "15"))
    sweeps = {
        "correlational": ("eta_w", [0.01, 0.02, 0.05, 0.1, 0.2]),
        "causal": ("causal_lr", [0.5, 1.0, 2.0, 4.0, 8.0]),
        "causal_lowvar": ("causal_lr", [1.0, 2.0, 4.0, 8.0, 16.0]),
    }
    res = {}
    for credit, (knob, vals) in sweeps.items():
        res[credit] = []
        for v in vals:
            acq = np.zeros(nf); ret = np.zeros(nf)
            kw = {knob: v}
            for s in range(nf):
                R, _ = run_sequence(s, "frozen", credit=credit, **kw)
                acq[s], ret[s] = R[1, 1], R[1, 0]
            res[credit].append((v, float(acq.mean()), float(ret.mean())))
            print(f"  {credit:14s} {knob}={v:<5}: new-task={acq.mean():.2f}  "
                  f"retention={ret.mean():.2f}", flush=True)
    with open(os.path.join(OUT, "continual_p3_frontier.json"), "w") as f:
        json.dump({"nf": nf, **res}, f, indent=2)
    np.savez(os.path.join(OUT, "continual_p3_frontier.npz"),
             **{f"{c}_pts": np.array([r[1:] for r in res[c]]) for c in res},
             **{f"{c}_knob": np.array([r[0] for r in res[c]]) for c in res})
    print("\nwrote continual_p3_frontier.{json,npz}", flush=True)


if __name__ == "__main__":
    if "frontier" in sys.argv[1:]:
        main_frontier()
    elif "p2" in sys.argv[1:]:
        main_p2()
    else:
        main()
