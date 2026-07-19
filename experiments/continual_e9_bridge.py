"""Track 3c bridge — does a *learned* conjunctive basis (E9) resolve continual-
learning interference? (E9 ↔ 3c)

3c/P4 showed capacity resolves interference, and that a single shared head with a
task-context input cut forgetting ~3× — but only to avg 0.66, because the context
projection was *frozen random*. The prediction: a **learned** (stimulus × context)
conjunction basis should push the single head toward the per-task ceiling. E9 grows
exactly that basis by a reward-free competitive-Hebbian rule.

This reuses E9's substrate and its emergent / frozen / wired bases, but runs the
**sequential reversal** continual protocol (train rule 0 = identity, then rule 1 =
reversal; E9's `trial_overlap` target is `x ^ rule`, so the two rules are the two
reversal tasks) and measures the continual-learning metric instead of E9's
switching accuracy. Same reward-phase readout (Line A) for all three bases; the
only difference is whether the conjunction basis is learned (emergent), random
(frozen), or hand-wired (wired).
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))
import e9_emergent_conjunction as e9  # noqa: E402
import ghca_stats as st  # noqa: E402

N = int(os.environ.get("STATS_N", "10"))
N_BLOCKS = int(os.environ.get("STATS_NBLOCKS", "20"))   # training blocks per task
BLOCK_LEN = 25
TASKS = [0, 1]                       # E9 rules: 0 = identity (x^0), 1 = reversal (x^1)
CHANCE = 1.0 / e9.A
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def evaluate(net, roles, rule, rng, n_blocks=4):
    hits = tot = 0
    for _ in range(n_blocks):
        e9.ignite_block(net, roles, rule)             # keep the context ring alive
        for _ in range(BLOCK_LEN):
            x = int(rng.integers(e9.K))
            r, _ = e9.trial_overlap(net, roles, x, rule, rng, learn=False)
            hits += r; tot += 1
    return hits / tot


def run_sequence(seed, kind):
    net, roles = e9.make(seed, kind)
    if kind == "emergent":
        e9.selforg(net, roles, np.random.default_rng(seed + 3))   # grow the basis
    rng = np.random.default_rng(seed + 7)
    R = np.zeros((len(TASKS), len(TASKS))); just = np.zeros(len(TASKS))
    for k in TASKS:
        for _ in range(N_BLOCKS):
            e9.ignite_block(net, roles, k)
            for _ in range(BLOCK_LEN):
                x = int(rng.integers(e9.K))
                e9.trial_overlap(net, roles, x, k, rng, learn=True)
        just[k] = evaluate(net, roles, k, rng)
        for j in TASKS:
            R[k, j] = evaluate(net, roles, j, rng)
    return R, just


def cl_metrics(R, just):
    T = R.shape[0]
    avg = float(R[-1].mean())
    bwt = float(np.mean([R[-1, k] - just[k] for k in range(T - 1)]))
    fwt = float(np.mean([R[k - 1, k] - CHANCE for k in range(1, T)]))
    return avg, bwt, fwt


def main():
    kinds = ["frozen", "emergent", "wired"]
    out = {}
    for kind in kinds:
        avg = np.zeros(N); bwt = np.zeros(N); fwt = np.zeros(N)
        Rsum = np.zeros((len(TASKS), len(TASKS)))
        for s in range(N):
            R, jt = run_sequence(s, kind)
            a, b, f = cl_metrics(R, jt)
            avg[s], bwt[s], fwt[s] = a, b, f
            Rsum += R
            print(f"    {kind} seed {s+1}/{N}: avg={a:.2f} bwt={b:.2f}", flush=True)
        out[kind] = {"avg_acc": avg, "bwt": bwt, "fwt": fwt, "Rmean": Rsum / N}
        print(f"=== {kind} (n={N}, chance={CHANCE:.2f}) ===", flush=True)
        for lab, arr in [("avg_acc", avg), ("backward_transfer", bwt)]:
            print("  " + st.fmt_row(st.summarise(lab, arr)), flush=True)
        print("  Rmean:\n" + "\n".join("    " + " ".join(f"{v:.2f}" for v in row)
                                       for row in out[kind]["Rmean"]), flush=True)
    np.savez(os.path.join(OUT, "continual_e9_bridge.npz"), n=N, chance=CHANCE,
             **{f"{k}_{m}": out[k][m] for k in out for m in ("avg_acc", "bwt", "fwt")},
             **{f"{k}_Rmean": out[k]["Rmean"] for k in out})
    with open(os.path.join(OUT, "continual_e9_bridge.json"), "w") as f:
        json.dump({"n": N, "chance": CHANCE,
                   "rows": {f"{k}:{m}": st.summarise(m, out[k][m])
                            for k in out for m in ("avg_acc", "bwt")}}, f, indent=2, default=float)
    print("\nwrote continual_e9_bridge.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
