"""
E5 - Executive control / task switching (options).

Shows a SLOW LOOP acting as an OPTION that gates fast stimulus->response
mappings, and tests rule reversal. See docs/learning_experiments.md section 5, E5.

Task
----
K=2 stimuli, A=2 actions. Two rules alternate in BLOCKS:
    Rule 0 (identity):  action = x         (stim0->act0, stim1->act1)
    Rule 1 (reversal):  action = 1 - x     (stim0->act1, stim1->act0)
So the correct action is  action = x XOR rule  -- it depends on BOTH the
stimulus AND the currently-active rule.

Mechanism (the "option")
------------------------
Each rule owns a SLOW directed-ring loop (the E2 working-memory mechanism: a ring
of length L sustains a rotating pulse indefinitely when tau < L). At the start of
each block the block's context cue is presented briefly, IGNITING that rule's ring
(and the other ring is reset/extinguished). The cue is then removed; the slow ring
PERSISTS for the whole block and supplies a standing context signal to the hidden
medium -- this persistence is the option holding the rule.

Hidden units receive BOTH a channel-biased stimulus input (a preferred stimulus)
AND a context-biased ring input (a preferred rule), so they behave as
(stimulus x rule) conjunction cells: the substrate's innate repertoire. Line A
(reward-driven conduction plasticity) learns the hidden->motor readout that routes
the four conjunction patterns to the two correct actions (the XOR). Nothing is
taught per-node; the only signal is a strict scalar reward r = 1[action==rule(x)],
from which the learner forms delta = r - V internally.

Predictions
-----------
* Accuracy dips at each block switch and RECOVERS; the switch cost SHRINKS across
  blocks as both context-conditioned mappings consolidate (options consolidate).
* The active rule is decodable from the slow ring's state (near-perfect while it
  persists).
* DISCRIMINATOR: ablating the slow loop (set ring tau >= L so it dies right after
  the cue, i.e. context is NOT held) abolishes flexible SWITCHING while leaving
  SINGLE-RULE performance intact.

Outputs
-------
docs/figures/e5_switching.png     : learning curve across blocks + switch-cost decay
docs/figures/e5_discriminator.png : intact vs ablated, switching vs single-rule
result/e5/e5_data.npz
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
DATADIR = os.path.join(ROOT, "result", "e5")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

# --- substrate / task constants ---
K = A = 2
N_S = 8                       # sensory nodes per stimulus channel
N_H = 120                     # hidden conjunction medium
N_M = 8                       # motor nodes per action channel
L_RING = 16                   # context ring length
TAU_SLOW = 12                 # ring tau < L  -> ring sustains  (the option holds)
TAU_DEAD = 18                 # ring tau >= L -> ring dies ~L steps after cue (ablation)
ACT = 2

SETTLE = 6                    # free steps before stimulus each trial
CUE = 3                       # stimulus clamp duration
WWIN = 6                      # response window
CTX_CUE = 6                   # context-cue duration at a block start


def build(seed=0, channel_bias=0.85, w_ring=1.0,
          w_sh=0.7, w_ch=0.7, w_hm=0.5, theta_h=3.5, theta_m=1.0,
          fanin_sh=4, fanin_hm=14, jitter=0.15):
    """Build the S + two context-rings + H + M graph.

    Node layout:
        sensory : K channels of N_S            -> [s0 ..)
        rings   : A rules, each a directed ring of L_RING nodes
        hidden  : N_H conjunction medium
        motor   : A channels of N_M

    Each hidden unit has a preferred stimulus (channel-biased sensory fan-in) and a
    preferred rule (it fans in from ALL nodes of that rule's ring, so the context
    signal is phase-invariant: while the ring is alive its ~`act` rotating-pulse
    nodes deliver a steady drive regardless of pulse phase). Thresholds are tuned
    (calibrated) to a HARD COINCIDENCE (AND) gate: stimulus-alone is subthreshold
    (~0 active), ring-alone is subthreshold (~0), and only their SUM fires the unit
    (~0.83 active). So during a block only the (current stimulus x current rule)
    conjunction subpopulation is active and the other three are silent -- no
    cross-talk. The alive ring is thus a hard OPTION that selects which routing
    subpopulation exists; H->M is the plastic readout (Line A) that maps each
    conjunction subpopulation to its rewarded action (the XOR). Because the option
    is required for hidden firing, a held (persistent) ring is what lets routing be
    conditioned on context across a block -- see the single-rule control, which
    re-cues every trial so it needs no persistence.
    """
    rng = np.random.default_rng(seed)
    s0 = 0
    r0 = s0 + K * N_S
    h0 = r0 + A * L_RING
    m0 = h0 + N_H
    N = m0 + A * N_M

    sensory = [np.arange(s0 + c * N_S, s0 + (c + 1) * N_S) for c in range(K)]
    rings = [np.arange(r0 + g * L_RING, r0 + (g + 1) * L_RING) for g in range(A)]
    hidden = np.arange(h0, h0 + N_H)
    motor = [np.arange(m0 + c * N_M, m0 + (c + 1) * N_M) for c in range(A)]
    allsen = np.concatenate(sensory)

    W = np.zeros((N, N))
    plastic = np.zeros((N, N), dtype=bool)
    theta = np.zeros(N)
    theta[allsen] = 1.0

    # directed rings (each rule's context loop): i -> i+1 around the ring
    for g in range(A):
        rg = rings[g]
        for a in range(L_RING):
            W[rg[(a + 1) % L_RING], rg[a]] = w_ring
        theta[rg] = 1.0

    # hidden units: preferred stimulus (channel-biased) + preferred rule (all of
    # that ring) -> (stimulus x rule) conjunction cells, gated by ring boost.
    pref_s = rng.integers(K, size=N_H)
    pref_g = rng.integers(A, size=N_H)
    for i, h in enumerate(hidden):
        sp = sensory[pref_s[i]]
        ns = int(round(channel_bias * fanin_sh))
        src_s = np.unique(np.concatenate([
            rng.choice(sp, min(ns, len(sp)), replace=False),
            rng.choice(allsen, max(fanin_sh - ns, 0), replace=False)]))
        W[h, src_s] = w_sh * (1 + jitter * rng.standard_normal(len(src_s)))
        rg = rings[pref_g[i]]                        # boost from ALL preferred-ring nodes
        W[h, rg] = w_ch * (1 + jitter * rng.standard_normal(L_RING))
        theta[h] = theta_h

    # hidden -> motor readout (plastic, Line A); starts weak + jittered
    for c in range(A):
        for m in motor[c]:
            src = rng.choice(hidden, min(fanin_hm, N_H), replace=False)
            W[m, src] = w_hm * (1 + jitter * rng.standard_normal(len(src)))
            plastic[m, src] = True
            theta[m] = theta_m

    roles = {"sensory": sensory, "rings": rings, "hidden": hidden, "motor": motor,
             "N": N, "K": K, "A": A, "pref_s": pref_s, "pref_g": pref_g}
    return np.clip(W, 0, None), plastic, roles, theta


def make(seed, ablate=False, eta_w=0.06, p_s=3e-3):
    W, plastic, roles, theta = build(seed=seed)
    tau_ring = TAU_DEAD if ablate else TAU_SLOW
    # Exploration (p_s) is confined to the hidden+motor medium: it supplies the
    # variance the reward-driven readout needs (as in E1), while leaving the
    # sensory cue and the context RINGS deterministic -- critically, masked p_s
    # never spuriously re-ignites an ablated (dead) ring, so the discriminator is
    # clean.
    ps_mask = np.zeros(roles["N"], dtype=bool)
    ps_mask[roles["hidden"]] = True
    for c in range(roles["A"]):
        ps_mask[roles["motor"][c]] = True
    net = GHLearner(W, plastic, roles, line="A", act=ACT, pas=tau_ring - ACT,
                    theta=theta, p_s=p_s, p_s_mask=ps_mask, eta_w=eta_w,
                    lam=0.8, w_max=4.0, seed=seed + 41)
    # per-node tau: rings at tau_ring; everything else is a gapless relay so the
    # cue/stimulus propagate promptly to motor within the response window.
    relay = np.concatenate([np.concatenate(roles["sensory"]), roles["hidden"]]
                           + roles["motor"])
    net.tau[relay] = ACT
    for g in range(roles["A"]):
        net.tau[roles["rings"][g]] = tau_ring
    return net, roles


def ignite_block(net, roles, rule):
    """Start a block: extinguish both rings, then ignite the active rule's ring.

    A SINGLE seed on the ring head launches a rotating pulse. With ring tau < L it
    sustains indefinitely (the option holds the rule); with tau >= L (ablation) it
    dies ~L steps later. A clamped multi-step ignition instead builds a smearing
    block that self-annihilates, so a single seed is used."""
    for g in range(roles["A"]):
        net.phi[roles["rings"][g]] = 0
    net.phi[roles["rings"][rule][0]] = 1


def ring_alive(net, roles, rule):
    return int(net.active_mask()[roles["rings"][rule]].sum() > 0)


def trial(net, roles, x, rule, rng, learn=True):
    net.reset_traces()
    # free settle (rings keep rotating on their own; stimulus/relay are rested)
    for _ in range(SETTLE):
        net.step_learn(None)
    feats = net.features()
    V = net.value()
    for _ in range(CUE):
        net.step_learn(net.sensory_drive(x))
    sc = np.zeros(roles["A"])
    for _ in range(WWIN):
        net.step_learn(None)
        sc += net.motor_scores()
    action = int(np.argmax(sc + 1e-6 * rng.standard_normal(len(sc)))) if sc.sum() > 0 else -1
    target = x ^ rule
    r = 1.0 if action == target else 0.0
    if learn:
        net.learn(r - V)
        net.update_critic(r - V, feats)
    return r, ring_alive(net, roles, rule)


def run_switching(seed, ablate=False, n_blocks=40, block_len=25):
    """Alternating-rule blocks; returns per-trial accuracy, switch flags, ring-alive."""
    net, roles = make(seed, ablate=ablate)
    rng = np.random.default_rng(seed + 5)
    acc, switch, alive = [], [], []
    for b in range(n_blocks):
        rule = b % 2
        ignite_block(net, roles, rule)
        for i in range(block_len):
            x = int(rng.integers(K))
            r, al = trial(net, roles, x, rule, rng, learn=True)
            acc.append(r)
            switch.append(1 if i == 0 else 0)
            alive.append(al)
    return np.array(acc), np.array(switch), np.array(alive), net, roles


def run_single_rule(seed, ablate=False, rule=0, n_trials=600):
    """Control: one fixed rule, ring RE-CUED EVERY TRIAL so persistence is not
    needed. This isolates the routing itself from the held-context requirement:
    because the ring is freshly ignited each trial, it is alive during the readout
    in BOTH the intact and ablated substrates, so basic single-rule routing should
    work in both. Contrast run_switching, which cues only at each block start and
    therefore depends on the ring PERSISTING across the block."""
    net, roles = make(seed, ablate=ablate)
    rng = np.random.default_rng(seed + 9)
    acc = []
    for t in range(n_trials):
        ignite_block(net, roles, rule)       # re-cue each trial: persistence irrelevant
        x = int(rng.integers(K))
        r, _ = trial(net, roles, x, rule, rng, learn=True)
        acc.append(r)
    return np.array(acc)


def switch_cost_by_block(acc, block_len, n_blocks, post=4):
    """Mean accuracy in the first `post` trials of each block (post-switch)."""
    a = acc.reshape(n_blocks, block_len)
    return a[:, :post].mean(1)


def moving_avg(x, w=25):
    return np.convolve(x, np.ones(w) / w, mode="valid")


def main():
    n_seeds = 5
    n_blocks, block_len = 40, 25

    print("E5: switching (intact) ...")
    acc_i = np.zeros((n_seeds, n_blocks * block_len))
    alive_i = np.zeros((n_seeds, n_blocks * block_len))
    for s in range(n_seeds):
        a, sw, al, net, roles = run_switching(s, ablate=False,
                                              n_blocks=n_blocks, block_len=block_len)
        acc_i[s], alive_i[s] = a, al
    print("E5: switching (ablated slow loop) ...")
    acc_a = np.zeros((n_seeds, n_blocks * block_len))
    alive_a = np.zeros((n_seeds, n_blocks * block_len))
    for s in range(n_seeds):
        a, sw, al, _, _ = run_switching(s, ablate=True,
                                        n_blocks=n_blocks, block_len=block_len)
        acc_a[s], alive_a[s] = a, al

    print("E5: single-rule (intact / ablated) ...")
    single_i = np.array([run_single_rule(s, ablate=False) for s in range(n_seeds)])
    single_a = np.array([run_single_rule(s, ablate=True) for s in range(n_seeds)])

    # switch cost per block (post-switch accuracy), first vs last third of training
    sc_i = np.array([switch_cost_by_block(acc_i[s], block_len, n_blocks) for s in range(n_seeds)])
    sc_a = np.array([switch_cost_by_block(acc_a[s], block_len, n_blocks) for s in range(n_seeds)])

    # headline scalars
    final_i = acc_i[:, -block_len * 6:].mean()
    final_a = acc_a[:, -block_len * 6:].mean()
    postswitch_early_i = sc_i[:, :6].mean()
    postswitch_late_i = sc_i[:, -6:].mean()
    single_final_i = single_i[:, -150:].mean()
    single_final_a = single_a[:, -150:].mean()
    rule_decode_i = alive_i.mean()          # fraction of trials ring holds the rule
    rule_decode_a = alive_a.mean()
    switch_final_a = acc_a[:, -block_len * 6:].mean()

    print(f"  switching final acc:  intact={final_i:.2f}  ablated={switch_final_a:.2f}")
    print(f"  post-switch acc:      early={postswitch_early_i:.2f}  late={postswitch_late_i:.2f} (intact)")
    print(f"  single-rule final:    intact={single_final_i:.2f}  ablated={single_final_a:.2f}")
    print(f"  ring holds rule (frac trials alive): intact={rule_decode_i:.2f}  ablated={rule_decode_a:.2f}")

    # ---- figure 1: learning curve across blocks + switch-cost decay ----
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    m = acc_i.mean(0)
    mv = moving_avg(m, 15)
    axes[0].plot(np.arange(len(mv)), mv, color="seagreen", label="intact (slow loop)")
    ma = moving_avg(acc_a.mean(0), 15)
    axes[0].plot(np.arange(len(ma)), ma, color="crimson", alpha=0.8, label="ablated (loop dies)")
    for b in range(1, n_blocks):
        axes[0].axvline(b * block_len, color="k", alpha=0.05)
    axes[0].axhline(0.5, ls="--", color="k", alpha=0.4)
    axes[0].set_xlabel("trial (15-trial moving average)")
    axes[0].set_ylabel("accuracy")
    axes[0].set_ylim(0, 1)
    axes[0].set_title("E5: rule-switching accuracy across blocks")
    axes[0].legend(fontsize=8)

    bl = np.arange(n_blocks)
    axes[1].plot(bl, sc_i.mean(0), "o-", color="seagreen", label="intact")
    axes[1].plot(bl, sc_a.mean(0), "o-", color="crimson", alpha=0.8, label="ablated")
    axes[1].axhline(0.5, ls="--", color="k", alpha=0.4, label="chance")
    axes[1].set_xlabel("block number")
    axes[1].set_ylabel(f"post-switch accuracy (first 4 trials)")
    axes[1].set_ylim(0, 1)
    axes[1].set_title("Switch cost: post-switch accuracy rises\nas options consolidate (intact only)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    p1 = os.path.join(FIGDIR, "e5_switching.png")
    fig.savefig(p1, dpi=110)
    print("wrote", p1)

    # ---- figure 2: discriminator ----
    fig, ax = plt.subplots(figsize=(7, 5))
    groups = ["switching\n(flexible)", "single-rule\n(fixed)"]
    intact = [final_i, single_final_i]
    ablated = [switch_final_a, single_final_a]
    xpos = np.arange(len(groups))
    w = 0.36
    ax.bar(xpos - w / 2, intact, w, color="seagreen", label="intact slow loop",
           yerr=[acc_i[:, -block_len*6:].mean(1).std()/np.sqrt(n_seeds),
                 single_i[:, -150:].mean(1).std()/np.sqrt(n_seeds)], capsize=4)
    ax.bar(xpos + w / 2, ablated, w, color="crimson", label="ablated slow loop",
           yerr=[acc_a[:, -block_len*6:].mean(1).std()/np.sqrt(n_seeds),
                 single_a[:, -150:].mean(1).std()/np.sqrt(n_seeds)], capsize=4)
    ax.axhline(0.5, ls="--", color="k", alpha=0.5, label="chance")
    ax.set_xticks(xpos)
    ax.set_xticklabels(groups)
    ax.set_ylabel("final accuracy")
    ax.set_ylim(0, 1)
    ax.set_title("E5 discriminator: the slow loop is needed for SWITCHING,\n"
                 "not for single-rule mapping")
    ax.legend(fontsize=9)
    fig.tight_layout()
    p2 = os.path.join(FIGDIR, "e5_discriminator.png")
    fig.savefig(p2, dpi=110)
    print("wrote", p2)

    np.savez(os.path.join(DATADIR, "e5_data.npz"),
             acc_intact=acc_i, acc_ablated=acc_a,
             alive_intact=alive_i, alive_ablated=alive_a,
             single_intact=single_i, single_ablated=single_a,
             switchcost_intact=sc_i, switchcost_ablated=sc_a,
             n_blocks=n_blocks, block_len=block_len)
    print("wrote", os.path.join(DATADIR, "e5_data.npz"))


if __name__ == "__main__":
    main()
