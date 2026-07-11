"""
E9 - Emergent conjunction cells (retiring the "afforded vs learned" caveat).

E1 and E5 WIRE stimulus/context selectivity into the hidden medium: E1 assigns
each hidden unit a preferred sensory channel (and draws 85% of its S->H fan-in
from it); E5 goes further and assigns each unit BOTH a preferred stimulus and a
preferred rule, so the (stimulus x rule) conjunction basis the XOR routing needs
is hand-built. Line A then only learns the H->M readout on top of a pre-existing
basis. The extensions audit (docs/extensions_review.md, caveat 1) flags this as
the central "afforded, not learned" limitation, and it is the oldest deferred
item in the programme (E1: "let the net DISCOVER selective hidden
representations"). See docs/next_steps.md, Track 1a.

This experiment removes the wiring. Hidden units start with UNBIASED input
fan-in -- drawn uniformly across BOTH sensory channels and BOTH context rings,
with no preferred stimulus and no preferred rule -- and a local, reward-free
competitive Hebbian rule on those input edges GROWS the conjunction basis from
the input statistics alone. The E5 AND-gate threshold is kept: it is what turns
coincidence-Hebbian into CONJUNCTION tuning (a unit fires only when its
strengthened stimulus AND its strengthened rule are both present). This is the
mechanism the tactile-coding review (npj/Microsyst. Nanoeng.,
doi:10.1038/s41378-025-01074-3) calls the move from INDEPENDENT population
coding (per-neuron tuning wired in) to CORRELATED population coding (tuning that
emerges from microscopic activity correlations, a la STDP).

Design
------
Two learning processes, both local, neither given task labels:
  1. UNSUPERVISED self-organisation (reward-free): present the four (stim x rule)
     combinations in random order; a competitive Hebbian rule + per-input-type
     weight normalisation + per-unit homeostatic thresholds differentiate the
     hidden units until each becomes selective for one conjunction, tiling all
     four. No reward, no labels -- only the input stream.
  2. REWARD-driven routing (Line A, reused verbatim from E1/E5): with the
     emergent input basis FROZEN, a strict scalar reward carves the H->M readout
     that maps the four conjunction subpopulations to the two correct actions
     (the XOR), exactly as in E5.

Conditions compared (same task, same reward phase):
  * EMERGENT  - unbiased init + unsupervised Hebbian self-organisation (this work)
  * WIRED     - E5's hand-built pref_s/pref_g basis (reference upper bound)
  * FROZEN    - unbiased init, NO Hebbian (random basis; the ablation that shows
                the self-organisation is load-bearing)

Predictions
-----------
* Conjunction selectivity of the hidden basis RISES over unsupervised training
  (emergence curve), from ~frozen-random up toward ~wired; all four combos tiled.
* Reward-driven switching on the EMERGENT basis ~ matches WIRED (~0.86, E5).
* On the FROZEN random basis, conjunction selectivity stays low and switching
  collapses -- the emergence, not the reward routing, supplies the basis.

Outputs
-------
docs/figures/e9_emergence.png     : conj-selectivity emergence curve + combo tiling
docs/figures/e9_switching.png     : switching accuracy, emergent vs wired vs frozen
result/e9/e9_data.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_learn import GHLearner
# reuse E5's substrate helpers (generic in net+roles) and its wired basis
from e5_executive import (build as build_wired, ignite_block, ring_alive, moving_avg,
                          K, A, N_S, N_H, N_M, L_RING, TAU_SLOW, TAU_DEAD, ACT,
                          SETTLE, CUE, WWIN)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e9")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

# self-organisation phase length / schedule
N_SELFORG = 500              # unsupervised combo presentations
SNAP_EVERY = 25              # snapshot the basis this often (emergence curve)
SO_SETTLE = 4                # free steps before the combo drive
SO_DRIVE = 8                 # steps the stimulus is cued (ring already alive)

# probe window for measuring hidden-unit conjunction tuning
PROBE_SETTLE, PROBE_DRIVE = 4, 8


# ----------------------------------------------------------------------------
# Unbiased graph: same layout as E5 but NO pref_s / pref_g, and the input edges
# (S->H, ring->H) are plastic (Hebbian). H->M stays the reward-plastic readout.
# ----------------------------------------------------------------------------

def build_unbiased(seed=0, w_ring=1.0, w_sh=0.29, w_ch=0.575, w_hm=0.5,
                   theta_h=3.5, theta_m=1.0, fanin_s=8, fanin_hm=14, jitter=0.15):
    """S + two context-rings + H + M, with UNBIASED hidden input fan-in.

    Each hidden unit fans in from a random subset of ALL sensory nodes (spanning
    both channels) and from ALL ring nodes of BOTH rings -- there is no preferred
    stimulus and no preferred rule. Input weights start weak/uniform+jittered; the
    Hebbian self-organiser concentrates each unit's weight onto one channel and
    one ring. The AND-gate threshold theta_h is inherited from E5 so that, once a
    unit's weight has concentrated, its preferred (stim x rule) combo is
    suprathreshold while either feature alone is not.
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
    allring = np.concatenate(rings)

    W = np.zeros((N, N))
    plastic_in = np.zeros((N, N), dtype=bool)     # S->H, ring->H (Hebbian)
    plastic_rd = np.zeros((N, N), dtype=bool)     # H->M (reward / Line A)
    theta = np.zeros(N)
    theta[allsen] = 1.0

    # directed context rings (each rule's loop)
    for g in range(A):
        rg = rings[g]
        for a in range(L_RING):
            W[rg[(a + 1) % L_RING], rg[a]] = w_ring
        theta[rg] = 1.0

    # hidden units: UNBIASED fan-in from both channels + both rings, plastic.
    for h in hidden:
        src_s = rng.choice(allsen, min(fanin_s, len(allsen)), replace=False)
        W[h, src_s] = w_sh * (1 + jitter * rng.standard_normal(len(src_s)))
        plastic_in[h, src_s] = True
        # fan in from ALL ring nodes of BOTH rings (phase-invariant context drive)
        W[h, allring] = w_ch * (1 + jitter * rng.standard_normal(len(allring)))
        plastic_in[h, allring] = True
        theta[h] = theta_h

    # hidden -> motor readout (reward-plastic, Line A); weak + jittered
    for c in range(A):
        for m in motor[c]:
            src = rng.choice(hidden, min(fanin_hm, N_H), replace=False)
            W[m, src] = w_hm * (1 + jitter * rng.standard_normal(len(src)))
            plastic_rd[m, src] = True
            theta[m] = theta_m

    W = np.clip(W, 0, None)
    roles = {"sensory": sensory, "rings": rings, "hidden": hidden, "motor": motor,
             "N": N, "K": K, "A": A, "allsen": allsen, "allring": allring}
    return W, plastic_in, plastic_rd, roles, theta


# ----------------------------------------------------------------------------
# The self-organising learner: reward-free competitive Hebbian on input edges.
# ----------------------------------------------------------------------------

class ConjLearner(GHLearner):
    """GHLearner + a local COMPETITIVE Hebbian rule on the hidden INPUT edges.

    Reward-free and label-free. For each (stim x rule) presentation the rule:
      1. runs the substrate and reads each hidden unit's analog input DRIVE to
         the current combo, plus the presynaptic activity profile of the combo;
      2. picks the top `win_frac` of hidden units by drive as WINNERS (a k-WTA /
         lateral-inhibition competition -- the same winner-take-all motif E4 uses
         for attention). Only winners learn this combo;
      3. strengthens each winner's input edges toward the combo's active nodes,
         then renormalises that winner's sensory- and ring-input weight blocks
         back to fixed totals (so weight concentrated on one channel/ring is
         weight taken from the others).
    Competition breaks the symmetry that plain firing-gated Hebbian could not
    (every unit fired for every combo, so nothing differentiated); winners for
    (s0,r0) concentrate onto s0+r0 and lose the other combos to other units,
    tiling all four. Thresholds are held fixed at the E5 AND-gate value so that,
    once weight has concentrated, a unit's preferred combo is suprathreshold
    while either feature alone is not.
    """

    def __init__(self, W, plastic_rd, roles, plastic_in,
                 eta_hebb=0.25, win_frac=0.25, consc_c=100.0, consc_rate=0.05, **kw):
        # Line A operates on the READOUT edges only (reward phase); input edges
        # are driven by the competitive rule below, not by reward.
        super().__init__(W, plastic_rd, roles, line="A", **kw)
        self.plastic_in = np.asarray(plastic_in, bool)
        self.hidden = roles["hidden"]
        self.rings = roles["rings"]
        self._insen = self._in_set(roles["allsen"])
        self._inring = self._in_set(roles["allring"])
        # per-unit input-edge index arrays (sensory block, ring block)
        self.s_edge = {h: np.where(self.plastic_in[h] & self._insen)[0] for h in self.hidden}
        self.r_edge = {h: np.where(self.plastic_in[h] & self._inring)[0] for h in self.hidden}
        # target block sums = the initial block sums (preserve the AND balance)
        self.s_target = {h: self.W[h, self.s_edge[h]].sum() for h in self.hidden}
        self.r_target = {h: self.W[h, self.r_edge[h]].sum() for h in self.hidden}
        self.eta_hebb = eta_hebb
        self.win_frac = win_frac
        # DeSieno "conscience": equalise win rates so no conjunction class is left
        # unclaimed (dead cluster). p_win tracks each unit's win frequency; units
        # winning more than their fair share (win_frac) are biased down.
        self.consc_c = consc_c
        self.consc_rate = consc_rate
        self.p_win = np.full(len(self.hidden), win_frac)

    def _in_set(self, idx):
        m = np.zeros(self.N, bool)
        m[idx] = True
        return m

    def selforg_present(self, roles, s, r, settle, drive_steps):
        """One competitive-Hebbian presentation of combo (s, r)."""
        ignite_block(self, roles, r)
        act_in = np.zeros(self.N)
        drive_h = np.zeros(len(self.hidden))
        Wsub = self.W[self.hidden]
        drv = self.sensory_drive(s)
        for k in range(settle + drive_steps):
            self.step(drv if k >= settle else None)
            a = self.active_mask().astype(float)
            act_in += a * (self._insen | self._inring).astype(float)
            drive_h += Wsub @ a
        # k-WTA competition among hidden units, with conscience bias
        n_win = max(1, int(self.win_frac * len(self.hidden)))
        score = drive_h - self.consc_c * (self.p_win - self.win_frac)
        win_idx = np.argsort(-score)[:n_win]
        winners = self.hidden[win_idx]
        won = np.zeros(len(self.hidden), bool)
        won[win_idx] = True
        self.p_win += self.consc_rate * (won.astype(float) - self.p_win)
        # Ring Hebbian target is RING-LEVEL (uniform within a ring), not per-node:
        # the context feature a unit selects is a ring identity (a rule), and its
        # weights within the chosen ring must stay uniform so the drive is
        # phase-invariant as the pulse rotates (E5's design). Per-node ring
        # concentration instead ties a unit to a phase band, so a short cue at a
        # random phase misses it. Sensory can stay per-node (the stimulus drives
        # all channel nodes uniformly, so within-channel weight stays uniform).
        ring_act = np.zeros(self.N)
        for g in range(len(self.rings)):
            rg = self.rings[g]
            ring_act[rg] = act_in[rg].mean()
        # Hebbian on winners toward the combo's active input profile, then
        # renormalise each winner's sensory / ring block to its target sum.
        for h in winners:
            se, re = self.s_edge[h], self.r_edge[h]
            if se.size:
                self.W[h, se] = np.clip(self.W[h, se] + self.eta_hebb * act_in[se], 0.0, None)
                tot = self.W[h, se].sum()
                if tot > 0:
                    self.W[h, se] *= self.s_target[h] / tot
            if re.size:
                self.W[h, re] = np.clip(self.W[h, re] + self.eta_hebb * ring_act[re], 0.0, None)
                tot = self.W[h, re].sum()
                if tot > 0:
                    self.W[h, re] *= self.r_target[h] / tot
        self.adj = self.W > 0


# ----------------------------------------------------------------------------
# make(): assemble a learner on a given basis kind.
# ----------------------------------------------------------------------------

def make(seed, kind, ablate=False, eta_w=0.06, p_s=3e-3):
    """kind in {'emergent','frozen'} use the unbiased graph; 'wired' uses E5's."""
    tau_ring = TAU_DEAD if ablate else TAU_SLOW
    if kind == "wired":
        W, plastic_rd, roles, theta = build_wired(seed=seed)
        # E5's roles lack the flat index arrays ConjLearner needs; not self-organised
        plastic_in = np.zeros_like(plastic_rd)
    else:
        W, plastic_in, plastic_rd, roles, theta = build_unbiased(seed=seed)

    ps_mask = np.zeros(roles["N"], dtype=bool)
    ps_mask[roles["hidden"]] = True
    for c in range(roles["A"]):
        ps_mask[roles["motor"][c]] = True

    common = dict(act=ACT, pas=tau_ring - ACT, theta=theta, p_s=p_s,
                  p_s_mask=ps_mask, eta_w=eta_w, lam=0.8, w_max=4.0, seed=seed + 41)
    if kind == "wired":
        net = GHLearner(W, plastic_rd, roles, line="A", **common)
    else:
        net = ConjLearner(W, plastic_rd, roles, plastic_in, **common)

    relay = np.concatenate([np.concatenate(roles["sensory"]), roles["hidden"]]
                           + roles["motor"])
    net.tau[relay] = ACT
    for g in range(roles["A"]):
        net.tau[roles["rings"][g]] = tau_ring
    return net, roles


# ----------------------------------------------------------------------------
# Unsupervised self-organisation phase (EMERGENT only).
# ----------------------------------------------------------------------------

def selforg(net, roles, rng, n_pres=N_SELFORG, snap_every=SNAP_EVERY):
    """Present (stim x rule) combos in random order, reward-free, applying the
    competitive Hebbian rule. Snapshot mean conjunction selectivity periodically.
    Returns (snap_steps, snap_conj)."""
    # initial snapshot (p=0) so the emergence curve shows the rise from baseline
    snap_steps, snap_conj = [0], [conj_scores(probe_responses(net, roles, rng)).mean()]
    combos = [(s, r) for s in range(K) for r in range(A)]
    for p in range(n_pres):
        # equal exposure: reshuffle the four combos each cycle
        if p % len(combos) == 0:
            order = rng.permutation(len(combos))
        s, r = combos[order[p % len(combos)]]
        net.selforg_present(roles, s, r, SO_SETTLE, SO_DRIVE)
        # snapshot finely over the first ~48 presentations (the conjunction
        # basis converges within ~3 exposures/combo) then coarsely to show it holds
        fine = (p + 1) <= 48 and (p + 1) % 3 == 0
        if fine or (p + 1) % snap_every == 0:
            R = probe_responses(net, roles, rng)
            snap_steps.append(p + 1)
            snap_conj.append(conj_scores(R).mean())
    return np.array(snap_steps), np.array(snap_conj)


# ----------------------------------------------------------------------------
# Selectivity probes.
# ----------------------------------------------------------------------------

def probe_responses(net, roles, rng):
    """Hidden-unit response to each of the 4 (stim x rule) combos.

    Returns R of shape (N_H, 4); column index = s*A + r. Response = fraction of
    the drive window the unit was active. Probing does not learn."""
    hidden = roles["hidden"]
    R = np.zeros((len(hidden), K * A))
    # snapshot phi so probing is side-effect free on the ongoing dynamics; also
    # silence spontaneous firing so the measured tuning is clean (driven only).
    phi_save, t_save, ps_save = net.phi.copy(), net.t, net.p_s
    net.p_s = 0.0
    for s in range(K):
        for r in range(A):
            for h in roles["rings"]:                # extinguish all rings
                net.phi[h] = 0
            net.phi[net.roles["hidden"]] = 0
            ignite_block(net, roles, r)
            for _ in range(PROBE_SETTLE):
                net.step(None)
            acc = np.zeros(len(hidden))
            for _ in range(PROBE_DRIVE):
                net.step(net.sensory_drive(s))
                acc += net.active_mask()[hidden].astype(float)
            R[:, s * A + r] = acc / PROBE_DRIVE
    net.phi, net.t, net.p_s = phi_save, t_save, ps_save
    return R


def conj_scores(R):
    """Per-unit conjunction score in [0,1].

    For each unit let p be its preferred combo (si,rj). A pure conjunction (AND)
    cell fires for (si,rj) but NOT for the combos that share only the stimulus
    (si, r_other) or only the rule (s_other, rj). score = (resp(p) - max(shared))
    / (resp(p)+eps): ~1 = pure conjunction; ~0 = single-feature (stim- or
    rule-only) tuning or unselective."""
    N = R.shape[0]
    sc = np.zeros(N)
    for h in range(N):
        p = int(np.argmax(R[h]))
        ps, pr = p // A, p % A
        shares_stim = ps * A + (1 - pr)
        shares_rule = (1 - ps) * A + pr
        top = R[h, p]
        shared = max(R[h, shares_stim], R[h, shares_rule])
        sc[h] = (top - shared) / (top + 1e-9)
    return sc


def combo_coverage(R):
    """Count of hidden units whose preferred combo is each of the 4."""
    pref = R.argmax(1)
    return np.array([(pref == c).sum() for c in range(K * A)])


# ----------------------------------------------------------------------------
# Reward-driven switching on a (frozen) basis -- E5 protocol.
# ----------------------------------------------------------------------------

def trial_overlap(net, roles, x, rule, rng, learn=True):
    """One trial, reading the action over an epoch that OVERLAPS the cue.

    Same as E5's trial except the motor is scored over the whole cue+response
    epoch (skipping the 1-step conduction delay) rather than only after the cue
    is removed. A perfect AND conjunction cell needs the stimulus present, so it
    goes silent the instant the cue ends; reading the action while the evidence
    is still on is the natural choice and is applied IDENTICALLY to all three
    bases, so the emergent-vs-wired-vs-frozen comparison is fair. (Because the
    readout epoch differs from E5's post-cue window, the absolute numbers here
    are not directly comparable to e5_results.md.)"""
    net.reset_traces()
    for _ in range(SETTLE):
        net.step_learn(None)
    feats = net.features()
    V = net.value()
    sc = np.zeros(roles["A"])
    for k in range(CUE + WWIN):
        net.step_learn(net.sensory_drive(x) if k < CUE else None)
        if k >= 1:                                   # skip the 1-step delay
            sc += net.motor_scores()
    action = int(np.argmax(sc + 1e-6 * rng.standard_normal(len(sc)))) if sc.sum() > 0 else -1
    target = x ^ rule
    r = 1.0 if action == target else 0.0
    if learn:
        net.learn(r - V)
        net.update_critic(r - V, feats)
    return r, ring_alive(net, roles, rule)


def run_switching(net, roles, seed, n_blocks=40, block_len=25):
    rng = np.random.default_rng(seed + 5)
    acc, alive = [], []
    for b in range(n_blocks):
        rule = b % 2
        ignite_block(net, roles, rule)
        for i in range(block_len):
            x = int(rng.integers(K))
            r, al = trial_overlap(net, roles, x, rule, rng, learn=True)
            acc.append(r)
            alive.append(al)
    return np.array(acc), np.array(alive)


# ----------------------------------------------------------------------------
# Main.
# ----------------------------------------------------------------------------

def main():
    n_seeds = 5
    n_blocks, block_len = 40, 25
    kinds = ["emergent", "wired", "frozen"]

    conj_pre = {k: [] for k in kinds}
    conj_post = {k: [] for k in kinds}
    coverage = {k: np.zeros(K * A) for k in kinds}
    emergence = []                       # (steps, conj) curves for emergent
    switch_acc = {k: np.zeros((n_seeds, n_blocks * block_len)) for k in kinds}

    for s in range(n_seeds):
        for k in kinds:
            net, roles = make(s, kind=k)
            rng = np.random.default_rng(s + 3)
            # basis selectivity BEFORE any self-organisation
            R0 = probe_responses(net, roles, rng)
            conj_pre[k].append(conj_scores(R0).mean())
            # self-organise the basis (emergent only)
            if k == "emergent":
                steps, curve = selforg(net, roles, rng)
                if s == 0:
                    emergence.append((steps, curve))
            R1 = probe_responses(net, roles, rng)
            conj_post[k].append(conj_scores(R1).mean())
            coverage[k] += combo_coverage(R1)
            # freeze the input basis: reward phase touches only H->M (Line A).
            # (ConjLearner is already line='A'; wired GHLearner too. The input
            # edges are simply never updated during switching -- no Hebbian calls.)
            acc, _ = run_switching(net, roles, s, n_blocks, block_len)
            switch_acc[k][s] = acc
        print(f"seed {s} done")

    # headline scalars
    def m(x):
        return float(np.mean(x))
    final = {k: switch_acc[k][:, -block_len * 6:].mean() for k in kinds}
    cpre = {k: m(conj_pre[k]) for k in kinds}
    cpost = {k: m(conj_post[k]) for k in kinds}

    print("\n=== E9 emergent conjunction cells ===")
    print(f"  conj selectivity (pre -> post self-org):")
    for k in kinds:
        print(f"    {k:9s}: {cpre[k]:.2f} -> {cpost[k]:.2f}")
    print(f"  combo coverage (units per (s,r) class, summed over seeds):")
    for k in kinds:
        print(f"    {k:9s}: {coverage[k].astype(int)}")
    print(f"  switching final accuracy (last 6 blocks):")
    for k in kinds:
        print(f"    {k:9s}: {final[k]:.2f}")

    # ---- figure 1: emergence curve + tiling ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    if emergence:
        steps, curve = emergence[0]
        axes[0].plot(steps, curve, "o-", color="seagreen", label="emergent (seed 0)")
    axes[0].axhline(cpost["wired"], ls="--", color="slategray", label="wired basis")
    axes[0].axhline(cpost["frozen"], ls=":", color="crimson", label="frozen (no Hebbian)")
    axes[0].set_xlabel("unsupervised presentation")
    axes[0].set_ylabel("mean conjunction selectivity")
    axes[0].set_ylim(0, 1)
    axes[0].set_title("E9: conjunction tuning self-organises\n(reward-free Hebbian)")
    axes[0].legend(fontsize=8)

    labels = [f"s{c//A}xr{c%A}" for c in range(K * A)]
    xp = np.arange(K * A)
    w = 0.27
    for i, k in enumerate(kinds):
        axes[1].bar(xp + (i - 1) * w, coverage[k] / n_seeds, w, label=k)
    axes[1].set_xticks(xp)
    axes[1].set_xticklabels(labels)
    axes[1].set_ylabel("hidden units (mean / seed)")
    axes[1].set_title("Combo tiling: units per conjunction class")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    p1 = os.path.join(FIGDIR, "e9_emergence.png")
    fig.savefig(p1, dpi=110)
    print("wrote", p1)

    # ---- figure 2: switching accuracy by basis ----
    fig, ax = plt.subplots(figsize=(7.5, 5))
    colors = {"emergent": "seagreen", "wired": "slategray", "frozen": "crimson"}
    for k in kinds:
        mv = moving_avg(switch_acc[k].mean(0), 15)
        ax.plot(np.arange(len(mv)), mv, color=colors[k], label=f"{k} (final {final[k]:.2f})")
    ax.axhline(0.5, ls="--", color="k", alpha=0.4, label="chance")
    ax.set_xlabel("trial (15-trial moving average)")
    ax.set_ylabel("switching accuracy")
    ax.set_ylim(0, 1)
    ax.set_title("E9: reward-driven switching on emergent vs wired vs frozen basis")
    ax.legend(fontsize=9)
    fig.tight_layout()
    p2 = os.path.join(FIGDIR, "e9_switching.png")
    fig.savefig(p2, dpi=110)
    print("wrote", p2)

    np.savez(os.path.join(DATADIR, "e9_data.npz"),
             conj_pre={k: np.array(conj_pre[k]) for k in kinds},
             conj_post={k: np.array(conj_post[k]) for k in kinds},
             coverage={k: coverage[k] for k in kinds},
             switch_acc={k: switch_acc[k] for k in kinds},
             emergence_steps=emergence[0][0] if emergence else np.array([]),
             emergence_curve=emergence[0][1] if emergence else np.array([]),
             n_blocks=n_blocks, block_len=block_len, allow_pickle=True)
    print("wrote", os.path.join(DATADIR, "e9_data.npz"))


if __name__ == "__main__":
    main()
