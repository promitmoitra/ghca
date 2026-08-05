"""Figures for the lattice arc: the transmission edge, the tonic trade, the identity bind.

Three figures, one per load-bearing claim. Each reads only committed `result/stats/*.json`,
so the numbers in the notes and the numbers in the plots cannot drift apart.

  lattice_transmission_edge.png — the action primitive. Transmission is a monotone STEP whose
      position is the learned interval, and the edge tracks the reward delay 1:1. The right
      panel is the anti-coincidence check: unpaired's edge is pinned near its random-delay
      distribution mean regardless of D, which at D=50 alone made it indistinguishable from
      trained. Only the sweep separates them — the reason a single delay would have misled.

  lattice_tonic_trade.png — why tonic drive has no window. Reentrant flooding and τ-learning
      collapse arrive BEFORE the monotonicity violation they were supposed to buy, and neither
      a band-pass peak nor an avoidance dip ever appears.

  lattice_identity_bind.png — why plastic identity failed structurally. Differentiation and
      function are anti-correlated across every rate and both set-points: wherever θ moves
      enough to build a gradient, propagation reach and τ-learning break.

Style follows the repo's other figure scripts: plain matplotlib, Agg, dpi=110, dotted grey
reference lines. Per-seed spreads are drawn wherever the source JSON carries them (house rule:
report spreads, not just means).
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SD = os.path.join(ROOT, "result", "stats")
FIG = os.path.join(ROOT, "docs", "figures")
os.makedirs(FIG, exist_ok=True)

C_TRAINED, C_UNTRAINED, C_UNPAIRED = "#1b6ca8", "0.55", "#c0392b"
C_GOOD, C_BAD, C_NEUTRAL = "#1b6ca8", "#c0392b", "#e08214"


def load(name):
    p = os.path.join(SD, name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def pick(*names):
    """Prefer the n=20 artifact when present, else the n=3 one; report which was used."""
    for n in names:
        d = load(n)
        if d is not None:
            return d, n
    return None, None


# ---------------------------------------------------------------------------------
def fig_transmission(d, src):
    rows, delays = d["rows"], d["delays"]
    # one axis per delay PLUS one for the summary -- an earlier version reused the
    # last delay's axis and silently destroyed that panel.
    fig, ax = plt.subplots(1, len(delays) + 1, figsize=(4.7 * (len(delays) + 1), 4.6))

    for i, D in enumerate(delays):
        a = ax[i]
        for mode, col, lab in (("trained", C_TRAINED, "trained"),
                               ("untrained", C_UNTRAINED, "untrained (random τ)"),
                               ("unpaired", C_UNPAIRED, "unpaired (no contingency)")):
            r = rows[f"D{D}_{mode}"]
            M = np.array(r["motor"])
            t = np.array(r["t"], float)
            m, sd = M.mean(0), M.std(0)
            a.plot(t, m, "-o", ms=3.5, color=col, label=lab)
            a.fill_between(t, m - sd, m + sd, color=col, alpha=0.18, lw=0)
        a.axvline(D, ls=":", color="0.35")
        a.annotate(f"reward at D={D}", xy=(D, 0.0), xycoords=("data", "axes fraction"),
                   xytext=(4, 6), textcoords="offset points", fontsize=8, color="0.25",
                   rotation=90, va="bottom")
        a.set_xlabel("probe time"); a.set_title(f"D = {D}")
        if i == 0:
            a.set_ylabel("motor spikes / cell\n(probe's causal contribution)")
            a.legend(fontsize=8, loc="upper left")

    a = ax[len(delays)]
    for mode, col, lab in (("trained", C_TRAINED, "trained"),
                           ("untrained", C_UNTRAINED, "untrained"),
                           ("unpaired", C_UNPAIRED, "unpaired")):
        cr = [rows[f"D{D}_{mode}"]["cross"] for D in delays]
        a.plot(delays, cr, "-o", ms=5, color=col, label=lab)
    a.plot(delays, delays, ls="--", color="0.3", lw=1, label="edge = D (identity)")
    # unpaired's random reward delays average here, which is why it collides with trained at D=50
    a.axhline(52, ls=":", color=C_UNPAIRED, lw=1)
    a.annotate("unpaired's edge sits near its random-delay\n"
               "mean (~52), not D. At D=50 alone it is\n"
               "indistinguishable from trained.",
               xy=(0.03, 0.97), xycoords="axes fraction", fontsize=7.5,
               color=C_UNPAIRED, va="top")
    a.set_xlabel("reward delay D"); a.set_ylabel("transmission edge (half-max crossing)")
    a.set_title("the edge tracks the learned interval")
    a.legend(fontsize=8, loc="lower right")

    fig.suptitle("The action is transmission: a step whose position is the learned interval"
                 f"   ({src}, n={d['n']}; bands = ±1 sd across seeds)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = os.path.join(FIG, "lattice_transmission_edge.png")
    fig.savefig(out, dpi=110); plt.close(fig)
    return out


# ---------------------------------------------------------------------------------
def fig_tonic(d, src):
    q = np.array(d["q_list"], float)
    x = np.arange(len(q))
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    g = lambda rule, k: np.array([d["rows"][f"{rule}_q{v:g}"][k] for v in q], float)
    gm = lambda rule, k: np.array([np.nanmean(d["rows"][f"{rule}_q{v:g}"][k]) for v in q])
    gs = lambda rule, k: np.array([np.nanstd(d["rows"][f"{rule}_q{v:g}"][k]) for v in q])

    a = ax[0]
    w2, sdw = gm("approach", "aligned_frac"), gs("approach", "aligned_frac")
    a.plot(x, w2, "-o", ms=4, color=C_GOOD, label="τ within ±2 (learning)")
    a.fill_between(x, w2 - sdw, w2 + sdw, color=C_GOOD, alpha=0.18, lw=0)
    re, sdr = gm("approach", "reentrant"), gs("approach", "reentrant")
    a.plot(x, re, "-s", ms=4, color=C_BAD, label="reentrant activity (flooding)")
    a.fill_between(x, re - sdr, re + sdr, color=C_BAD, alpha=0.18, lw=0)
    drops = np.array([d["rows"][f"approach_q{v:g}"]["max_drop"] for v in q], float)
    # first APPRECIABLE violation, not the first numerically-nonzero one: at the latter
    # (q=3e-7) learning is still 0.98, so labelling it "already halved" would be false.
    iv = int(np.argmax(drops >= 0.01)) if (drops >= 0.01).any() else len(q) - 1
    a.axvline(iv, ls=":", color="0.35")
    a.annotate(f"monotonicity violation reaches +0.01 here\n"
               f"(learning already down to {w2[iv]:.2f}; flooding {re[iv]:.2f})",
               xy=(iv, 0.5), xytext=(8, 0), textcoords="offset points",
               fontsize=8, color="0.25")
    a.set_xticks(x); a.set_xticklabels([("0" if v == 0 else f"{v:g}") for v in q],
                                       rotation=45, fontsize=8)
    a.set_xlabel("tonic drive q"); a.set_ylabel("fraction")
    a.set_title("the cost arrives before the benefit"); a.legend(fontsize=8)

    a = ax[1]
    for rule, col in (("approach", C_GOOD), ("avoid-margin", C_NEUTRAL)):
        dr = np.array([d["rows"][f"{rule}_q{v:g}"]["max_drop"] for v in q], float)
        a.plot(x, dr, "-o", ms=4, color=col, label=f"{rule}: monotonicity drop")
    a.axhline(0, color="k", lw=0.8)
    a.set_xticks(x); a.set_xticklabels([("0" if v == 0 else f"{v:g}") for v in q],
                                       rotation=45, fontsize=8)
    a.set_xlabel("tonic drive q"); a.set_ylabel("max decrease along the curve")
    a.set_title("the violation it buys is marginal (≤ +0.05)"); a.legend(fontsize=8)

    a = ax[2]
    pk = np.array([d["rows"][f"approach_q{v:g}"]["peakiness"] for v in q], float)
    dp = np.array([d["rows"][f"avoid-margin_q{v:g}"]["dip"] for v in q], float)
    a.plot(x, np.nan_to_num(pk, nan=0.0), "-o", ms=4, color=C_GOOD,
           label="band-pass peakiness (need > 1)")
    a.axhline(1.0, ls="--", color=C_GOOD, lw=1)
    nan_dip = np.isnan(dp)
    a.plot(x[nan_dip], np.zeros(nan_dip.sum()), "x", ms=9, color=C_BAD,
           label="avoidance dip: absent (nan) at every q")
    a.set_ylim(-0.15, 1.4)
    a.set_xticks(x); a.set_xticklabels([("0" if v == 0 else f"{v:g}") for v in q],
                                       rotation=45, fontsize=8)
    a.set_xlabel("tonic drive q")
    a.set_title("neither stalled thread unlocks"); a.legend(fontsize=8, loc="upper left")

    fig.suptitle("Tonic drive: no window. Flooding and learning-collapse precede the "
                 f"monotonicity violation they were meant to buy   ({src}, n={d['n']})",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = os.path.join(FIG, "lattice_tonic_trade.png")
    fig.savefig(out, dpi=110); plt.close(fig)
    return out


# ---------------------------------------------------------------------------------
def fig_identity(mv, kp, srcs):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    sets = [("a* = 0.08  (θ free to move)", mv, "-o", C_BAD),
            ("a* = 0.16  (function preserved)", kp, "-s", C_GOOD)]

    a = ax[0]
    for lab, d, mk, col in sets:
        r = np.array(d["r_list"], float)
        cc = np.array([np.nanmean(d["rows"][f"r{v:g}"]["corr_theta_x"]) for v in r])
        sd = np.array([np.nanstd(d["rows"][f"r{v:g}"]["corr_theta_x"]) for v in r])
        xx = np.arange(len(r))
        a.plot(xx, cc, mk[1], ms=4, ls="-", color=col, label=lab)
        a.fill_between(xx, cc - sd, cc + sd, color=col, alpha=0.18, lw=0)
    a.axhline(0, color="k", lw=0.8)
    a.set_xticks(np.arange(len(r))); a.set_xticklabels([f"{v:g}" for v in r], fontsize=8)
    a.set_xlabel("identity rate r"); a.set_ylabel("corr(θ, x)   (more negative = gradient)")
    a.set_title("a gradient does form"); a.legend(fontsize=8, loc="lower left")

    a = ax[1]
    for lab, d, mk, col in sets:
        r = np.array(d["r_list"], float); xx = np.arange(len(r))
        rc = np.array([np.nanmean(d["rows"][f"r{v:g}"]["reach"]) for v in r])
        sd = np.array([np.nanstd(d["rows"][f"r{v:g}"]["reach"]) for v in r])
        a.plot(xx, rc, mk[1], ms=4, ls="-", color=col, label=lab)
        a.fill_between(xx, rc - sd, rc + sd, color=col, alpha=0.18, lw=0)
    a.set_xticks(np.arange(len(r))); a.set_xticklabels([f"{v:g}" for v in r], fontsize=8)
    a.set_xlabel("identity rate r"); a.set_ylabel("wave propagation reach")
    a.set_title("...but function is lost where it does"); a.legend(fontsize=8, loc="lower left")

    a = ax[2]
    for lab, d, mk, col in sets:
        r = np.array(d["r_list"], float)
        cc = np.array([abs(np.nanmean(d["rows"][f"r{v:g}"]["corr_theta_x"])) for v in r])
        rc = np.array([np.nanmean(d["rows"][f"r{v:g}"]["reach"]) for v in r])
        a.plot(cc, rc, mk[1], ms=6, ls="-", color=col, alpha=0.85, label=lab)
        for j, v in enumerate(r):
            a.annotate(f"r={v:g}", (cc[j], rc[j]), xytext=(4, 4),
                       textcoords="offset points", fontsize=7, color=col)
    a.set_xlabel("|corr(θ, x)|   differentiation →")
    a.set_ylabel("propagation reach   function →")
    a.set_title("the bind: no point is high on both axes")
    a.legend(fontsize=8, loc="upper right")
    a.margins(x=0.12, y=0.12)

    fig.suptitle("Plastic identity: differentiation and function are anti-correlated at every "
                 f"rate and both set-points   ({', '.join(srcs)}, n={mv['n']})", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = os.path.join(FIG, "lattice_identity_bind.png")
    fig.savefig(out, dpi=110); plt.close(fig)
    return out


def main():
    made = []
    sm, src = pick("lattice_sensorimotor_n20.json", "lattice_sensorimotor.json")
    if sm:
        made.append(fig_transmission(sm, src))
    else:
        print("skip transmission: no sensorimotor artifact", flush=True)

    to, src = pick("lattice_tonic_n5.json", "lattice_tonic.json")
    if to:
        made.append(fig_tonic(to, src))
    else:
        print("skip tonic: no artifact", flush=True)

    mv, s1 = pick("lattice_identity_n20_move.json", "lattice_identity_move.json")
    kp, s2 = pick("lattice_identity_n20_keep.json", "lattice_identity_keep.json")
    if mv and kp:
        made.append(fig_identity(mv, kp, [s1, s2]))
    else:
        print("skip identity: need both set-point artifacts", flush=True)

    for m in made:
        print("wrote", os.path.relpath(m, ROOT), flush=True)


if __name__ == "__main__":
    main()
