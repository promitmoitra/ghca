"""E6 mechanism animation -- three demon meters read one frozen medium.

The climax visual for the scrollytelling page (beat 09) and the dynamic
counterpart to docs/figures/e6_horde.png. E6's point is that memory / attention
/ executive are EMERGENT CATEGORIES, not modules: one unchanging Greenberg-
Hastings substrate carries three motifs at once --

  - two MEMORY rings   (a rotating pulse that sometimes dies mid-epoch),
  - an ATTENTION chain  (two waves collide; one captures the centre),
  - an EXECUTIVE ring   (a slow pulse sweeping toward the next context switch),

-- and three identical linear GVF demons read the SAME whole-network active mask,
each answering a different question. Left panel = the substrate (nodes coloured
by excitable state, dark-native palette so it reads on the page's black bg).
Right = three meters; each demon's value (colour) tracks its ground truth (grey)
over a sliding window with a "now" marker.

Reuses the real E6 machinery from experiments/e6_horde.py (build_substrate, the
schedule, features/td_continuing/mc_forecast/mc_return) -- only the phase stream
is additionally captured, for colouring.

Output: docs/figures/e6_horde.gif  (+ a still poster e6_horde_poster.png)
Run:    python3 experiments/e6_animation.py
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation, gridspec

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from ghca_net import Network
import e6_horde as e6

FIGDIR = os.path.join(ROOT, "docs", "figures")
os.makedirs(FIGDIR, exist_ok=True)

# dark-native state palette (this panel sits on black, unlike the light-field
# GIFs): rested = dim (off), active = bright red (on), refractory = fading tail.
RED, BLUE, GREEN = "#d4171c", "#4c8fd0", "#3fb37a"
INK, BG, PANEL = "#ececee", "#0b0b0d", "#141418"
_D_REST = np.array([0.20, 0.20, 0.24, 1.0])
_D_ACT = np.array([0.86, 0.11, 0.13, 1.0])
_D_RHEAD = np.array([0.52, 0.05, 0.09, 1.0])
_D_RTAIL = np.array([0.16, 0.13, 0.16, 1.0])


def state_colors(phi_row, act, tau):
    phi_row = np.asarray(phi_row)
    tau = np.broadcast_to(np.asarray(tau), phi_row.shape)
    rgba = np.tile(_D_REST, (len(phi_row), 1))
    active = (phi_row >= 1) & (phi_row <= act)
    refr = phi_row > act
    rgba[active] = _D_ACT
    if refr.any():
        span = np.maximum(tau[refr] - act, 1)
        frac = ((phi_row[refr] - act) / span)[:, None]
        rgba[refr] = (1 - frac) * _D_RHEAD + frac * _D_RTAIL
    return rgba


def run_capture(seed=0, T=2000):
    """Replicate e6.run_world but also capture net.phi each step (for colouring)."""
    W, roles, tau = e6.build_substrate(seed)
    net = Network(W, act=e6.ACT, pas=0, theta=np.ones(W.shape[0]), p_s=0.0, seed=seed)
    net.tau = tau.copy()
    net.phi[:] = 0
    rng = np.random.default_rng(seed + 3)
    N = roles["N"]
    mem, att, exe = roles["mem"], roles["att"], roles["exe"]
    net.phi[exe[0]] = 1

    phi = np.zeros((T, N), np.int64)
    mask = np.zeros((T, N), np.float32)
    mem_alive = np.zeros(T); att_winner = np.full(T, -1)
    att_forecast = np.zeros(T, bool); exe_switch = np.zeros(T)
    cur_stim = 0; winner = -1; collide_t = 0; s0 = s1 = -1
    for t in range(T):
        if t % e6.T_MEM == 0:
            for rg in mem:
                net.phi[rg] = 0
            cur_stim = int(rng.integers(2))
            dying = rng.random() < e6.P_DISRUPT
            net.tau[mem[cur_stim]] = e6.TAU_MEM_DIE if dying else e6.TAU_MEM_LIVE
            net.phi[mem[cur_stim][0]] = 1
        if t % e6.T_ATT == 0:
            net.phi[att] = 0
            bias = int(rng.integers(-3, 4)); j0, j1 = rng.integers(0, 2, size=2)
            s0 = t + max(0, -bias) + int(j0); s1 = t + max(0, bias) + int(j1)
            winner = 0 if s0 < s1 else (1 if s1 < s0 else int(rng.integers(2)))
            collide_t = max(s0, s1) + e6.CENTER
        if t % e6.T_SWITCH == 0 and t > 0:
            exe_switch[t] = 1.0
        drive = np.zeros(N, bool)
        if t == s0:
            drive[att[0]] = True
        if t == s1:
            drive[att[-1]] = True
        net.step(drive if drive.any() else None)
        phi[t] = net.phi
        mask[t] = net.active_mask().astype(np.float32)
        mem_alive[t] = 1.0 if net.active_mask()[mem[cur_stim]].any() else 0.0
        att_winner[t] = winner
        att_forecast[t] = (t >= max(s0, s1)) and (t < collide_t)
    sig = {"mem_alive": mem_alive, "att_winner": att_winner,
           "att_forecast": att_forecast, "exe_switch": exe_switch, "roles": roles}
    return phi, mask, sig, net.tau.copy(), roles


def node_positions(roles):
    """Custom 2-D layout for the four regions (indexed by node id)."""
    N = roles["N"]
    pos = np.zeros((N, 2))
    for k, rg in enumerate(roles["mem"]):        # two memory rings, upper-left
        cx, cy, r = -1.7, (1.35 - 1.5 * k), 0.62
        ang = 2 * np.pi * np.arange(len(rg)) / len(rg) + np.pi / 2
        pos[rg, 0] = cx + r * np.cos(ang); pos[rg, 1] = cy + r * np.sin(ang)
    att = roles["att"]                            # attention chain, bottom
    pos[att, 0] = np.linspace(-0.55, 1.35, len(att)); pos[att, 1] = -1.55
    exe = roles["exe"]                            # executive ring, right
    ang = 2 * np.pi * np.arange(len(exe)) / len(exe) + np.pi / 2
    pos[exe, 0] = 1.35 + 0.8 * np.cos(ang); pos[exe, 1] = 0.55 + 0.8 * np.sin(ang)
    return pos


def _norm(a):
    a = np.asarray(a, float)
    lo, hi = np.nanpercentile(a, 2), np.nanpercentile(a, 98)
    return np.clip((a - lo) / (hi - lo + 1e-9), 0, 1)


def main():
    phi, mask, sig, tau, roles = run_capture(seed=0, T=2000)
    X = e6.features(mask)
    split = len(mask) // 2
    w_mem = e6.td_continuing(X[:split], sig["mem_alive"][:split], gamma=0.9)
    w_exe = e6.td_continuing(X[:split], sig["exe_switch"][:split], gamma=0.8)
    w_att = e6.mc_forecast(X[:split], sig["att_forecast"][:split], sig["att_winner"][:split])
    V_mem = _norm(X @ w_mem); G_mem = _norm(e6.mc_return(sig["mem_alive"], 0.9))
    V_exe = _norm(X @ w_exe); G_exe = _norm(e6.mc_return(sig["exe_switch"], 0.8))
    V_att = np.clip(X @ w_att, -0.2, 1.2)

    pos = node_positions(roles)
    t0, t1, stride = 240, 400, 2
    frames = list(range(t0, t1, stride))
    WIN = 130

    fig = plt.figure(figsize=(9.6, 4.9), facecolor=BG)
    gs = gridspec.GridSpec(3, 2, width_ratios=[1.35, 1.0], height_ratios=[1, 1, 1],
                           wspace=0.18, hspace=0.5, left=0.02, right=0.975,
                           top=0.88, bottom=0.1)
    axS = fig.add_subplot(gs[:, 0]); axS.set_facecolor(BG)
    axS.set_xlim(-2.5, 2.35); axS.set_ylim(-2.2, 2.2); axS.set_aspect("equal"); axS.axis("off")
    scat = axS.scatter(pos[:, 0], pos[:, 1], s=64,
                       c=state_colors(phi[t0], e6.ACT, tau), edgecolors="none")
    axS.text(-1.7, 2.15, "MEMORY", color=BLUE, ha="center", fontsize=9, fontweight="bold")
    axS.text(0.4, -1.95, "ATTENTION", color=RED, ha="center", fontsize=9, fontweight="bold")
    axS.text(1.35, 1.65, "CONTROL", color=GREEN, ha="center", fontsize=9, fontweight="bold")

    meters = [("still remembering?", V_mem, G_mem, BLUE),
              ("who is winning?", V_att, None, RED),
              ("switch coming?", V_exe, G_exe, GREEN)]
    vlines, glines, nowlines = [], [], []
    for i, (label, V, G, col) in enumerate(meters):
        ax = fig.add_subplot(gs[i, 1]); ax.set_facecolor(PANEL)
        ax.set_xlim(0, WIN); ax.set_ylim(-0.15, 1.42); ax.set_xticks([]); ax.set_yticks([])
        for spn in ax.spines.values():
            spn.set_color("#2a2a30")
        ax.text(0.03, 0.86, label, transform=ax.transAxes, color=INK, fontsize=9.5,
                fontweight="bold")
        (gl,) = ax.plot([], [], color="#6a6a72", lw=1.4)
        (vl,) = ax.plot([], [], color=col, lw=2.2)
        nl = ax.axvline(WIN - 1, color=INK, lw=0.8, alpha=0.5)
        if label == "who is winning?":
            ax.axhline(0.5, ls=":", color="#6a6a72", lw=1)
        vlines.append(vl); glines.append(gl); nowlines.append(nl)

    fig.suptitle("E6 — three questions, one medium", color=INK, fontsize=13,
                 fontweight="bold", y=0.975)

    def update(t):
        arts = [scat]
        scat.set_facecolors(state_colors(phi[t], e6.ACT, tau))
        lo = max(0, t - WIN)
        xs = np.arange(t - lo)
        for i, (label, V, G, col) in enumerate(meters):
            vlines[i].set_data(xs, V[lo:t])
            glines[i].set_data(xs, G[lo:t]) if G is not None else glines[i].set_data([], [])
            nowlines[i].set_xdata([t - lo - 1, t - lo - 1])
            arts += [vlines[i], glines[i], nowlines[i]]
        return arts

    ani = animation.FuncAnimation(fig, update, frames=frames, blit=False)
    out = os.path.join(FIGDIR, "e6_horde.gif")
    ani.save(out, writer=animation.PillowWriter(fps=12), savefig_kwargs={"facecolor": BG})
    plt.close(fig)
    raw = os.path.getsize(out) / 1024

    # shrink for the web: downscale + adaptive 96-colour palette + optimize
    from PIL import Image, ImageSequence
    im = Image.open(out)
    TARGET_W = 680

    def prep(f):
        rgb = f.copy().convert("RGB")
        if rgb.width > TARGET_W:
            rgb = rgb.resize((TARGET_W, round(rgb.height * TARGET_W / rgb.width)), Image.LANCZOS)
        return rgb.quantize(colors=96, method=Image.FASTOCTREE)
    fr = [prep(f) for f in ImageSequence.Iterator(im)]
    fr[0].save(out, save_all=True, append_images=fr[1:],
               duration=im.info.get("duration", 83), loop=0, optimize=True)
    fr[0].convert("RGB").save(os.path.join(FIGDIR, "e6_horde_poster.png"))
    kb = os.path.getsize(out) / 1024
    print(f"wrote {out} ({raw:.0f} KB raw -> {kb:.0f} KB optimised, {len(frames)} frames)")


if __name__ == "__main__":
    main()
