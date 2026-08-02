"""Clean role-based sprite sheets for the plasticity beats of the scrollytelling page.

Renders the closed-loop plasticity substrate by its actual roles --
  input (sensory, left)  ->  recurrent hidden ring (centre)  ->  output (motor, right)
-- so "the signal routes to a target" is legible, instead of cramming 198 nodes
into a featureless square. No baked-in text (the page's HTML caption narrates);
dark-native palette (dim = rest, bright red = active). Motor nodes carry a
short persistence trace so the output channel that wins stays visibly lit.

    plasticity.png -- 3 substrates side by side: unlearned / weight-only / tri-axis
    retention.png  -- 2 substrates: weight-only (forgets) / tri-axis (retains)

Outputs -> docs/scroll/assets/.  Run from a checkout with the substrate code.
"""
import os, sys, io, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "experiments"))
from ghca_plasticity import make_closed_loop_learner

OUT = os.path.join(ROOT, "docs", "scroll", "assets")
os.makedirs(OUT, exist_ok=True)

BG = "#0b0b0d"
_REST = np.array([0.20, 0.20, 0.24]); _ACT = np.array([0.86, 0.11, 0.13])
_RH = np.array([0.52, 0.05, 0.09]); _RT = np.array([0.16, 0.13, 0.16])


def _state_colors(phi, act, tau):
    phi = np.asarray(phi); tau = np.broadcast_to(np.asarray(tau), phi.shape)
    c = np.tile(_REST, (len(phi), 1))
    a = (phi >= 1) & (phi <= act); r = phi > act
    c[a] = _ACT
    if r.any():
        span = np.maximum(tau[r] - act, 1); f = ((phi[r] - act) / span)[:, None]
        c[r] = (1 - f) * _RH + f * _RT
    return c


def _layout(roles):
    """Positions + sizes indexed by node id: sensory left, hidden ring, motor right."""
    N = roles["N"]; pos = np.zeros((N, 2)); sz = np.full(N, 10.0)
    for ch, arr in enumerate(roles["sensory"]):
        y0 = 0.62 if ch == 0 else -0.62
        pos[arr, 0] = -1.82 + 0.15 * (np.arange(len(arr)) % 2)
        pos[arr, 1] = y0 + np.linspace(-0.30, 0.30, len(arr)); sz[arr] = 60
    hid = roles["hidden"]
    ang = 2 * np.pi * np.arange(len(hid)) / len(hid) + np.pi / 2
    pos[hid, 0] = 0.92 * np.cos(ang); pos[hid, 1] = 0.92 * np.sin(ang); sz[hid] = 9
    for ch, arr in enumerate(roles["motor"]):
        y0 = 0.62 if ch == 0 else -0.62
        pos[arr, 0] = 1.82 - 0.15 * (np.arange(len(arr)) % 2)
        pos[arr, 1] = y0 + np.linspace(-0.30, 0.30, len(arr)); sz[arr] = 88
    return pos, sz


def _train(net, rng, n, reversal=False):
    for _ in range(n):
        x = int(rng.integers(net.roles["K"])); tgt = (1 - x) if reversal else x
        net.reset_traces()
        for _ in range(4):
            net.step_learn(net.sensory_drive(x))
        sc = sum(net.motor_scores() for _ in range(6))
        net.learn((1.0 if np.argmax(sc) == tgt else 0.0) - 0.5)


def _rollout(net, x_test, T, drive_steps=9):
    """Drive input x_test; per frame colour the substrate by GH state, and light the
    OUTPUT channel the readout is currently routing to (the winner glows, the loser
    stays dim) -- that decision, not raw motor firing, is the 'routes to target' signal."""
    roles = net.roles
    cum = np.zeros(roles["A"])                    # cumulative per-channel motor score
    frames = []
    for t in range(T):
        net.step_learn(net.sensory_drive(x_test) if t < drive_steps else None)
        col = _state_colors(net.phi, net.act, net.tau)
        cum += np.asarray(net.motor_scores(), float)
        # light the output channel the readout is routing to (winner glows).
        b = cum / cum.max() if cum.max() > 1e-9 else np.zeros_like(cum)
        for ch, arr in enumerate(roles["motor"]):
            g = float(b[ch]) ** 1.6
            col[arr] = (1 - g) * _REST + g * _ACT
        frames.append(col)
    return frames


def _panel(ax, pos, sz, col):
    ax.set_facecolor(BG)
    ax.scatter(pos[:, 0], pos[:, 1], s=sz, c=col, edgecolors="none")
    ax.set_xlim(-2.15, 2.15); ax.set_ylim(-1.25, 1.25)
    ax.set_aspect("equal"); ax.axis("off")


def _fig_to_img(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", facecolor=BG, dpi=95)
    plt.close(fig); buf.seek(0); return Image.open(buf).convert("RGB")


def _sheet(frames_img, name, cols=8):
    n = len(frames_img); fw, fh = frames_img[0].size; rows = (n + cols - 1) // cols
    sheet = Image.new("RGB", (cols * fw, rows * fh), (11, 11, 13))
    for i, im in enumerate(frames_img):
        r, c = divmod(i, cols); sheet.paste(im, (c * fw, r * fh))
    # web-shrink: quantise
    q = sheet.quantize(colors=64, method=Image.FASTOCTREE).convert("RGB")
    png = os.path.join(OUT, f"{name}.png"); q.save(png, optimize=True)
    json.dump({"frame_w": fw, "frame_h": fh, "cols": cols, "count": n},
              open(os.path.join(OUT, f"{name}.json"), "w"))
    print(f"{name}: {n} frames {fw}x{fh} grid {cols}x{rows} -> {os.path.getsize(png)//1024} KB")


def render_plasticity(T=44, cols=8):
    seed = 42
    nets = [make_closed_loop_learner(axes=(), seed=seed),
            make_closed_loop_learner(axes=("w",), seed=seed),
            make_closed_loop_learner(axes=("tau", "theta", "w"), seed=seed)]
    for i, net in enumerate(nets):
        if i > 0:
            _train(net, np.random.default_rng(seed), 130)
    pos, sz = _layout(nets[0].roles)
    rolls = [_rollout(net, 0, T) for net in nets]
    imgs = []
    for t in range(T):
        fig, axs = plt.subplots(1, 3, figsize=(9.6, 2.7), facecolor=BG)
        fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005, wspace=0.03)
        for k in range(3):
            _panel(axs[k], pos, sz, rolls[k][t])
        imgs.append(_fig_to_img(fig))
    _sheet(imgs, "plasticity", cols=cols)


def render_retention(T=44, cols=8):
    seed = 7
    w_net = make_closed_loop_learner(axes=("w",), seed=seed)
    m_net = make_closed_loop_learner(axes=("tau", "theta", "w"), seed=seed)
    for net in (w_net, m_net):
        rng = np.random.default_rng(seed)
        _train(net, rng, 130, reversal=False)   # Task A
        _train(net, rng, 130, reversal=True)    # Task B
    pos, sz = _layout(w_net.roles)
    rw, rm = _rollout(w_net, 0, T), _rollout(m_net, 0, T)   # re-test Task A cue
    imgs = []
    for t in range(T):
        fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.7), facecolor=BG)
        fig.subplots_adjust(left=0.008, right=0.992, top=0.995, bottom=0.005, wspace=0.04)
        _panel(axs[0], pos, sz, rw[t]); _panel(axs[1], pos, sz, rm[t])
        imgs.append(_fig_to_img(fig))
    _sheet(imgs, "retention", cols=cols)


if __name__ == "__main__":
    render_plasticity()
    render_retention()
