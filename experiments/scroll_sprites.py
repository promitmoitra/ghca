"""Sprite sheets for the scroll-scrubbed beats of the scrollytelling page.

The scrollytelling walkthrough (docs/scroll/index.html on the deploy-viz-page
branch) drives three animations *frame-by-frame off the scroll position* rather
than as autoplay loops. Each needs a frame sequence, not a sealed GIF, so this
script renders each as ONE sprite-sheet PNG (a grid of frames) plus the frame
geometry the page's <canvas> scrubber reads from data-* attributes:

    emergence.png  -- the raw spiral emerging from noise            (beat 02)
    triptych.png   -- one threshold, three fates: die/organise/boil (beat 03)
    necessity.png  -- intact vs do(theta)-ablated, side by side     (beat 09-> C6)

Frames are rendered directly from state_colors (no matplotlib chrome) and
nearest-neighbour upscaled to keep the lattice crisp. Outputs land in
docs/scroll/assets/ (curated onto the deploy-viz-page branch for the site).

Run from a checkout that has the substrate code (this `main` branch):
    python3 experiments/scroll_sprites.py
"""
import os
import sys
import json
import numpy as np
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from ghca_net import Network, lattice2d
from ghca_net_viz import state_colors
import e7_spiral_option as sp
import e7_learning as e7

OUT = os.path.join(ROOT, "docs", "scroll", "assets")
os.makedirs(OUT, exist_ok=True)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _to_img(colors_hw4, scale):
    """(H,W,4) float RGBA in 0..1 -> nearest-upscaled PIL RGB image."""
    arr = np.clip(colors_hw4[..., :3] * 255.0, 0, 255).astype(np.uint8)
    img = Image.fromarray(arr, mode="RGB")
    return img.resize((img.width * scale, img.height * scale), Image.NEAREST)


def _sheet(frames, name, cols=8):
    """Tile equal-size PIL frames into a grid sprite sheet + a JSON of geometry."""
    n = len(frames)
    fw, fh = frames[0].size
    rows = (n + cols - 1) // cols
    sheet = Image.new("RGB", (cols * fw, rows * fh), (10, 10, 12))
    for i, fr in enumerate(frames):
        r, c = divmod(i, cols)
        sheet.paste(fr, (c * fw, r * fh))
    png = os.path.join(OUT, f"{name}.png")
    sheet.save(png, optimize=True)
    with open(os.path.join(OUT, f"{name}.json"), "w") as fh_:
        json.dump({"frame_w": fw, "frame_h": fh, "cols": cols, "count": n}, fh_)
    print(f"{name}: {n} frames {fw}x{fh}  grid {cols}x{rows}  "
          f"-> {png} ({os.path.getsize(png) // 1024} KB)")


# ---------------------------------------------------------------------------
# beat 02 -- spiral emergence (the ghca_net_viz __main__ demo operating point)
# ---------------------------------------------------------------------------

def render_emergence(scale=7, stride=2, warmup=40, cols=8):
    net = Network(lattice2d(44, r=2), act=6, pas=8, theta=4.0, p_s=5e-3, seed=3)
    net.seed_random(0.15, 0.15)
    roll = net.run(150, record=True)["phi"]
    L = 44
    frames = [_to_img(state_colors(roll[t], net.act, net.tau).reshape(L, L, 4), scale)
              for t in range(warmup, roll.shape[0], stride)]
    _sheet(frames, "emergence", cols=cols)


# ---------------------------------------------------------------------------
# beat 03 -- E0 threshold triptych: one knob (theta), three fates
# same substrate point (r=2, a=6, tau=14), p_s=0 so nothing is noise-driven;
# only theta differs:  theta=6 dies / theta=4 organises / theta=2 saturates.
# ---------------------------------------------------------------------------

def render_triptych(scale=4, stride=2, T=200, cols=8):
    L, GAP = 44, 3

    def roll(theta, seed=3, dens=0.30):
        net = Network(lattice2d(L, r=2, periodic=True), act=6, pas=8,
                      theta=theta, p_s=0.0, seed=seed)
        net.seed_random(dens, dens)
        return net.run(T, record=True)["phi"], net

    die, net = roll(6.0)      # subcritical -> dies
    org, _ = roll(4.0)        # organised band -> spirals
    sat, _ = roll(2.0)        # supercritical -> saturated / turbulent
    tau = net.tau
    cols_px = 3 * L + 2 * GAP
    white = np.array([1.0, 1.0, 1.0, 1.0])
    frames = []
    for t in range(0, T, stride):
        c = np.zeros((L, cols_px, 4))
        c[:, :L] = state_colors(die[t], 6, tau).reshape(L, L, 4)
        c[:, L:L + GAP] = white
        c[:, L + GAP:2 * L + GAP] = state_colors(org[t], 6, tau).reshape(L, L, 4)
        c[:, 2 * L + GAP:2 * L + 2 * GAP] = white
        c[:, 2 * L + 2 * GAP:] = state_colors(sat[t], 6, tau).reshape(L, L, 4)
        frames.append(_to_img(c, scale))
    _sheet(frames, "triptych", cols=cols)


# ---------------------------------------------------------------------------
# beat 09 (causal test) -- intact vs do(theta)-ablated, side by side (C6)
# ---------------------------------------------------------------------------

def render_necessity(scale=5, stride=2, t_run=110, cols=8):
    L, GAP, CHIR = sp.L, 3, +1

    def one(ablate):
        net = e7.make_spiral(0, ablate=ablate)
        rng = np.random.default_rng(0)
        sp.seed_spiral(net, CHIR, jitter=0.02, rng=rng)
        phi = np.zeros((t_run, L * L), np.int64)
        for t in range(t_run):
            net.step(None)
            phi[t] = net.phi
        return phi, net.tau

    intact, tau = one(False)
    ablated, _ = one(True)
    cols_px = 2 * L + GAP
    white = np.array([1.0, 1.0, 1.0, 1.0])
    frames = []
    for t in range(0, t_run, stride):
        c = np.zeros((L, cols_px, 4))
        c[:, :L] = state_colors(intact[t], sp.ACT, tau).reshape(L, L, 4)
        c[:, L:L + GAP] = white
        c[:, L + GAP:] = state_colors(ablated[t], sp.ACT, tau).reshape(L, L, 4)
        frames.append(_to_img(c, scale))
    _sheet(frames, "necessity", cols=cols)


def main():
    render_emergence()
    render_triptych()
    render_necessity()


if __name__ == "__main__":
    main()
