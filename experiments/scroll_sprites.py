"""Sprite sheets for the scroll-scrubbed beats of the scrollytelling page.

The scrollytelling walkthrough (docs/scroll/index.html on the deploy-viz-page
branch) drives three animations *frame-by-frame off the scroll position* rather
than as autoplay loops. Each needs a frame sequence, not a sealed GIF, so this
script renders each as ONE sprite-sheet PNG (a grid of frames) plus the frame
geometry the page's <canvas> scrubber reads from data-* attributes:

    emergence.png  -- the raw spiral emerging from noise            (beat 02)
    triptych.png   -- one threshold, three fates: die/organise/boil (beat 03)
    necessity.png  -- intact vs do(theta)-ablated, side by side     (beat 09-> C6)

Frames are rendered directly from state_colors, upscaled, and annotated with crisp
panel headers, status badges, and step counters so the visual comparison is unmistakable.

Outputs land in docs/scroll/assets/.
"""

import os
import sys
import json
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))

from ghca_net import Network, lattice2d
from ghca_net_viz import state_colors
import e7_spiral_option as sp
import e7_learning as e7

OUT = os.path.join(ROOT, "docs", "scroll", "assets")
os.makedirs(OUT, exist_ok=True)


def _get_font(size=11):
    try:
        return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", size)
    except Exception:
        return ImageFont.load_default()


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
# beat 02 -- spiral emergence
# ---------------------------------------------------------------------------

def render_emergence(scale=7, stride=2, warmup=40, cols=8):
    net = Network(lattice2d(44, r=2), act=6, pas=8, theta=4.0, p_s=5e-3, seed=3)
    net.seed_random(0.15, 0.15)
    roll = net.run(150, record=True)["phi"]
    L = 44
    font = _get_font(12)

    frames = []
    total_t = roll.shape[0] - warmup
    for idx, t in enumerate(range(warmup, roll.shape[0], stride)):
        grid_img = _to_img(state_colors(roll[t], net.act, net.tau).reshape(L, L, 4), scale)
        gw, gh = grid_img.size

        hdr_h = 32
        ftr_h = 24
        out = Image.new("RGB", (gw, gh + hdr_h + ftr_h), (16, 16, 20))
        out.paste(grid_img, (0, hdr_h))

        draw = ImageDraw.Draw(out)
        draw.rectangle([0, 0, gw, hdr_h - 1], fill=(24, 24, 30))
        draw.text((10, 8), "ORGANIZED BAND (θ = 4.0)", fill=(230, 230, 240), font=font)
        draw.text((gw - 90, 8), f"STEP {t:03d}", fill=(212, 23, 28), font=font)

        draw.rectangle([0, gh + hdr_h, gw, gh + hdr_h + ftr_h], fill=(12, 12, 16))
        active_cnt = np.sum((roll[t] > 0) & (roll[t] <= net.act))
        pct = (active_cnt / (L * L)) * 100
        draw.text((10, gh + hdr_h + 5), f"Active: {active_cnt} ({pct:.1f}%) | Self-Sustaining Spiral", fill=(160, 160, 170), font=_get_font(10))

        frames.append(out)

    _sheet(frames, "emergence", cols=cols)


# ---------------------------------------------------------------------------
# beat 03 -- E0 threshold triptych: one knob (theta), three fates
# ---------------------------------------------------------------------------

def render_triptych(scale=4, stride=2, T=200, cols=8):
    L, GAP = 44, 4

    def roll(theta, seed=3, dens=0.30):
        net = Network(lattice2d(L, r=2, periodic=True), act=6, pas=8,
                      theta=theta, p_s=0.0, seed=seed)
        net.seed_random(dens, dens)
        return net.run(T, record=True)["phi"], net

    die, net = roll(6.0)      # subcritical -> dies
    org, _ = roll(4.0)        # organised band -> spirals
    sat, _ = roll(2.0)        # supercritical -> saturated / turbulent
    tau = net.tau

    pw = L * scale
    gw = 3 * pw + 2 * GAP * scale
    gh = L * scale
    hdr_h = 36
    ftr_h = 28

    font_hdr = _get_font(11)
    font_badge = _get_font(10)

    white = np.array([1.0, 1.0, 1.0, 1.0])
    frames = []

    for t in range(0, T, stride):
        c = np.zeros((L, 3 * L + 2 * GAP, 4))
        c[:, :L] = state_colors(die[t], 6, tau).reshape(L, L, 4)
        c[:, L:L + GAP] = white
        c[:, L + GAP:2 * L + GAP] = state_colors(org[t], 6, tau).reshape(L, L, 4)
        c[:, 2 * L + GAP:2 * L + 2 * GAP] = white
        c[:, 2 * L + 2 * GAP:] = state_colors(sat[t], 6, tau).reshape(L, L, 4)

        grid_img = _to_img(c, scale)
        out = Image.new("RGB", (gw, gh + hdr_h + ftr_h), (14, 14, 18))
        out.paste(grid_img, (0, hdr_h))

        draw = ImageDraw.Draw(out)
        draw.rectangle([0, 0, gw, hdr_h - 1], fill=(22, 22, 28))

        # Panel 1 Header
        draw.text((10, 10), "1. SUBCRITICAL (θ=6.0)", fill=(220, 100, 100), font=font_hdr)
        # Panel 2 Header
        draw.text((pw + GAP * scale + 10, 10), "2. CRITICAL (θ=4.0)", fill=(100, 220, 120), font=font_hdr)
        # Panel 3 Header
        draw.text((2 * pw + 2 * GAP * scale + 10, 10), "3. SUPERCRITICAL (θ=2.0)", fill=(240, 200, 80), font=font_hdr)

        # Bottom Badges
        draw.rectangle([0, gh + hdr_h, gw, gh + hdr_h + ftr_h], fill=(10, 10, 14))

        act_die = np.sum((die[t] > 0) & (die[t] <= 6))
        act_org = np.sum((org[t] > 0) & (org[t] <= 6))
        act_sat = np.sum((sat[t] > 0) & (sat[t] <= 6))

        txt_die = "DIES OUT" if act_die == 0 else f"FADING ({act_die})"
        col_die = (255, 80, 80) if act_die == 0 else (200, 120, 120)

        draw.text((10, gh + hdr_h + 7), f"✖ {txt_die}", fill=col_die, font=font_badge)
        draw.text((pw + GAP * scale + 10, gh + hdr_h + 7), f"✔ SPIRAL ({act_org})", fill=(100, 230, 120), font=font_badge)
        draw.text((2 * pw + 2 * GAP * scale + 10, gh + hdr_h + 7), f"⚠ BOIL/TURBULENT ({act_sat})", fill=(240, 200, 80), font=font_badge)

        frames.append(out)

    _sheet(frames, "triptych", cols=cols)


# ---------------------------------------------------------------------------
# beat 09 (causal test) -- intact vs do(theta)-ablated, side by side (C6)
# ---------------------------------------------------------------------------

def render_necessity(scale=5, stride=2, t_run=110, cols=8):
    L, GAP, CHIR = sp.L, 4, +1

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

    pw = L * scale
    gw = 2 * pw + GAP * scale
    gh = L * scale
    hdr_h = 36
    ftr_h = 28

    font_hdr = _get_font(11)
    font_badge = _get_font(10)

    white = np.array([1.0, 1.0, 1.0, 1.0])
    frames = []

    for t in range(0, t_run, stride):
        c = np.zeros((L, 2 * L + GAP, 4))
        c[:, :L] = state_colors(intact[t], sp.ACT, tau).reshape(L, L, 4)
        c[:, L:L + GAP] = white
        c[:, L + GAP:] = state_colors(ablated[t], sp.ACT, tau).reshape(L, L, 4)

        grid_img = _to_img(c, scale)
        out = Image.new("RGB", (gw, gh + hdr_h + ftr_h), (14, 14, 18))
        out.paste(grid_img, (0, hdr_h))

        draw = ImageDraw.Draw(out)
        draw.rectangle([0, 0, gw, hdr_h - 1], fill=(22, 22, 28))

        draw.text((10, 10), "LEFT: INTACT (θ = 4.0)", fill=(100, 220, 120), font=font_hdr)
        draw.text((pw + GAP * scale + 10, 10), "RIGHT: do(θ → 6.0) ABLATED", fill=(240, 100, 100), font=font_hdr)

        draw.rectangle([0, gh + hdr_h, gw, gh + hdr_h + ftr_h], fill=(10, 10, 14))

        act_intact = np.sum((intact[t] > 0) & (intact[t] <= sp.ACT))
        act_ablated = np.sum((ablated[t] > 0) & (ablated[t] <= sp.ACT))

        draw.text((10, gh + hdr_h + 7), f"✔ SPIRAL SUSTAINED ({act_intact})", fill=(100, 230, 120), font=font_badge)

        if act_ablated == 0 or t > 25:
            draw.text((pw + GAP * scale + 10, gh + hdr_h + 7), "⚠ SPIRAL DESTROYED (Extinguished)", fill=(255, 80, 80), font=font_badge)
        else:
            draw.text((pw + GAP * scale + 10, gh + hdr_h + 7), f"⚠ COLLAPSING... ({act_ablated})", fill=(240, 180, 80), font=font_badge)

        frames.append(out)

    _sheet(frames, "necessity", cols=cols)


def main():
    print("Rendering annotated sprite sheets for Beats 02, 03, 09...")
    render_emergence()
    render_triptych()
    render_necessity()
    print("Done!")


if __name__ == "__main__":
    main()
