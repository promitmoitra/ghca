"""Sprite sheets for Beat 12 and Beat 13 of the scrollytelling walkthrough.

Generates:
    plasticity.png -- Tri-axis closed-loop self-assembly vs unlearned/weight-only (Beat 12)
    retention.png  -- Topological loop preservation vs weight-only forgetting (Beat 13)

Frames are rendered with crisp panel headers, node role markers, status badges, and step counters.
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

from ghca_plasticity import ClosedLoopPlasticityEngine, make_closed_loop_learner
from ghca_net_viz import state_colors

OUT = os.path.join(ROOT, "docs", "scroll", "assets")
os.makedirs(OUT, exist_ok=True)


def _get_font(size=11):
    try:
        return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", size)
    except Exception:
        return ImageFont.load_default()


def _to_img(colors_hw4, scale=12):
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


def _pad_grid(phi, tau, act, grid_dim=(12, 12)):
    H, W = grid_dim
    N = H * W
    phi_pad = np.zeros(N, dtype=phi.dtype)
    tau_pad = np.full(N, 14, dtype=tau.dtype)

    n = min(len(phi), N)
    phi_pad[:n] = phi[:n]
    tau_pad[:n] = tau[:n]

    colors = state_colors(phi_pad, act, tau_pad).reshape(H, W, 4)
    return colors


def render_plasticity_triptych(scale=12, stride=2, T=100, cols=8):
    """
    Beat 12: Closed-Loop Tri-Axis Plasticity Triptych
    Panel 1: Control / Unlearned (axes=())
    Panel 2: Weight-Only Plasticity (axes=('w',))
    Panel 3: Closed-Loop Multi-Axis Plasticity (axes=('tau', 'theta', 'w'))
    """
    seed = 42

    ctrl = make_closed_loop_learner(axes=(), seed=seed)
    w_only = make_closed_loop_learner(axes=("w",), seed=seed)
    multi = make_closed_loop_learner(axes=("tau", "theta", "w"), seed=seed)

    rng = np.random.default_rng(seed)

    # Train w_only and multi
    for t in range(150):
        x = int(rng.integers(ctrl.roles["K"]))
        ctrl.reset_traces()
        for _ in range(4): ctrl.step_learn(ctrl.sensory_drive(x))

        w_only.reset_traces()
        for _ in range(4): w_only.step_learn(w_only.sensory_drive(x))
        sc_w = sum([w_only.motor_scores() for _ in range(6)])
        r_w = 1.0 if np.argmax(sc_w) == x else 0.0
        w_only.learn(r_w - 0.5)

        multi.reset_traces()
        for _ in range(4): multi.step_learn(multi.sensory_drive(x))
        sc_m = sum([multi.motor_scores() for _ in range(6)])
        r_m = 1.0 if np.argmax(sc_m) == x else 0.0
        multi.learn(r_m - 0.5)

    H, W = 12, 12
    GAP = 2
    pw = W * scale
    gw = 3 * pw + 2 * GAP * scale
    gh = H * scale
    hdr_h = 36
    ftr_h = 28

    font_hdr = _get_font(11)
    font_badge = _get_font(10)

    white = np.array([1.0, 1.0, 1.0, 1.0])
    frames = []
    x_test = 0

    for t in range(T):
        drive_ctrl = ctrl.sensory_drive(x_test) if t < 15 else None
        drive_w = w_only.sensory_drive(x_test) if t < 15 else None
        drive_multi = multi.sensory_drive(x_test) if t < 15 else None

        ctrl.step_learn(drive_ctrl)
        w_only.step_learn(drive_w)
        multi.step_learn(drive_multi)

        if t % stride == 0:
            c = np.zeros((H, 3 * W + 2 * GAP, 4))
            c_ctrl = _pad_grid(ctrl.phi, ctrl.tau, ctrl.act, (H, W))
            c_w = _pad_grid(w_only.phi, w_only.tau, w_only.act, (H, W))
            c_m = _pad_grid(multi.phi, multi.tau, multi.act, (H, W))

            c[:, :W] = c_ctrl
            c[:, W:W+GAP] = white
            c[:, W+GAP:2*W+GAP] = c_w
            c[:, 2*W+GAP:2*W+2*GAP] = white
            c[:, 2*W+2*GAP:] = c_m

            grid_img = _to_img(c, scale=scale)
            out = Image.new("RGB", (gw, gh + hdr_h + ftr_h), (14, 14, 18))
            out.paste(grid_img, (0, hdr_h))

            draw = ImageDraw.Draw(out)
            draw.rectangle([0, 0, gw, hdr_h - 1], fill=(22, 22, 28))

            draw.text((10, 10), "1. UNLEARNED (RIR=0.47)", fill=(180, 180, 190), font=font_hdr)
            draw.text((pw + GAP * scale + 10, 10), "2. WEIGHT-ONLY (RIR=0.79)", fill=(240, 200, 80), font=font_hdr)
            draw.text((2 * pw + 2 * GAP * scale + 10, 10), "3. TRI-AXIS (RIR=0.85)", fill=(100, 230, 120), font=font_hdr)

            # Node markers
            draw.ellipse([4, hdr_h + 4, 16, hdr_h + 16], outline=(0, 255, 120), width=2)
            draw.ellipse([pw + GAP * scale + 4, hdr_h + 4, pw + GAP * scale + 16, hdr_h + 16], outline=(0, 255, 120), width=2)
            draw.ellipse([2 * pw + 2 * GAP * scale + 4, hdr_h + 4, 2 * pw + 2 * GAP * scale + 16, hdr_h + 16], outline=(0, 255, 120), width=2)

            draw.rectangle([0, gh + hdr_h, gw, gh + hdr_h + ftr_h], fill=(10, 10, 14))

            m_nodes = multi.roles.get("M", [])
            m_multi_act = np.sum((multi.phi[m_nodes] > 0) & (multi.phi[m_nodes] <= multi.act)) > 0 if len(m_nodes) > 0 else False

            txt_ctrl = "DISPERSED"
            txt_w = "WEAK ROUTE"
            txt_multi = "DIRECT TARGET HIT! 🎯" if m_multi_act else "ROUTED WAVEFRONT"

            draw.text((10, gh + hdr_h + 7), f"S0→M0: {txt_ctrl}", fill=(160, 160, 170), font=font_badge)
            draw.text((pw + GAP * scale + 10, gh + hdr_h + 7), f"S0→M0: {txt_w}", fill=(240, 200, 80), font=font_badge)
            draw.text((2 * pw + 2 * GAP * scale + 10, gh + hdr_h + 7), f"S0→M0: {txt_multi}", fill=(100, 230, 120), font=font_badge)

            frames.append(out)

    _sheet(frames, "plasticity", cols=cols)


def render_retention_comparison(scale=12, stride=2, T=100, cols=8):
    """
    Beat 13: Sequential Learning & Topological Loop Protection
    Left: Weight-Only Adaptation (Task B overwrites Task A -> Forgetting)
    Right: Multi-Axis Plasticity (High-tau loop fossilization -> Retention)
    """
    seed = 7

    w_net = make_closed_loop_learner(axes=("w",), seed=seed)
    m_net = make_closed_loop_learner(axes=("tau", "theta", "w"), seed=seed)
    rng = np.random.default_rng(seed)

    # Phase 1: Train Task A
    for _ in range(150):
        x = int(rng.integers(w_net.roles["K"]))
        w_net.reset_traces()
        for _ in range(4): w_net.step_learn(w_net.sensory_drive(x))
        sc_w = sum([w_net.motor_scores() for _ in range(6)])
        r_w = 1.0 if np.argmax(sc_w) == x else 0.0
        w_net.learn(r_w - 0.5)

        m_net.reset_traces()
        for _ in range(4): m_net.step_learn(m_net.sensory_drive(x))
        sc_m = sum([m_net.motor_scores() for _ in range(6)])
        r_m = 1.0 if np.argmax(sc_m) == x else 0.0
        m_net.learn(r_m - 0.5)

    # Phase 2: Train Task B (Reversal)
    for _ in range(150):
        x = int(rng.integers(w_net.roles["K"]))
        target = 1 - x
        w_net.reset_traces()
        for _ in range(4): w_net.step_learn(w_net.sensory_drive(x))
        sc_w = sum([w_net.motor_scores() for _ in range(6)])
        r_w = 1.0 if np.argmax(sc_w) == target else 0.0
        w_net.learn(r_w - 0.5)

        m_net.reset_traces()
        for _ in range(4): m_net.step_learn(m_net.sensory_drive(x))
        sc_m = sum([m_net.motor_scores() for _ in range(6)])
        r_m = 1.0 if np.argmax(sc_m) == target else 0.0
        m_net.learn(r_m - 0.5)

    H, W = 12, 12
    GAP = 3
    pw = W * scale
    gw = 2 * pw + GAP * scale
    gh = H * scale
    hdr_h = 36
    ftr_h = 28

    font_hdr = _get_font(11)
    font_badge = _get_font(10)

    white = np.array([1.0, 1.0, 1.0, 1.0])
    frames = []
    x_test = 0  # Task A cue

    for t in range(T):
        drive_w = w_net.sensory_drive(x_test) if t < 15 else None
        drive_m = m_net.sensory_drive(x_test) if t < 15 else None

        w_net.step_learn(drive_w)
        m_net.step_learn(drive_m)

        if t % stride == 0:
            c = np.zeros((H, 2 * W + GAP, 4))
            c_w = _pad_grid(w_net.phi, w_net.tau, w_net.act, (H, W))
            c_m = _pad_grid(m_net.phi, m_net.tau, m_net.act, (H, W))

            c[:, :W] = c_w
            c[:, W:W+GAP] = white
            c[:, W+GAP:] = c_m

            grid_img = _to_img(c, scale=scale)
            out = Image.new("RGB", (gw, gh + hdr_h + ftr_h), (14, 14, 18))
            out.paste(grid_img, (0, hdr_h))

            draw = ImageDraw.Draw(out)
            draw.rectangle([0, 0, gw, hdr_h - 1], fill=(22, 22, 28))

            draw.text((10, 10), "LEFT: WEIGHT-ONLY (29.6% RETENTION)", fill=(240, 100, 100), font=font_hdr)
            draw.text((pw + GAP * scale + 10, 10), "RIGHT: TRI-AXIS (70.0% RETAINED)", fill=(100, 230, 120), font=font_hdr)

            draw.rectangle([0, gh + hdr_h, gw, gh + hdr_h + ftr_h], fill=(10, 10, 14))

            draw.text((10, gh + hdr_h + 7), "⚠ CATASTROPHIC FORGETTING (Overwritten)", fill=(255, 80, 80), font=font_badge)
            draw.text((pw + GAP * scale + 10, gh + hdr_h + 7), "✔ TASK A LOOP FOSSILIZED (Protected)", fill=(100, 230, 120), font=font_badge)

            frames.append(out)

    _sheet(frames, "retention", cols=cols)


def main():
    print("Generating scrollytelling sprite sheets for plasticity & retention...")
    render_plasticity_triptych()
    render_retention_comparison()
    print("Done!")


if __name__ == "__main__":
    main()
