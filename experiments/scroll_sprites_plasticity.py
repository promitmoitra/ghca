"""Sprite sheets for Beat 10 and Beat 11 of the scrollytelling walkthrough.

Generates:
    plasticity.png -- Tri-axis closed-loop self-assembly vs unlearned/weight-only (Beat 10)
    retention.png  -- Topological loop preservation vs weight-only forgetting (Beat 11)

Outputs land in docs/scroll/assets/.
"""

import os
import sys
import json
import numpy as np
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))

from ghca_plasticity import ClosedLoopPlasticityEngine, make_closed_loop_learner
from ghca_net_viz import state_colors

OUT = os.path.join(ROOT, "docs", "scroll", "assets")
os.makedirs(OUT, exist_ok=True)


def _to_img(colors_hw4, scale=6):
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


def render_plasticity_triptych(scale=6, stride=2, T=100, cols=8):
    """
    Beat 10: Closed-Loop Tri-Axis Plasticity Triptych
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
    cols_px = 3 * W + 2 * GAP
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
            c = np.zeros((H, cols_px, 4))
            
            c_ctrl = _pad_grid(ctrl.phi, ctrl.tau, ctrl.act, (H, W))
            c_w = _pad_grid(w_only.phi, w_only.tau, w_only.act, (H, W))
            c_m = _pad_grid(multi.phi, multi.tau, multi.act, (H, W))

            c[:, :W] = c_ctrl
            c[:, W:W+GAP] = white
            c[:, W+GAP:2*W+GAP] = c_w
            c[:, 2*W+GAP:2*W+2*GAP] = white
            c[:, 2*W+2*GAP:] = c_m

            frames.append(_to_img(c, scale=scale))

    _sheet(frames, "plasticity", cols=cols)


def render_retention_comparison(scale=6, stride=2, T=100, cols=8):
    """
    Beat 11: Sequential Learning & Topological Loop Protection
    Left: Weight-Only Adaptation (Task B overwrites Task A -> Forgetting)
    Right: Multi-Axis Plasticity (High-tau loop fossilization -> Retention)
    """
    seed = 7
    
    w_net = make_closed_loop_learner(axes=("w",), seed=seed)
    m_net = make_closed_loop_learner(axes=("tau", "theta", "w"), seed=seed)
    rng = np.random.default_rng(seed)

    # Phase 1: Train Task A (Identity 0->0, 1->1)
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

    # Phase 2: Train Task B (Reversal 0->1, 1->0)
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

    # Test Task A retention
    H, W = 12, 12
    GAP = 2
    cols_px = 2 * W + GAP
    white = np.array([1.0, 1.0, 1.0, 1.0])

    frames = []
    x_test = 0  # Task A cue
    for t in range(T):
        drive_w = w_net.sensory_drive(x_test) if t < 15 else None
        drive_m = m_net.sensory_drive(x_test) if t < 15 else None

        w_net.step_learn(drive_w)
        m_net.step_learn(drive_m)

        if t % stride == 0:
            c = np.zeros((H, cols_px, 4))
            
            c_w = _pad_grid(w_net.phi, w_net.tau, w_net.act, (H, W))
            c_m = _pad_grid(m_net.phi, m_net.tau, m_net.act, (H, W))

            c[:, :W] = c_w
            c[:, W:W+GAP] = white
            c[:, W+GAP:] = c_m

            frames.append(_to_img(c, scale=scale))

    _sheet(frames, "retention", cols=cols)


def main():
    print("Generating scrollytelling sprite sheets for plasticity & retention...")
    render_plasticity_triptych()
    render_retention_comparison()
    print("Done!")


if __name__ == "__main__":
    main()
