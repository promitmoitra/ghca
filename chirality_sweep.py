"""
Chirality-controlled density sweep for resting-cell-free GHCA pairs.

Tests whether chirality alignment between adjacent cores affects the lower
nucleation threshold (d_crit) or the upper extinction transition (d → 1).

Three seeding modes compared:
  mixed  — draw uniformly from all persistent configs (CW + CCW)
  cw     — draw only CW (chirality=+1) configs
  ccw    — draw only CCW (chirality=-1) configs

Usage:
    python chirality_sweep.py
    python chirality_sweep.py --pair 3 13 --L 80 --trials 60
"""

import numpy as np
import argparse
import os

from ghca_cluster import (step_vec, run_vec, make_lattice_persistent,
                           sweep_density, critical_density,
                           upper_critical_density, decode_all,
                           chirality_batch)


def sweep_chirality(pair, L, densities, n_trials, T, data_path, rng,
                    chirality_modes=(None, +1, -1)):
    """Run density sweep for each chirality mode; return dict of results."""
    import ghca_core as cacore

    act, pas = pair
    fname = data_path + 'act_config_ids_states-({0:02d},{1:02d})_core-04.npy'.format(
        act, pas)
    persistent_ids = np.load(fname)
    core_len = 2

    # Verify the pair is resting-cell-free
    configs = decode_all(persistent_ids, act + pas + 1)
    if np.any(configs == 0):
        raise ValueError("Pair {} has resting cells — upper extinction not expected".format(pair))

    results = {}
    mode_labels = {None: 'mixed', +1: 'cw', -1: 'ccw', 0: 'chirality0'}

    for chi_mode in chirality_modes:
        label = mode_labels[chi_mode]
        probs, lifetimes = [], []
        for density in densities:
            survived, lts = 0, []
            for _ in range(n_trials):
                lat = make_lattice_persistent(
                    L, core_len, density, pair, persistent_ids, rng, cacore,
                    chirality=chi_mode)
                sv, lt = run_vec(lat, act, pas, T)
                survived += sv
                lts.append(lt)
            probs.append(survived / n_trials)
            lifetimes.append(float(np.mean(lts)))

        probs = np.array(probs)
        lifetimes = np.array(lifetimes)
        dc_lo = critical_density(densities, probs)
        dc_hi = upper_critical_density(densities, probs)
        results[label] = dict(probs=probs, lifetimes=lifetimes,
                               dc_lo=dc_lo, dc_hi=dc_hi)
        print("  {:7s}  d_crit_lo={:s}  d_crit_hi={:s}".format(
            label,
            '{:.4f}'.format(dc_lo) if dc_lo is not None else 'N/A',
            '{:.4f}'.format(dc_hi) if dc_hi is not None else 'N/A'))

    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pair', type=int, nargs=2, default=[3, 13],
                        metavar=('ACT', 'PAS'))
    parser.add_argument('--L', type=int, default=80)
    parser.add_argument('--T', type=int, default=2000)
    parser.add_argument('--trials', type=int, default=60)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--data_path', default='./result/')
    parser.add_argument('--out', default='./result/chirality_sweep.npz')
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()

    pair = tuple(args.pair)
    act, pas = pair
    rng = np.random.default_rng(args.seed)

    # Dense grid around d=0..1, finer near the boundaries
    densities = np.array([
        0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20,
        0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.93,
        0.95, 0.97, 0.98, 0.99, 1.00])

    print("Chirality sweep for pair ({},{}) tau0={}".format(act, pas, act+pas))
    print("L={} T={} trials={}\n".format(args.L, args.T, args.trials))

    results = sweep_chirality(
        pair, args.L, densities, args.trials, args.T,
        args.data_path, rng, chirality_modes=(None, +1, -1))

    np.savez(args.out, densities=densities, pair=np.array(pair),
             **{'{}_probs'.format(k): v['probs'] for k, v in results.items()},
             **{'{}_lifetimes'.format(k): v['lifetimes'] for k, v in results.items()})
    print('\nResults saved to {}'.format(args.out))

    print("\n=== Summary ===")
    print("{:>8s}  {:>10s}  {:>10s}".format("mode", "d_crit_lo", "d_crit_hi"))
    for label, r in results.items():
        lo = '{:.4f}'.format(r['dc_lo']) if r['dc_lo'] is not None else 'N/A'
        hi = '{:.4f}'.format(r['dc_hi']) if r['dc_hi'] is not None else 'N/A'
        print("{:>8s}  {:>10s}  {:>10s}".format(label, lo, hi))

    if args.plot:
        _plot(pair, densities, results, args)


def _plot(pair, densities, results, args):
    import matplotlib.pyplot as plt

    colors = {'mixed': 'black', 'cw': 'tab:blue', 'ccw': 'tab:red'}
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for label, r in results.items():
        col = colors.get(label, 'gray')
        axes[0].plot(densities, r['probs'], '-o', color=col, label=label, markersize=4)
        axes[1].plot(densities, r['lifetimes'], '-o', color=col, label=label, markersize=4)
        if r['dc_lo']:
            axes[0].axvline(r['dc_lo'], color=col, linestyle=':', alpha=0.6, lw=1)
        if r['dc_hi']:
            axes[0].axvline(r['dc_hi'], color=col, linestyle='--', alpha=0.6, lw=1)

    axes[0].axhline(0.5, color='k', linestyle='--', lw=0.8, alpha=0.5)
    axes[0].set_xlabel('Seeding density')
    axes[0].set_ylabel('P(activity persists)')
    axes[0].set_title('Chirality-controlled seeding: ({},{})  L={}'.format(
        pair[0], pair[1], args.L))
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel('Seeding density')
    axes[1].set_ylabel('Mean lifetime (steps)')
    axes[1].set_title('Mean lifetime vs density')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = args.out.replace('.npz', '.png')
    plt.savefig(out_png, dpi=150)
    print('Plot saved to {}'.format(out_png))
    plt.show()


if __name__ == '__main__':
    main()
