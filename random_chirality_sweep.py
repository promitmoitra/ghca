"""
Random-seeding chirality sweep — regime 2 collective nucleation.

Tests whether imposing a chirality bias on random initial configs shifts the
lower nucleation threshold d_crit for collective wave formation.

Four seeding modes compared:
  uniform   — draw uniformly from all n_states^4 configs (baseline)
  no_rest   — reject configs with any resting (state-0) cell
  cw        — rejection-sample for CW configs (W=+1)
  ccw       — rejection-sample for CCW configs (W=-1)

The uniform/no_rest comparison isolates the effect of resting-cell removal.
The cw/ccw comparison isolates pure handedness at equal density.

Target pair: (4,12)  tau0=16  n_states=17
  - d_crit_random ≈ 0.02 (low enough that the transition is in [0,0.5])
  - resting-cell-free pairs have n_configs=0 for (4,16) and (4,17) so we
    avoid those; (4,12) has many persistent configs and a non-trivial threshold.

Usage:
    python random_chirality_sweep.py
    python random_chirality_sweep.py --pair 4 12 --L 80 --trials 60
    python random_chirality_sweep.py --pair 6 17 --L 80 --trials 60
"""

import numpy as np
import argparse

from ghca_cluster import (step_vec, run_vec, make_lattice_random,
                           make_lattice_random_chirality,
                           critical_density, decode_all, chirality_batch)


def sweep_random_chirality(pair, L, densities, n_trials, T, rng):
    """Run density sweep for each random-seeding mode; return dict of results."""
    act, pas = pair
    tau0 = act + pas

    modes = {
        'uniform': dict(chirality=None, no_resting=False),
        'no_rest': dict(chirality=None, no_resting=True),
        'cw':      dict(chirality=+1,   no_resting=False),
        'ccw':     dict(chirality=-1,   no_resting=False),
    }

    results = {}
    for label, opts in modes.items():
        probs, lifetimes = [], []
        for density in densities:
            survived, lts = 0, []
            for _ in range(n_trials):
                if opts['chirality'] is None and not opts['no_resting']:
                    lat = make_lattice_random(L, 2, density, pair, rng)
                else:
                    lat = make_lattice_random_chirality(
                        L, 2, density, pair, rng,
                        chirality=opts['chirality'],
                        no_resting=opts['no_resting'])
                sv, lt = run_vec(lat, act, pas, T)
                survived += sv
                lts.append(lt)
            probs.append(survived / n_trials)
            lifetimes.append(float(np.mean(lts)))

        probs = np.array(probs)
        lifetimes = np.array(lifetimes)
        dc = critical_density(densities, probs)
        results[label] = dict(probs=probs, lifetimes=lifetimes, dc=dc)
        print("  {:8s}  d_crit={}".format(
            label,
            '{:.4f}'.format(dc) if dc is not None else 'N/A'))

    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pair', type=int, nargs=2, default=[4, 12],
                        metavar=('ACT', 'PAS'))
    parser.add_argument('--L', type=int, default=80)
    parser.add_argument('--T', type=int, default=2000)
    parser.add_argument('--trials', type=int, default=60)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--out', default='./result/random_chirality_sweep.npz')
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()

    pair = tuple(args.pair)
    act, pas = pair
    rng = np.random.default_rng(args.seed)

    # Focused on low-density regime where collective nucleation happens
    densities = np.array([
        0.005, 0.008, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05,
        0.07, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50])

    print("Random-seeding chirality sweep for pair ({},{}) tau0={}".format(
        act, pas, act + pas))
    print("L={} T={} trials={}\n".format(args.L, args.T, args.trials))

    # Report acceptance rates for CW/CCW rejection sampling
    _report_accept_rates(pair, rng_seed=args.seed + 1)

    results = sweep_random_chirality(pair, args.L, densities, args.trials,
                                     args.T, rng)

    np.savez(args.out,
             densities=densities,
             pair=np.array(pair),
             **{'{}_probs'.format(k): v['probs'] for k, v in results.items()},
             **{'{}_lifetimes'.format(k): v['lifetimes'] for k, v in results.items()})
    print('\nResults saved to {}'.format(args.out))

    print("\n=== Summary ===")
    print("{:>8s}  {:>10s}".format("mode", "d_crit"))
    for label, r in results.items():
        dc = '{:.4f}'.format(r['dc']) if r['dc'] is not None else 'N/A'
        print("{:>8s}  {:>10s}".format(label, dc))

    _print_significance(results, densities)

    if args.plot:
        _plot(pair, densities, results, args)


def _report_accept_rates(pair, n_sample=50000, rng_seed=99):
    """Profile the chirality distribution of random configs for this pair."""
    act, pas = pair
    tau0 = act + pas
    n_states = tau0 + 1
    rng2 = np.random.default_rng(rng_seed)
    sample = rng2.integers(0, n_states, size=(n_sample, 4)).astype(np.int8)
    chi = chirality_batch(sample, tau0)
    no_zero = ~np.any(sample == 0, axis=1)
    n_cw  = int(np.sum(chi == +1))
    n_ccw = int(np.sum(chi == -1))
    n_noz = int(np.sum(no_zero))
    print("  Random config profile (n={:,}):".format(n_sample))
    print("    CW={:.3f}  CCW={:.3f}  no-resting={:.3f}".format(
        n_cw/n_sample, n_ccw/n_sample, n_noz/n_sample))
    print("    Expected draws/slot: CW={:.1f}  CCW={:.1f}  no-rest={:.1f}\n".format(
        n_sample/max(n_cw, 1), n_sample/max(n_ccw, 1), n_sample/max(n_noz, 1)))


def _print_significance(results, densities):
    """Compare CW vs CCW and uniform vs no_rest at each density."""
    print("\n=== CW vs CCW comparison ===")
    print("{:>6s}  {:>6s}  {:>6s}  {:>8s}".format(
        "d", "P_cw", "P_ccw", "delta"))
    cw_p  = results.get('cw',  {}).get('probs', None)
    ccw_p = results.get('ccw', {}).get('probs', None)
    if cw_p is not None and ccw_p is not None:
        for d, p1, p2 in zip(densities, cw_p, ccw_p):
            print("{:>6.3f}  {:>6.3f}  {:>6.3f}  {:>+8.3f}".format(
                d, p1, p2, p1 - p2))

    print("\n=== uniform vs no_rest comparison ===")
    print("{:>6s}  {:>6s}  {:>6s}  {:>8s}".format(
        "d", "P_uni", "P_nor", "delta"))
    uni_p = results.get('uniform', {}).get('probs', None)
    nor_p = results.get('no_rest', {}).get('probs', None)
    if uni_p is not None and nor_p is not None:
        for d, p1, p2 in zip(densities, uni_p, nor_p):
            print("{:>6.3f}  {:>6.3f}  {:>6.3f}  {:>+8.3f}".format(
                d, p1, p2, p1 - p2))


def _plot(pair, densities, results, args):
    import matplotlib.pyplot as plt

    colors = {'uniform': 'black', 'no_rest': 'gray',
              'cw': 'tab:blue', 'ccw': 'tab:red'}
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for label, r in results.items():
        col = colors.get(label, 'green')
        axes[0].plot(densities, r['probs'], '-o', color=col,
                     label=label, markersize=4)
        axes[1].plot(densities, r['lifetimes'], '-o', color=col,
                     label=label, markersize=4)
        if r['dc']:
            axes[0].axvline(r['dc'], color=col, linestyle=':', alpha=0.6, lw=1)

    axes[0].axhline(0.5, color='k', linestyle='--', lw=0.8, alpha=0.5)
    axes[0].set_xlabel('Seeding density (random)')
    axes[0].set_ylabel('P(activity persists)')
    axes[0].set_title('Random-chirality seeding: ({},{}) L={}'.format(
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
