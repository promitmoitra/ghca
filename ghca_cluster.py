"""
Multi-core embedding experiment for Greenberg-Hastings CA.

Addresses: at what seeding density (clustering level) does activity persist
when cores are embedded in a large resting lattice?

Three regimes / transitions are distinguished:

  1. Pre-selected persistent cores — individual cores are spiral-wave sources
     that survive indefinitely in absorbing boundaries at any non-zero density,
     UNLESS the persistent config set is resting-cell-free (see regime 3).

  2. Uniformly random cores (lower nucleation threshold) — individually
     non-persistent cores can *collectively nucleate* stable activity once the
     seeding density exceeds a critical threshold.  This threshold scales with
     act/tau0: lower active fraction → higher required density.

  3. Upper extinction threshold (resting-cell starvation) — observed for
     parameter pairs whose entire persistent config set is resting-cell-free.
     Five such pairs exist in the tested 1..17 × 1..17 grid:

       Pair    tau0  n_configs  n_multisets  structure
       (2,10)   12       8           1       pure phase-wave
       (3,14)   17       8           1       pure phase-wave
       (3,13)   16      40           5       resting-cell-free (mixed)
       (4,16)   20     120          15       resting-cell-free (mixed)
       (4,17)   21      40           5       resting-cell-free (mixed)

     Pure phase-wave pairs ((2,10) and (3,14)) have all configs drawn from the
     single multiset {act, 2·act, 3·act, 4·act} — uniform inter-state spacing
     of `act`, no resting cells anywhere.  At full tiling density (d=1) the
     lattice contains zero resting cells so the GHCA infection mechanism
     (resting → active) can never fire; activity collapses in a single step.
     Mixed resting-cell-free pairs share the same extinction mechanism but
     spread probability mass over several distinct multisets.

     The transition sharpens as density approaches 1: P(persist) drops from
     ~1 at d≈0.90 to 0 at d=1.00, tracking the vanishing supply of resting
     cells from the few remaining empty slots.

Parameter space notes (act in 1..17, pas in 1..17):
  - Pairs with act/tau0 < 1/4 tend to have n_configs=0 (no persistent cores
    exist for isolated 2×2 patches — activity requires cooperation).
  - For pairs near the "extinction boundary" (n_configs just above zero) the
    config set often collapses to a single multiset (n_multisets=1).
  - All five resting-cell-free pairs sit at or near that boundary.

Usage:
    python ghca_cluster.py                     # default sweep
    python ghca_cluster.py --L 80 --T 2000 --trials 40
    python ghca_cluster.py --plot              # requires matplotlib
"""

import numpy as np
import os
import argparse


# ---------------------------------------------------------------------------
# Vectorised CA step (Von Neumann / r=1 neighbourhood — matches ghca_main.py)
# ---------------------------------------------------------------------------

def step_vec(p, act, tau0, periodic=False):
    """Single GHCA step on a 2-D array, matching the original Python-loop rules."""
    p = p.astype(np.int16)
    L = p.shape[0]
    # 4-connected (cardinal) neighbours only — r=1 Von Neumann
    if periodic:
        nbrs = np.stack(
            [np.roll(np.roll(p, dr, 0), dc, 1)
             for dr in (-1, 0, 1) for dc in (-1, 0, 1)
             if abs(dr) + abs(dc) == 1],
            axis=0)
    else:
        pad = np.pad(p, 1, mode='constant', constant_values=0)
        nbrs = np.stack(
            [pad[1+dr:1+dr+L, 1+dc:1+dc+L]
             for dr in (-1, 0, 1) for dc in (-1, 0, 1)
             if abs(dr) + abs(dc) == 1],
            axis=0)
    has_active_nbr = np.any((nbrs >= 1) & (nbrs <= act), axis=0)
    q = np.where(p == tau0, 0, p)
    q = np.where((p > 0) & (p < tau0), p + 1, q)
    q = np.where((p == 0) & has_active_nbr, 1, q)
    return q.astype(np.int8)


def run_vec(lattice, act, pas, T, periodic=False):
    """Run up to T steps; return (survived, lifetime)."""
    tau0 = act + pas
    p = lattice.copy()
    for t in range(T):
        if not np.any((p > 0) & (p <= act)):
            if not np.any(p > act):
                return False, t
        p = step_vec(p, act, tau0, periodic)
    return bool(np.any((p > 0) & (p <= act))), T


# ---------------------------------------------------------------------------
# Lattice builders
# ---------------------------------------------------------------------------

def _decode(config_id, n_states, core_size, cacore):
    s = cacore.base_conv(int(config_id), n_states)
    s = s.rjust(core_size, '0')
    return cacore.str_to_state(s, n_states)


def make_lattice_persistent(L, core_len, density, pair, persistent_ids, rng, cacore):
    """Seed from the pre-computed persistent config set."""
    act, pas = pair
    n_states = act + pas + 1
    core_size = core_len ** 2
    n_grid = L // core_len
    n_slots = n_grid ** 2
    n_cores = max(1, int(round(density * n_slots)))

    all_slots = [(r * core_len, c * core_len)
                 for r in range(n_grid) for c in range(n_grid)]
    lat = np.zeros((L, L), dtype=np.int8)
    for i in rng.choice(n_slots, min(n_cores, n_slots), replace=False):
        r, c = all_slots[i]
        cid = int(rng.choice(persistent_ids))
        state = _decode(cid, n_states, core_size, cacore)
        lat[r:r+core_len, c:c+core_len] = state
    return lat


def make_lattice_random(L, core_len, density, pair, rng):
    """Seed with uniformly random configs (not filtered to persistent set)."""
    act, pas = pair
    n_states = act + pas + 1
    n_grid = L // core_len
    n_slots = n_grid ** 2
    n_cores = max(1, int(round(density * n_slots)))

    all_slots = [(r * core_len, c * core_len)
                 for r in range(n_grid) for c in range(n_grid)]
    lat = np.zeros((L, L), dtype=np.int8)
    for i in rng.choice(n_slots, min(n_cores, n_slots), replace=False):
        r, c = all_slots[i]
        lat[r:r+core_len, c:c+core_len] = rng.integers(
            0, n_states, size=(core_len, core_len)).astype(np.int8)
    return lat


# ---------------------------------------------------------------------------
# Sweep utilities
# ---------------------------------------------------------------------------

def sweep_density(pair, L, densities, mode, core_len, n_trials, T,
                  data_path, rng, persistent_ids=None, cacore=None):
    """
    Estimate P(persist) and mean lifetime across seeding densities.

    mode : 'random'     — uniform random cores
           'persistent' — draw from the pre-computed persistent set
    """
    act, pas = pair

    persist_probs, mean_lifetimes = [], []
    for density in densities:
        survived, lifetimes = 0, []
        for _ in range(n_trials):
            if mode == 'persistent':
                lat = make_lattice_persistent(
                    L, core_len, density, pair, persistent_ids, rng, cacore)
            else:
                lat = make_lattice_random(L, core_len, density, pair, rng)
            sv, lt = run_vec(lat, act, pas, T)
            survived += sv
            lifetimes.append(lt)
        persist_probs.append(survived / n_trials)
        mean_lifetimes.append(float(np.mean(lifetimes)))

    return np.array(persist_probs), np.array(mean_lifetimes)


def critical_density(densities, probs, threshold=0.5):
    """Linearly interpolated density where P(persist) crosses threshold (upward)."""
    for i in range(len(probs) - 1):
        if probs[i] < threshold <= probs[i + 1]:
            d0, d1 = densities[i], densities[i + 1]
            p0, p1 = probs[i], probs[i + 1]
            return d0 + (threshold - p0) / (p1 - p0) * (d1 - d0)
    return None


def upper_critical_density(densities, probs, threshold=0.5):
    """Linearly interpolated density where P(persist) crosses threshold (downward)."""
    for i in range(len(probs) - 1):
        if probs[i] >= threshold > probs[i + 1]:
            d0, d1 = densities[i], densities[i + 1]
            p0, p1 = probs[i], probs[i + 1]
            return d0 + (threshold - p0) / (p1 - p0) * (d1 - d0)
    return None


def decode_all(ids, n_states, core_size=4):
    """Vectorised base-n decoder: convert integer config IDs to cell-state arrays.

    Returns array of shape (len(ids), core_size) with dtype int8.
    ~10,000× faster than the ghca_core pure-Python loop for large id arrays.
    """
    ids = np.array(ids, dtype=np.int64)
    out = np.zeros((len(ids), core_size), dtype=np.int8)
    tmp = ids.copy()
    for pos in range(core_size - 1, -1, -1):
        out[:, pos] = tmp % n_states
        tmp //= n_states
    return out


def config_resting_cell_fraction(pair, data_path, core_size=4):
    """Fraction of persistent configs for `pair` that have at least one resting cell.

    Zero means resting-cell-free — necessary (and nearly sufficient) for the
    upper extinction transition to appear.
    """
    act, pas = pair
    n_states = act + pas + 1
    ids = np.load(
        data_path + 'act_config_ids_states-({0:02d},{1:02d})_core-04.npy'.format(act, pas))
    if len(ids) == 0:
        return float('nan')
    configs = decode_all(ids, n_states, core_size)
    has_resting = np.any(configs == 0, axis=1)
    return float(has_resting.mean())


def analyze_config_structure(pair, data_path, core_size=4):
    """Classify the persistent config set for a (act, pas) pair.

    Returns a dict with keys:
      n_configs      — number of persistent configs
      resting_frac   — fraction that contain at least one state-0 cell
      n_multisets    — number of distinct sorted-state multisets
      is_phase_wave  — True if all configs equal {act, 2*act, 3*act, 4*act}
      resting_free   — True if resting_frac == 0 and n_configs > 0
    """
    act, pas = pair
    n_states = act + pas + 1
    tau0 = act + pas
    ids = np.load(
        data_path + 'act_config_ids_states-({0:02d},{1:02d})_core-04.npy'.format(act, pas))
    n = len(ids)
    if n == 0:
        return dict(n_configs=0, resting_frac=float('nan'),
                    n_multisets=0, is_phase_wave=False, resting_free=False)

    configs = decode_all(ids, n_states, core_size)
    has_resting = np.any(configs == 0, axis=1)
    resting_frac = float(has_resting.mean())

    configs_sorted = np.sort(configs, axis=1)
    n_multisets = len(np.unique(configs_sorted, axis=0))

    pw = np.array([act, 2 * act, 3 * act, 4 * act], dtype=np.int8)
    valid_pw = bool((pw >= 1).all() and (pw <= tau0).all())
    if valid_pw:
        is_phase_wave = bool(np.all(configs_sorted == pw))
    else:
        is_phase_wave = False

    return dict(
        n_configs=n,
        resting_frac=resting_frac,
        n_multisets=n_multisets,
        is_phase_wave=is_phase_wave,
        resting_free=(resting_frac == 0.0),
    )


def nn_distance_cells(density, L, core_len):
    """Estimated mean nearest-neighbour distance in cells at a given density."""
    n_grid = L // core_len
    n_cores = density * n_grid ** 2
    if n_cores <= 0:
        return float('inf')
    return core_len * n_grid / (2 * np.sqrt(n_cores))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--L', type=int, default=80)
    parser.add_argument('--T', type=int, default=2000)
    parser.add_argument('--trials', type=int, default=40)
    parser.add_argument('--n_densities', type=int, default=15)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--data_path', default='./result/')
    parser.add_argument('--out', default='./result/cluster_sweep.npz')
    parser.add_argument('--plot', action='store_true',
                        help='Generate plots (requires matplotlib)')
    args = parser.parse_args()

    import ghca_core as cacore  # imported here so ghca_cluster is importable without it

    rng = np.random.default_rng(args.seed)
    core_len = 2
    core_size = core_len ** 2

    # Long-refractory pairs are where the threshold is most interesting
    pairs = [
        (1, 1), (2, 2), (4, 4), (8, 8),  # symmetric — should persist easily
        (6, 17), (4, 12), (3, 14), (2, 10), (1, 7),  # long refractory
        (17, 6), (12, 4), (7, 1),  # short refractory — should persist easily
    ]
    pairs = [p for p in pairs if os.path.exists(
        args.data_path + 'act_config_ids_states-({0:02d},{1:02d})_core-04.npy'.format(*p))]

    # Coarse sweep to find the interesting density range, then fine around it
    coarse = np.array([0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.20, 0.35, 0.50, 1.0])

    all_results = {}
    print(f"{'pair':12s} {'tau0':>5s} {'d_crit':>8s} {'nn_crit':>9s}")

    for pair in pairs:
        act, pas = pair
        fname = args.data_path + 'act_config_ids_states-({0:02d},{1:02d})_core-04.npy'.format(
            act, pas)
        persistent_ids = np.load(fname)

        probs, lifetimes = sweep_density(
            pair, args.L, coarse, mode='random',
            core_len=core_len, n_trials=args.trials, T=args.T,
            data_path=args.data_path, rng=rng,
            persistent_ids=persistent_ids, cacore=cacore)

        dc = critical_density(coarse, probs)
        nn = nn_distance_cells(dc, args.L, core_len) if dc else float('nan')
        tau0 = act + pas

        all_results[str(pair)] = dict(
            probs=probs, lifetimes=lifetimes, dc=dc, nn=nn)

        dc_str = f'{dc:.4f}' if dc else '<0.005'
        nn_str = f'{nn:.1f}' if dc else '>20'
        print(f'{str(pair):12s} {tau0:5d} {dc_str:>8s} {nn_str:>9s} cells')

    np.savez(args.out,
             densities=coarse,
             pairs=np.array(pairs),
             **{f'probs_{p[0]}_{p[1]}': all_results[str(p)]['probs'] for p in pairs},
             **{f'lifetimes_{p[0]}_{p[1]}': all_results[str(p)]['lifetimes'] for p in pairs},
             **{f'dc_{p[0]}_{p[1]}': np.array([all_results[str(p)]['dc'] or np.nan])
                for p in pairs})
    print(f'\nResults saved to {args.out}')

    if args.plot:
        _make_plots(pairs, coarse, all_results, args)


def _make_plots(pairs, densities, all_results, args):
    import matplotlib.pyplot as plt

    # Sort pairs by tau0 for colour ordering
    pairs_sorted = sorted(pairs, key=lambda p: p[0]+p[1])
    cmap = plt.cm.viridis(np.linspace(0, 1, len(pairs_sorted)))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for pair, col in zip(pairs_sorted, cmap):
        r = all_results[str(pair)]
        label = f"({pair[0]},{pair[1]}) τ={pair[0]+pair[1]}"
        axes[0].plot(densities, r['probs'], '-o', color=col, label=label)
        axes[1].plot(densities, r['lifetimes'], '-o', color=col, label=label)
        if r['dc']:
            axes[0].axvline(r['dc'], color=col, linestyle=':', alpha=0.5, linewidth=0.8)

    axes[0].axhline(0.5, color='k', linestyle='--', linewidth=0.8, alpha=0.6)
    axes[0].set_xlabel('Seeding density (fraction of core slots filled)')
    axes[0].set_ylabel('P(activity persists at T)')
    axes[0].set_title(f'Critical seeding density  (L={args.L}, T={args.T})')
    axes[0].legend(fontsize=7, ncol=2)
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel('Seeding density')
    axes[1].set_ylabel('Mean activity lifetime (steps)')
    axes[1].set_title('Mean lifetime vs seeding density')
    axes[1].legend(fontsize=7, ncol=2)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = args.out.replace('.npz', '.png')
    plt.savefig(out_png, dpi=150)
    print(f'Plot saved to {out_png}')
    plt.show()


if __name__ == '__main__':
    main()
