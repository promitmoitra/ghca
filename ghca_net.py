"""
Greenberg-Hastings excitable dynamics on a graph.

This generalises the lattice `Population` in ``ghca_main.py`` to an arbitrary
directed weighted graph, with per-node timescales, weighted-threshold
excitation, spontaneous firing and a homeostatic threshold. It implements the
substrate specified in ``docs/learning_experiments.md`` (sections 1 and 7,
steps 1-4). No learning is implemented here; this is the medium the E0
characterisation runs on and that later experiments make plastic.

State encoding (matches ``ghca_main.Population``):
    phi = 0            -> rest / quiescent (excitable)
    phi in [1, a]      -> active / excited (excites neighbours)
    phi in [a+1, tau_i]-> refractory (counts down, cannot be excited)
advancing past tau_i wraps back to 0.
"""

import numpy as np


# ----------------------------------------------------------------------------
# Topology builders. Each returns a dense (N, N) weight matrix W where
# W[i, j] > 0 means "j is an input neighbour of i" (j can excite i).
# ----------------------------------------------------------------------------

def ring(N, k=2, weight=1.0):
    """1D ring; each node connects to its k nearest neighbours on each side."""
    W = np.zeros((N, N))
    for i in range(N):
        for d in range(1, k + 1):
            W[i, (i - d) % N] = weight
            W[i, (i + d) % N] = weight
    return W


def lattice2d(L, r=1, periodic=True, weight=1.0):
    """2D lattice (torus if periodic) with Euclidean-radius-r neighbourhood.

    Uses the same neighbourhood convention as ``ghca_main`` (offsets with
    ``sqrt(dx^2+dy^2) <= r``), so r=1 gives the 4 von Neumann neighbours.
    Node index is row-major: node = x * L + y.
    """
    offs = [(dx, dy) for dx in range(-r, r + 1) for dy in range(-r, r + 1)
            if (dx or dy) and np.hypot(dx, dy) <= r]
    N = L * L
    W = np.zeros((N, N))
    for x in range(L):
        for y in range(L):
            i = x * L + y
            for dx, dy in offs:
                nx, ny = x + dx, y + dy
                if periodic:
                    nx, ny = nx % L, ny % L
                elif not (0 <= nx < L and 0 <= ny < L):
                    continue
                W[i, nx * L + ny] = weight
    return W


def smallworld(N, k=6, beta=0.1, weight=1.0, seed=0):
    """Directed Watts-Strogatz-style ring lattice with rewired shortcuts."""
    rng = np.random.default_rng(seed)
    W = ring(N, k=k // 2, weight=weight)
    for i in range(N):
        for j in np.where(W[i] > 0)[0]:
            if rng.random() < beta:
                W[i, j] = 0.0
                choices = np.setdiff1d(np.arange(N), np.append(np.where(W[i] > 0)[0], i))
                if len(choices):
                    W[i, rng.choice(choices)] = weight
    return W


def rgg(N, radius=0.12, weight=1.0, seed=0):
    """Random geometric graph on the unit square (undirected, symmetric W)."""
    rng = np.random.default_rng(seed)
    pos = rng.random((N, 2))
    d = np.sqrt(((pos[:, None, :] - pos[None, :, :]) ** 2).sum(-1))
    W = ((d <= radius) & (d > 0)).astype(float) * weight
    return W


# ----------------------------------------------------------------------------
# The excitable network.
# ----------------------------------------------------------------------------

class Network:
    """Greenberg-Hastings excitable automaton on a weighted graph.

    Parameters
    ----------
    W : (N, N) ndarray
        Input weight matrix; W[i, j] is the weight of edge j -> i.
    act : int
        Excited (active) duration `a`, shared across nodes.
    pas : int
        Initial refractory duration; sets tau_i = act + pas per node.
    theta : float
        Initial excitation threshold (weighted active-neighbour sum required).
    p_s : float
        Per-step spontaneous firing probability for a rested node.
    rho_star : float or None
        Homeostatic target active fraction. If None, threshold is static.
    eta_theta : float
        Homeostatic threshold adaptation rate.
    theta_bounds : (float, float)
        Clip range for the (global) homeostatic threshold.
    seed : int
        RNG seed for spontaneous firing and initialisation.
    """

    def __init__(self, W, act=2, pas=8, theta=2.0, p_s=5e-3,
                 rho_star=None, eta_theta=1e-3, theta_bounds=(0.5, 20.0),
                 seed=0, p_s_mask=None):
        self.W = np.asarray(W, dtype=float)
        self.N = self.W.shape[0]
        self.act = int(act)
        self.tau = np.full(self.N, int(act) + int(pas), dtype=np.int64)  # per-node
        self.theta = np.full(self.N, float(theta)) if np.isscalar(theta) \
            else np.asarray(theta, dtype=float).copy()               # per-node
        self.p_s = float(p_s)
        # optional per-node spontaneous-firing mask: if given, only masked nodes
        # fire spontaneously (lets an experiment confine exploration to a medium
        # while leaving driven / context nodes deterministic). None = all nodes.
        self.p_s_mask = None if p_s_mask is None else np.asarray(p_s_mask, bool)
        self.rho_star = rho_star
        self.eta_theta = float(eta_theta)
        self.theta_bounds = theta_bounds
        self.rng = np.random.default_rng(seed)
        self.phi = np.zeros(self.N, dtype=np.int64)
        self.t = 0

    # -- initialisation helpers ---------------------------------------------
    def seed_random(self, frac_active=0.05, frac_refractory=0.10):
        """Random initial condition: some active, some refractory, rest rested."""
        self.phi[:] = 0
        idx = self.rng.permutation(self.N)
        na = int(frac_active * self.N)
        nr = int(frac_refractory * self.N)
        self.phi[idx[:na]] = self.rng.integers(1, self.act + 1, size=na)
        self.phi[idx[na:na + nr]] = self.rng.integers(self.act + 1,
                                                      self.tau[idx[na:na + nr]] + 1)
        return self

    def seed_point(self, node=0):
        """Single active seed (for wave/persistence probing)."""
        self.phi[:] = 0
        self.phi[node] = 1
        return self

    # -- dynamics ------------------------------------------------------------
    def active_mask(self):
        return (self.phi >= 1) & (self.phi <= self.act)

    def step(self, drive=None):
        """Advance one time step.

        `drive` : optional bool/0-1 array of nodes to force-excite this step
        (used by later experiments to inject cues). Returns the active mask
        *after* the update.
        """
        active = self.active_mask()
        rested = (self.phi == 0)

        # weighted active-neighbour input to each node
        inp = self.W @ active.astype(float)
        excite = rested & (inp >= self.theta)

        # spontaneous firing among rested nodes (optionally confined by a mask)
        if self.p_s > 0:
            spont = rested & (self.rng.random(self.N) < self.p_s)
            if self.p_s_mask is not None:
                spont &= self.p_s_mask
            excite |= spont

        # forced external drive (cues), if any
        if drive is not None:
            excite |= (rested & np.asarray(drive, dtype=bool))

        # advance the clock: 1..tau-1 -> +1 ; at/after tau -> wrap to 0
        moving = (self.phi >= 1)
        self.phi[moving] += 1
        self.phi[self.phi > self.tau] = 0  # wrap (uses per-node tau)

        # apply new excitations
        self.phi[excite] = 1

        # homeostatic threshold update
        if self.rho_star is not None:
            rho = self.active_mask().mean()
            self.theta = np.clip(
                self.theta + self.eta_theta * (rho - self.rho_star),
                *self.theta_bounds)

        self.t += 1
        return self.active_mask()

    def run(self, T, drive_fn=None, record=True):
        """Run T steps. Returns dict of recorded observables.

        `drive_fn(t)` -> optional per-step drive array.
        """
        A = np.empty(T)          # active fraction
        R = np.empty(T)          # phase coherence (Kuramoto order parameter)
        phis = np.empty((T, self.N), dtype=np.int64) if record else None
        for k in range(T):
            drive = drive_fn(k) if drive_fn is not None else None
            self.step(drive)
            A[k] = self.active_mask().mean()
            R[k] = self.coherence()
            if record:
                phis[k] = self.phi
        return {"A": A, "R": R, "phi": phis}

    # -- observables (order parameters / regime probes) ---------------------
    def coherence(self):
        """Kuramoto order parameter over cyclic phases mapped to angles."""
        ang = 2 * np.pi * self.phi / self.tau
        return np.abs(np.exp(1j * ang).mean())


def dominant_period(A, min_period=2):
    """Estimate the dominant oscillation period of a scalar time series via FFT.

    Returns np.inf if there is no clear spectral peak away from DC.
    """
    A = np.asarray(A, dtype=float)
    A = A - A.mean()
    if np.allclose(A, 0):
        return np.inf
    spec = np.abs(np.fft.rfft(A)) ** 2
    freqs = np.fft.rfftfreq(len(A))
    valid = freqs > (1.0 / (len(A)))          # ignore DC / ultra-low
    valid &= (freqs > 0) & (1.0 / np.where(freqs > 0, freqs, np.inf) >= min_period)
    if not valid.any() or spec[valid].max() <= 0:
        return np.inf
    f_peak = freqs[valid][np.argmax(spec[valid])]
    return 1.0 / f_peak if f_peak > 0 else np.inf


def count_phase_singularities(phi_grid, act, tau):
    """Count spiral cores (phase singularities) on a 2D periodic lattice.

    Sums |topological charge| around each 2x2 plaquette using the cyclic phase.
    Defined for lattice2d only; returns an integer count of defect cores.
    """
    L = phi_grid.shape[0]
    ang = 2 * np.pi * phi_grid / tau

    def wrap(d):
        return (d + np.pi) % (2 * np.pi) - np.pi

    count = 0
    for x in range(L):
        for y in range(L):
            a = ang[x, y]
            b = ang[(x + 1) % L, y]
            c = ang[(x + 1) % L, (y + 1) % L]
            d = ang[x, (y + 1) % L]
            circ = wrap(b - a) + wrap(c - b) + wrap(d - c) + wrap(a - d)
            if abs(circ) > np.pi:            # nonzero winding
                count += 1
    return count


if __name__ == "__main__":
    # smoke test: a small torus should sustain activity at low threshold
    net = Network(lattice2d(24, r=1), act=2, pas=8, theta=1.0, p_s=1e-3, seed=1)
    net.seed_random()
    out = net.run(200, record=False)
    print("mean active fraction:", round(out["A"][100:].mean(), 3),
          "| mean coherence:", round(out["R"][100:].mean(), 3),
          "| dominant period:", round(dominant_period(out["A"][50:]), 2))
