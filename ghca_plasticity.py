"""
Closed-Loop Multi-Axis Substrate Plasticity Engine for GHCA.

Provides simultaneous local adaptation across three physical substrate axes:
1. Intrinsic Timescale Plasticity (Axis Tau): Refractory/recovery duration adaptation.
2. Homeostatic & Reward Threshold Plasticity (Axis Theta): Firing sensitivity adaptation.
3. Structural Hebbian Routing Plasticity (Axis W): Conduction weight adaptation.

Enables substrate-level credit assignment, allowing excitable media to self-organize
and retain task representations without relying on external task heads or trained readouts.
"""

import numpy as np
from ghca_net import Network
from ghca_learn import GHLearner, layered_graph


class ClosedLoopPlasticityEngine(GHLearner):
    """
    Substrate-level learner incorporating multi-axis local plasticity (tau, theta, W).

    Parameters
    ----------
    W : (N, N) ndarray
        Weight matrix.
    plastic : (N, N) ndarray
        Boolean mask of plastic edges.
    roles : dict
        Node role definitions (sensory, hidden, motor, etc.).
    axes : tuple or list or set
        Subset of ('tau', 'theta', 'w') indicating active plasticity axes.
    eta_w : float
        Learning rate for weight adaptation (Axis W).
    eta_tau : float
        Learning rate for timescale adaptation (Axis Tau).
    eta_theta : float
        Learning rate for threshold adaptation (Axis Theta).
    gamma_theta : float
        Running average factor for activity trace tracking.
    target_activity : float
        Target homeostatic firing rate (rho_star / p_star).
    theta_min, theta_max : float
        Min/max bounds for node thresholds.
    w_max : float
        Max weight bound for plastic edges.
    tau_min, tau_max : float
        Min/max bounds for per-node refractory/recovery timescales.
    seed : int
        RNG seed for deterministic initialization and perturbations.
    """

    def __init__(self, W, plastic, roles,
                 axes=("tau", "theta", "w"),
                 eta_w=0.02, eta_tau=0.05, eta_theta=0.02,
                 gamma_theta=0.05, target_activity=0.10,
                 theta_min=1.0, theta_max=10.0,
                 w_max=4.0, tau_min=2, tau_max=25,
                 lam=0.9, gamma=0.95, alpha_v=0.1,
                 seed=0, **kwargs):
        
        # Standardize axes set (case-insensitive)
        self.axes = {str(a).lower() for a in axes}

        # Initialize base GHLearner
        line_str = ""
        if "w" in self.axes:
            line_str += "A"
        if "tau" in self.axes:
            line_str += "B"

        super().__init__(W, plastic, roles, line=line_str if line_str else "AB",
                         eta_w=eta_w, eta_tau=eta_tau, lam=lam, gamma=gamma,
                         alpha_v=alpha_v, w_max=w_max, tau_min=tau_min,
                         tau_max=tau_max, seed=seed, **kwargs)

        self.eta_theta = float(eta_theta)
        self.gamma_theta = float(gamma_theta)
        self.target_activity = float(target_activity)
        self.theta_bounds = (float(theta_min), float(theta_max))

        # Per-node activity trace tracking for homeostatic threshold plasticity
        self.activity_trace = np.full(self.N, self.target_activity, dtype=float)
        self.theta_base = self.theta.copy()

    def step_learn(self, drive=None):
        """Advance one time step and update multi-axis local eligibility traces."""
        active_prev = self.active_mask()
        self._prev = self.phi.copy()
        
        # Advance network dynamics step
        self.step(drive)
        fired = (self.phi == 1) & (self._prev == 0)

        # Update per-node activity trace
        self.activity_trace = (1.0 - self.gamma_theta) * self.activity_trace + \
                              self.gamma_theta * fired.astype(float)

        # Homeostatic threshold adaptation (continuous background)
        if "theta" in self.axes:
            d_theta_homeo = self.eta_theta * (self.activity_trace - self.target_activity)
            self.theta = np.clip(self.theta + d_theta_homeo, *self.theta_bounds)

        # Axis W: edge eligibility (pre = active neighbour j, post = i just fired)
        if "w" in self.axes:
            self.E *= self.lam
            fi, aj = np.where(fired)[0], np.where(active_prev)[0]
            if fi.size and aj.size:
                self.E[np.ix_(fi, aj)] += 1.0

        # Axis Tau: period eligibility (observed inter-fire interval vs current tau)
        if "tau" in self.axes:
            fi = np.where(fired)[0]
            for i in fi:
                if self.last_fire[i] >= 0:
                    self.tau_elig[i] = (self.t - self.last_fire[i]) - self.tau[i]
                self.last_fire[i] = self.t

        return fired

    def learn(self, delta, delta_b=None):
        """Apply scalar reward RPE delta across all active plasticity axes."""
        db = delta if delta_b is None else delta_b

        # 1. Axis W (Structural Routing Plasticity)
        if "w" in self.axes:
            dW = self.eta_w * delta * self.E
            if "tau" in self.axes:
                # Metaplastic protection: longer timescale tau stabilizes synaptic weights
                tau_prot = 1.0 / (1.0 + 0.30 * (self.tau[:, None] - self.tau_min))
                dW *= tau_prot
            self.W[self.plastic] = np.clip(
                self.W[self.plastic] + dW[self.plastic],
                0.0, self.w_max
            )
            self.adj = self.W > 0

        # 2. Axis Tau (Intrinsic Timescale Plasticity)
        if "tau" in self.axes:
            if self.tau_shared:
                self.tau_scalar = float(np.clip(
                    self.tau_scalar + self.eta_tau * db * self._eps_s,
                    self.tau_min, self.tau_max
                ))
                self.tau[self.tau_mask] = self.tau_scalar
            elif self.tau_sigma > 0:
                self.tau_base = np.clip(
                    self.tau_base + self.eta_tau * db * self._eps,
                    self.tau_min, self.tau_max
                )
                self.tau = self.tau_base.copy()
            else:
                # Reward-driven timescale expansion on active pathway nodes
                d_tau = self.eta_tau * (delta * self.activity_trace + 0.05 * self.tau_elig)
                self.tau_base = np.clip(
                    self.tau_base + d_tau,
                    self.tau_min, self.tau_max
                )
                self.tau = self.tau_base.copy()

        # 3. Axis Theta (Reward-Modulated Threshold Sensitivity)
        if "theta" in self.axes:
            # Reward decreases threshold for nodes active during the trial
            d_theta_reward = -self.eta_theta * delta * self.activity_trace
            self.theta = np.clip(self.theta + d_theta_reward, *self.theta_bounds)

    def reset_traces(self):
        """Reset eligibility traces and activity memory between trials."""
        super().reset_traces()
        self.activity_trace[:] = self.target_activity


def compute_readout_independence_ratio(fixed_acc, trained_acc, eps=1e-6):
    """
    Compute Readout Independence Ratio (RIR).

    RIR = fixed_readout_accuracy / trained_readout_accuracy
    An RIR close to 1.0 (or > 0.85) indicates that substrate dynamics alone
    performing spatial/temporal credit assignment account for the task performance.
    """
    fixed_acc = float(fixed_acc)
    trained_acc = float(trained_acc)
    if trained_acc < eps:
        return 0.0
    return float(fixed_acc / (trained_acc + eps))


def make_closed_loop_learner(axes=("tau", "theta", "w"), seed=0, cfg=None):
    """
    Factory helper building a layered graph and initializing a ClosedLoopPlasticityEngine.
    """
    if cfg is None:
        cfg = dict(act=3, pas=5, theta=4.0, p_s=3e-3, eta_w=0.05, eta_tau=0.01,
                   eta_theta=0.01, w_hm=0.6, w_hh=0.25)

    W, plastic, roles = layered_graph(
        seed=seed,
        w_hm=cfg.get("w_hm", 0.6),
        w_hh=cfg.get("w_hh", 0.25)
    )

    return ClosedLoopPlasticityEngine(
        W, plastic, roles,
        axes=axes,
        act=cfg.get("act", 3),
        pas=cfg.get("pas", 5),
        theta=cfg.get("theta", 4.0),
        p_s=cfg.get("p_s", 3e-3),
        eta_w=cfg.get("eta_w", 0.05),
        eta_tau=cfg.get("eta_tau", 0.12),
        eta_theta=cfg.get("eta_theta", 0.02),
        seed=seed
    )
