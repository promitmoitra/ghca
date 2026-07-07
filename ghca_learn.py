"""
Reward-modulated learning on a Greenberg-Hastings network.

Implements the learning framework of docs/learning_experiments.md section 2 on
top of the ``ghca_net.Network`` substrate:

  * strict scalar reward from the environment -> internal TD error delta,
  * an order-parameter critic (value read from the medium's collective state,
    updated on a faster timescale than plasticity),
  * two parallel plasticity lines sharing the same delta broadcast:
      Line A - conduction weights w_ij with edge eligibility traces,
      Line B - per-node timescales tau_i with period eligibility.

The learner is deliberately generic; individual experiments (E1..E5) supply the
graph, the node roles, and the trial protocol.
"""

import numpy as np
from ghca_net import Network


# ----------------------------------------------------------------------------
# Layered graph builder shared by the conditioning experiments.
# ----------------------------------------------------------------------------

def layered_graph(n_s=12, n_h=150, n_m=12, K=2, A=2,
                  fanin_sh=8, fanin_hm=14, hh_k=6, hh_beta=0.2,
                  w_sh=1.0, w_hm=0.4, w_hh=0.25, channel_bias=0.85,
                  init_jitter=0.15, seed=0):
    """Build a S -> H -> M graph with a recurrent hidden medium.

    Layout of node indices:
        sensory : K channels of n_s nodes   -> [0, K*n_s)
        hidden  : n_h recurrent nodes        -> [K*n_s, K*n_s + n_h)
        motor   : A channels of n_m nodes    -> trailing block

    Hidden nodes are assigned a preferred sensory channel and draw a fraction
    `channel_bias` of their S->H inputs from it, so the two stimuli evoke
    partly-separable hidden patterns (the substrate's innate repertoire that
    reward later binds to actions). The plastic H->M readout starts weak and
    jittered so the winning channel is not pre-determined and there is variance
    for exploration.

    Returns (W, plastic, roles); `plastic` marks the S->H and H->M edges Line A
    may modify.
    """
    rng = np.random.default_rng(seed)
    n_sen, n_mot = K * n_s, A * n_m
    N = n_sen + n_h + n_mot
    s0, h0, m0 = 0, n_sen, n_sen + n_h
    W = np.zeros((N, N))
    plastic = np.zeros((N, N), dtype=bool)

    sensory = [np.arange(s0 + c * n_s, s0 + (c + 1) * n_s) for c in range(K)]
    hidden = np.arange(h0, h0 + n_h)
    motor = [np.arange(m0 + c * n_m, m0 + (c + 1) * n_m) for c in range(A)]
    all_sensory = np.concatenate(sensory)
    hidden_pref = rng.integers(K, size=n_h)          # preferred channel per hidden

    # S -> H : channel-biased sampling (plastic)
    for idx, h in enumerate(hidden):
        pref = sensory[hidden_pref[idx]]
        n_pref = int(round(channel_bias * fanin_sh))
        src = np.concatenate([
            rng.choice(pref, size=min(n_pref, len(pref)), replace=False),
            rng.choice(all_sensory, size=max(fanin_sh - n_pref, 0), replace=False)])
        src = np.unique(src)
        W[h, src] = w_sh * (1.0 + init_jitter * rng.standard_normal(len(src)))
        plastic[h, src] = True

    # H -> H : recurrent small-world medium (fixed reservoir)
    for a, h in enumerate(hidden):
        for d in range(1, hh_k // 2 + 1):
            for nb in (hidden[(a - d) % n_h], hidden[(a + d) % n_h]):
                if rng.random() < hh_beta:
                    nb = rng.choice(hidden)
                if nb != h:
                    W[h, nb] = w_hh

    # H -> M : weak jittered readout (plastic)
    for c in range(A):
        for m in motor[c]:
            src = rng.choice(hidden, size=min(fanin_hm, n_h), replace=False)
            W[m, src] = w_hm * (1.0 + init_jitter * rng.standard_normal(len(src)))
            plastic[m, src] = True

    W = np.clip(W, 0.0, None)
    roles = {"sensory": sensory, "hidden": hidden, "motor": motor,
             "s0": s0, "h0": h0, "m0": m0, "N": N, "K": K, "A": A,
             "n_s": n_s, "n_m": n_m, "hidden_pref": hidden_pref}
    return W, plastic, roles


# ----------------------------------------------------------------------------
# The learner.
# ----------------------------------------------------------------------------

class GHLearner(Network):
    """GH network with reward-modulated conduction / timescale plasticity."""

    def __init__(self, W, plastic, roles, line="AB",
                 eta_w=0.02, eta_tau=0.05, lam=0.9, gamma=0.95, alpha_v=0.1,
                 w_max=4.0, tau_min=2, tau_max=23, tau_sigma=0.0, tau_mask=None,
                 tau_shared=False, **kw):
        super().__init__(W, **kw)
        self.plastic = np.asarray(plastic, dtype=bool)
        self.adj = self.W > 0
        self.roles = roles
        self.line = line                      # 'A', 'B' or 'AB'
        self.eta_w, self.eta_tau = eta_w, eta_tau
        self.lam, self.gamma, self.alpha_v = lam, gamma, alpha_v
        self.w_max = w_max
        self.tau_min, self.tau_max = tau_min, tau_max
        self.tau = self.tau.astype(float)     # timescales become continuous
        self.tau_base = self.tau.copy()       # perturbation baseline (Line B)
        self.tau_sigma = tau_sigma            # >0 enables tau exploration
        self.tau_mask = np.ones(self.N, bool) if tau_mask is None \
            else np.asarray(tau_mask, bool)
        self.tau_shared = tau_shared          # one shared timescale per group
        self.tau_scalar = float(self.tau[self.tau_mask].mean())
        self._eps = np.zeros(self.N)
        self._eps_s = 0.0

        self.E = np.zeros_like(self.W)         # edge eligibility (Line A)
        self.tau_elig = np.zeros(self.N)       # period eligibility (Line B)
        self.last_fire = np.full(self.N, -1.0)
        self.v_w = np.zeros(3)                 # critic weights on [A, R, 1]
        self._prev = self.phi.copy()

    def perturb_tau(self):
        """Draw a fresh timescale perturbation for a trial (node-perturbation
        exploration for Line B). Call at the start of each trial. In shared
        mode a single scalar perturbs the whole group (one regional timescale),
        which avoids the weakest-link credit-assignment problem of independent
        per-node perturbation (see e2_results.md)."""
        explore = self.tau_sigma > 0 and "B" in self.line
        if self.tau_shared:
            self._eps_s = self.tau_sigma * np.random.standard_normal() if explore else 0.0
            self.tau = self.tau_base.copy()
            self.tau[self.tau_mask] = np.clip(self.tau_scalar + self._eps_s,
                                              self.tau_min, self.tau_max)
        else:
            self._eps = self.tau_sigma * np.random.standard_normal(self.N) * self.tau_mask \
                if explore else np.zeros(self.N)
            self.tau = np.clip(self.tau_base + self._eps, self.tau_min, self.tau_max)

    # -- one learning step ---------------------------------------------------
    def step_learn(self, drive=None):
        active_prev = self.active_mask()
        self._prev = self.phi.copy()
        self.step(drive)
        fired = (self.phi == 1) & (self._prev == 0)

        # Line A edge eligibility: pre = active neighbour j, post = i just fired
        self.E *= self.lam
        fi, aj = np.where(fired)[0], np.where(active_prev)[0]
        if fi.size and aj.size:
            self.E[np.ix_(fi, aj)] += 1.0

        # Line B period eligibility: observed inter-fire interval vs current tau
        for i in fi:
            if self.last_fire[i] >= 0:
                self.tau_elig[i] = (self.t - self.last_fire[i]) - self.tau[i]
            self.last_fire[i] = self.t
        return fired

    # -- critic (order-parameter value) --------------------------------------
    def features(self):
        return np.array([self.active_mask().mean(), self.coherence(), 1.0])

    def value(self):
        return float(self.v_w @ self.features())

    # -- apply the reward-driven updates -------------------------------------
    def learn(self, delta, delta_b=None):
        """Apply reward-modulated plasticity. `delta` gates Line A (conduction);
        `delta_b` gates Line B (timescale) -- pass a separate value for FACTORED
        credit (each line credited by the error component it controls; see
        e3_factored_credit.py). If `delta_b` is None the same `delta` gates both
        (the original shared-scalar rule)."""
        db = delta if delta_b is None else delta_b
        if "A" in self.line:
            dW = self.eta_w * delta * self.E
            self.W[self.plastic] = np.clip(self.W[self.plastic] + dW[self.plastic],
                                           0.0, self.w_max)
            self.adj = self.W > 0
        if "B" in self.line:
            if self.tau_shared:
                # reinforce the shared-timescale perturbation that earned reward
                self.tau_scalar = float(np.clip(
                    self.tau_scalar + self.eta_tau * db * self._eps_s,
                    self.tau_min, self.tau_max))
                self.tau[self.tau_mask] = self.tau_scalar
            elif self.tau_sigma > 0:
                # per-node perturbation: reinforce the timescale perturbation
                # that earned reward. Reward is only earned when a loop
                # *sustains* (tau below the loop transit time), so this drives
                # tau toward the sustaining regime -- but hits a weakest-link
                # problem (the loop dies at the slowest node; see e2_results.md).
                self.tau_base = np.clip(self.tau_base + self.eta_tau * db * self._eps,
                                        self.tau_min, self.tau_max)
                self.tau = self.tau_base.copy()
            else:
                # resonance rule: drive tau toward the observed re-fire interval.
                # NOTE: this targets tau = loop period, i.e. the marginal (death)
                # boundary; it maintains an already-sustaining loop but does not
                # by itself find the sustaining regime (see e2_results.md).
                self.tau_base = np.clip(self.tau_base + self.eta_tau * db * self.tau_elig,
                                        self.tau_min, self.tau_max)
                self.tau = self.tau_base.copy()

    def update_critic(self, delta, feats):
        self.v_w += self.alpha_v * delta * feats

    def reset_traces(self):
        self.E[:] = 0.0
        self.tau_elig[:] = 0.0
        self.last_fire[:] = -1.0

    # -- I/O helpers ---------------------------------------------------------
    def sensory_drive(self, x):
        d = np.zeros(self.N, dtype=bool)
        d[self.roles["sensory"][x]] = True
        return d

    def motor_scores(self):
        act = self.active_mask()
        return np.array([act[m].sum() for m in self.roles["motor"]], dtype=float)
