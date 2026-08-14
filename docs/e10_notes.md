# E10 (WIP) — Emergent timescale hierarchy: handoff notes

**Status: paused after a negative mechanism result. Not merged to main.**
This branch (`claude/e10-timescale-hierarchy`) records an exploration of
[`next_steps.md`](next_steps.md) **Track 4a** — *do learned per-node timescales
`τ` self-organise into a fast/slow hierarchy with cross-frequency coupling?* — and
the finding that the **existing Line B rule cannot do it**. Kept so the next
attempt starts from the diagnosis, not from scratch. Evidence:
[`experiments/e10_diagnostics.py`](../experiments/e10_diagnostics.py).

## The goal

With per-node `τ` plastic (Line B), drive a recurrent medium with a two-rhythm
input (fast period `P_f`, slow `P_s`) and ask whether the learned `τ` distribution
becomes **bimodal** — a fast cluster near `P_f` and a slow cluster near `P_s` —
i.e. a timescale hierarchy self-organises, ideally with theta–gamma-style
cross-frequency coupling. This would make E8.5's *two hand-set timescales*
emergent, and is the deepest genuinely-new phenomenon the substrate could show.

## The finding (why it's paused)

**The existing Line B resonance rule `τ ← τ + η·(interfire − τ)` is structurally
incapable of building a fast/slow hierarchy.** It can only ratchet `τ` **upward**,
toward *multiples* of a drive period, never down to a fundamental. Three
diagnostics, most decisive last (all in `e10_diagnostics.py`):

| diagnostic | setup | result |
|---|---|---|
| **1. naive** | recurrent pool, two-rhythm drive | two-rhythm ≡ fast-only; **no** cluster at `P_f`; the slow rhythm is invisible (fast drive fires every node first) |
| **2. competitive** | + E9-style k-WTA/conscience channel specialisation (`e10_proto2`, not kept) | still no `P_f` cluster; fast channel (active 33% of the time) swamps the slow (8%) |
| **3. ideal (decisive)** | isolated nodes (no recurrence), each hand-fed a **single** rhythm | fast group (fed `P_f=6`) → `τ≈25`; slow group (fed `P_s=24`) → `τ` rails at ceiling |

The mechanism of the failure is clear from diagnostic 3: once a node's `τ` exceeds
the fundamental period, it becomes refractory through the next pulse, so its
*observed* inter-fire interval is a larger multiple of the period, which pulls `τ`
up further — a one-way ratchet toward the ceiling. A small-`τ` "fast" population
therefore cannot form from any reasonable initialisation. This is the concrete
form of the risk `next_steps.md` flagged ("the `perturb_tau` mechanism is coarse").

## Correction to the earlier de-risking claim

I had said E9 "de-risks 4a" by supplying reusable representation-learning
machinery. That is **only half right**: E9's competitive k-WTA + conscience does
solve the **grouping** problem (which nodes handle which scale), and is reusable.
But 4a's actual blocker is the **`τ`-value learning rule**, which is orthogonal to
E9 and is the broken part. So 4a is a *mechanism-design* task, not a
*reuse-E9* task. (Merging E9 first was still fine — its grouping machinery is
genuinely useful here — but the effort estimate was wrong.)

## Proposed path for the next attempt

Replace the inter-fire resonance rule with a **bidirectional rule that tunes `τ`
to the dominant period of a node's *input*** (not its own firing) — measurable
locally as the interval between input-drive peaks the node receives, and able to
move `τ` **down** as well as up. Combine with:

- **E9's competition** (k-WTA + DeSieno conscience) to specialise each node to one
  rhythm *channel* (fast vs slow), so a fast-specialised node's input period is
  `P_f` and its `τ` is pulled to `P_f` — restoring the reusable half of E9.
- **Balanced channel activity** — give the fast and slow channels comparable
  activity mass (e.g. equal drive amplitude × duty cycle) so competition is not
  swamped by the more frequent rhythm (diagnostic 2's failure).

Validation targets: `τ` histogram unimodal → **bimodal** at (`P_f`, `P_s`);
control (single-rhythm drive) stays unimodal; then measure phase–amplitude
coupling between the slow and fast populations, and a functional benefit
(hierarchical network tracks a two-timescale signal better than a `τ`-homogeneous
one). Effort: real but bounded mechanism design; medium risk it needs iteration.

## Files on this branch

- `experiments/e10_diagnostics.py` — reproduces diagnostics 1 and 3 (the decisive
  ratchet result). Diagnostic 2's competitive prototype was exploratory and not
  kept; its conclusion (fast channel swamps slow) is folded into the notes above.
- `docs/e10_notes.md` — this file.

Nothing here is a validated result; it is a signposted dead end plus a concrete
next step.
