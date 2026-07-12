"""Scoring harness + reproduce entry point.

`score(method_fn)` runs a user's causal method over the registry and reports, per
scenario, whether it matches the ground truth. `method_fn=None` runs the built-in
reference implementations (executable C1 scenarios); substrate scenarios (C2–C7)
carry their ground truth as data and are reproduced by running their C-script.
"""

from dataclasses import dataclass, field
from typing import Any, Optional

from .registry import REGISTRY


@dataclass
class ScoreRow:
    scenario: str
    task: str
    quantity: str
    expected: Any
    got: Any
    ok: Optional[bool]        # None = carried / not scored in-process
    note: str = ""


@dataclass
class ScoreReport:
    rows: list = field(default_factory=list)

    @property
    def scored(self):
        return [r for r in self.rows if r.ok is not None]

    @property
    def passed(self):
        return sum(1 for r in self.scored if r.ok)

    @property
    def total(self):
        return len(self.scored)

    @property
    def ok(self):
        return self.total > 0 and self.passed == self.total

    def __str__(self):
        lines = []
        for r in self.rows:
            if r.ok is None:
                mark = "·"
                val = r.note
            else:
                mark = "OK " if r.ok else "XX "
                val = f"expected {r.expected}, got {r.got}"
            lines.append(f"  [{mark}] {r.scenario:16s} {r.task} {r.quantity:14s} {val}")
        lines.append(f"\n  {self.passed}/{self.total} scored checks passed"
                     f"  ({sum(1 for r in self.rows if r.ok is None)} substrate scenarios carried)")
        return "\n".join(lines)


def _tol(expected):
    try:
        return max(0.02, 0.05 * abs(float(expected)))
    except (TypeError, ValueError):
        return 0.02


def score(method_fn=None, scenarios=REGISTRY, tasks=None):
    """Score `method_fn` (or the reference impls) against the ground truth.

    `method_fn(scenario) -> dict` of measured quantities keyed like the
    scenario's ground-truth scalars (+ optional 'verdict'). `None` runs the
    executable reference scenarios; non-executable scenarios are reported carried.
    """
    report = ScoreReport()
    for name in sorted(scenarios):
        sc = scenarios[name]
        if tasks and sc.task not in tasks:
            continue
        gt = sc.ground_truth
        if method_fn is None:
            if not sc.executable:
                report.rows.append(ScoreRow(sc.name, sc.task, "(carried)", "-", "-",
                                            None, f"run {sc.script}"))
                continue
            got = sc.run()
        else:
            got = method_fn(sc)
        if gt.verdict is not None and "verdict" in got:
            report.rows.append(ScoreRow(sc.name, sc.task, "verdict", gt.verdict,
                                        got["verdict"], got["verdict"] == gt.verdict))
        for q, exp in gt.scalars.items():
            if q in got:
                ok = abs(float(got[q]) - float(exp)) <= _tol(exp)
                report.rows.append(ScoreRow(sc.name, sc.task, q, exp, round(float(got[q]), 4), ok))
    return report


def reproduce(mode="report"):
    """Print ground truth + the model's own reference results.

    mode 'report' (human), 'regression' (exit nonzero on any scored FAIL).
    Returns the ScoreReport.
    """
    print("=== ghca_testbed reproduce ===\n")
    print("C1 scenarios executed in-process (Definition-1 + Theorem-1):")
    rep = score(method_fn=None, tasks=("T1",))
    print(rep)
    print("\nSubstrate scenarios (C2–C7) — ground truth carried; reproduce by script:")
    for name in sorted(REGISTRY):
        sc = REGISTRY[name]
        if sc.task == "T1":
            continue
        gt = sc.ground_truth
        print(f"  {sc.track} {sc.name:12s} [{sc.task}]  {sc.script}")
        for q, v in gt.scalars.items():
            print(f"       {q:18s} = {v}")
    return rep


if __name__ == "__main__":
    import sys
    m = sys.argv[1] if len(sys.argv) > 1 else "report"
    rep = reproduce(mode=m)
    if m == "regression" and not rep.ok:
        sys.exit(1)
