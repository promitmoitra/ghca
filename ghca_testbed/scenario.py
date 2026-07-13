"""Scenario / ground-truth dataclasses.

A Scenario is a named, seeded ground-truth SCM. It OWNS its ground truth (data,
not code that re-derives it), so the harness and any user method read the answer
rather than recompute it. `kind` selects how it is realised:
  'scm'    — pure SCM, executable in-process (the C1 graphs);
  'cnet' / 'e3' / 'spiral' — substrate scenarios; ground truth is carried as
             verified data and full reproduction runs the referenced C-script.
"""

from dataclasses import dataclass, field
from typing import Any, Callable, Optional


@dataclass(frozen=True)
class GraphSpec:
    edges: tuple                       # tuple of (src, dst)
    S: tuple                           # observed conditioning set
    expected_verdict: str             # 'epiphenomenal' | 'causal'


@dataclass(frozen=True)
class GroundTruth:
    task: str                          # 'T1'..'T5'
    verdict: Optional[str] = None      # for T1
    scalars: dict = field(default_factory=dict)     # named expected values
    ordering: tuple = ()               # robust qualitative claims (strings)
    exact: bool = False                # T1 categorical/analytic == exact-portable
    notes: str = ""


@dataclass(frozen=True)
class Scenario:
    name: str
    track: str                         # 'C1'..'C7'
    task: str                          # 'T1'..'T5'
    kind: str                          # 'scm' | 'cnet' | 'e3' | 'spiral'
    ground_truth: GroundTruth
    graph: Optional[GraphSpec] = None
    seed: int = 0
    params: dict = field(default_factory=dict)
    run: Optional[Callable[..., Any]] = None   # in-process runner (scm only)
    script: str = ""                   # experiments/*.py for substrate reproduction

    @property
    def executable(self) -> bool:
        return self.run is not None
