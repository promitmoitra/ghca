# Girth decides whether dwelling attractors coexist — once Lemma E is factored out

**Status:** theory. The discriminator that
[`coherent_core.md`](coherent_core.md) named in its honest scope and did not
run. A prediction is recorded, **falsified at exactly one cell**, and then
repaired into a form certified 80/80.

Experiment: [`experiments/girth_parity.py`](../experiments/girth_parity.py) ·
Archive: `result/topology/girth_parity.npz`

Setting: arbitrary graphs, `θ = 1`, states 0 (receptive), `1..τa` (active),
`τa+1..S` (passive), `S = τa + τp`, `B = S + 1`. Exhaustive and deterministic —
no RNG anywhere in this doc.

---

## What was open

[`coherent_core.md`](coherent_core.md) proved Theorem C (the dwell-free
attractor set is exactly the static predicate `C`) and Lemma E (when `C` is
empty), but left one claim **certified, not proven**:

> every *live* attractor is dwell-free when `S+1 ≥ 4`

— certified only on square lattices, where 4 is the girth, and resting on the
single anomalous cell (1,1). The honest scope named the fix: a lattice of
different girth **and** parity. Nine graphs, chosen so that girth and parity
vary independently. Every property below is **measured**, not assumed.

| graph | n | girth | bipartite | cycle lengths |
| :--- | ---: | ---: | :---: | :--- |
| C3 ring | 3 | 3 | no | {3} |
| C4 ring | 4 | 4 | yes | {4} |
| C5 ring | 5 | 5 | no | {5} |
| C6 ring | 6 | 6 | yes | {6} |
| C7 ring | 7 | 7 | no | {7} |
| C8 ring | 8 | 8 | yes | {8} |
| square 3x3 | 9 | 4 | yes | {4, 6, 8} |
| triangular 3x3 | 9 | 3 | no | {3, 4, 5, 6, 7, 8, 9} |
| honeycomb 2-hex | 10 | 6 | yes | {6, 10} |

Rings are the clean instrument: a ring's girth is its length and nothing else
varies. The triangular patch supplies girth 3 with a full odd-and-even cycle
spectrum; the honeycomb supplies girth 6 with the same parity as the square
lattice.

## The prediction, and its falsification

> **H_g** — a live attractor dwells ⟺ `B < girth`.

**H_g is FALSIFIED: 79/80 rows, one failure** — square 3×3 at `(τa, τp) = (1,3)`,
where `B = 5 ≥ girth 4` yet all 28 live attractors dwell.

The failure is diagnostic, not fatal. At (1,3) the coherent set is **empty** for
Lemma E's arithmetic reason: `τa = 1` forces the lag sum around a cycle to equal
the cycle length exactly, and no multiple of `B = 5` lies in `{4, 6, 8}`. By
Theorem C, no dwell-free attractor can exist at all — girth notwithstanding — so
any live attractor must dwell.

## The repair

> **H_g′** — a dwelling live attractor exists ⟺ live attractors exist **and**
> (`B < girth` **or** `C` is empty).

**Certified 80/80.** It decomposes into one triviality and two substantive
halves, each separately checked:

| claim | statement | violations | rows in scope |
| :---: | :--- | ---: | ---: |
| **A** (trivial, from Theorem C) | `C` empty and some live attractor exists ⟹ *every* live attractor dwells | 0 | — |
| **B** (substantive) | `C` non-empty and `B ≥ girth` ⟹ **no** live attractor dwells | 0 | 39 |
| **C** (substantive) | `C` non-empty and `B < girth` ⟹ **some** live attractor dwells | 0 | 21 |

So the girth reading is **rescued, correctly conditioned**. There are *two*
obstructions, and the square lattice never separated them at the cells
previously tested:

- **Lemma E's arithmetic** decides whether dwell-free attractors can exist at all;
- **given that they can, girth** decides whether dwelling ones coexist.

`coherent_core.md`'s verdict that "the obstruction is arithmetic, not metric"
was therefore half the story, and is corrected in place there.

## The discriminating comparison

Girth does real work, and one pair isolates it with everything else held fixed —
same cell `(2,1)`, same `B = 4`, both bipartite:

| graph | girth | cell | B | B vs girth | live | dwelling |
| :--- | ---: | :---: | ---: | :---: | ---: | ---: |
| square 3x3 | 4 | (2, 1) | 4 | B ≥ girth | 12,049 | 0 |
| triangular 3x3 | 3 | (2, 1) | 4 | B ≥ girth | 24,907 | 0 |
| honeycomb 2-hex | 6 | (2, 1) | 4 | B < girth | 20,401 | 10 |

Same dynamics, same parameters, different girth, opposite answer. The
triangular lattice makes the point from the other side: girth 3 means
`B ≥ girth` **always** (since `S ≥ 2`), and it has **zero** dwelling attractors
at every affordable cell. The anomaly that made (1,1) special on the square
lattice simply does not exist there.

## Lemma E, extended off the square lattice

`coherent_core.py` proved Lemma E for any graph but only ever checked it on
square lattices — both bipartite, girth 4. It is re-checked here on
**non-bipartite** graphs (C3, C5, C7, triangular) and at girths 3, 5, 6, 7, 8:

- **0 necessity violations over 80 rows**, 39 of them non-bipartite;
- the converse held on all 80 as well — still **reported, not claimed**.

Theorem C and Theorem Z are also re-asserted inside every census call, so both
now carry empirical support on rings, triangular and honeycomb graphs rather
than square lattices alone.

## Full census

| graph | girth | (τa, τp) | B | \|C\| | live | dwell-free | dwelling | H_g | H_g' |
| :--- | ---: | :---: | ---: | ---: | ---: | ---: | ---: | :---: | :---: |
| C3 ring | 3 | (1, 1) | 3 | 6 | 2 | 2 | 0 | ✓ | ✓ |
| C3 ring | 3 | (2, 1) | 4 | 36 | 9 | 9 | 0 | ✓ | ✓ |
| C3 ring | 3 | (1, 2) | 4 | 0 | 0 | 0 | 0 | ✓ | ✓ |
| C3 ring | 3 | (2, 2) | 5 | 30 | 6 | 6 | 0 | ✓ | ✓ |
| C3 ring | 3 | (3, 1) | 5 | 90 | 18 | 18 | 0 | ✓ | ✓ |
| C3 ring | 3 | (1, 3) | 5 | 0 | 0 | 0 | 0 | ✓ | ✓ |
| C3 ring | 3 | (3, 2) | 6 | 102 | 17 | 17 | 0 | ✓ | ✓ |
| C3 ring | 3 | (2, 3) | 6 | 12 | 2 | 2 | 0 | ✓ | ✓ |
| C3 ring | 3 | (4, 1) | 6 | 174 | 29 | 29 | 0 | ✓ | ✓ |
| C3 ring | 3 | (1, 4) | 6 | 0 | 0 | 0 | 0 | ✓ | ✓ |
| C3 ring | 3 | (3, 3) | 7 | 84 | 12 | 12 | 0 | ✓ | ✓ |
| C4 ring | 4 | (1, 1) | 3 | 0 | 2 | 0 | 2 | ✓ | ✓ |
| C4 ring | 4 | (2, 1) | 4 | 84 | 21 | 21 | 0 | ✓ | ✓ |
| C4 ring | 4 | (1, 2) | 4 | 8 | 2 | 2 | 0 | ✓ | ✓ |
| C4 ring | 4 | (2, 2) | 5 | 40 | 8 | 8 | 0 | ✓ | ✓ |
| C4 ring | 4 | (3, 1) | 5 | 330 | 66 | 66 | 0 | ✓ | ✓ |
| C4 ring | 4 | (1, 3) | 5 | 0 | 0 | 0 | 0 | ✓ | ✓ |
| C4 ring | 4 | (3, 2) | 6 | 378 | 63 | 63 | 0 | ✓ | ✓ |
| C4 ring | 4 | (2, 3) | 6 | 72 | 12 | 12 | 0 | ✓ | ✓ |
| C4 ring | 4 | (4, 1) | 6 | 846 | 141 | 141 | 0 | ✓ | ✓ |
| C4 ring | 4 | (1, 4) | 6 | 0 | 0 | 0 | 0 | ✓ | ✓ |
| C4 ring | 4 | (3, 3) | 7 | 224 | 32 | 32 | 0 | ✓ | ✓ |
| C5 ring | 5 | (1, 1) | 3 | 0 | 2 | 0 | 2 | ✓ | ✓ |
| C5 ring | 5 | (2, 1) | 4 | 220 | 57 | 55 | 2 | ✓ | ✓ |
| C5 ring | 5 | (1, 2) | 4 | 0 | 2 | 0 | 2 | ✓ | ✓ |
| C5 ring | 5 | (2, 2) | 5 | 20 | 4 | 4 | 0 | ✓ | ✓ |
| C5 ring | 5 | (3, 1) | 5 | 1,320 | 264 | 264 | 0 | ✓ | ✓ |
| C5 ring | 5 | (1, 3) | 5 | 10 | 2 | 2 | 0 | ✓ | ✓ |
| C5 ring | 5 | (3, 2) | 6 | 1,350 | 225 | 225 | 0 | ✓ | ✓ |
| C5 ring | 5 | (2, 3) | 6 | 60 | 10 | 10 | 0 | ✓ | ✓ |
| C5 ring | 5 | (4, 1) | 6 | 4,410 | 735 | 735 | 0 | ✓ | ✓ |
| C5 ring | 5 | (1, 4) | 6 | 0 | 0 | 0 | 0 | ✓ | ✓ |
| C5 ring | 5 | (3, 3) | 7 | 280 | 40 | 40 | 0 | ✓ | ✓ |
| C6 ring | 6 | (1, 1) | 3 | 6 | 4 | 2 | 2 | ✓ | ✓ |
| C6 ring | 6 | (2, 1) | 4 | 748 | 189 | 187 | 2 | ✓ | ✓ |
| C6 ring | 6 | (1, 2) | 4 | 0 | 2 | 0 | 2 | ✓ | ✓ |
| C6 ring | 6 | (2, 2) | 5 | 150 | 32 | 30 | 2 | ✓ | ✓ |
| C6 ring | 6 | (3, 1) | 5 | 5,800 | 1,162 | 1,160 | 2 | ✓ | ✓ |
| C6 ring | 6 | (1, 3) | 5 | 0 | 2 | 0 | 2 | ✓ | ✓ |
| C6 ring | 6 | (3, 2) | 6 | 6,336 | 1,056 | 1,056 | 0 | ✓ | ✓ |
| C6 ring | 6 | (2, 3) | 6 | 24 | 4 | 4 | 0 | ✓ | ✓ |
| C6 ring | 6 | (4, 1) | 6 | 24,066 | 4,011 | 4,011 | 0 | ✓ | ✓ |
| C6 ring | 6 | (1, 4) | 6 | 12 | 2 | 2 | 0 | ✓ | ✓ |
| C6 ring | 6 | (3, 3) | 7 | 1,344 | 192 | 192 | 0 | ✓ | ✓ |
| C7 ring | 7 | (1, 1) | 3 | 0 | 4 | 0 | 4 | ✓ | ✓ |
| C7 ring | 7 | (2, 1) | 4 | 2,212 | 555 | 553 | 2 | ✓ | ✓ |
| C7 ring | 7 | (1, 2) | 4 | 0 | 2 | 0 | 2 | ✓ | ✓ |
| C7 ring | 7 | (2, 2) | 5 | 350 | 72 | 70 | 2 | ✓ | ✓ |
| C7 ring | 7 | (3, 1) | 5 | 24,570 | 4,916 | 4,914 | 2 | ✓ | ✓ |
| C7 ring | 7 | (1, 3) | 5 | 0 | 2 | 0 | 2 | ✓ | ✓ |
| C7 ring | 7 | (3, 2) | 6 | 27,804 | 4,636 | 4,634 | 2 | ✓ | ✓ |
| C7 ring | 7 | (2, 3) | 6 | 252 | 44 | 42 | 2 | ✓ | ✓ |
| C7 ring | 7 | (4, 1) | 6 | 129,318 | 21,555 | 21,553 | 2 | ✓ | ✓ |
| C7 ring | 7 | (1, 4) | 6 | 0 | 2 | 0 | 2 | ✓ | ✓ |
| C7 ring | 7 | (3, 3) | 7 | 5,530 | 790 | 790 | 0 | ✓ | ✓ |
| C8 ring | 8 | (1, 1) | 3 | 0 | 6 | 0 | 6 | ✓ | ✓ |
| C8 ring | 8 | (2, 1) | 4 | 6,500 | 1,627 | 1,625 | 2 | ✓ | ✓ |
| C8 ring | 8 | (1, 2) | 4 | 8 | 4 | 2 | 2 | ✓ | ✓ |
| C8 ring | 8 | (2, 2) | 5 | 360 | 74 | 72 | 2 | ✓ | ✓ |
| C8 ring | 8 | (3, 1) | 5 | 103,530 | 20,708 | 20,706 | 2 | ✓ | ✓ |
| C8 ring | 8 | (1, 3) | 5 | 0 | 2 | 0 | 2 | ✓ | ✓ |
| C8 ring | 8 | (3, 2) | 6 | 117,042 | 19,509 | 19,507 | 2 | ✓ | ✓ |
| C8 ring | 8 | (2, 3) | 6 | 840 | 142 | 140 | 2 | ✓ | ✓ |
| C8 ring | 8 | (4, 1) | 6 | 693,630 | 115,607 | 115,605 | 2 | ✓ | ✓ |
| C8 ring | 8 | (1, 4) | 6 | 0 | 2 | 0 | 2 | ✓ | ✓ |
| square 3x3 | 4 | (1, 1) | 3 | 84 | 62 | 28 | 34 | ✓ | ✓ |
| square 3x3 | 4 | (2, 1) | 4 | 48,196 | 12,049 | 12,049 | 0 | ✓ | ✓ |
| square 3x3 | 4 | (1, 2) | 4 | 136 | 34 | 34 | 0 | ✓ | ✓ |
| square 3x3 | 4 | (2, 2) | 5 | 47,520 | 9,504 | 9,504 | 0 | ✓ | ✓ |
| square 3x3 | 4 | (3, 1) | 5 | 763,040 | 152,608 | 152,608 | 0 | ✓ | ✓ |
| square 3x3 | 4 | (1, 3) | 5 | 0 | 28 | 0 | 28 | **✗** | ✓ |
| triangular 3x3 | 3 | (1, 1) | 3 | 1,602 | 534 | 534 | 0 | ✓ | ✓ |
| triangular 3x3 | 3 | (2, 1) | 4 | 99,628 | 24,907 | 24,907 | 0 | ✓ | ✓ |
| triangular 3x3 | 3 | (1, 2) | 4 | 2,624 | 656 | 656 | 0 | ✓ | ✓ |
| triangular 3x3 | 3 | (2, 2) | 5 | 284,390 | 56,878 | 56,878 | 0 | ✓ | ✓ |
| triangular 3x3 | 3 | (3, 1) | 5 | 1,152,780 | 230,556 | 230,556 | 0 | ✓ | ✓ |
| triangular 3x3 | 3 | (1, 3) | 5 | 2,480 | 496 | 496 | 0 | ✓ | ✓ |
| honeycomb 2-hex | 6 | (1, 1) | 3 | 42 | 24 | 14 | 10 | ✓ | ✓ |
| honeycomb 2-hex | 6 | (2, 1) | 4 | 81,564 | 20,401 | 20,391 | 10 | ✓ | ✓ |
| honeycomb 2-hex | 6 | (1, 2) | 4 | 0 | 10 | 0 | 10 | ✓ | ✓ |

## Substrate / analysis boundary

- **Substrate:** the step map, its cycles, and dwell events.
- **Analysis:** `C`, girth, bipartiteness and the cycle spectrum. These are
  constructs imposed on the enumeration; the substrate does not compute them.

## Honest scope

- **H_g′ is certified, not proven.** It is an exhaustive check over 80
  (graph, cell) rows on nine graphs of at most 10 nodes. No argument here
  transfers it to a family.
- **Small graphs.** The exhaustive budget is 4,000,000 configurations, so the
  honeycomb reaches only `B ≤ 4` and the triangular patch `B ≤ 5`. The
  honeycomb's `B < girth` half is therefore tested at `B = 3, 4` only — enough
  to discriminate against the square lattice at the same cell, not enough to
  see it switch over at `B = 6`.
- **The falsified prediction is asserted as falsified**, at exactly the named
  cell, so a future change that silently "fixes" H_g trips the assertion rather
  than quietly rewriting history.
- **Rings are not lattices.** Six of the nine graphs are cycles, chosen because
  they vary girth cleanly. The 2-D evidence is three patches.
- **`θ = 1`** throughout, as everywhere else in this thread.

## Reproduce

```bash
python3 experiments/girth_parity.py
```

Exhaustive, deterministic, RNG-free; asserts every claim above — including the
H_g failure set — and rewrites its archive. Runtime a few minutes.
