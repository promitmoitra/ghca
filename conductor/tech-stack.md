# Technology Stack — GHCA

## Core Environment
- **Language:** Python 3.10+
- **Numerical Computation:** NumPy (using `default_rng`), SciPy
- **Graph & Topology:** NetworkX
- **Data Analysis & ML:** scikit-learn (used in evaluation suites), PyTorch (optional GPU acceleration)
- **Visualization & Animation:** Matplotlib, ImageIO / Pillow (GIF generation), MkDocs (`mkdocs-material` for docs site)

## Testing & Quality Assurance
- **Test Runner:** `pytest`
- **Integrity Harness:** `reproduce_all.py` (automated verification of E0–E9, C0–C7, test suites)
- **RNG Auditor:** `.claude/skills/experiment-review/review_helper.py audit-rng`
