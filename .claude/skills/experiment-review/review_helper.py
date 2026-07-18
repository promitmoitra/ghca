#!/usr/bin/env python3
# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""CLI helper script for the Greenwich-Hastings Causal Analysis (GHCA) dual-track workflow.

Provides commands to:
1. audit-rng: Scan the codebase for global NumPy RNG usage (preventing the perturb_tau bug).
2. scaffold-review: Create a structured core or extensions review template.
3. scaffold-plan: Create a structured planning/roadmap template.
"""

import argparse
import os
import re
import sys


def find_python_files(root_dir):
  """Finds all core python files and experiment scripts in the repo."""
  python_files = []
  
  # Scan root directory for ghca_*.py files and test files
  for f in os.listdir(root_dir):
    if f.endswith('.py') and (f.startswith('ghca_') or f.startswith('test_')):
      python_files.append(os.path.join(root_dir, f))
      
  # Scan experiments/ directory
  exp_dir = os.path.join(root_dir, 'experiments')
  if os.path.exists(exp_dir):
    for f in os.listdir(exp_dir):
      if f.endswith('.py'):
        python_files.append(os.path.join(exp_dir, f))
        
  return sorted(python_files)


def audit_rng_command(args):
  """Scans Python files for forbidden global random usage."""
  print(">>> Scanning codebase for global RNG usage...", file=sys.stderr)
  
  root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
  py_files = find_python_files(root_dir)
  
  global_rng_patterns = [
    # Global numpy random calls (excluding default_rng and Generator)
    (re.compile(r'\bnp\.random\.(?!default_rng\b|Generator\b)[a-zA-Z0-9_]+\b'), "Forbidden np.random.<method> call"),
    # Direct use of standard library random module's global functions (excluding np.random)
    (re.compile(r'(?<!np\.)\brandom\.(choice|randint|uniform|shuffle|random|seed)\b'), "Forbidden global random.<method> call"),
  ]
  
  violations = []
  
  for file_path in py_files:
    rel_path = os.path.relpath(file_path, root_dir)
    try:
      with open(file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    except OSError as e:
      print(f"Warning: Could not read {rel_path}: {e}", file=sys.stderr)
      continue
      
    in_docstring = False
    docstring_delim = None
    
    for line_idx, line in enumerate(lines, 1):
      # Handle docstring tracking
      stripped = line.strip()
      
      # Basic docstring toggle
      if '"""' in stripped:
        if not in_docstring:
          in_docstring = True
          docstring_delim = '"""'
        elif docstring_delim == '"""':
          in_docstring = False
          docstring_delim = None
      elif "'''" in stripped:
        if not in_docstring:
          in_docstring = True
          docstring_delim = "'''"
        elif docstring_delim == "'''":
          in_docstring = False
          docstring_delim = None
          
      if in_docstring:
        continue
        
      # Remove inline comment for check
      code_part = line.split('#')[0]
      
      for pattern, desc in global_rng_patterns:
        matches = pattern.finditer(code_part)
        for m in matches:
          # Double check if it is part of variable names or imports
          match_str = m.group(0)
          violations.append({
            'file': rel_path,
            'line': line_idx,
            'content': line.strip(),
            'match': match_str,
            'desc': desc
          })
          
  if violations:
    print("\n❌ CRITICAL AUDIT FAILURE: Found global RNG usage!", file=sys.stderr)
    print("All experiments and substrate code MUST thread seeded Generator instances (using default_rng) explicitly.", file=sys.stderr)
    print("Never use the global NumPy RNG to avoid the 'perturb_tau' bug.\n", file=sys.stderr)
    for v in violations:
      print(f"  [{v['file']}:{v['line']}] - {v['desc']}: '{v['match']}'", file=sys.stderr)
      print(f"    Code: {v['content']}\n", file=sys.stderr)
    sys.exit(1)
  else:
    print("\n✅ AUDIT SUCCESS: No global RNG usage found! All files conform to the seeding guidelines.", file=sys.stderr)
    sys.exit(0)


def scaffold_review_command(args):
  """Generates a structured review document template."""
  review_type = args.review_type
  output_path = args.output
  
  os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
  
  if review_type == 'core':
    title = "Core Series Review — Integrity & Overreach Audit (E0–E6, C0–C4)"
    coverage_desc = "independent audit, E0–E6 / C0–C4"
    table_content = """| Series | Data reproduces? | Verdict |
|--------|:----------------:|---------|
| E0 substrate characterisation | [Yes/No] | [Pending/Supported] |
| E1 conditioning | [Yes/No] | [Pending/Supported] |
| E2 working memory | [Yes/No] | [Pending/Supported] |
| E3 timed response | [Yes/No] | [Pending/Supported/Overstated - see E3 deep-dive] |
| E4 attention | [Yes/No] | [Pending/Supported] |
| E5 executive control | [Yes/No] | [Pending/Supported] |
| E6 emergent categories | [Yes/No] | [Pending/Supported] |
| C0 instrumentation | [Yes/No] | [Pending/Supported] |
| C1 certificates | [Yes/No] | [Pending/Supported] |
| C2 fat-handed `do(W)` | [Yes/No] | [Pending/Supported] |
| C3 `do(θ)` well-posed | [Yes/No] | [Pending/Supported] |
| C4 outcome-relativity | [Yes/No] | [Pending/Supported] |"""
    deep_dive_section = """## E3 deep-dive — honest numbers, overstated framing

[Analyze the variance, bimodal spreads, or hand-picked operating points in the core series. Recall the E3 lesson: joint success identity AND timing must be verified per-seed, not hidden under a clean mean.]"""
  else:
    title = "Extensions Self-Audit — E7, C5–C7, E8.x, E9.x"
    coverage_desc = "self-audit, E7 / C5–C7 / E8.x / E9.x"
    table_content = """| experiment | reproduces? | headline (re-run == doc) |
|---|:--:|---|
| E7 Phase A (spiral mechanism) | [Yes/No] | [Details...] |
| E7 Phase B (rotation = rule) | [Yes/No] | [Details...] |
| C5 (`do(χ)` fat-handed) | [Yes/No] | [Details...] |
| C6 (`do(θ_χ)`, necessity) | [Yes/No] | [Details...] |
| C7 (outcome-relativity) | [Yes/No] | [Details...] |
| E8.0–8.4 (prediction) | [Yes/No] | [Details...] |
| E8.5 (nested) | [Yes/No] | [Details...] |
| E8.6 (order-preserving) | [Yes/No] | [Details...] |
| E8.7 (conditional) | [Yes/No] | [Details...] |
| E9 (emergent conjunction) | [Yes/No] | [Details...] |"""
    deep_dive_section = """## Substrate/Analysis Boundary Caveats

[Analyze computed readouts/features vs emergent circuits. Ensure we make it clear that 'the substrate does X' translates to 'the substrate's dynamics make X linearly decodable', rather than 'the substrate learned X end-to-end' when hand-built features/readouts are used.]"""

  template = f"""# {title}

*An integrity/overreach audit of the work in this repository ({coverage_desc}), checking for the failure modes an AI-assisted programme is prone to: fabricated numbers, claims not backed by code or data, invented citations, and grandiose theoretical framing disconnected from what the code actually does.*

Reviewer pass date: [YYYY-MM-DD] (produced on the independent-review branch).

## Bottom line

[Skeptical summary of the audit. Outline the bottom-line truth of the claims, any residual risk of ambitious framing resting on small hand-built toy simulations, and specific material exceptions.]

## Method

For each series, the auditor did three checks:
1. **Number match** — loaded saved data and compared the saved arrays against specific quantitative claims in the doc.
2. **Citation check** — verified external references actually exist and are accurately characterized.
3. **Prose/altitude** — audited the result docs for overreach or unhedged claims.

## Per-series verdicts

{table_content}

## Citations — check for hallucinations

[Verify all external references cited in the docs. Provide titles, authors, dates, and confirm their contents match the claims.]

{deep_dive_section}

## What is genuinely good (honesty markers)

- [List sections where the docs explicitly detail limitations, hand-chosen operating points, small n, or failure conditions.]

## Not fully verified

- [List files, mathematical proofs, or conditions that were not exhaustively checked in this audit.]
"""

  try:
    with open(output_path, 'w', encoding='utf-8') as f:
      f.write(template)
    print(f"Success! Review scaffold written to: {output_path}")
  except OSError as e:
    print(f"Error writing review scaffold to {output_path}: {e}", file=sys.stderr)
    sys.exit(1)


def scaffold_plan_command(args):
  """Generates a structured planning/roadmap document template."""
  output_path = args.output
  
  os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
  
  template = """# Next Steps — Planning / Roadmap

*Candidate directions after recent experiment series and audits. This is a menu, not a commitment: each option lists what it is, why it matters, what it would take, effort/risk, and how it connects to existing work and to the honesty caveats the review passes surfaced. Nothing here is built yet.*

## Where the project stands

- **E-series (learning).** [Recap status of E-series experiments...]
- **C-series (causality).** [Recap status of C-series causal analysis...]
- **Reviews.** [Recap most recent independent and self-audits, detailing any reproducibility bugs resolved.]

## The three tensions that should steer the roadmap

*The audits converge on three honest limitations. Good next steps retire one of these.*

1. **Afforded vs learned.** [Details on readouts/features vs emergent circuits.]
2. **Narrow evidence.** [Details on small n, single substrate, hand-chosen operating points.]
3. **Illustrative, not proof.** [Details on toy simulations vs rigorous math/SCMs.]

---

## Track 1 — Close the substrate/analysis gap

### Track 1a. [Title]
- **What.** [Explain what will be built...]
- **Why.** [Explain how it retires active tensions and the value of this build...]
- **Effort/Risk.** [Estimate effort and risk factors...]
- **Connection.** [Connect with existing experiments and files...]

---

## Track 2 — Expand evidence base

### Track 2a. [Title]
- **What.** [Explain what will be built...]
- **Why.** [Explain why this matters...]
- **Effort/Risk.** [Estimate effort...]
- **Connection.** [Explain connections...]
"""

  try:
    with open(output_path, 'w', encoding='utf-8') as f:
      f.write(template)
    print(f"Success! Planning scaffold written to: {output_path}")
  except OSError as e:
    print(f"Error writing planning scaffold to {output_path}: {e}", file=sys.stderr)
    sys.exit(1)


def main():
  parser = argparse.ArgumentParser(description='GHCA Dual-Track Experiment-Review CLI Helper')
  subparsers = parser.add_subparsers(dest='command', required=True)
  
  # --- Subcommand: audit-rng ---
  subparsers.add_parser(
    'audit-rng',
    help='Scan codebase for forbidden global random usage'
  )
  
  # --- Subcommand: scaffold-review ---
  p_review = subparsers.add_parser(
    'scaffold-review',
    help='Scaffold a fresh review audit document'
  )
  p_review.add_argument(
    '--type',
    required=True,
    dest='review_type',
    choices=['core', 'extensions'],
    help='Type of review (core: E0-E6/C0-C4, extensions: E7/C5-C7/E8+)'
  )
  p_review.add_argument(
    '--output',
    required=True,
    help='Path to write the scaffold file (e.g. docs/extensions_review.md)'
  )
  
  # --- Subcommand: scaffold-plan ---
  p_plan = subparsers.add_parser(
    'scaffold-plan',
    help='Scaffold a fresh roadmap/planning document'
  )
  p_plan.add_argument(
    '--output',
    required=True,
    help='Path to write the scaffold file (e.g. docs/next_steps.md)'
  )
  
  args = parser.parse_args()
  
  if args.command == 'audit-rng':
    audit_rng_command(args)
  elif args.command == 'scaffold-review':
    scaffold_review_command(args)
  elif args.command == 'scaffold-plan':
    scaffold_plan_command(args)
  else:
    print(f"Unknown command: {args.command}", file=sys.stderr)
    sys.exit(1)


if __name__ == '__main__':
  main()
