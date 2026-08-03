# Draft: collaboration enquiry — Komoto et al., *ACS Nano* 2026 (autonomous nanopore)

> **Status: internal draft, not published.** Not in the MkDocs nav. Placeholders
> `[[ ]]` need filling before sending. Rationale and the technical case are in
> [`aria_positioning.md`](aria_positioning.md) (Addendum).

**To:** Makusu Tsutsui — `makusu32@sanken.osaka-u.ac.jp` (corresponding author, SANKEN,
The University of Osaka)
**Cc (optional, second round):** Yuki Komoto (first author, SANKEN); Denis Garoli
(IIT Genova / Modena) — see notes below on whether to cc at all.

---

## The email

**Subject:** Spike-sequence structure in your autonomous nanopore — a reanalysis
proposal (acsnano.6c08258)

Dear Professor Tsutsui,

I read your recent *ACS Nano* paper on the autonomous solid-state nanopore
(doi:10.1021/acsnano.6c08258) with a lot of interest, and I am writing with a specific
reanalysis proposal that I think could add to the result. I have no nanopore
background — my work is computational, on learning in excitable media — which is
precisely why one thing in the paper stood out to me.

You make two claims that sit in productive tension. The pore is **stateful**: as you
put it, "ionic spikes already retain a short-lived history, as their amplitude and
duration reflect the preceding chemical evolution of the pore." But the classification
treats each spike **independently** — per-event waveform features into XGBoost. If the
statefulness is real, then the *sequence* — the interspike-interval structure, and how
each spike's features depend on the ones before it — should carry analyte information
that a per-event classifier cannot see. Your own measurements already point this way:
dwell time and firing rate both order by nucleotide, and those are properties of the
spike train rather than of any single spike.

The concrete proposal is to test that on your existing recordings, with a stated kill
condition:

> **Hypothesis.** A sequence model over the spike train (interspike-interval structure
> and history-conditioned features) improves four-nucleotide discrimination over the
> per-spike baseline.
> **Kill condition.** If it does not beat per-spike XGBoost under identical
> cross-validation, the pore's statefulness carries no additional decodable analyte
> information, and that is a clean negative result worth knowing.

Given F_m = 0.68 for four nucleotides and 0.37 for seven amino acids, and your own note
that "improved waveform models will be required," there seems to be headroom worth
testing. This needs no new experiments on your side.

**The ask, in ascending order of effort for you:**

1. **Minimal —** an event table for the 50 mM nucleotide conditions: one row per
   detected spike, with timestamp, amplitude, dwell time, and the nucleotide label.
   That alone is enough to test the interspike-interval part of the hypothesis.
2. **Fuller —** the raw ionic-current traces for those conditions, which would also
   let the waveform and the sequence structure be modelled jointly.

**Terms, to be clear up front:** it is your data and your call. I would publish nothing
without your agreement, and I would be equally happy to simply run the analysis and
send you the results and code to use as you see fit, with no authorship expectation. If
it becomes a joint output, the authorship arrangement is entirely yours to set.

My code and results are open by default — `[[ repo URL ]]` — so anything I did with the
data would be fully reproducible on your side.

If this is not of interest, or the data cannot be shared, no reply is needed and I
will not follow up.

With thanks and best wishes,
`[[ Name ]]`
`[[ email ]]` · independent researcher

---

## Notes on sending (not part of the email)

**Why this framing.** The ask is small, the hypothesis is falsifiable with a stated
kill condition, and the offer is explicitly *no-strings* — a careful group is far more
likely to share data with someone who pre-commits to a negative result than with
someone who sounds like they are mining for a paper. The "no reply needed" close
removes any social cost to ignoring it.

**Why option 1 matters.** An event table is a spreadsheet export from an analysis they
have already done. Raw traces may involve institutional data-sharing friction, large
files, or reluctance. Offering the small version first makes yes much cheaper, and the
ISI hypothesis — the core of the claim — is testable from the event table alone.

**On cc'ing.** Recommend emailing Tsutsui alone first. Komoto (first author) did the
analysis and is the person who would actually produce the export, but going to the
corresponding author first is the convention and avoids appearing to route around him.
Add Komoto on the second exchange if it progresses. Garoli (IIT) is a separate group
on the paper — only relevant if the collaboration broadens.

**What to have ready if they say yes.** A one-page analysis plan, pre-registered:
baseline reproduction (their per-spike XGBoost first, to establish we can match their
numbers before claiming to beat them), then the sequence model, with the CV protocol
fixed in advance. Reproducing their baseline first is both good practice and the
fastest way to earn trust.

**Realistic expectation.** Cold data requests to busy groups usually go unanswered.
This costs one email and gates a whole line of work, so it is worth sending, but the
plan should not depend on a reply. If there is no response in ~3 weeks, the fallbacks
are public nanopore spike datasets from other groups, or the synthetic route (simulate
an excitable-pore model, which the repo's substrate already essentially is).
