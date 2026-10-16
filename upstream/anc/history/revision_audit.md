> HISTORICAL AUDIT — superseded by version 4 for equation (3.5).
> The current mathematical status is in ../revision_audit_v4.md and
> ../proofs/critical_norm_bridge.md. This file records the earlier audit only.

# Revision 3: Broadhurst Paris 2017 citation and overlap audit

**Manuscript:** *A quadratic critical-value conjecture for the fifth Bessel moment*  
**Author:** Jonas Matuzas  
**Revision date:** 15 September 2026  
**Outcome:** citation and contextual revisions; all mathematical formulas, proofs, programs, and eight output records retained.

## Decision

The supplied lecture should be cited directly. It was absent from the revision-2 bibliography, which cited Broadhurst's 2016 paper instead. Revision 3 adds the exact Paris talk as reference [5], retains the 2016 reference, and gives pinpoint citations in the introduction.

## Source and overlap

D. Broadhurst, *L-series from Feynman diagrams with up to 22 loops*, Workshop on Multi-Loop Calculations (methods and applications), Paris, 7 June 2017. [Official slides](https://multi-loop-2017.sciencesconf.org/data/program/Broadhurst.pdf); [official program](https://multi-loop-2017.sciencesconf.org/program.html).

The relevant PDF pages are 4 (moment notation), 4–5 and 9 (weighted five-Bessel examples), and 13 and 15 (even-power/twisted-series counterparts and matrix relations). The title page credits D. P. Roberts for the joint work. The explicit formula of Conjecture 1.1 was not found in this 18-page document. This is a comparison with the supplied source, not an exhaustive novelty search.

The introduction now identifies the target as `M(0,5,0)` and disclaims novelty for the general Bessel–L-value framework. The displayed odd-power, level-15 examples are distinguished from the present `dt`, level-60 conjecture. The source's empirical announcements are not promoted to proofs of this note's conjectures.

The Paris title is not silently replaced by *Critical L-values for products of up to 20 Bessel functions*, the differently titled lecture listed as [Bro17] in [Lim–Tu–Yu](https://arxiv.org/pdf/2609.13576v1). No claim that those are the same document is made.

## A notation correction worth making

Section 3 now explains that the even/odd subscripts in Lim–Tu–Yu's determinant notation refer to the powers of `I_0`, while all entries here contain even powers of `t`. This follows directly from their equation (42), not from an inference based on the names of the determinants. It prevents confusion between the two uses of parity.

## What remains unchanged

The leading formula, all 9 numbered mathematical displays, the two proved modular identities, Chuang's attributed determinant corollary, and every question mark are unchanged. In particular, no proof of the Bessel comparison or of the additional twisted symmetric-square critical-value comparison is asserted. The abstract, author/email block, preamble, and ChatGPT-6-pro disclosure are also unchanged.

The review checks the new attribution and guards against unintended mathematical changes. It is not a new independent analytic proof audit of every theorem used in revision 2. The earlier detailed six-point review remains in `revision_audit_v2.md`.

## Executed checks

`logs/preservation_and_citation_check_v3.json` records byte-for-byte preservation of the scientific source files. Both reproduction routes were run:

```sh
python -B anc/verify.py
python -B anc/verify_primary.py
```

The default route uses only Python's standard library: nine jobs passed and all eight freshly regenerated JSON records matched the delivered records. The optional route uses SymPy and mpmath: all thirteen exact reconstruction jobs passed. Current logs are `logs/v3_default_replay.txt` and `logs/v3_primary_replay.txt`.

The reproduced conservative absolute-discrepancy bounds remain below `10^-358` for the individual integral and below `10^-352` for each determinant comparison. The intervals contain zero; no exact equality follows from these bounds. No additional digits or newly proved identity is claimed.

The manuscript compiled twice to four pages without warnings, unresolved references, or overfull/underfull boxes. All four pages were rendered and visually checked. All eleven bibliography entries are defined and cited, and the new PDF hyperlink points to the exact supplied URL. A focused source diff is in `logs/manuscript_v2_to_v3.diff`.

The final source archive was checked for CRC and manifest integrity, extracted into a fresh directory, compiled there to a byte-identical PDF, and replayed using both verification routes. The archive contains the source, proofs, programs, and outputs, but no public article PDFs, font files, cached bytecode, or nested archives. No arXiv submission has been made.
