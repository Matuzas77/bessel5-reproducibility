# What this release does and does not establish

All the work was done by AI. Jonas Matuzas is the responsible human author. No independent human review or proof-assistant formalization is asserted here.

## Mathematical status in the source

| Manuscript statement | Status asserted by the published v4 manuscript |
|---|---|
| Equation (1.2), the proposed evaluation of the pure fifth Bessel moment | Conjectural |
| Equation (2.2), modular norm and coefficient-twist identities | Modular propositions with supplied mathematical arguments and finite witnesses |
| Equation (3.2), relation between the two determinants | Corollary of the cited period formulas, with basis and orientation conventions in the ancillary notes |
| Equations (3.3) and the Bessel equality in (3.4) | Conjectural Bessel-to-modular comparisons |
| Equation (3.5), twisted symmetric-square critical-value/Petersson-norm identity | Modular proposition with supplied analytic argument and finite witnesses |

The code checks the finite computational pieces. Passing them is not independent validation of every analytic step in those arguments.

## Fresh numerical computation

The standard-library replay recomputes the source periods through hypergeometric series, recurrence relations and directed arithmetic, and computes the modular periods with separate decimal and binary implementations. It encloses rounding and truncation errors according to the supplied algorithms and analytic bounds. The odd Bessel row uses the inherited pairing and orientation inputs; it is not a completely separate direct quadrature measurement.

The primary runner rebuilds exact modular and Manin-symbol data, not numerical Bessel quadrature. Those original suites constitute a fresh replay of the **supplied computational pipeline**, not a new independent derivation of all its analytic inputs.

Software v1.1.0 additionally provides a **separate direct quadrature notebook**. It evaluates K0(t)^5 directly and the eta-defined Mellin integral separately, with no saved results or upstream numerical routines as inputs. It changes precision, coordinates and partitions and rejects wrong phase/normalization controls. This new route supplies empirical convergence evidence, not a rigorous global error enclosure. Its two sides share mpmath, and its AI implementation/review is not independent human peer review. [Detailed methods](INDEPENDENT_CHECK.md).

## Remaining limits

- No finite residual bound proves an exact equality.
- Analytic continuation, Stokes arguments, all-index identities and source theorems remain written mathematics; the programs do not formalize them.
- The package does not establish novelty or priority against all prior literature.
- A green workflow means the programmed checks passed for that commit and environment.
- The original paper's source version and claims remain unchanged by this packaging release.

See [upstream verification notes](../upstream/anc/verification_notes.md), [original provenance](../upstream/anc/PROVENANCE.md), and the [full mathematical arguments](../upstream/anc/proofs/).

