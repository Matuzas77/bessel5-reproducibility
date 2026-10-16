# Direct-check notebook release v1.1.0

Date: 16 October 2026. **All the work was done by AI.**

- Adds direct numerical K0 quadrature, separate eta-Mellin quadrature, convergence checks and negative controls; no archived numerical answers are read.
- Publishes a self-contained, executed notebook with an Open in Colab link; ordinary CPU execution, no Drive or credentials in code.
- Records actual 40/60/80-digit runs and a clean-kernel notebook execution. The 80-digit cross-comparison passes the predeclared relative threshold 1e-70; empirical agreement does not prove the conjecture.
- Adds clean-kernel notebook CI and source-parity validation. Keeps all 65 upstream files untouched and all original certificate/reconstruction suites.
- No new scientific DOI, manuscript revision or external peer review is implied. The original v1.0.0 release remains available.

See [methods, limitations and evidence](INDEPENDENT_CHECK.md).

# Reproducibility packaging release v1.0.0 (historical)

Date: 16 October 2026.

All the work was done by AI. Jonas Matuzas is responsible for this release and its scientific claims.

This release preserves the 65 files in `Bessel5_Quadratic_Conjecture_v4_Source.zip` byte-for-byte and adds a reproducible execution environment, real execution logs, CI, documentation, citation metadata and deterministic packaging. It does not alter the manuscript, proof files, verifiers or saved scientific output records.

The scientific source is **Zenodo v4, 16 September 2026**, DOI **10.5281/zenodo.22793737**. The GitHub tag **v1.0.0** refers only to this reproducibility packaging. It is not manuscript revision 5, an independently peer-reviewed edition, or a new DOI deposit.

The source's `submission_metadata.txt` retains an older Zenodo URL (`22773594`) and the proposed arXiv category `math.NT` because upstream bytes are preserved. These are historical submission notes, not current indexing metadata: the reproduced source is Zenodo `22793737`, and the current arXiv record is classified `math.GM`. Use the DOI and record in `provenance/source.json` to identify the exact archive reproduced here. arXiv version numbering is also separate from Zenodo version numbering.

The dated local evidence records all three upstream suites. CI replays them separately on GitHub; consult the actual workflow result before describing a hosted run as passed. This release performs no PDF typesetting or direct Bessel quadrature.

