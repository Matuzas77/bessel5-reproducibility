# Provenance and attribution

Revision 3 for Jonas Matuzas, 15 September 2026. The proof and numerical provenance below was established in revision 2; the current revision adds the citation audit recorded in `revision_audit.md`. The four-page article includes two proved modular identities, a specialization of Chuang's determinant relation, and explicitly conjectural Bessel-to-modular comparisons. The author requested attribution to OpenAI's ChatGPT-6-pro for obtaining the conjecture and developing the research and verification programs. This is a process disclosure, not a citation to an invented model technical paper or a claim of external mathematical review.

## Reused numerical engines

The following source engines come from `BESSEL5_CONTACT_CORRECTED_COMPARISON_20260913.zip`, SHA-256:

```
54d3fe64251fb676879dd9947950eb07c0a26b1e727e1319f54f9261180b86d4
```

That archive is not nested or required at runtime. Its source/target certificate was rerun during the original preparation and its actual log remains labelled historical. The current article's own numerical replay was run again during this revision.

The original bytes are retained for:

- `engine/fixed_interval.py`: `bf804bfea110067c89c7bd502c1a6f70faf8185da652401c645acae0df874619`.
- `engine/certify_source.py`: `01ef2111aa03a3373ffe4f135752f7715a4a7ee05b0666f59331c5dc4a41c4a0`.
- `engine/certify_independent.py`: `1cd1ba3ddea9a46e42d2775d7d76291361a10c9e613cb1c99675e34b61567ed7`.
- `engine/certify_eta_independent.py`: `bd6433cf79567adc546c36617a6a0dc8b75226c790c617e8e4288b3a2b3c2312`.

The publication-specific `certify_lvalue.py`, `check_paper.py`, and `audit_algebra.py` independently define the eta critical value and compare the stated candidates. No stored source or target decimal initializes a calculation.

## Recovered modular proofs

The complete norm and twist proofs were recovered from these actual archives:

- `BESSEL5_PETERSSON_NORM_AND_RIGOROUS_CANDIDATE_20260905.zip`: SHA-256 `52549d117b2252733e57aeef43f8ad73f502449a172255d8ddb86510439b03f6`.
- `BESSEL5_EXACT_TWIST_PERIOD_REDUCTION_20260906.zip`: SHA-256 `9af0c6b1466acac2137340c49d6209aa6108aa8f0e37e90201d56ef96a837836`.

The full norm archive replay and the seven exact twist jobs were rerun during the revision-2 audit. Relevant proof text and exact source programs are included without nesting either archive. The written proofs were edited to remove obsolete numerical sections and make current file locations explicit; the mathematical norm and twist arguments were retained and checked.

`engine/verify_norm.py` adapts the independently written logarithmic-eta/Fraction norm verifier only in its output location and removal of an elapsed-time field. `engine/verify_twist_symbols.py` similarly adapts the independent ceiling-continued-fraction verifier. It reads just the raw 33-term witness in `twist_witness.json`; it does not load the primary reduced matrices. `primary/` retains the exact reconstruction programs, with temporary-workspace orchestration in `verify_primary.py`.

`engine/verify_modular_bridge.py` was added for this audit. It checks that the norm proof, twist proof and numerical eta calculation use the same form and embedding; it verifies the depleted coefficient twist through the level-300 Sturm bound, the relevant good-Hecke orthogonality data, the critical-value coefficient algebra, and the constants in Chuang's cofactor calculation. It explicitly records that no twisted symmetric-square special value was evaluated.

All eight distributed default JSON records were freshly regenerated. A replay compares them only after independently recomputing them in an empty workspace.

## Article format and author block

The class, font packages, A4 geometry, link colors and theorem conventions were copied from the recovered source `Shifted_Poles_and_Chamber_Cancellation_v11.tex`. In the revision-2 audit the public first page of arXiv:2609.00290 was inspected directly: it prints only Jonas Matuzas and jonas.matuzas@gmail.com, without an affiliation. The present author block matches that public block. No institutional affiliation is inferred from a template or personal history. No mathematical dependency on the Witten-zeta paper is asserted.

## Scholarly credit and boundaries

Lim--Tu--Yu supply the arithmetic modularity theorem, the determinant notation and the existing equation (44) symmetric-square conjecture. Ping-Hsun Chuang is the sole author of the Ramanujan Journal article; the appendix is in collaboration with Jeng-Daw Yu. Chuang's period pairing and determinant formulas supply the even/odd cofactor identity. Bailey--Borwein--Broadhurst--Glasser and Zhou supply the source Bessel identities. The added Broadhurst, Evans, Fresan--Sabbah--Yu and Zhou references explain the lineage and distinguish the weighted level-15 sector.

The older exact modular norm and twist proofs alone do not establish the special-value relation involving `L(chi_-4 Sym^2 f,2)`. Version 4 adds the independent analytic bridge in `proofs/critical_norm_bridge.md`, proving this purely modular equation (3.5). The Bessel question marks are retained. No Bessel equality, exhaustive priority result, independent human review, or proof-assistant formalization is claimed. No public article PDF or font file is redistributed.

## Revision-3 preservation check

The input is `Bessel5_Quadratic_Conjecture_Audited_v2_Arxiv.zip`. That historical
revision changed attribution, contextual explanation and parity terminology;
it did not alter the theorem proofs, verification code or eight numerical/
exact JSON output records. The new citation uses the exact supplied Paris
slides title, not the differently titled 2017 talk cited as [Bro17] by
Lim--Tu--Yu. Public PDFs and font files are not included. Version-4 replay logs are in `logs/`; historical v3 logs are not redistributed.


## Version 4, 16 September 2026

The uploaded modular critical-norm bridge bundle supplies the new analytic
argument and two finite reconstructions. Original archive SHA-256: `d5db74b4b9ed48e2a6597aee34cc16b267d1304983f447f786b794d9380442e4`.
The complete proof is edited in proofs/critical_norm_bridge.md. The independent
script changes only its output filename for publication use. The three
primary reconstruction scripts are retained verbatim. A new boundary/CRT
audit is supplied as verify_audit.py. Historical revision audits are under
history/ and are not current claim-status documents.
