# Sources, inherited inputs, and limits

Standard inputs retained from the 5 and 6 September 2026 proofs; reviewed for this revision on 15 September 2026.

## Mathematical source ledger

1. W. Stein, *Modular Forms: A Computational Approach*, “Eisenstein Series and Bernoulli Numbers”, especially Theorems 5.8–5.10. Used for the standard modularity and coefficients of Eisenstein series. Our norm coefficient and finite witness are computed here, not taken from this source.
   `https://wstein.org/books/modform/modform/eisenstein.html`

2. A. Agashe and W. Stein, *Generating the Hecke algebra*. Used for the finite Sturm/Hecke method; the standard characteristic-zero coefficient-pairing argument is applied to the cusp character space. The bound 36 is computed here from weight 3 and index 144.
   `https://wstein.org/papers/generating_hecke/generating_hecke/`

3. W. Stein, *Modular Forms: A Computational Approach*, “Dimension Formulas”. Standard dimension framework. The particular full-space dimension 30 is derived in the proof from valence and Riemann–Roch, with all cusp and torsion data checked exactly.
   `https://wstein.org/books/modform/modform/dimension_formulas.html`

4. NIST Digital Library of Mathematical Functions, 10.31.2, 10.32.9 and 5.11.2 with 5.11(ii). The logarithmic K0 series, positive integral representation, and first-neglected-term bound for the positive-real digamma expansion support the retained Bessel enclosure proof.
   `https://dlmf.nist.gov/10.31`
   `https://dlmf.nist.gov/10.32`
   `https://dlmf.nist.gov/5.11`

5. Classical Dedekind-eta transformation and cusp-order formulas, Petersson adjointness, and regular weight-two Eisenstein continuation. All conventions actually used are printed in `norm_proof.md`. The exact Fricke coefficient is computed by rational Dedekind sums and cyclotomic arithmetic, and full proportionality is proved by the new full-space eigenline calculation.

6. P.-H. Chuang (with an appendix in collaboration with J.-D. Yu), *On the periods of twisted moments of the Kloosterman connection*, arXiv:2306.15216 and published in *The Ramanujan Journal*. Background for Bessel period relations only. It is not used as a proof of any of the three proposed individual Bessel-to-L-value formulas.
   `https://arxiv.org/abs/2306.15216`



# Sources and logical dependencies

## Primary mathematical sources

1. William Stein, *Modular Forms: A Computational Approach*, chapter “General Modular Symbols.”
   https://wstein.org/books/modform/modform/modular_symbols.html
   Used for the Manin-symbol presentation with character (Theorem 1.25 and the character subsection), the adjugate polynomial action, the critical-value integration normalization, and Hecke compatibility (Theorem 1.42). The specialized rational witness for level 60 and conductor 5 is printed and independently reconstructed in this package. The source does not assert the fifth-Bessel identity.

2. William Stein, “Bounding the number of generators,” in the note on generating the Hecke algebra.
   https://wstein.org/papers/generating_hecke/generating_hecke/node3.html
   Used for the finite Sturm bound in characteristic zero. With character, the same bound follows by raising a hypothetical nonzero difference to the finite order of its character and applying the trivial-character vanishing bound. The particular indices 144,720 and bounds 36,180 are calculated here.

3. David Savitt, *An elementary proof of Newman's eta-quotient theorem*, arXiv:2507.16225, Theorem 1.1 and Remark 1.2.
   https://arxiv.org/html/2507.16225v1
   Supplies the general integer-weight eta-quotient criterion including the nontrivial Nebentypus character, and the cusp positivity condition. The squareclasses in the present case are explicitly verified as 15 times rational squares, with odd weight giving character−15.

   Jeremy Rouse and John J. Webb, *On spaces of modular forms spanned by eta-quotients*, arXiv:1311.1460.
   https://arxiv.org/abs/1311.1460
   Used here for the cusp-order formula on the compact modular curve. The nontrivial-character criterion is supplied by Savitt/Newman above. Every specialized congruence, character squareclass, and positive cusp order is explicit in the proof or generated data.

4. NIST DLMF, Section23.18, modular transformations of Dedekind's eta function.
   https://dlmf.nist.gov/23.18
   Used for the imaginary-axis inversion in the numerical check and as the classical eta transformation input. The cusp criterion and theta identity have separate sources and finite checks.

5. NIST DLMF, Section10.32, modified Bessel integral representations.
   https://dlmf.nist.gov/10.32
   Context for the directly evaluated target and prior attempts. No Bessel transformation from this source is being claimed to prove the candidate.

