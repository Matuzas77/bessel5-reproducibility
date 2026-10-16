# Direct Bessel quadrature and a separate modular calculation

**All the work was done by AI.** Jonas Matuzas is responsible for the claims. This is new AI-executed numerical evidence, not independent human peer review or a proof of the conjecture.

[Run in Colab](https://colab.research.google.com/github/Matuzas77/bessel5-reproducibility/blob/main/notebooks/Bessel5_independent_check.ipynb) · [Executed notebook](../notebooks/Bessel5_independent_check.ipynb) · [Actual 40/60/80-digit report](../results/2026-10-16-independent/run.json)

## Independence and normalization

The two new modules use only the mathematical definitions in equations (1.1), (1.2) and (2.1) of the [v4 manuscript](../upstream/main.tex). They do not import upstream verification code, saved numerical values, coefficient tables, interval certificates or the opposite side's result. `bessel.py` computes the unweighted integral with measure **dt**. `lvalue.py` computes the classical value **L(f,2)=4π²∫y f(iy)dy**, then forms the proposed expression from the **complex square L²**. The numerical paths meet only in `compare.py`.

Both modules use mpmath 1.3.0 arithmetic and tanh-sinh quadrature. They are alternative routes to the published certificate pipeline, not independently developed arithmetic libraries. The additional reviewer is also AI. Its [review record](../results/2026-10-16-independent/review.json) includes a separate literal Euler-product comparison with the eta implementation and checks of the measure, phase, normalization and endpoint inequalities.

## Bessel calculation

Split at t=1. The primary run uses t=exp(-x) on the lower part and truncates at X=(p+35)log(10). The alternate uses t=u^16 with Jacobian 16u^15, avoiding a finite lower cutoff. The upper part directly quadratures K0(t)^5 with different partitions and cutoffs in the two variants.

The integral representation K0(t)=∫exp(-t cosh u)du and cosh u≥1+u²/2 give K0(t)≤sqrt(π/(2t))exp(-t). Consequently the upper omitted tail is at most `(π/(2T))^(5/2) exp(-5T)/5`. The convergent small-t series gives K0(t)≤exp(1/4)(1−log t) for 0<t≤1. Integrating its fifth power yields the explicit lower-tail polynomial in the code. These analytic inequalities are evaluated in ordinary floating point; their displayed values do not certify rounding or interior quadrature error. [DLMF integral representation](https://dlmf.nist.gov/10.32.E9), [DLMF small-argument series](https://dlmf.nist.gov/10.31.E2).

Arithmetic begins with p+20 digits. Upper segments [a,b] are scaled by exp(5a) and use max(25,p+20−floor(5a/log(10))) working digits, then are rescaled at full precision. This retains an intended absolute guard budget and reduces runtime. It is a numerical strategy, not a rounding-error proof. The JSON records every segment's precision and cutoff.

## Modular calculation

Dedekind eta on the imaginary axis is computed from its Euler pentagonal expansion. Its S transformation maps x<1 to 1/x, keeping q≤exp(-2π). No precomputed Fourier table is used. Fresh formal integer Euler products through q^7 verify the specified complex embedding before integration. [DLMF eta transformation](https://dlmf.nist.gov/23.18).

The primary Mellin integral uses a partition of the full positive y-axis. The alternate uses u=log(y) with different partitions, over −12≤u≤8. Its omitted tails decay exponentially beyond the requested precision, but the program does not supply a rigorous global enclosure. The saved eta-series truncation bound applies to that series only, not the whole integral.

## Actual result and predeclared checks

The 16 October 2026 run computed:

```text
A0 ≈ 135.26830258086883759422627964619220742030588935942352678469351371045888711773849
```

| Requested precision | Bessel / modular route | Relative cross-discrepancy | Required threshold |
|---|---|---:|---:|
|40|logarithmic / direct y|4.35e-41|1e-30|
|60|power substitution / log y|3.59e-60|1e-50|
|80|logarithmic / direct y|9.73e-81|1e-70|

Discrepancies include decimal serialization and are **not upper bounds on the true mathematical error**. Precision/coordinate convergence passes at the lower run's threshold. Three deliberate errors are detected: conjugating the embedding changes the answer by about 17.44%; omitting the 4π² Mellin factor by 99.94%; replacing L² by |L|² by 5.33%. A synthetic relative 1e-12 perturbation also fails the correct comparison tolerance. These thresholds are set in code before reading values.

The [clean-kernel notebook report](../results/2026-10-16-independent/notebook-run.json) repeats the default 40/60-digit calculation from an empty working directory. The notebook embeds the three module implementations; its source hashes and parity are checked by `tools/check_notebook.py`. Saved outputs are evidence, never calculation inputs.

## Run locally or regenerate the notebook

After installing `requirements.lock` as in the main README:

```sh
.venv/bin/python tools/run_independent.py --digits 40,60,80 --output results/current/direct-run.json
.venv/bin/python -m pip install -r requirements-notebook.txt
.venv/bin/python tools/check_notebook.py --execute --output results/current/executed.ipynb --report results/current/notebook-run.json
```

Use a new output path to preserve previous standalone results. Browser execution uses the host's Python runtime, recorded in its output, with pinned mpmath. Local and CI runs pin Python 3.12.12. Colab may require the user's normal Google sign-in and imposes its own runtime limits.

If the module code changes, regenerate with `tools/build_notebook.py`, execute it again, and retain fresh evidence. Regeneration clears saved outputs. The numerical environment is pinned with wheel hashes; notebook-authoring tools are pinned by version in `requirements-notebook.txt`.

The outstanding mathematical task is an analytic proof of the Bessel-to-modular identity. Increasing precision alone cannot supply that proof.
