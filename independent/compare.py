"""Compare freshly computed decimal strings without lowering precision.

All the work was done by AI. These are empirical tests, not interval proofs.
"""
import mpmath


def assess_checks(runs):
    """Apply thresholds fixed in advance: requested digits minus ten.

    Each run compares separate computations. Adjacent runs test sensitivity to
    precision and coordinates. Negative controls use an incorrect phase and an
    omitted Mellin normalization; a perturbed result tests the comparison itself.
    """
    if len(runs) < 2:
        raise ValueError("At least two precision/coordinate runs are required")
    precisions = [run["dps"] for run in runs]
    if precisions != sorted(set(precisions)) or min(precisions) < 30:
        raise ValueError("Use increasing, distinct precisions >=30")
    mp = mpmath.mp.clone()
    mp.dps = max(run["dps"] for run in runs) + 30
    rows = []
    previous = None
    for run in runs:
        dps = run["dps"]
        a = mp.mpf(run["bessel"]["value"])
        rhs = mp.mpf(run["modular"]["rhs"])
        tolerance = mp.mpf(10) ** (-(dps - 10))
        error = abs(a-rhs)/abs(a)
        row = {"dps": dps, "relative_discrepancy": mp.nstr(error, 12),
               "threshold": mp.nstr(tolerance, 8), "passed": bool(error < tolerance)}
        if previous is not None:
            low, old_a, old_rhs = previous
            convergence_tolerance = mp.mpf(10) ** (-(low - 10))
            da, dr = abs(a-old_a)/abs(a), abs(rhs-old_rhs)/abs(rhs)
            row["convergence"] = {
                "previous_dps": low, "bessel_relative_change": mp.nstr(da, 12),
                "rhs_relative_change": mp.nstr(dr, 12),
                "threshold": mp.nstr(convergence_tolerance, 8),
                "passed": bool(max(da, dr) < convergence_tolerance)}
        rows.append(row)
        previous = (dps, a, rhs)
    final = runs[-1]
    L = mp.mpc(final["modular"]["real"], final["modular"]["imag"])
    r = mp.sqrt(5)
    def expression(z):
        return mp.pi/4*((125-13*r)*mp.re(z**2)+(41*r-125)*mp.im(z**2))
    controls = {
        "conjugate_complex_embedding": expression(mp.conj(L)),
        "omit_mellin_4pi_squared": expression(L/(4*mp.pi**2)),
        "replace_L_squared_by_abs_L_squared": mp.pi/4*(125-13*r)*abs(L)**2,
    }
    negative = {}
    for name, value in controls.items():
        error = abs(a-value)/abs(a)
        negative[name] = {"value": mp.nstr(value, 24),
                          "relative_discrepancy": mp.nstr(error, 12),
                          "rejected": bool(error > mp.mpf("0.001"))}
    # Meaningful sensitivity test: injecting a 1e-12 relative change must fail.
    perturbation_detected = bool(abs(a*(1+mp.mpf("1e-12"))-rhs)/abs(a) > tolerance)
    passed = (all(x["passed"] and x.get("convergence", {}).get("passed", True) for x in rows)
              and all(x["rejected"] for x in negative.values()) and perturbation_detected)
    return {"passed": passed, "comparisons": rows, "negative_controls": negative,
            "relative_1e_minus_12_perturbation_detected": perturbation_detected,
            "rigorous_enclosure": False,
            "interpretation": "Empirical numerical agreement; not a proof or human peer review.",
            "disclosure": "All the work was done by AI."}
