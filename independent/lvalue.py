"""Fresh numerical Mellin integration for the Bessel 5 modular expression.

All the work was done by AI.

Only the eta-quotient mathematical definition in upstream/main.tex is used.
No archived numerical values, upstream routines, or coefficient tables are
imported. This is numerical convergence evidence, not a rigorous enclosure or
a proof of the Bessel conjecture.
"""

from __future__ import annotations

import json
import math
import time

import mpmath


def _eta_imaginary(mp, x, terms: int):
    """Dedekind eta(i*x), x>0, using its S transformation and pentagonal sum.

    eta(i*x) = exp(-pi*x/12) prod_{n>=1}(1-exp(-2*pi*n*x)).
    Euler's pentagonal theorem evaluates the product as
    1 + sum_{n>=1} (-1)^n (q^(n(3n-1)/2)+q^(n(3n+1)/2)).
    Modular transformation keeps q <= exp(-2*pi); all square roots are positive.
    """
    if x <= 0:
        raise ValueError("eta argument must be positive")
    z = x if x >= 1 else 1 / x
    q = mp.exp(-2 * mp.pi * z)
    product = mp.mpf(1)
    for n in range(1, terms + 1):
        pair = q ** (n * (3 * n - 1) // 2) + q ** (n * (3 * n + 1) // 2)
        product += pair if n % 2 == 0 else -pair
    eta = mp.exp(-mp.pi * z / 12) * product
    return eta if x >= 1 else eta / mp.sqrt(x)


def _embedding_check(mp, r, omega, epsilon):
    """Fresh integer Euler-product coefficients through q^7, not a data table."""
    degree = 7
    definitions = (
        (1, {2: 1, 6: 2, 10: 4, 30: -1}),
        (3, {2: 4, 6: -1, 10: 1, 30: 2}),
        (3, {2: -1, 6: 4, 10: 2, 30: 1}),
        (5, {2: 2, 6: 1, 10: -1, 30: 4}),
    )
    rows = []
    for shift, powers in definitions:
        row = [1] + [0] * degree
        for d, exponent in powers.items():
            for step in range(d, degree + 1, d):
                if exponent == -1:
                    # Multiplication by the exact formal geometric series.
                    for k in range(step, degree + 1):
                        row[k] += row[k-step]
                else:
                    for _ in range(exponent):
                        for k in range(degree, step-1, -1):
                            row[k] -= row[k-step]
        rows.append([0]*shift + row[:degree+1-shift])
    weights = (mp.mpc(1), omega, 3*epsilon, 3*epsilon*mp.conj(omega))
    expected = {1: mp.mpc(1), 3: r-2j, 5: r*(1+2j), 7: mp.mpc(0, -8)}
    errors = {}
    values = {}
    for n, target in expected.items():
        value = sum(weight*row[n] for weight, row in zip(weights, rows))
        error = abs(value-target)
        errors[str(n)] = mp.nstr(error, 8)
        values[str(n)] = {"real": mp.nstr(mp.re(value), 18), "imag": mp.nstr(mp.im(value), 18)}
        if error > mp.mpf(10)**(-mp.dps+5):
            raise ArithmeticError(f"Unexpected q^{n} coefficient / complex embedding")
    return {"passed": True, "method": "Fresh formal integer Euler products through q^7",
            "coefficients": values, "absolute_errors": errors}


def compute_lvalue(dps: int, variant: str = "primary") -> dict:
    """Return L(f,2) and the proposed Bessel expression as decimal strings.

    primary: integrate y*f(iy) directly on a partition of (0,infinity).
    shifted: integrate after y=exp(u), so the integrand is exp(2u)*f(i*exp(u)),
             on a differently partitioned finite interval, with rapidly
             decaying numerical tails truncated outside it.

    Values have dps significant digits; arithmetic uses 20 guard digits.
    Neither mpmath's adaptive integration nor the truncated tails supplies
    a rigorous interval certificate. Compare independent runs/partitions.
    """
    if not isinstance(dps, int) or not 20 <= dps <= 300:
        raise ValueError("dps must be an integer between 20 and 300")
    if variant not in ("primary", "shifted"):
        raise ValueError("variant must be 'primary' or 'shifted'")
    start = time.perf_counter()
    mp = mpmath.mp.clone()
    mp.dps = dps + 20
    terms = math.ceil(math.sqrt((mp.dps + 10) * math.log(10) / (3 * math.pi))) + 2
    r = mp.sqrt(5)
    omega = (1 - 2j) / r
    epsilon = (5 + 4 * r + 1j * (2 * r - 10)) / 15
    embedding_check = _embedding_check(mp, r, omega, epsilon)
    calls = 0

    def f(y):
        nonlocal calls
        calls += 1
        e2, e6, e10, e30 = [_eta_imaginary(mp, d * y, terms) for d in (2, 6, 10, 30)]
        p1 = e2 * e6**2 * e10**4 / e30
        p2 = e2**4 * e10 * e30**2 / e6
        p3 = e6**4 * e10**2 * e30 / e2
        p4 = e2**2 * e6 * e30**4 / e10
        return p1 + omega * p2 + 3 * epsilon * p3 + 3 * epsilon * mp.conj(omega) * p4

    if variant == "primary":
        partition = [mp.mpf(0), mp.mpf(1)/60, mp.mpf(1)/30, mp.mpf(1)/10,
                     mp.mpf(1)/6, mp.mpf(1)/2, mp.mpf(1), mp.mpf(3), mp.inf]
        value = 4 * mp.pi**2 * mp.quad(lambda y: y * f(y), partition)
        variable = "y; full positive real axis"
        tails = "No explicit tail cutoff; mpmath handles infinite endpoints numerically."
    else:
        # The slowest small-y eta quotient has exponential exp(-pi/(30*y))
        # times a power y^-3. At y=exp(-12), its exponential is below exp(-17000).
        # At y=exp(8), the large-y leading q term decays below exp(-18000).
        partition = [mp.mpf(k) for k in (-12, -6, -4, -3, -2, -1, 0, 1, 3, 8)]
        value = 4 * mp.pi**2 * mp.quad(lambda u: mp.exp(2*u) * f(mp.exp(u)), partition)
        variable = "u=log(y); y in [exp(-12),exp(8)]"
        tails = ("Numerical tails omitted outside [exp(-12),exp(8)]; eta asymptotic "
                 "decay is far below requested precision, but no interval bound is claimed.")
    square = value**2
    rhs = mp.pi / 4 * ((125 - 13*r) * mp.re(square) + (41*r - 125) * mp.im(square))
    # A conservative magnitude bound for the omitted pentagonal terms at z>=1:
    # 2*q^a/(1-q), a=(N+1)(3(N+1)-1)/2, since exponents increase by >=1.
    n = terms + 1
    qmax = mp.exp(-2 * mp.pi)
    omitted_product_bound = 2*qmax**(n*(3*n-1)//2)/(1-qmax)
    return {
        "real": mp.nstr(mp.re(value), dps),
        "imag": mp.nstr(mp.im(value), dps),
        "rhs": mp.nstr(rhs, dps),
        "settings": {
            "dps": dps, "working_dps": mp.dps, "variant": variant,
            "quadrature": "mpmath tanh-sinh", "variable": variable,
            "partition": [mp.nstr(x, 18) for x in partition],
            "eta_method": "S transformation and fresh Euler pentagonal sum",
            "eta_terms_each_side": terms,
            "eta_product_truncation_bound": mp.nstr(omitted_product_bound, 8),
            "tails": tails,
            "integrand_evaluations": calls,
            "mpmath_version": mpmath.__version__,
            "rigorous_enclosure": False,
            "source": "Mathematical eta-quotient definition in upstream/main.tex only",
            "embedding_check": embedding_check,
        },
        "seconds": round(time.perf_counter()-start, 6),
        "disclosure": "All the work was done by AI.",
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--variant", choices=("primary", "shifted"), default="primary")
    args = parser.parse_args()
    print(json.dumps(compute_lvalue(args.dps, args.variant), indent=2))
