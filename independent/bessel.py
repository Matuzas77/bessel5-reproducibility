"""Fresh quadrature of A0 = integral_0^infinity K0(t)^5 dt.

All the work was done by AI.

Only the defining Bessel integral and mpmath's special functions enter this
calculation. No archived Bessel values, certificates, modular constants, or
conjectured right-hand side are read. Results are arbitrary-precision numerical
approximations, not interval-certified enclosures or proofs of the conjecture.
"""
from __future__ import annotations

import time
import mpmath


def compute_bessel(dps: int, variant: str = "primary") -> dict:
    """Compute the fifth K0 moment without using the conjectured answer.

    ``primary`` maps 0<t<1 to t=exp(-x), truncating at x=X.
    ``alternate`` maps 0<t<1 to t=u**16, 0<u<1, and changes the upper
    integration partitions. Both use mpmath's tanh-sinh quadrature.

    The returned value has extra display digits for convergence comparisons;
    the requested precision is not an assertion of certified accuracy.
    Endpoint truncation bounds are analytic inequalities evaluated in ordinary
    arbitrary-precision floating point. The quadrature error is heuristic.
    """
    if isinstance(dps, bool) or not isinstance(dps, int) or not 20 <= dps <= 200:
        raise ValueError("dps must be an integer from 20 through 200")
    if variant not in {"primary", "alternate"}:
        raise ValueError("variant must be 'primary' or 'alternate'")

    started = time.perf_counter()
    # A private context avoids changing a notebook's global mpmath precision.
    ctx = mpmath.mp.clone()
    guard = 20
    ctx.dps = dps + guard
    stringify = lambda value: ctx.nstr(value, dps + 8)
    zero = ctx.mpf(0)
    one = ctx.mpf(1)
    low_error = zero

    if variant == "primary":
        X = (dps + 35) * ctx.log(10)
        low_nodes = [ctx.mpf(v) for v in [0, 1, 4, 10, 25, 60] if v < X] + [X]

        def low_integrand(x):
            t = ctx.exp(-x)
            return t * ctx.besselk(0, t) ** 5

        low, low_error = ctx.quad(low_integrand, low_nodes, error=True)
        # For 0<t<=1, z=t**2/4 and the convergent K0/I0 series give
        # K0(t) <= exp(z)*(log(2/t)-EulerGamma+z)
        #        <= exp(1/4)*(1-log(t)).
        # Integrating this fifth power over 0<t<exp(-X) yields:
        y = X + 1
        polynomial = y**5 + 5*y**4 + 20*y**3 + 60*y**2 + 120*y + 120
        lower_tail = ctx.exp(ctx.mpf(5)/4 - X) * polynomial
        low_settings = {
            "substitution": "t = exp(-x)",
            "x_partitions": [stringify(v) for v in low_nodes],
            "small_t_cutoff": stringify(ctx.exp(-X)),
            "omitted_lower_tail_bound_formula":
                "exp(5/4-X)*((X+1)^5+5*(X+1)^4+20*(X+1)^3+60*(X+1)^2+120*(X+1)+120)",
        }
    else:
        low_nodes = [zero, ctx.mpf("0.25"), ctx.mpf("0.5"), ctx.mpf("0.75"), one]

        def low_integrand(u):
            if not u:
                return zero  # The limiting value is zero.
            return 16 * u**15 * ctx.besselk(0, u**16) ** 5

        low, low_error = ctx.quad(low_integrand, low_nodes, error=True)
        lower_tail = zero  # No finite lower-end truncation in this variant.
        low_settings = {
            "substitution": "t = u^16",
            "u_partitions": [stringify(v) for v in low_nodes],
            "small_t_cutoff": "0",
            "omitted_lower_tail_bound_formula":
                "0 (no explicit lower-end truncation; quadrature error remains numerical)",
        }

    # The rigorous real-variable inequality K0(t)<=sqrt(pi/(2*t))*exp(-t)
    # follows from K0(t)=int_0^infinity exp(-t*cosh(u))du and
    # cosh(u)>=1+u^2/2. Consequently the omitted upper tail is at most
    # (pi/(2*T))^(5/2)*exp(-5*T)/5. We evaluate that expression numerically.
    tail_extra = 12 if variant == "primary" else 17
    T = (dps + tail_extra) * ctx.log(10) / 5
    candidates = [1, 2, 5, 10, 15, 20, 30, 45, 65, 90, 120]
    if variant == "alternate":
        candidates = [1, "1.5", 3, 7, 12, 18, 26, 38, 55, 75, 100, 125]
    upper_nodes = [ctx.mpf(v) for v in candidates if ctx.mpf(v) < T] + [T]
    high = zero
    high_error = zero
    segment_details = []
    for a, b in zip(upper_nodes[:-1], upper_nodes[1:]):
        # Scale away exp(-5*a). This permits lower local working precision
        # while preserving an absolute-error budget of approximately dps+20
        # digits. This is a numerical efficiency choice, NOT an error proof.
        local_digits = max(25, dps + guard - int(5*a/ctx.log(10)))
        with ctx.workdps(local_digits):
            normalization = ctx.exp(5*a)
            normalized, error = ctx.quad(
                lambda t: normalization * ctx.besselk(0, t)**5,
                [a, b], error=True,
            )
        scale = ctx.exp(-5*a)
        high += normalized * scale
        high_error += abs(error) * scale
        segment_details.append({
            "interval": [stringify(a), stringify(b)],
            "working_dps": local_digits,
        })

    upper_tail = (ctx.pi/(2*T))**(ctx.mpf(5)/2) * ctx.exp(-5*T) / 5
    total = low + high
    return {
        "value": stringify(total),
        "requested_dps": dps,
        "working_dps": dps + guard,
        "variant": variant,
        "backend": "mpmath " + mpmath.__version__,
        "method": "Direct K0 evaluation and tanh-sinh quadrature",
        "integral": "integral_0^infinity K0(t)^5 dt (no factor t)",
        "low_integral": stringify(low),
        "high_integral": stringify(high),
        "low_settings": low_settings,
        "upper_t_cutoff": stringify(T),
        "upper_segments": segment_details,
        "omitted_lower_tail_bound_estimate": stringify(lower_tail),
        "omitted_upper_tail_bound_estimate": stringify(upper_tail),
        "omitted_tail_upper_bound_estimate": stringify(lower_tail + upper_tail),
        "upper_tail_bound_formula": "(pi/(2*T))^(5/2)*exp(-5*T)/5",
        "quadrature_error_heuristic": stringify(abs(low_error) + high_error),
        "error_status":
            "Analytic tail inequalities evaluated in floating point; quadrature and roundoff are not rigorously enclosed. Check convergence across precision and substitutions.",
        "elapsed_seconds": round(time.perf_counter() - started, 6),
        "ai_disclosure": "All the work was done by AI.",
    }
