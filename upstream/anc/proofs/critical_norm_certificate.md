# Finite certificate for the modular critical-norm bridge

The theorem is proved by the analytic kernel argument in critical_norm_bridge.md plus these exact modular identities. The scripts do not formalize that analytic argument. The original Bessel/eta equality remains unproved.

## Normalizations

- Weight 3; level 60; character chi_-15; Petersson integral not divided by index144.
- q=exp(2 pi i tau); the original embedding has a3=sqrt5-2i, a7=-8i, a11=4i sqrt5.
- The critical Poincare kernel is analytically continued to sigma=0. It is not a raw absolutely convergent critical m-sum.
- The twisted symmetric-square L-function has conductor960 and the stated local factors at2,3,5.

## Eta cusp forms

The columns in the following exponents are the divisors

`1,2,3,4,5,6,10,12,15,20,30,60`.

```json
[
  [
    3,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3,
    0,
    0,
    0
  ],
  [
    2,
    0,
    0,
    1,
    0,
    0,
    0,
    0,
    2,
    0,
    0,
    1
  ],
  [
    2,
    0,
    0,
    0,
    0,
    1,
    1,
    0,
    2,
    0,
    0,
    0
  ],
  [
    1,
    0,
    0,
    2,
    0,
    0,
    0,
    0,
    1,
    0,
    0,
    2
  ],
  [
    1,
    0,
    0,
    1,
    0,
    1,
    1,
    0,
    1,
    0,
    0,
    1
  ],
  [
    1,
    0,
    0,
    0,
    0,
    2,
    2,
    0,
    1,
    0,
    0,
    0
  ],
  [
    0,
    2,
    0,
    0,
    0,
    0,
    0,
    1,
    0,
    1,
    2,
    0
  ],
  [
    0,
    1,
    0,
    0,
    0,
    0,
    0,
    2,
    0,
    2,
    1,
    0
  ],
  [
    0,
    0,
    3,
    0,
    3,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  ],
  [
    0,
    0,
    2,
    0,
    2,
    0,
    0,
    1,
    0,
    1,
    0,
    0
  ],
  [
    0,
    0,
    0,
    3,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3
  ],
  [
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    3,
    0,
    3,
    0,
    0
  ],
  [
    -1,
    3,
    3,
    0,
    0,
    0,
    1,
    0,
    0,
    0,
    0,
    0
  ],
  [
    -1,
    2,
    3,
    0,
    0,
    1,
    0,
    0,
    0,
    0,
    1,
    0
  ],
  [
    -1,
    1,
    4,
    0,
    1,
    1,
    0,
    0,
    0,
    0,
    0,
    0
  ],
  [
    -1,
    1,
    3,
    0,
    0,
    1,
    0,
    1,
    0,
    1,
    0,
    0
  ],
  [
    -1,
    0,
    3,
    1,
    0,
    0,
    0,
    0,
    0,
    3,
    0,
    0
  ],
  [
    2,
    -1,
    0,
    0,
    0,
    2,
    0,
    0,
    2,
    0,
    1,
    0
  ]
]
```

All22 quotients (these18 plus the four P products) have the correct character, both congruences modulo24 and strictly positive order at every cusp. Both implementations check these conditions.

## Explicit critical-kernel identity

In the above basis the exact coefficient vector is

```json
[
  "0",
  "-15/8",
  "269/128",
  "-5/2",
  "127/32",
  "-235/64",
  "-4",
  "35/8",
  "29/128",
  "11/16",
  "5/4",
  "1/6",
  "37/64",
  "15/64",
  "33/128",
  "-3/16",
  "13/48",
  "87/128"
]
```

It agrees with the finite class-number formula through q^200 in the primary calculation and q^240 in the independent calculation. The equality theorem uses the analytically established modularity and the Sturm bound36, not merely the longer prefix.

## Hecke witness

`K - H_f = (T7^2+64) R`. The four-product coefficients of H_f are

```json
["49/144", "-5/48", "5/16", "1/16"]
```

The18 coefficients of R are

```json
[
  "1/192",
  "-43/55296",
  "17/27648",
  "-619/27648",
  "-871/27648",
  "-11/384",
  "-269/3456",
  "1/6",
  "-1/6144",
  "-1/3072",
  "0",
  "0",
  "49/2880",
  "-13/8640",
  "1/72",
  "-1/160",
  "0",
  "0"
]
```

This identity is checked directly through the Sturm bound using eta coefficients through q^1764. The source data for it are the Euler products, not a stored period or approximate Hecke matrix.

## Spectral coefficient and exact scalar

The q-coefficients of U3^j H_f for j=0,1,2,3 are

```json
["49/144", "-19/144", "-91/144", "-119/144"]
```

They give
\[
 \delta_f=\frac{(7-i\sqrt5)(7-i)}{576}.
\]
The independent program verifies the full Galois decomposition of H_f in the four normalized eigenforms. Their good T7/T11 eigenpairs are distinct, so Petersson orthogonality applies.

With
\[
 R(3)=\frac{16(7+i\sqrt5)(7+i)}{675},
 \qquad \zeta^{(30)}(2)=\frac{8\pi^2}{75},
\]
the exact scalar check is
\[
 \frac{128\bar\delta_f}{75R(3)}=\frac18.
\]
This is the factor in `L(chi_-4 Sym^2 f,3)=pi^4/8*<f,f>_60`.
The conductor-960 functional equation then gives
`L(chi_-4 Sym^2 f,2)=sqrt15*pi^2*<f,f>_60`.

## Critical-limit audit

A compact-support Fourier transform alone does not justify dropping the infinite hyperbolic sum. The full proof `critical_norm_bridge.md` first proves that every arithmetic Dirichlet factor is regular at sigma=0. The potentially principal discriminants are Delta=60u^2 and Delta=0. In those cases the complete local factor at2 vanishes at zero: it is the integral of chi_-4(r) times the character of the odd part of a norm `v^2-15u^2`. That odd part is1mod4, so the average iszero. It cancels the zeta pole before the Fourier cutoff is used. Uniform summable estimates then justify the limit. Omitting this argument would leave a gap.

The generic Gauss/CRT formulas are also checked by20 direct cyclotomic-ring calculations; no numerical roots of unity are used. Their full records are in `../outputs/critical_norm_primary.json`.

## Scope

No numerical Bessel or eta period, numerical L-value, PSLQ or numerical rank is used. No original Bessel-to-L theorem or Domb-to-period identity is proved. Neither an integral motive correspondence nor a new source projector is claimed.
