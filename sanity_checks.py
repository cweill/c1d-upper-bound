"""Sanity checks for a C_1d step-function construction that do NOT rely on the knot formula.

1. Continuous evaluation: (f*f)(t) computed exactly in rationals at arbitrary t (including
   random off-knot points) from the overlap-length formula
       (f*f)(t) = sum_{i,j} v_i v_j * |[a_i, a_i+w] ∩ [t - a_j - w, t - a_j]|,
   confirming that no off-knot t exceeds the knot maximum (i.e. f*f really is piecewise
   linear with the claimed knots, and the reported max is the true sup).
2. The sup over t is attained: report the maximizing knot(s) and check |f*f| there equals the
   claimed ratio * (int f)^2 exactly.
3. Definition checks: support exactly [-1/4, 1/4]; sign convention (|.| inside the max);
   compare with the |max f*f| variant (a different constant, e.g. EinsteinArena problem 4).

Usage: python3 sanity_checks.py construction.json [n_random_t]
"""

import json
import random
import sys
from fractions import Fraction

vals = json.load(open(sys.argv[1]))["values"]
n_rand = int(sys.argv[2]) if len(sys.argv) > 2 else 300
v = [Fraction(x) for x in vals]
n = len(v)
w = Fraction(1, 2 * n)
a = [Fraction(-1, 4) + i * w for i in range(n)]  # left edges of cells
S = sum(v)
I = S * w  # int f


def overlap(lo1, hi1, lo2, hi2):
    return max(Fraction(0), min(hi1, hi2) - max(lo1, lo2))


# exact knot values c_m = sum_{i+j=m} v_i v_j, (f*f)(-1/2 + (m+1) w) = w c_m.
# Large n: exact integer polynomial squaring via FLINT (pure-Python O(n^2) takes hours at 50k+).
_den = max(x.denominator for x in v)
_ints = [x.numerator * (_den // x.denominator) for x in v]
try:
    import flint

    c = [Fraction(int(x), _den * _den) for x in (flint.fmpz_poly(_ints) ** 2).coeffs()]
    c += [Fraction(0)] * (2 * n - 1 - len(c))
except ImportError:
    if n > 4000:
        sys.exit("python-flint is required for n > 4000 (pip install python-flint)")
    c = [
        sum(v[i] * v[m - i] for i in range(max(0, m - n + 1), min(m, n - 1) + 1))
        for m in range(2 * n - 1)
    ]
knot_max = max(abs(x) for x in c) * w
ratio = knot_max / (I * I)
argmax = [m for m, x in enumerate(c) if abs(x) * w == knot_max]


def ff(t):
    """Exact (f*f)(t) = sum_j v_j * int_{a_j}^{a_j+w} f(t - y) dy, via cell overlaps."""
    total = Fraction(0)
    # f(t - y) for y in cell j is f on [t - a_j - w, t - a_j]; only cells i overlapping matter
    for j in range(n):
        lo, hi = t - a[j] - w, t - a[j]
        i0 = max(0, int((lo - a[0]) / w) - 1)
        i1 = min(n - 1, int((hi - a[0]) / w) + 1)
        for i in range(i0, i1 + 1):
            ov = overlap(a[i], a[i] + w, lo, hi)
            if ov:
                total += v[i] * v[j] * ov
    return total


print(f"n = {n}, ratio (knot formula) = {float(ratio):.16f}")
print(f"support: [{float(a[0])}, {float(a[-1] + w)}]  (must be [-0.25, 0.25])")
print(
    f"max |f*f| attained at {len(argmax)} knot(s); first at t = {float(Fraction(-1, 2) + (argmax[0] + 1) * w):.6f}"
)

# 1. exact continuous evaluation at the argmax knot, at neighbours/midpoints, and at random t
t_star = Fraction(-1, 2) + (argmax[0] + 1) * w
assert abs(ff(t_star)) == knot_max, (
    "continuous formula disagrees with knot formula at argmax"
)
rng = random.Random(0)
worst = Fraction(0)
checked = 0
test_ts = [t_star + w / 2, t_star - w / 2, t_star + w / 3]
test_ts += [Fraction(rng.randrange(-(10**9), 10**9), 2 * 10**9) for _ in range(n_rand)]
for t in test_ts:
    val = abs(ff(t))
    worst = max(worst, val)
    checked += 1
    assert val <= knot_max, (
        f"off-knot value {float(val)} exceeds knot max at t={float(t)}"
    )
print(
    f"continuous check: {checked} off-knot t values, max |f*f| = {float(worst / (I * I)):.16f} <= knot max  OK"
)

# 2./3. the other variant, for contrast
var = abs(max(c)) * w / (I * I)
print(
    f"|max f*f| / (int f)^2 = {float(var):.16f}   (different constant; abs outside the max)"
)
print(
    f"max f*f attained on the {'positive' if max(c) >= -min(c) else 'negative'} side; "
    f"min f*f / max|f*f| = {float(min(c) / max(abs(x) for x in c)):.4f}"
)
