"""Standalone exact verifier for an upper bound on Tao's constant C_1d (Python stdlib only).

C_1d is the largest C with  max_t |(f*f)(t)| >= C (int f)^2  for all real f supported on
[-1/4, 1/4]; any f gives C_1d <= max_t |f*f| / (int f)^2.

For a step function with n equal cells of height v_i on [-1/4, 1/4] (cell width w = 1/(2n)),
f*f is continuous and piecewise linear with knots at t = -1/2 + k w, so the maximum of |f*f| is
attained at a knot, where (f*f)(-1/2 + (m+1) w) = w * sum_{i+j=m} v_i v_j. Hence

    max_t |f*f| / (int f)^2 = 2n * max_m |sum_{i+j=m} v_i v_j| / (sum_i v_i)^2.

Every value in the JSON is read as the exact binary rational of its IEEE double.

Usage: python3 verify_standalone.py construction.json [--bound DECIMAL] [--exact P/Q] [--sha256 HEX]

Input is checked strictly: a JSON object with the single key "values", a non-empty list of JSON
numbers (no strings, booleans, NaN or Infinity) with nonzero sum. Optional claims are checked and
the exit status is 1 if any fails: --bound (the printed decimal must be >= the exact value, i.e. a
true upper bound), --exact (must equal the exact value), --sha256 (hash of the file's bytes).
"""

import argparse
import decimal
import hashlib
import json
import math
import sys
from fractions import Fraction

ap = argparse.ArgumentParser()
ap.add_argument("construction")
ap.add_argument("--bound", help="claimed upper bound as a decimal string")
ap.add_argument("--exact", help="claimed exact value P/Q")
ap.add_argument("--sha256", help="claimed SHA-256 of the construction file")
args = ap.parse_args()


def refuse(msg):
    print("REFUSED:", msg)
    sys.exit(1)


TABLE_ENTRY = Fraction(
    38534420343892614563110610373646557184000,
    26432274393314964174848236787785831382961,
)  # SCiEL 2026 (800 cells), last merged table entry, = 1.4578548849219886...
OPEN_PR_372 = Fraction(
    21267647941510998263770007610303168512,
    14626277862557668614837367632277156521,
)  # open PR #372 (51200 cells), best known as of 2026-11-16, = 1.4540711000681048...

raw = open(args.construction, "rb").read()
if args.sha256 and hashlib.sha256(raw).hexdigest() != args.sha256.lower():
    refuse(f"sha256 {hashlib.sha256(raw).hexdigest()} != claimed {args.sha256}")


def _no_const(c):
    raise ValueError(f"non-finite JSON constant {c}")


try:
    doc = json.loads(raw, parse_constant=_no_const)
except ValueError as e:
    refuse(f"invalid JSON: {e}")
if not isinstance(doc, dict) or set(doc) != {"values"}:
    refuse("top level must be an object with exactly one key, 'values'")
vals = doc["values"]
if not isinstance(vals, list) or not vals:
    refuse("'values' must be a non-empty list")
if not all(type(x) in (int, float) and math.isfinite(x) for x in vals):
    refuse(
        "every height must be a finite JSON number (no strings, booleans, NaN, Infinity)"
    )
v = [Fraction(x) for x in vals]
n = len(v)
# Exact integers: scale every value by the common denominator D (a power of two), u_i = D v_i.
D = 1
for x in v:
    D = D * x.denominator // math.gcd(D, x.denominator)
u = [int(x * D) for x in v]
assert all(Fraction(a, D) == x for a, x in zip(u, v))
S = sum(u)
if S == 0:
    refuse("sum of heights is 0 (ratio undefined)")
# c_m = sum_{i+j=m} u_i u_j, exactly, by Kronecker substitution with Python big integers
# (pure-Python O(n^2) would take hours at 100k cells). Signed values: split u = P - N, then
# c = P*P + N*N - 2 P*N, each a product of polynomials with nonnegative coefficients. The slot
# width B exceeds the bit length of every coefficient of these products, so slots never carry.
P = [a if a > 0 else 0 for a in u]
N = [-a if a < 0 else 0 for a in u]
W = (
    2 * max(abs(a) for a in u).bit_length() + n.bit_length() + 16
) // 8  # bytes per slot


def pack(c):
    return int.from_bytes(b"".join(x.to_bytes(W, "little") for x in c), "little")


def unpack(z):
    b = z.to_bytes(W * (2 * n) + 8, "little")
    return [int.from_bytes(b[k * W : (k + 1) * W], "little") for k in range(2 * n - 1)]


pP, pN = pack(P), pack(N)
c = [
    a + b - 2 * d for a, b, d in zip(unpack(pP * pP), unpack(pN * pN), unpack(pP * pN))
]
best = max(abs(x) for x in c)
if n <= 3000:  # small n: also check against the direct O(n^2) sum
    assert c == [
        sum(u[i] * u[m - i] for i in range(max(0, m - n + 1), min(m, n - 1) + 1))
        for m in range(2 * n - 1)
    ]
# 2n * max|sum v_i v_j| / (sum v)^2 = 2n * max|c| / S^2  (the D^2 factors cancel)
r = Fraction(2 * n * best, S * S)
print("n =", n)
print("C_1d <=", f"{r.numerator}/{r.denominator}")
print("       =", f"{float(r):.16f}")
# 17-significant-digit decimal rounded toward +infinity: safe to print as an upper bound
with decimal.localcontext() as ctx:
    ctx.prec, ctx.rounding = 17, decimal.ROUND_CEILING
    up = decimal.Decimal(r.numerator) / decimal.Decimal(r.denominator)
print("print as (rounded up):", up)
print(
    "below merged table entry 1.4578548849219886:",
    r < TABLE_ENTRY,
    f"(difference {float(TABLE_ENTRY - r):+.6e})",
)
print(
    "below open PR #372      1.4540711000681048:",
    r < OPEN_PR_372,
    f"(difference {float(OPEN_PR_372 - r):+.6e})",
)
if args.exact is not None:
    try:
        claimed = Fraction(args.exact)
    except (ValueError, ZeroDivisionError):
        refuse(f"--exact {args.exact!r} is not a fraction")
    if claimed != r:
        refuse(f"exact value {r} != claimed {claimed}")
if args.bound is not None:
    try:
        claimed = Fraction(args.bound)
    except ValueError:
        refuse(f"--bound {args.bound!r} is not a decimal")
    if claimed < r:
        refuse(
            f"claimed bound {args.bound} is below the exact value by {float(r - claimed):.3e}"
        )
    print(
        f"claimed bound {args.bound} holds (exceeds the exact value by {float(claimed - r):.3e})"
    )
print("ACCEPTED")
