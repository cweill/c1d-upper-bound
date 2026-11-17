"""Independent exact check of a 1d construction: R = 2n*max_m|sum_{i+j=m} v_i v_j|/(sum v)^2.
No code shared with constants/1d/verify.py: exact float->int scaling, big-int Kronecker substitution
with signed split (P-N)^2 = P^2 - 2PN + N^2, plus a float FFT sanity check."""

import json, sys, math
from fractions import Fraction
import numpy as np

v = json.load(open(sys.argv[1]))["values"]
n = len(v)
assert all(isinstance(x, (int, float)) and math.isfinite(x) for x in v)
# exact: every float is m*2^e; scale all to integers with a common power of two
fr = [Fraction(x) for x in v]
den = 1
for f in fr:
    den = den * f.denominator // math.gcd(den, f.denominator)
ints = [int(f * den) for f in fr]
assert all(Fraction(a, den) == f for a, f in zip(ints, fr))
P = [a if a > 0 else 0 for a in ints]
N = [-a if a < 0 else 0 for a in ints]
B = (
    2 * max(map(abs, ints)).bit_length() + n.bit_length() + 2
)  # slot width > any coefficient


def pack(c):
    return int.from_bytes(
        b"".join(x.to_bytes((B + 7) // 8, "little") for x in c), "little"
    )


def unpack(z, m):
    w = (B + 7) // 8
    bs = z.to_bytes(w * m + 8, "little")
    return [int.from_bytes(bs[i * w : (i + 1) * w], "little") for i in range(m)]


pP, pN = pack(P), pack(N)
m = 2 * n - 1
pp, nn, pn = unpack(pP * pP, m), unpack(pN * pN, m), unpack(pP * pN, m)
conv = [a + b - 2 * c for a, b, c in zip(pp, nn, pn)]
S = sum(ints)
M = max(abs(c) for c in conv)
R = Fraction(2 * n * M, S * S)
# float FFT sanity
a = np.array(v, float) / max(map(abs, v))
L = 1 << (2 * n).bit_length()
cf = np.fft.irfft(np.fft.rfft(a, L) ** 2, L)[:m]
Rf = 2 * n * np.max(np.abs(cf)) / a.sum() ** 2
print(
    json.dumps(
        {
            "n": n,
            "neg_cells": sum(x < 0 for x in ints),
            "zero_cells": ints.count(0),
            "argmax_knot": max(range(m), key=lambda i: abs(conv[i])),
            "R_float": float(R),
            "R_exact": f"{R.numerator}/{R.denominator}",
            "R_fft": float(Rf),
        }
    )
)
