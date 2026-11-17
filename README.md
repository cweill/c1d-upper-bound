# Upper bound C_1d ≤ 1.4513064006879361 (signed third autocorrelation constant)

A certificate for constant **1d** in Terence Tao's
[optimization constants table](https://teorth.github.io/optimizationproblems/constants/1d.html):
C_1d is the largest C with max_t |(f*f)(t)| ≥ C (∫f)² for every real f supported on [-1/4, 1/4],
so every such f gives an upper bound.

The file `c1d_n102400_1.4513064007.json` holds a step function with **n = 102400** equal cells on
[-1/4, 1/4]. Its key `values` lists the signed cell heights. Each height is an IEEE binary64
number, read as its exact binary rational.

- Exact value: 2n · max_m |Σ_{i+j=m} v_i v_j| / (Σ v_i)² =
  `513655524515252183593642627316297728000 / 353926313748615382700541734286037793169`
  = 1.4513064006879361…
- SHA-256 of the JSON: `391cc726acbab3a7f4216dedab2ec5dd74cb797dc12189966dac11a0dee31e2d`

## Verify

```
python3 verify_standalone.py c1d_n102400_1.4513064007.json \
    --bound 1.4513064006879361 \
    --exact 513655524515252183593642627316297728000/353926313748615382700541734286037793169 \
    --sha256 391cc726acbab3a7f4216dedab2ec5dd74cb797dc12189966dac11a0dee31e2d   # stdlib only, ~10 s
python3 controls.py c1d_n102400_1.4513064007.json 1.4513064006879361 \
    513655524515252183593642627316297728000/353926313748615382700541734286037793169 \
    391cc726acbab3a7f4216dedab2ec5dd74cb797dc12189966dac11a0dee31e2d       # stdlib only, ~2 min
pip install python-flint numpy
python3 sanity_checks.py c1d_n102400_1.4513064007.json 6               # exact off-knot checks
python3 independent_check.py c1d_n102400_1.4513064007.json             # second exact implementation
```

- `verify_standalone.py` uses only the Python standard library. f*f is piecewise linear, with
  knots on the cell grid, so its maximum is attained at a knot. The script computes every knot
  value exactly: it scales the heights to integers, then squares the polynomial by Kronecker
  substitution with big integers. For n ≤ 3000 it also checks the result against the direct sum.
- `sanity_checks.py` evaluates (f*f)(t) exactly from cell overlaps at off-knot points and at the
  maximizing knot, without using the knot formula. It also checks that the support is exactly
  [-1/4, 1/4] and that the absolute value is taken inside the max.
- `independent_check.py` is a separate exact implementation. It splits the heights into
  positive and negative parts before the big-integer substitution, and adds a float FFT check.

Output of the three checks:

```
n = 102400
C_1d <= 513655524515252183593642627316297728000/353926313748615382700541734286037793169
       = 1.4513064006879361
print as (rounded up): 1.4513064006879361
below merged table entry 1.4578548849219886: True (difference +6.548484e-03)
below open PR #372      1.4540711000681048: True (difference +2.764699e-03)
claimed bound 1.4513064006879361 holds (exceeds the exact value by 1.203e-18)
ACCEPTED

n = 102400, ratio (knot formula) = 1.4513064006879361
support: [-0.25, 0.25]  (must be [-0.25, 0.25])
max |f*f| attained at 1 knot(s); first at t = -0.352100
continuous check: 9 off-knot t values, max |f*f| = 1.4513063935908062 <= knot max  OK
|max f*f| / (int f)^2 = 1.4513064006879361   (different constant; abs outside the max)
max f*f attained on the positive side; min f*f / max|f*f| = -1.0000

{"n": 102400, "neg_cells": 32150, "zero_cells": 0, "argmax_knot": 30289, "R_float": 1.4513064006879361, "R_exact": "513655524515252183593642627316297728000/353926313748615382700541734286037793169", "R_fft": 1.451306400687936}
```

## Rounding and controls

The printed bound 1.4513064006879361 is **above** the exact value, by 1.2e-18, so it holds as
printed. `verify_standalone.py --bound` checks this, and the script also prints the value rounded
toward +infinity.

`controls.py` shows that the verifier refuses altered certificates and false claims. Each of 13
red controls must be refused: perturbed heights, a height given as a string, `true` or NaN, an
extra or missing key, a removed cell, a zero sum, a bound 1e-15 too low, the exact value
truncated, a wrong hash, and the reversed cells presented under the original hash. Each of 4
positive controls must be accepted: the certificate itself, and the reversed, negated and doubled
heights, which have the same value.

```
red   ok   (refused)  one height scaled by 1+1e-12: REFUSED: exact value 16436976784488069874996564074121527296/11325642039955692968231811878054700605 != claimed 
red   ok   (refused)  first height scaled by 1-1e-9: REFUSED: exact value 513655524511314296948243744043544576000/353926313747829249075990680923389757329 != claime
red   ok   (refused)  a height as a string: REFUSED: every height must be a finite JSON number (no strings, booleans, NaN, Infinity)
red   ok   (refused)  a height as true: REFUSED: every height must be a finite JSON number (no strings, booleans, NaN, Infinity)
red   ok   (refused)  a height as NaN: REFUSED: invalid JSON: non-finite JSON constant NaN
red   ok   (refused)  an extra key: REFUSED: top level must be an object with exactly one key, 'values'
red   ok   (refused)  'values' missing: REFUSED: top level must be an object with exactly one key, 'values'
red   ok   (refused)  last cell removed: REFUSED: exact value 124575466081318752082036684179446626833/70514399478255773290643693602171397405 != claimed
red   ok   (refused)  sum of heights 0: REFUSED: sum of heights is 0 (ratio undefined)
red   ok   (refused)  bound 1e-15 below the claim: REFUSED: claimed bound 1.451306400687935 is below the exact value by 1.099e-15
red   ok   (refused)  bound = exact value truncated to 16 s.f.: REFUSED: claimed bound 1.451306400687936e+00 is below the exact value by 9.880e-17
red   ok   (refused)  wrong sha256: REFUSED: sha256 391cc726acbab3a7f4216dedab2ec5dd74cb797dc12189966dac11a0dee31e2d != claimed 000000000000000000
red   ok   (refused)  cells reversed, with the original sha256: REFUSED: sha256 16d5f172ff04c64f8efd1f378cb9aef50bf91c6a9040a3fe538d90262436be5f != claimed 391cc726acbab3a7f4
pos   ok   (accepted)  the certificate with all claims: ACCEPTED
pos   ok   (accepted)  cells reversed (same value): ACCEPTED
pos   ok   (accepted)  every height negated (same value): ACCEPTED
pos   ok   (accepted)  every height doubled (scale invariant): ACCEPTED
13 red, 4 positive controls; 0 failed
```

## Credit

Derived from the public certificates of open pull requests #372 and #371 on
teorth/optimizationproblems. Found and prepared with AI assistance. No human mathematician has
checked this value, and no independent party has replayed it.
