"""Red and positive controls for verify_standalone.py: show that it refuses altered certificates.

Usage: python3 controls.py construction.json BOUND P/Q SHA256
Each red control alters the certificate (or a claim) and must be REFUSED; each positive control
must be ACCEPTED. Exits 1 if any control behaves otherwise. Python standard library only.
"""

import json
import os
import subprocess
import sys
import tempfile
from fractions import Fraction

path, bound, exact, sha = sys.argv[1:5]
HERE = os.path.dirname(os.path.abspath(__file__))
vals = json.load(open(path))["values"]
n = len(vals)
r = Fraction(exact)

# a heavy cell near the middle of the support (any cell works for the perturbation test)
i_mid = max(range(n // 2 - 50, n // 2 + 50), key=lambda i: abs(vals[i]))


def variant(name, doc=None, text=None):
    d = tempfile.mkdtemp()
    p = os.path.join(d, name + ".json")
    with open(p, "w") as f:
        f.write(text if text is not None else json.dumps(doc))
    return p


def run(p, *claims):
    out = subprocess.run(
        [sys.executable, "-I", os.path.join(HERE, "verify_standalone.py"), p, *claims],
        capture_output=True,
        text=True,
    )
    last = (out.stdout.strip().splitlines() or ["(no output)"])[-1]
    return out.returncode == 0 and last == "ACCEPTED", last


def bump(i, rel):
    w = list(vals)
    w[i] = w[i] * (1 + rel)
    return {"values": w}


low_dec = f"{float(Fraction(bound) - Fraction(10**-15)):.17g}"
full = ["--bound", bound, "--exact", exact]
red = [
    ("one height scaled by 1+1e-12", variant("bump", bump(i_mid, 1e-12)), full),
    ("first height scaled by 1-1e-9", variant("bump0", bump(0, -1e-9)), full),
    (
        "a height as a string",
        variant("str", {"values": [str(vals[0])] + vals[1:]}),
        full,
    ),
    ("a height as true", variant("bool", {"values": [True] + vals[1:]}), full),
    ("a height as NaN", None, full),
    ("an extra key", variant("extra", {"values": vals, "note": "x"}), full),
    ("'values' missing", variant("missing", {"heights": vals}), full),
    ("last cell removed", variant("short", {"values": vals[:-1]}), full),
    ("sum of heights 0", variant("zero", {"values": [1.0, -1.0]}), []),
    ("bound 1e-15 below the claim", path, ["--bound", low_dec]),
    ("bound = exact value truncated to 16 s.f.", path, ["--bound", f"{float(r):.15e}"]),
    ("wrong sha256", path, full + ["--sha256", "0" * 64]),
    (
        "cells reversed, with the original sha256",
        variant("rev", {"values": vals[::-1]}),
        full + ["--sha256", sha],
    ),
]
nan_text = json.dumps({"values": vals}).replace(json.dumps(vals[0]), "NaN", 1)
red[4] = ("a height as NaN", variant("nan", text=nan_text), full)
pos = [
    ("the certificate with all claims", path, full + ["--sha256", sha]),
    ("cells reversed (same value)", variant("rev2", {"values": vals[::-1]}), full),
    (
        "every height negated (same value)",
        variant("neg", {"values": [-x for x in vals]}),
        full,
    ),
    (
        "every height doubled (scale invariant)",
        variant("dbl", {"values": [2 * x for x in vals]}),
        full,
    ),
]
bad = 0
for name, p, claims in red:
    ok, last = run(p, *claims)
    bad += ok
    print(
        f"red   {'FAIL (accepted)' if ok else 'ok   (refused)'}  {name}: {last[:110]}"
    )
for name, p, claims in pos:
    ok, last = run(p, *claims)
    bad += not ok
    print(
        f"pos   {'ok   (accepted)' if ok else 'FAIL (refused)'}  {name}: {last[:110]}"
    )
print(f"{len(red)} red, {len(pos)} positive controls; {bad} failed")
sys.exit(1 if bad else 0)
