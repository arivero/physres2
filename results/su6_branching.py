"""
Branching rules for SU(6) representations under subgroup chains.

Part 1: SU(6) -> SU(5) x U(1)
Part 2: SU(6) -> SU(3) x SU(3) x U(1)
Part 3: Dimension checks
Part 4: Interesting features
"""

from itertools import combinations
from collections import defaultdict


# ─────────────────────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────────────────────

def comb(n, k):
    if k < 0 or k > n:
        return 0
    from math import factorial
    return factorial(n) // (factorial(k) * factorial(n - k))


def sym_dim(n, k):
    """Dimension of Sym^k(n): C(n+k-1, k)."""
    return comb(n + k - 1, k)


def su_rep_dim_from_partition(partition, n):
    """
    Dimension of SU(n) representation from a partition (Young diagram row lengths).
    Uses the hook-length formula: prod_{boxes} (n + content) / hook_length
    where content of box (i,j) is j-i (0-indexed).
    """
    if not partition:
        return 1
    # Pad partition to have n rows
    rows = list(partition) + [0] * (n - len(partition))
    # Build the Young diagram boxes
    boxes = []
    for i, r in enumerate(rows):
        for j in range(r):
            boxes.append((i, j))
    # Compute numerator: product of (n + content) for each box
    num = 1
    for (i, j) in boxes:
        num *= (n + j - i)
    # Compute denominator: product of hook lengths
    denom = 1
    for (i, j) in boxes:
        # hook length = arm + leg + 1
        arm = rows[i] - j - 1  # boxes to the right
        leg = sum(1 for ii in range(i + 1, n) if rows[ii] > j)  # boxes below
        denom *= (arm + leg + 1)
    return num // denom


# ─────────────────────────────────────────────────────────────────────────────
# Part 1: SU(6) -> SU(5) x U(1)
# ─────────────────────────────────────────────────────────────────────────────
#
# Convention: indices 1..6.  Indices 1..5 transform under SU(5); index 6 is
# the "singlet" direction.  The U(1) generator is
#
#   T = diag(1,1,1,1,1,-5) / sqrt(30)     (traceless, Tr T^2 = 1/2)
#
# so each index i carries U(1) charge
#   q_i = 1  for i in {1,2,3,4,5}
#   q_6 = -5
#
# For a multi-index tensor the U(1) charge is the sum of individual charges.
#
# For the conjugate representation the charges are negated:
#   q_i = -1  for i in {1,2,3,4,5}
#   q_6 = +5
#
# SU(5) representations that arise:
#   k indices from {1..5} -> ∧^k(5)  with dimension C(5,k)
#
# The standard SU(5) rep labels are
#   ∧^0 = 1, ∧^1 = 5, ∧^2 = 10, ∧^3 = 10̄, ∧^4 = 5̄, ∧^5 = 1

CHARGE = {i: 1 for i in range(1, 6)}
CHARGE[6] = -5

SU5_ANTISYM = {0: '1', 1: '5', 2: '10', 3: '10̄', 4: '5̄', 5: '1'}


def antisym_su6_to_su5xU1(k):
    """
    Decompose ∧^k(6) under SU(6) -> SU(5) x U(1).

    A basis vector is a sorted k-tuple of indices from {1..6}.
    We split: how many indices are from {1..5} (call it m) vs from {6} (0 or 1).

    Returns: list of (SU5_rep_label, dim, U1_charge)
    """
    results = []
    # m indices from {1..5}, (k-m) copies of index 6 (but index 6 can appear at most once)
    for m in range(k + 1):
        n6 = k - m   # number of 6's
        if n6 > 1:
            continue  # antisymmetric -> 6 can appear at most once
        if n6 < 0:
            continue
        # SU(5) rep: ∧^m(5)
        su5_dim = comb(5, m)
        if su5_dim == 0:
            continue
        su5_label = SU5_ANTISYM.get(m, f'∧^{m}(5)')
        # U(1) charge: m * 1 + n6 * (-5)
        q = m * 1 + n6 * (-5)
        results.append((su5_label, su5_dim, q))
    return results


def sym2_su6_to_su5xU1():
    """
    Decompose Sym²(6) under SU(6) -> SU(5) x U(1).

    Basis: symmetric pairs (i,j) with i<=j from {1..6}.
    Split by how many indices are 6.
      - 0 sixes: both from {1..5} -> Sym²(5) = 15
      - 1 six:   one from {1..5}, one is 6 -> 5 (the SU(5) vector)
      - 2 sixes: both are 6 -> 1 (SU(5) singlet)
    """
    results = []
    # 0 sixes: Sym^2(5), dim = C(6,2) = 15
    results.append(('15', 15, 2 * 1))         # q = 1+1 = 2
    # 1 six: 5 of SU(5), dim = 5
    results.append(('5', 5, 1 + (-5)))         # q = 1 + (-5) = -4
    # 2 sixes: singlet, dim = 1
    results.append(('1', 1, (-5) + (-5)))      # q = -10
    return results


def sym3_su6_to_su5xU1():
    """
    Decompose Sym³(6) under SU(6) -> SU(5) x U(1).

    Basis: symmetric triples (i,j,k) with i<=j<=k from {1..6}.
    Split by number of sixes (0,1,2,3).
      - 0 sixes: all from {1..5} -> Sym³(5), dim = C(7,3) = 35
      - 1 six: two from {1..5}, one is 6 -> Sym²(5) on the SU(5) part = 15
      - 2 sixes: one from {1..5}, two are 6 -> 5 of SU(5)
      - 3 sixes: all 6 -> 1
    """
    results = []
    # 0 sixes: Sym^3(5), dim = C(7,3) = 35
    results.append(('35(SU5)', sym_dim(5, 3), 3 * 1))         # q = 3
    # 1 six: Sym^2(5), dim = 15
    results.append(('15(SU5)', sym_dim(5, 2), 2 * 1 + (-5)))   # q = 2-5 = -3
    # 2 sixes: 5 of SU(5), dim = 5
    results.append(('5', 5, 1 + 2 * (-5)))                      # q = 1-10 = -9
    # 3 sixes: singlet
    results.append(('1', 1, 3 * (-5)))                          # q = -15
    return results


def adjoint_su6_to_su5xU1():
    """
    Decompose the adjoint 35 of SU(6) under SU(6) -> SU(5) x U(1).

    The adjoint of SU(6) = 6 ⊗ 6̄ - 1 (traceless).
    Under SU(5) x U(1):
      6 -> 5(+1) + 1(-5)
      6̄ -> 5̄(-1) + 1(+5)

    6 ⊗ 6̄ = [5(+1) + 1(-5)] ⊗ [5̄(-1) + 1(+5)]
           = (5 ⊗ 5̄)(0) + (5)(+6) + (5̄)(-6) + (1)(0)
           = [24(0) + 1(0)] + 5(+6) + 5̄(-6) + 1(0)

    Subtract the trace (the overall SU(6) singlet, which is 1(0)):
    35 = 24(0) + 5(+6) + 5̄(-6) + 1(0)

    U(1) charges are in units where T = diag(1,1,1,1,1,-5)/sqrt(30).
    For the adjoint T_ij = q_i - q_j.
    """
    results = [
        ('24', 24, 0),    # SU(5) adjoint + singlet combined -> just 24 here
        ('1', 1, 0),      # diagonal SU(5) singlet (the 1 in SU(5) adjoint decomp is inside 24+1=25-1=24; let's redo)
        ('5', 5, 6),
        ('5̄', 5, -6),
    ]
    # Let me redo this more carefully.
    # 5 ⊗ 5̄ of SU(5) = 24 + 1 (adjoint + singlet)
    # So 6 ⊗ 6̄ = (24+1)(0) + 5(6) + 5̄(-6) + 1(0)  [the last 1(0) from 1(-5)⊗1(+5)]
    # Total = 24+1+5+5+1 = 36
    # Subtract 1 singlet (the trace): 35 = 24(0) + 1(0) + 5(6) + 5̄(-6)
    # But wait: we have two 1(0) terms from 5⊗5̄=24+1 and from 1⊗1=1.
    # The SU(6) trace is the overall singlet 1(0). After subtracting it:
    # 35 = 24(0) + 5(6) + 5̄(-6) + 1(0)
    # That last 1(0): from 5⊗5̄ we get 24+1, and separately 1⊗1=1. Both have q=0.
    # Subtracting one 1(0) leaves one 1(0) + 24(0).
    # Total dim check: 24 + 1 + 5 + 5 = 35. ✓
    results = [
        ('24', 24, 0),
        ('1', 1, 0),
        ('5', 5, 6),
        ('5̄', 5, -6),
    ]
    return results


def bar6_to_su5xU1():
    """6̄ -> 5̄(-1) + 1(+5)"""
    return [('5̄', 5, -1), ('1', 1, 5)]


# ─────────────────────────────────────────────────────────────────────────────
# Part 2: SU(6) -> SU(3)_A x SU(3)_B x U(1)
# ─────────────────────────────────────────────────────────────────────────────
#
# Split indices: A = {1,2,3}, B = {4,5,6}.
# U(1) charge: q_i = +1 for i in A, q_i = -1 for i in B  (normalized to ±1)
# (Traceless generator: diag(1,1,1,-1,-1,-1)/sqrt(6))
#
# SU(3) antisymmetric representations:
#   ∧^0 = 1, ∧^1 = 3, ∧^2 = 3̄, ∧^3 = 1
#   (since ∧^2(3) ≅ 3̄ as SU(3) reps, ∧^3(3) = det = singlet)
#
# For Sym representations of SU(3):
#   Sym^0 = 1, Sym^1 = 3, Sym^2 = 6, Sym^3 = 10

SU3_ANTISYM = {0: '1', 1: '3', 2: '3̄', 3: '1'}
SU3_SYM = {0: '1', 1: '3', 2: '6', 3: '10'}


def antisym_su6_to_su3xsu3xU1(k):
    """
    Decompose ∧^k(6) under SU(6) -> SU(3)_A x SU(3)_B x U(1).

    Basis: sorted k-tuples from {1..6}.
    Split: m indices from A={1,2,3}, n=(k-m) from B={4,5,6}.
    Constraint: m <= 3, n <= 3.

    SU(3)_A rep: ∧^m(3)
    SU(3)_B rep: ∧^n(3)
    U(1) charge: m*(+1) + n*(-1) = m - n
    """
    results = []
    for m in range(min(k, 3) + 1):
        n = k - m
        if n < 0 or n > 3:
            continue
        dim_A = comb(3, m)
        dim_B = comb(3, n)
        if dim_A == 0 or dim_B == 0:
            continue
        label_A = SU3_ANTISYM[m]
        label_B = SU3_ANTISYM[n]
        q = m - n
        results.append((label_A, label_B, dim_A * dim_B, q))
    return results


def sym2_su6_to_su3xsu3xU1():
    """
    Decompose Sym²(6) under SU(6) -> SU(3)_A x SU(3)_B x U(1).

    Basis: symmetric pairs (i,j) with i<=j from {1..6}.
    For pair (i,j) with m of them from A and n=2-m from B:
      - m=2, n=0: Sym²(A) = 6_A, dim=6
      - m=1, n=1: (A tensor B) with one from each = 3_A ⊗ 3_B, dim=9
                  (This is NOT Sym² of each separately; it's the mixed term)
      - m=0, n=2: Sym²(B) = 6_B, dim=6

    For the SU(3) x SU(3) content:
      Sym²(6): the 21 basis vectors split as
        - both from A: Sym²(A) = 6, dim 6
        - one from A, one from B: 3_A ⊗ 3_B = 9, dim 9 (the (3,3) of SU(3)xSU(3))
        - both from B: Sym²(B) = 6, dim 6
    """
    results = [
        ('6', '1', 6, 2),    # Sym^2(A), q = 2
        ('3', '3', 9, 0),    # mixed, q = 0
        ('1', '6', 6, -2),   # Sym^2(B), q = -2
    ]
    return results


def sym3_su6_to_su3xsu3xU1():
    """
    Decompose Sym³(6) under SU(6) -> SU(3)_A x SU(3)_B x U(1).

    Basis: sorted triples (i<=j<=k) from {1..6}.
    Split by m from A, n=3-m from B:
      - (m,n) = (3,0): Sym³(A) = 10_A, dim=10
      - (m,n) = (2,1): Sym²(A) ⊗ B = 6_A ⊗ 3_B, dim=18
      - (m,n) = (1,2): A ⊗ Sym²(B) = 3_A ⊗ 6_B, dim=18
      - (m,n) = (0,3): Sym³(B) = 10_B, dim=10
    """
    results = [
        ('10', '1', 10, 3),   # Sym^3(A), q=3
        ('6', '3', 18, 1),    # Sym^2(A) x B, q=1
        ('3', '6', 18, -1),   # A x Sym^2(B), q=-1
        ('1', '10', 10, -3),  # Sym^3(B), q=-3
    ]
    return results


def adjoint_su6_to_su3xsu3xU1():
    """
    Decompose adj(SU(6)) = 35 under SU(6) -> SU(3)_A x SU(3)_B x U(1).

    Adj(SU(6)) = 6 ⊗ 6̄ - 1 (traceless hermitian matrices).
    6 = (3_A, 1)(+1) + (1, 3_B)(-1)
    6̄ = (3̄_A, 1)(-1) + (1, 3̄_B)(+1)

    6 ⊗ 6̄ = [(3,1)(+1) + (1,3)(-1)] ⊗ [(3̄,1)(-1) + (1,3̄)(+1)]
           = (3⊗3̄, 1)(0) + (3, 3̄)(+2) + (3̄, 3)(-2) + (1, 3⊗3̄)(0)
           = [(8+1, 1)(0)] + (3, 3̄)(+2) + (3̄, 3)(-2) + [(1, 8+1)(0)]
           = (8,1)(0) + (1,1)(0) + (3,3̄)(2) + (3̄,3)(-2) + (1,8)(0) + (1,1)(0)
    dim check: 8+1+9+9+8+1 = 36. Subtract trace (overall 1(0)):
    35 = (8,1)(0) + (3,3̄)(+2) + (3̄,3)(-2) + (1,8)(0) + (1,1)(0)
    dim check: 8+9+9+8+1 = 35 ✓

    The (1,1)(0) is one of the two original (1,1) terms; one is subtracted.
    """
    results = [
        ('8', '1',  8, 0),
        ('1', '8',  8, 0),
        ('1', '1',  1, 0),   # the U(1) generator itself
        ('3', '3̄',  9, 2),
        ('3̄', '3',  9, -2),
    ]
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Collect and print results
# ─────────────────────────────────────────────────────────────────────────────

def check_dim(label, expected, components):
    total = sum(d for (_, d, *_) in components)
    status = "OK" if total == expected else f"FAIL (got {total})"
    return status


def check_dim2(label, expected, components):
    total = sum(c[2] for c in components)
    status = "OK" if total == expected else f"FAIL (got {total})"
    return status


def format_su5_table(components):
    """Format a list of (su5_label, dim, q) as a markdown table."""
    lines = []
    lines.append("| SU(5) rep | dim | U(1) charge |")
    lines.append("|-----------|-----|-------------|")
    total = 0
    for (label, dim, q) in components:
        lines.append(f"| {label} | {dim} | {q:+d} |")
        total += dim
    lines.append(f"| **Total** | **{total}** | |")
    return "\n".join(lines)


def format_su3su3_table(components):
    """Format a list of (su3A_label, su3B_label, dim, q) as a markdown table."""
    lines = []
    lines.append("| SU(3)_A rep | SU(3)_B rep | dim | U(1) charge |")
    lines.append("|-------------|-------------|-----|-------------|")
    total = 0
    for (lA, lB, dim, q) in components:
        lines.append(f"| {lA} | {lB} | {dim} | {q:+d} |")
        total += dim
    lines.append(f"| | **Total** | **{total}** | |")
    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Main computation
# ─────────────────────────────────────────────────────────────────────────────

# SU(6) rep info: (name, dimension)
SU6_REPS = [
    ('6',  6),
    ('6̄',  6),
    ('15', 15),
    ('20', 20),
    ('21', 21),
    ('35', 35),
    ('56', 56),
]

# Part 1 decompositions
p1 = {
    '6':  antisym_su6_to_su5xU1(1),
    '6̄':  bar6_to_su5xU1(),
    '15': antisym_su6_to_su5xU1(2),
    '20': antisym_su6_to_su5xU1(3),
    '21': sym2_su6_to_su5xU1(),
    '35': adjoint_su6_to_su5xU1(),
    '56': sym3_su6_to_su5xU1(),
}

# Part 2 decompositions
p2 = {
    '6':  antisym_su6_to_su3xsu3xU1(1),
    '6̄':  [(lA, lB, d, -q) for (lA, lB, d, q) in antisym_su6_to_su3xsu3xU1(1)],
    '15': antisym_su6_to_su3xsu3xU1(2),
    '20': antisym_su6_to_su3xsu3xU1(3),
    '21': sym2_su6_to_su3xsu3xU1(),
    '35': adjoint_su6_to_su3xsu3xU1(),
    '56': sym3_su6_to_su3xsu3xU1(),
}

# ─────────────────────────────────────────────────────────────────────────────
# Build markdown output
# ─────────────────────────────────────────────────────────────────────────────

md_lines = []

md_lines.append("# SU(6) Branching Rules")
md_lines.append("")
md_lines.append("Branching rules for selected SU(6) representations under two subgroup chains.")
md_lines.append("All dimension checks are included. U(1) charges use the conventions stated in each section.")
md_lines.append("")

# ────────────────────────────── PART 1 ────────────────────────────────────────

md_lines.append("## Part 1: SU(6) → SU(5) × U(1)")
md_lines.append("")
md_lines.append("**Convention.** Split index 6 from {1,…,5}. The U(1) generator is")
md_lines.append("")
md_lines.append("    T = diag(1, 1, 1, 1, 1, −5) / √30")
md_lines.append("")
md_lines.append("so each fundamental index carries charge +1 (indices 1–5) or −5 (index 6).")
md_lines.append("For a multi-index state the charge is the sum of individual charges.")
md_lines.append("Conjugate representations have all charges negated.")
md_lines.append("")
md_lines.append("SU(5) antisymmetric tensor representations appearing below:")
md_lines.append("")
md_lines.append("| k | ∧^k(5) | SU(5) label | dim |")
md_lines.append("|---|--------|-------------|-----|")
for k in range(6):
    md_lines.append(f"| {k} | ∧^{k}(5) | {SU5_ANTISYM[k]} | {comb(5,k)} |")
md_lines.append("")

rep_descriptions = {
    '6':  '∧¹(6) — fundamental',
    '6̄':  '∧⁵(6) — anti-fundamental',
    '15': '∧²(6) — antisymmetric 2-tensor',
    '20': '∧³(6) — antisymmetric 3-tensor',
    '21': 'Sym²(6) — symmetric 2-tensor',
    '35': 'adjoint (6 ⊗ 6̄ − 1)',
    '56': 'Sym³(6) — symmetric 3-tensor',
}

for name, dim in SU6_REPS:
    comps = p1[name]
    total = sum(d for (_, d, *_) in comps)
    status = "✓" if total == dim else f"FAIL (sum={total})"
    md_lines.append(f"### SU(6) rep **{name}** ({rep_descriptions[name]}, dim {dim})  —  dimension check: {status}")
    md_lines.append("")
    md_lines.append(format_su5_table(comps))
    md_lines.append("")

    # Singlets
    singlets = [(q) for (label, d, q) in comps if label == '1']
    if singlets:
        md_lines.append(f"Contains SU(5) singlet(s) at U(1) charge(s): {singlets}.")
    else:
        md_lines.append("Contains no SU(5) singlets.")
    md_lines.append("")

# ────────────────────────────── PART 2 ────────────────────────────────────────

md_lines.append("## Part 2: SU(6) → SU(3)_A × SU(3)_B × U(1)")
md_lines.append("")
md_lines.append("**Convention.** Split indices {1,2,3} → SU(3)_A, {4,5,6} → SU(3)_B.")
md_lines.append("The U(1) generator is")
md_lines.append("")
md_lines.append("    T = diag(1, 1, 1, −1, −1, −1) / √6")
md_lines.append("")
md_lines.append("so A-type indices carry charge +1 and B-type indices carry charge −1.")
md_lines.append("")
md_lines.append("SU(3) representations that appear:")
md_lines.append("")
md_lines.append("| label | description | dim |")
md_lines.append("|-------|-------------|-----|")
md_lines.append("| 1 | singlet | 1 |")
md_lines.append("| 3 | fundamental | 3 |")
md_lines.append("| 3̄ | anti-fundamental | 3 |")
md_lines.append("| 6 | symmetric 2-tensor | 6 |")
md_lines.append("| 8 | adjoint | 8 |")
md_lines.append("| 10 | symmetric 3-tensor | 10 |")
md_lines.append("")

for name, dim in SU6_REPS:
    comps = p2[name]
    total = sum(c[2] for c in comps)
    status = "✓" if total == dim else f"FAIL (sum={total})"
    md_lines.append(f"### SU(6) rep **{name}** ({rep_descriptions[name]}, dim {dim})  —  dimension check: {status}")
    md_lines.append("")
    md_lines.append(format_su3su3_table(comps))
    md_lines.append("")

    singlets = [q for (lA, lB, d, q) in comps if lA == '1' and lB == '1']
    if singlets:
        md_lines.append(f"Contains SU(3)_A × SU(3)_B singlet(s) at U(1) charge(s): {singlets}.")
    else:
        md_lines.append("Contains no SU(3)_A × SU(3)_B singlets.")
    md_lines.append("")

# ────────────────────────────── PART 3 ────────────────────────────────────────

md_lines.append("## Part 3: Dimension Check Summary")
md_lines.append("")
md_lines.append("| SU(6) rep | dim | Part 1 sum | Part 1 | Part 2 sum | Part 2 |")
md_lines.append("|-----------|-----|------------|--------|------------|--------|")
for name, dim in SU6_REPS:
    s1 = sum(d for (_, d, *_) in p1[name])
    s2 = sum(c[2] for c in p2[name])
    ok1 = "OK" if s1 == dim else f"FAIL({s1})"
    ok2 = "OK" if s2 == dim else f"FAIL({s2})"
    md_lines.append(f"| {name} | {dim} | {s1} | {ok1} | {s2} | {ok2} |")
md_lines.append("")

# ────────────────────────────── PART 4 ────────────────────────────────────────

md_lines.append("## Part 4: Notable Features")
md_lines.append("")

md_lines.append("### SU(6) → SU(5) × U(1)")
md_lines.append("")
md_lines.append("| SU(6) rep | SU(5) reps present | Contains singlet? | U(1) charge spectrum |")
md_lines.append("|-----------|-------------------|-------------------|----------------------|")
for name, dim in SU6_REPS:
    comps = p1[name]
    su5_reps = ", ".join(label for (label, d, q) in comps)
    singlet_q = [q for (label, d, q) in comps if label == '1']
    has_singlet = f"yes, q = {singlet_q}" if singlet_q else "no"
    charges = sorted(set(q for (_, d, q) in comps))
    md_lines.append(f"| {name} | {su5_reps} | {has_singlet} | {charges} |")
md_lines.append("")

md_lines.append("### SU(6) → SU(3)_A × SU(3)_B × U(1)")
md_lines.append("")
md_lines.append("| SU(6) rep | (A,B) reps present | Contains (1,1) singlet? | U(1) charge spectrum |")
md_lines.append("|-----------|-------------------|------------------------|----------------------|")
for name, dim in SU6_REPS:
    comps = p2[name]
    ab_reps = ", ".join(f"({lA},{lB})" for (lA, lB, d, q) in comps)
    singlet_q = [q for (lA, lB, d, q) in comps if lA == '1' and lB == '1']
    has_singlet = f"yes, q = {singlet_q}" if singlet_q else "no"
    charges = sorted(set(q for (_, _, d, q) in comps))
    md_lines.append(f"| {name} | {ab_reps} | {has_singlet} | {charges} |")
md_lines.append("")

md_lines.append("### Physical interpretation notes")
md_lines.append("")
md_lines.append("**6 = (5)(+1) + (1)(−5):** The fundamental splits into the SU(5) quintet"
                " plus one singlet with large negative U(1) charge.")
md_lines.append("")
md_lines.append("**15 = (10)(+2) + (5)(−4):** The antisymmetric 2-tensor of SU(6) branches"
                " to the SU(5) decuplet at charge +2 and the quintet at charge −4."
                " No SU(5) singlet. This is the representation relevant for antisymmetric"
                " diquark operators.")
md_lines.append("")
md_lines.append("**20 = (10)(+3) + (10̄)(−3) :** The antisymmetric 3-tensor is self-conjugate"
                " as an SU(6) branching set: the SU(5) 10 and its conjugate 10̄ appear"
                " with opposite charges. No SU(5) singlet.")
md_lines.append("")
md_lines.append("**21 = (15)(+2) + (5)(−4) + (1)(−10):** The symmetric 2-tensor contains"
                " one SU(5) singlet at charge −10 (high magnitude). The 15 of SU(5)"
                " (symmetric 2-tensor of SU(5)) appears at charge +2.")
md_lines.append("")
md_lines.append("**35 = (24)(0) + (1)(0) + (5)(+6) + (5̄)(−6):** The adjoint decomposes"
                " into the SU(5) adjoint (24) and one singlet, both at zero U(1) charge,"
                " plus the quintet/antiquintet pair at ±6. The zero-charge singlet is the"
                " U(1) generator itself.")
md_lines.append("")
md_lines.append("**56 = (35(SU5))(+3) + (15(SU5))(−3) + (5)(−9) + (1)(−15):**"
                " The symmetric cubic contains one SU(5) singlet at charge −15. The"
                " SU(5) representations 35 and 15 here are the symmetric 3- and 2-tensors"
                " of SU(5) (not the same as ∧²(5)=10 or ∧³(5)=10̄).")
md_lines.append("")
md_lines.append("**SU(3)×SU(3) split of 20:** The antisymmetric 3-tensor decomposes as"
                " (1,1)+( 3,3̄)+(3̄,3)+(1,1) at charges ±3 and ±1 — wait, let us be precise.")
md_lines.append("")

# Print precise 20 decomp
comps20 = p2['20']
for (lA, lB, d, q) in comps20:
    md_lines.append(f"- ({lA}, {lB}), dim={d}, q={q:+d}")
md_lines.append("")
md_lines.append("The two (1,1) singlets in the 20 under SU(3)×SU(3) are the totally"
                " antisymmetric tensors ε_{abc} and ε_{def} acting as scalars on each"
                " SU(3) factor — the volume forms.")
md_lines.append("")

md_lines.append("**SU(3)×SU(3) split of 35 (adjoint):** The decomposition"
                " (8,1)+(1,8)+(1,1)+(3,3̄)+(3̄,3) shows the two SU(3) sub-algebras"
                " as (8,1) and (1,8), the U(1) as (1,1), and the off-diagonal"
                " coset directions as (3,3̄)+(3̄,3). This is the standard pattern"
                " for a symmetric SU(3)×SU(3)×U(1) embedding.")
md_lines.append("")

md_lines.append("---")
md_lines.append("")
md_lines.append("*All branching rules derived by explicit index splitting."
                " Dimension checks: sum of component dimensions equals SU(6) rep dimension in all cases.*")

output = "\n".join(md_lines)
print(output)

# Save
with open("/home/codexssh/phys3/results/su6_branching.md", "w") as f:
    f.write(output)

print("\n\n[Script complete. Results written to su6_branching.md]")
