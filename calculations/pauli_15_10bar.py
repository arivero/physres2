#!/usr/bin/env python3
"""
Pauli exclusion analysis for diquarks in SU(5)_flavor representations.

Problem: diquarks (two-quark composites) with J=0, color 3-bar.
Which SU(5)_flavor representation do they fill: 15 (symmetric) or 10-bar (antisymmetric)?

SU(5) flavor: quarks q = (u, d, s, c, b) in the fundamental 5.
5 ⊗ 5 = 15_S ⊕ 10-bar_A

Under SU(3)_color × SU(2)_L, the 5 of SU(5) → (3,1) ⊕ (1,2).
But here SU(5) is FLAVOR, not GUT. All five quarks are color triplets.
"""

import itertools
from collections import defaultdict

# ============================================================
# Task 1: Explicit decomposition of 15 and 10-bar
# ============================================================

flavors = ['u', 'd', 's', 'c', 'b']
N = len(flavors)

print("=" * 70)
print("TASK 1: Decompose 5 ⊗ 5 = 15_S ⊕ 10-bar_A explicitly")
print("=" * 70)

# Symmetric product: 15
# Basis: { e_a ⊗ e_b + e_b ⊗ e_a } for a <= b  (unnormalized)
# Diagonal: e_a ⊗ e_a  (5 elements)
# Off-diagonal: e_a ⊗ e_b + e_b ⊗ e_a for a < b  (C(5,2) = 10 elements)
# Total: 5 + 10 = 15 ✓

sym_diag = []
sym_offdiag = []
for i, a in enumerate(flavors):
    for j, b in enumerate(flavors):
        if i == j:
            sym_diag.append((a, b))
        elif i < j:
            sym_offdiag.append((a, b))

print(f"\n15 (symmetric tensor D^{{ab}} = D^{{ba}}):")
print(f"  Diagonal (same-flavor diquarks): {len(sym_diag)} states")
for pair in sym_diag:
    print(f"    ({pair[0]}{pair[1]})")

print(f"  Off-diagonal (mixed-flavor diquarks): {len(sym_offdiag)} states")
for pair in sym_offdiag:
    print(f"    ({pair[0]}{pair[1]} + {pair[1]}{pair[0]})/√2")

print(f"  Total: {len(sym_diag) + len(sym_offdiag)} = 15 ✓")

# Antisymmetric product: 10-bar
# Basis: { e_a ⊗ e_b - e_b ⊗ e_a } for a < b
# Total: C(5,2) = 10
antisym = []
for i, a in enumerate(flavors):
    for j, b in enumerate(flavors):
        if i < j:
            antisym.append((a, b))

print(f"\n10-bar (antisymmetric tensor D^{{[ab]}} = -D^{{[ba]}}):")
print(f"  Mixed-flavor diquarks only: {len(antisym)} states")
for pair in antisym:
    print(f"    ({pair[0]}{pair[1]} - {pair[1]}{pair[0]})/√2")
print(f"  Total: {len(antisym)} = 10 ✓")

print(f"\nDimension check: 15 + 10 = {15 + 10} = N(N+1)/2 + N(N-1)/2 = {N*(N+1)//2} + {N*(N-1)//2} = {N*N} = 5² ✓")

# ============================================================
# Task 2: Color-spin-flavor factorization for J=0, color 3-bar
# ============================================================

print("\n" + "=" * 70)
print("TASK 2: Color-spin-flavor factorization")
print("=" * 70)

print("""
For a two-fermion state, the TOTAL wavefunction must be antisymmetric
under exchange of the two identical fermions (Pauli exclusion principle).

Ψ_total = ψ_spatial × ψ_color × ψ_spin × ψ_flavor

Exchange symmetry of each factor:
  Spatial (S-wave, L=0): SYMMETRIC  → sign = +1
  Color (3̄ channel):    ANTISYMMETRIC → sign = -1
     (3 ⊗ 3 = 6_S ⊕ 3̄_A; the 3̄ is antisymmetric)
  Spin (J=0 singlet):   ANTISYMMETRIC → sign = -1
     (↑↓ - ↓↑)/√2

Total symmetry of (spatial × color × spin):
  (+1) × (-1) × (-1) = +1  (SYMMETRIC)

For total to be antisymmetric:
  flavor wavefunction must be ANTISYMMETRIC → sign = -1

Therefore: J=0, S-wave, color-3̄ diquarks require flavor in the 10-bar (antisymmetric).

The 15 (symmetric flavor) is FORBIDDEN by Pauli for identical fermions
at J=0, L=0, color 3̄.
""")

# Systematic sign table
print("Sign table under particle exchange:")
print(f"  {'Factor':<20} {'Channel':<25} {'Symmetry':<15} {'Sign':>5}")
print(f"  {'-'*20} {'-'*25} {'-'*15} {'-'*5}")
factors = [
    ("Spatial", "S-wave (L=0)", "symmetric", +1),
    ("Color", "3̄ (antisym)", "antisymmetric", -1),
    ("Spin", "singlet (J=0)", "antisymmetric", -1),
]
product = 1
for name, channel, sym, sign in factors:
    print(f"  {name:<20} {channel:<25} {sym:<15} {sign:>+2}")
    product *= sign

print(f"  {'-'*20} {'-'*25} {'-'*15} {'-'*5}")
print(f"  {'Product':<20} {'':<25} {'':15} {product:>+2}")
print(f"\n  Pauli requires total = -1")
print(f"  Product without flavor = {product:+d}")
print(f"  Therefore flavor must have sign = {-1 * product // abs(product):+d} → {'ANTISYMMETRIC (10-bar)' if product > 0 else 'SYMMETRIC (15)'}")

# ============================================================
# Task 3: Mixed-flavor composites and the "different species" loophole
# ============================================================

print("\n" + "=" * 70)
print("TASK 3: Resolution via mixed-flavor composites?")
print("=" * 70)

print("""
Key question: Does Pauli exclusion apply to diquarks made of DIFFERENT
quark flavors, like (u,d) or (s,c)?

ANALYSIS:

Pauli exclusion applies to IDENTICAL fermions. Two quarks of different
flavor (e.g., u and d) are distinguishable particles. The full
wavefunction need NOT be antisymmetric under their exchange.

However, there is a subtlety: when we organize states into SU(5)_flavor
multiplets, the representation theory FORCES a definite exchange symmetry
on the flavor indices, regardless of whether the quarks are "identical."

Consider the 15 of SU(5)_flavor:
""")

print("States in the 15:")
print("\n  DIAGONAL (same-flavor — Pauli DOES apply):")
for pair in sym_diag:
    pauli_ok = False  # Same-flavor, symmetric → Pauli violated at J=0, L=0
    print(f"    D^{{{pair[0]}{pair[1]}}} = |{pair[0]}{pair[1]}⟩  "
          f"[Pauli: VIOLATED for J=0, L=0, color 3̄]")

print("\n  OFF-DIAGONAL (different-flavor — Pauli does NOT directly apply):")
for pair in sym_offdiag:
    print(f"    D^{{{pair[0]}{pair[1]}}} = (|{pair[0]}{pair[1]}⟩ + |{pair[1]}{pair[0]}⟩)/√2  "
          f"[Pauli: not directly applicable]")

print(f"""
VERDICT on the mixed-flavor loophole:

Strictly speaking, Pauli exclusion only constrains the 5 diagonal entries
(uu, dd, ss, cc, bb). The 10 off-diagonal entries of the 15 are NOT
individually forbidden by Pauli.

BUT: The 15 is an IRREDUCIBLE representation of SU(5). You cannot
consistently keep only the off-diagonal part — that is not a representation.
The off-diagonal elements of the 15 and the 10-bar share the same flavor
pairs (a,b) with a ≠ b; they differ only in symmetrization.

If we restrict to a ≠ b only:
  - From the 15: 10 states (symmetric combinations)
  - From the 10-bar: 10 states (antisymmetric combinations)
  These span the full 20-dimensional space of off-diagonal 5⊗5.
  Neither subset alone forms an irrep of SU(5).

So the "just use off-diagonal entries of the 15" approach:
  1. Breaks SU(5) flavor symmetry (not a valid irrep)
  2. Is equivalent to choosing a SPECIFIC LINEAR COMBINATION of 15 and 10-bar
  3. Cannot be physical if the theory has exact SU(5)_flavor symmetry

HOWEVER: If SU(5)_flavor is BROKEN (as it is in nature — quarks have
different masses), this loophole could work in principle. You would need
the 5 diagonal diquarks to be removed from the spectrum by some mechanism
(e.g., they are too heavy, or selection rules forbid them).
""")

# ============================================================
# Task 4: P-wave (L=1) resolution
# ============================================================

print("=" * 70)
print("TASK 4: Resolution via nonzero orbital angular momentum")
print("=" * 70)

print("""
If the diquark has L=1 (P-wave), the spatial wavefunction is ANTISYMMETRIC
under exchange.

Revised sign table:
""")

factors_pwave = [
    ("Spatial", "P-wave (L=1)", "antisymmetric", -1),
    ("Color", "3̄ (antisym)", "antisymmetric", -1),
    ("Spin", "singlet (J=0)", "antisymmetric", -1),
]
product_pwave = 1
print(f"  {'Factor':<20} {'Channel':<25} {'Symmetry':<15} {'Sign':>5}")
print(f"  {'-'*20} {'-'*25} {'-'*15} {'-'*5}")
for name, channel, sym, sign in factors_pwave:
    print(f"  {name:<20} {channel:<25} {sym:<15} {sign:>+2}")
    product_pwave *= sign

print(f"  {'-'*20} {'-'*25} {'-'*15} {'-'*5}")
print(f"  {'Product':<20} {'':<25} {'':15} {product_pwave:>+2}")
print(f"\n  Product without flavor = {product_pwave:+d}")
print(f"  Flavor must have sign = {-1 * product_pwave // abs(product_pwave):+d}")
print(f"  → {'SYMMETRIC (15)' if product_pwave < 0 else 'ANTISYMMETRIC (10-bar)'}")

print("""
YES: P-wave diquarks can be in the 15 (symmetric flavor)!

But wait — if L=1 and S=0 (spin singlet), then J = L = 1.
This is NOT a J=0 diquark anymore. It has J^P = 1^-.

Can we get J=0 with L=1? We need J = |L - S| to |L + S|.
  - S=0, L=1: J=1 only
  - S=1, L=1: J=0, 1, 2

For J=0 with L=1, we need S=1 (spin triplet, SYMMETRIC):

Revised for S=1, L=1, J=0:
""")

factors_pwave_triplet = [
    ("Spatial", "P-wave (L=1)", "antisymmetric", -1),
    ("Color", "3̄ (antisym)", "antisymmetric", -1),
    ("Spin", "triplet (S=1)", "symmetric", +1),
]
product_pt = 1
print(f"  {'Factor':<20} {'Channel':<25} {'Symmetry':<15} {'Sign':>5}")
print(f"  {'-'*20} {'-'*25} {'-'*15} {'-'*5}")
for name, channel, sym, sign in factors_pwave_triplet:
    print(f"  {name:<20} {channel:<25} {sym:<15} {sign:>+2}")
    product_pt *= sign

print(f"  {'-'*20} {'-'*25} {'-'*15} {'-'*5}")
print(f"  {'Product':<20} {'':<25} {'':15} {product_pt:>+2}")
print(f"\n  Product without flavor = {product_pt:+d}")
print(f"  Flavor sign needed = {-1 * product_pt // abs(product_pt):+d}")
print(f"  → {'SYMMETRIC (15)' if product_pt < 0 else 'ANTISYMMETRIC (10-bar)'}")

print("""
Result: J=0 with L=1, S=1 → flavor must be ANTISYMMETRIC → 10-bar again!

Let's check ALL combinations systematically:
""")

print(f"  {'L':>3} {'S':>3} {'J values':>12} {'Spatial':>12} {'Spin':>12} "
      f"{'Color':>12} {'Prod(no flav)':>15} {'Flavor needed':>15}")
print(f"  {'-'*3} {'-'*3} {'-'*12} {'-'*12} {'-'*12} {'-'*12} {'-'*15} {'-'*15}")

for L in range(3):
    for S in [0, 1]:
        spatial_sign = (-1)**L
        spin_sign = (-1)**(S+1)  # antisym for S=0 (singlet), sym for S=1 (triplet)
        color_sign = -1  # always antisym for 3-bar
        product = spatial_sign * spin_sign * color_sign
        flavor_needed = -1 * product // abs(product)
        J_vals = list(range(abs(L-S), L+S+1))
        J_str = ",".join(str(j) for j in J_vals)
        sp_str = "sym" if spatial_sign > 0 else "antisym"
        sn_str = "sym" if spin_sign > 0 else "antisym"
        fl_str = "15(sym)" if flavor_needed > 0 else "10-bar(antisym)"

        print(f"  {L:>3} {S:>3} {J_str:>12} {sp_str:>12} {sn_str:>12} "
              f"{'antisym':>12} {product:>+15} {fl_str:>15}")

print("""
PATTERN: The 15 (symmetric flavor) is allowed for:
  - L=0, S=1 → J=1 (vector diquark, J^P = 1^+)
  - L=1, S=0 → J=1 (pseudovector diquark, J^P = 1^-)
  - L=2, S=1 → J=1,2,3 (tensor diquarks)
  - Generally: L+S = odd → 15; L+S = even → 10-bar

For J=0 diquarks specifically:
  - J=0 requires L=S, and:
    - L=S=0: even sum → 10-bar  (J^P = 0^+)
    - L=S=1: even sum → 10-bar  (J^P = 0^-)
    - L=S=2: even sum → 10-bar  (J^P = 0^+)
  J=0 ALWAYS requires 10-bar. The 15 is NEVER accessible for J=0.

This is actually a theorem: J=0 means L=S, so L+S=2S (always even),
so flavor is always antisymmetric for J=0, color-3̄ diquarks.
""")

# ============================================================
# Task 5: SUSY resolution — chiral superfield composites
# ============================================================

print("=" * 70)
print("TASK 5: Resolution via SUSY (chiral superfield composites)")
print("=" * 70)

print("""
In SUSY QCD (SQCD), the fundamental fields are chiral superfields:
  Q^a = (q^a, ψ^a, F^a)
where q^a is a scalar squark, ψ^a is a Weyl fermion (quark), F^a is auxiliary.

The composite "diquark" (baryon-like) operator in SU(3)_color SQCD:
  B^{ab} = ε_{ijk} Q^{ia} Q^{jb}    (summing over color i,j,k)

Key insight: the ε_{ijk} makes the COLOR contraction antisymmetric.
This is the ε-tensor of SU(3) color, giving a color singlet.

Now, Q^{ia} Q^{jb} involves the PRODUCT of two chiral superfields.
In the SUSY context:

1. If Q are BOSONIC chiral superfields (as in the holomorphic description
   of the moduli space), there is NO Pauli exclusion. The superfields
   commute, and B^{ab} = ε_{ijk} Q^{ia} Q^{jb} is automatically
   ANTISYMMETRIC in (a,b) due to the ε_{ijk}:

   B^{ab} = ε_{ijk} Q^{ia} Q^{jb} = -ε_{jik} Q^{jb} Q^{ia} = -B^{ba}

   Wait — let's be more careful.

   ε_{ijk} Q^{ia} Q^{jb}: swap a↔b means swap the flavor indices:
   ε_{ijk} Q^{ib} Q^{ja}

   Now swap the dummy color indices i↔j:
   ε_{jik} Q^{ja} Q^{ib} = -ε_{ijk} Q^{ja} Q^{ib}

   For COMMUTING (bosonic) superfields: Q^{ja} Q^{ib} = Q^{ib} Q^{ja}
   So: ε_{ijk} Q^{ib} Q^{ja} = -ε_{ijk} Q^{ja} Q^{ib} = -ε_{ijk} Q^{ib} Q^{ja}

   Hmm, let's redo this cleanly.

CAREFUL ANALYSIS of B^{ab} = ε_{ijk} Q^{ia} Q^{jb}:

For N_c = 3 colors, N_f = 5 flavors:

Baryon operator (for N_c = 3):
  B^{abc} = ε_{ijk} Q^{ia} Q^{jb} Q^{kc}

This is the standard SQCD baryon for N_c = 3. It has THREE flavor indices.
The ε_{ijk} antisymmetrizes the color. For commuting superfields:

  B^{abc} = ε_{ijk} Q^{ia} Q^{jb} Q^{kc}

Under permutation of any two flavor indices, say a↔b:
  ε_{ijk} Q^{ib} Q^{ja} Q^{kc}

Swap dummy indices i↔j:
  ε_{jik} Q^{jb} Q^{ia} Q^{kc} = -ε_{ijk} Q^{jb} Q^{ia} Q^{kc}

Since Q's commute: Q^{jb} Q^{ia} = Q^{ia} Q^{jb}
  = -ε_{ijk} Q^{ia} Q^{jb} Q^{kc} = -B^{abc}

So B^{abc} is ANTISYMMETRIC in all flavor pairs → B ∈ ∧³(5) = 10-bar.

For N_c = 2 (or diquarks in a different context):
  B^{ab} = ε_{ij} Q^{ia} Q^{jb}  (for N_c = 2)

Swap a↔b: ε_{ij} Q^{ib} Q^{ja}
Swap i↔j: ε_{ji} Q^{ja} Q^{ib} = -ε_{ij} Q^{ja} Q^{ib}
Commuting: = -ε_{ij} Q^{ia} Q^{jb} = -B^{ab}

Again ANTISYMMETRIC → 10-bar of SU(5)_flavor.

CRITICAL POINT: Even for BOSONIC (commuting) chiral superfields, the
ε-tensor antisymmetrization on color FORCES antisymmetry on flavor.
This is PURELY from the epsilon tensor — no Pauli principle invoked!

For N_c = 3 SQCD baryons: B ∈ ∧³(5) = 10-bar   (dim 10)
For diquarks in color 3̄: D ∈ ∧²(5) = 10-bar    (dim 10)

The 15 (symmetric) requires the COLOR contraction to be in the
symmetric channel: 3 ⊗ 3 → 6. But 6 is NOT color-singlet or color-3̄.
For SU(3)_color: 3 ⊗ 3 = 3̄ ⊕ 6. Only the 3̄ can form a color-singlet
baryon (with one more quark). The 6 cannot.

SUSY RESOLUTION VERDICT:
In standard SQCD, the baryon/diquark composite superfields are
AUTOMATICALLY in the antisymmetric representation. This is not from
Pauli but from the ε-tensor color contraction. The 15 is associated
with the color-6 channel, which has no role in forming color singlets.

HOWEVER: The sBootstrap framework is NOT standard SQCD. The composites
could be identified differently — as entries in a meson matrix M^a_b
(adjoint of SU(5)) rather than as baryons. Mesons M^a_b = Q^a Q̄_b
have NO antisymmetrization requirement since a and b are in different
representations (5 and 5̄).
""")

# ============================================================
# Task 6: Representation theory checks
# ============================================================

print("=" * 70)
print("TASK 6: Representation theory consistency checks")
print("=" * 70)

print(f"""
For SU(N) with N = {N}:

1. Tensor products of fundamentals:
   5 ⊗ 5 = 15_S ⊕ 10-bar_A
   dim: {N*(N+1)//2} + {N*(N-1)//2} = {N*N} = {N}² ✓

2. Meson matrix (fundamental × antifundamental):
   5 ⊗ 5̄ = 24 ⊕ 1   (adjoint + singlet)
   dim: {N**2 - 1} + 1 = {N**2} = {N}² ✓

3. Claimed decomposition "24 = 15 ⊕ 10̄ − 1":
   This is NOT a valid group theory statement. Let's check what was meant.

   24 ≠ 15 ⊕ 10̄ as representations (they are in different tensor spaces).
   - 24 lives in 5 ⊗ 5̄ (mixed indices, M^a_b)
   - 15 and 10̄ live in 5 ⊗ 5 (both upper indices, D^{{ab}})

   Dimensionally: 24 = 15 + 10 - 1? → 24 = 24 ✓ (numerically works)

   But this is a COINCIDENCE of dimensions, not a decomposition.
   The three representations live in different tensor spaces.

4. Correct decompositions:
   Symmetric square:    S²(5) = 15     (dim {N*(N+1)//2})
   Antisymmetric square: ∧²(5) = 10    (dim {N*(N-1)//2})
   Note: ∧²(5) = 10 (not 10̄). More precisely:
     ∧²(5) is the 10 of SU(5)
     10̄ is its conjugate

   For the meson matrix:
     5 ⊗ 5̄ = 24 ⊕ 1

5. Dimension check of SU(5) representations:
""")

# List SU(5) representations and dimensions
reps = [
    ("1", 1, "singlet"),
    ("5", 5, "fundamental"),
    ("5̄", 5, "anti-fundamental"),
    ("10", 10, "∧²(5)"),
    ("10̄", 10, "∧²(5̄)"),
    ("15", 15, "S²(5)"),
    ("15̄", 15, "S²(5̄)"),
    ("24", 24, "adjoint"),
    ("35", 35, "∧³(5̄) [= 10̄ of 'baryon' type for SU(5)]"),
    ("40", 40, "mixed symmetry"),
    ("45", 45, "∧²(5)⊗5 traceless"),
    ("50", 50, "S³(5)"),
]

for name, dim, desc in reps:
    print(f"  {name:>4}: dim = {dim:>3}   ({desc})")

print(f"""
6. How 15 and 10̄ relate to the adjoint 24:

   The adjoint 24 is the traceless part of 5 ⊗ 5̄.
   The 15 is the symmetric part of 5 ⊗ 5.
   The 10 (or 10̄) is the antisymmetric part of 5 ⊗ 5.

   These are different objects. The coincidence 15 + 10 - 1 = 24
   reflects the identity:
     N(N+1)/2 + N(N-1)/2 - 1 = N² - 1
   which is trivially N²/2 + N²/2 - 1 = N² - 1.

   There IS a deeper connection via the REAL representation of SU(N):
   The N(N-1)/2-dimensional antisymmetric tensor = the fundamental of SO(N).
   And dim(adjoint of SU(N)) = dim(S²(N)) + dim(∧²(N)) - 1.
   But this is arithmetic, not a rep-theoretic decomposition.
""")

# ============================================================
# Task 7: Summary — when can diquarks fill the 15?
# ============================================================

print("=" * 70)
print("TASK 7: Summary — Under what conditions can diquarks fill the 15?")
print("=" * 70)

print("""
QUESTION: Under what conditions (if any) can the J=0 diquark spectrum
fill the 15 of SU(5)_flavor instead of the 10̄?

SHORT ANSWER: For standard J=0 diquarks, NEVER. But there are loopholes.

DETAILED ANALYSIS OF LOOPHOLES:

╔══════════════════════════════════════════════════════════════════════╗
║ Loophole 1: Nonzero orbital angular momentum (L ≠ 0)              ║
╠══════════════════════════════════════════════════════════════════════╣
║ Status: CLOSED for J=0                                             ║
║                                                                    ║
║ J=0 requires L=S. Then L+S = 2S (always even), so the flavor      ║
║ wavefunction is always antisymmetric for color-3̄ diquarks.        ║
║ The 15 requires L+S = odd, which gives J ≥ 1.                     ║
║                                                                    ║
║ Physical cost: Must abandon J=0. The lightest 15-type diquarks     ║
║ have J=1 (either S-wave spin-triplet or P-wave spin-singlet).      ║
╚══════════════════════════════════════════════════════════════════════╝

╔══════════════════════════════════════════════════════════════════════╗
║ Loophole 2: Different color channel (6 instead of 3̄)             ║
╠══════════════════════════════════════════════════════════════════════╣
║ Status: OPEN but costly                                            ║
║                                                                    ║
║ If diquarks are in the color-symmetric 6 channel:                  ║
║   Color: symmetric (+1)                                            ║
║   Spin singlet: antisymmetric (-1)                                 ║
║   Spatial S-wave: symmetric (+1)                                   ║
║   Product = -1 → flavor must be SYMMETRIC → 15 ✓                  ║
║                                                                    ║
║ Physical cost: Color-6 diquarks cannot form color-singlet baryons  ║
║ with a single antiquark. They would need exotic color structures.   ║
║ In standard QCD, these are NOT the diquarks relevant for baryons.  ║
╚══════════════════════════════════════════════════════════════════════╝

╔══════════════════════════════════════════════════════════════════════╗
║ Loophole 3: Broken SU(5)_flavor (off-diagonal entries only)       ║
╠══════════════════════════════════════════════════════════════════════╣
║ Status: OPEN but breaks SU(5)                                      ║
║                                                                    ║
║ For mixed-flavor diquarks (a ≠ b), Pauli doesn't directly apply.  ║
║ The 10 off-diagonal entries of the 15 could in principle exist.    ║
║                                                                    ║
║ Physical cost: These 10 states don't form an SU(5) irrep.          ║
║ Requires SU(5) flavor breaking to remove the 5 diagonal states.    ║
║ Conceptually ad hoc.                                               ║
╚══════════════════════════════════════════════════════════════════════╝

╔══════════════════════════════════════════════════════════════════════╗
║ Loophole 4: SUSY chiral superfield reinterpretation               ║
╠══════════════════════════════════════════════════════════════════════╣
║ Status: PARTIALLY OPEN                                             ║
║                                                                    ║
║ In SQCD, the composite B^{ab} = ε_{ij} Q^{ia} Q^{jb} is          ║
║ antisymmetric in flavor due to the ε-tensor, NOT due to Pauli.     ║
║ This gives 10̄ (or 10) automatically.                              ║
║                                                                    ║
║ To get the 15, one would need the symmetric color contraction:     ║
║   B̃^{ab} = δ_{(ij)} Q^{ia} Q^{jb} → color 6                     ║
║ This is symmetric in flavor, giving the 15. But it's color-6.     ║
║                                                                    ║
║ HOWEVER: if the composites are MESONS (M^a_b = Q^a Q̄_b) rather   ║
║ than diquarks, the Pauli/epsilon issue is completely irrelevant.   ║
║ Mesons live in 5 ⊗ 5̄ = 24 ⊕ 1, not in 5 ⊗ 5.                   ║
║ The 15 and 10̄ simply don't appear in the meson sector.            ║
║                                                                    ║
║ Physical cost for diquarks: need color-6 channel (exotic).         ║
║ For mesons: the 15 is irrelevant (mesons are in the adjoint).      ║
╚══════════════════════════════════════════════════════════════════════╝

╔══════════════════════════════════════════════════════════════════════╗
║ Loophole 5: Non-relativistic vs relativistic treatment             ║
╠══════════════════════════════════════════════════════════════════════╣
║ Status: CLOSED                                                     ║
║                                                                    ║
║ The Pauli exclusion principle is exact in QFT (spin-statistics     ║
║ theorem). It's not an artifact of the non-relativistic limit.      ║
║ The sign analysis using (spatial)(color)(spin)(flavor) is valid     ║
║ in the non-relativistic limit, and the conclusion extends to the   ║
║ fully relativistic case via the spin-statistics theorem.            ║
╚══════════════════════════════════════════════════════════════════════╝

╔══════════════════════════════════════════════════════════════════════╗
║ Loophole 6: Composite at different spacetime points                ║
╠══════════════════════════════════════════════════════════════════════╣
║ Status: OPEN (exotic)                                              ║
║                                                                    ║
║ Pauli exclusion in the strict sense applies when two fermions are  ║
║ at the same spacetime point (same quantum numbers including        ║
║ position). Extended composites with structure (like a string or    ║
║ flux tube) don't have both quarks at the same point.               ║
║                                                                    ║
║ In practice, this is encoded in the partial-wave expansion (the L  ║
║ quantum number). Non-pointlike composites naturally have L ≥ 1     ║
║ contributions, which can access the 15. But L ≥ 1 diquarks are    ║
║ heavier and have higher spin (back to Loophole 1).                 ║
╚══════════════════════════════════════════════════════════════════════╝

CONCLUSION:

There is NO way to put J=0, color-3̄ diquarks in the 15 of SU(5)_flavor.
The theorem is: for J=0 diquarks, L=S always, so L+S is even, so flavor
is antisymmetric → 10̄ obligatory.

The 15 requires EITHER:
  (a) J ≥ 1 (abandon scalar diquarks), OR
  (b) Color-6 channel (exotic, no standard baryons), OR
  (c) Broken SU(5) flavor (keep only off-diagonal part, ad hoc), OR
  (d) Reinterpreting the composites as something other than diquarks
      (e.g., Seiberg dual mesons in the adjoint of SU(5), which is a
      completely different tensor structure).

The most natural exit is (d): the physical scalars that fill the 15
are NOT diquarks but rather objects in a dual description where the
representation constraint arises differently.
""")

# ============================================================
# Write results to markdown
# ============================================================

output = """# Pauli Exclusion and the 15 vs 10-bar of SU(5) Flavor

## Setup

Five quark flavors (u, d, s, c, b) form the fundamental **5** of SU(5)_flavor.
The tensor product decomposes as:

    5 ⊗ 5 = 15_S ⊕ 10-bar_A

where 15 is the symmetric product (diquarks with symmetric flavor) and 10-bar is
the antisymmetric product.

## 1. Explicit decomposition

**15 (symmetric)**: 15 states total
- 5 diagonal (same-flavor): uu, dd, ss, cc, bb
- 10 off-diagonal (mixed-flavor): (ud+du)/√2, (us+su)/√2, ..., (cb+bc)/√2

**10-bar (antisymmetric)**: 10 states total
- 10 mixed-flavor only: (ud−du)/√2, (us−su)/√2, ..., (cb−bc)/√2

Check: 15 + 10 = 25 = 5²  ✓

## 2. Color-spin-flavor factorization

For J=0, S-wave, color-3-bar diquarks:

| Factor  | Channel       | Exchange symmetry | Sign |
|---------|---------------|-------------------|------|
| Spatial | S-wave (L=0)  | symmetric         | +1   |
| Color   | 3-bar         | antisymmetric     | −1   |
| Spin    | singlet (J=0) | antisymmetric     | −1   |
| **Product** |           |                   | **+1** |

Pauli requires total antisymmetry (−1), so flavor must be **antisymmetric → 10-bar**.

## 3. Mixed-flavor loophole

For diquarks of different flavors (e.g., u,d), Pauli exclusion does not directly
apply since they are distinguishable particles. The 10 off-diagonal entries of the
15 are not individually forbidden.

**But**: these 10 states do not form an irrep of SU(5). Keeping only off-diagonal
entries of the 15 breaks the representation. This loophole requires explicit SU(5)
flavor breaking and a mechanism to remove the 5 diagonal (same-flavor) states.

## 4. Orbital angular momentum

Systematic scan of all (L, S) combinations for color-3-bar diquarks:

| L | S | J values | L+S | Flavor symmetry | Representation |
|---|---|----------|-----|-----------------|----------------|
| 0 | 0 | 0        | 0   | antisymmetric   | 10-bar         |
| 0 | 1 | 1        | 1   | symmetric       | 15             |
| 1 | 0 | 1        | 1   | symmetric       | 15             |
| 1 | 1 | 0,1,2    | 2   | antisymmetric   | 10-bar         |
| 2 | 0 | 2        | 2   | antisymmetric   | 10-bar         |
| 2 | 1 | 1,2,3    | 3   | symmetric       | 15             |

**Theorem**: For J=0 diquarks, L must equal S (to get J = |L−S| = 0), so L+S = 2S
is always even, and flavor is always antisymmetric. J=0, color-3-bar diquarks
are **always** in the 10-bar. The 15 requires J ≥ 1.

## 5. SUSY chiral superfield composites

In SQCD, the baryon-like composite:

    B^{ab} = ε_{ij} Q^{ia} Q^{jb}   (N_c = 2)
    B^{abc} = ε_{ijk} Q^{ia} Q^{jb} Q^{kc}   (N_c = 3)

is automatically antisymmetric in flavor indices, due to the ε-tensor on color,
even though the superfields Q commute (no Pauli). The antisymmetry is algebraic,
not statistical.

To access the 15 (symmetric flavor), one would need the symmetric color
contraction (color 6 channel), but this cannot form color singlets.

**Meson alternative**: Mesons M^a_b = Q^a Q-bar_b live in 5 ⊗ 5-bar = 24 ⊕ 1.
The 15 vs 10-bar distinction is irrelevant for mesons — they live in a different
tensor space entirely.

## 6. Representation theory check

    5 ⊗ 5 = 15 ⊕ 10-bar:   15 + 10 = 25 = 5²   ✓
    5 ⊗ 5-bar = 24 ⊕ 1:    24 + 1 = 25 = 5²     ✓

The numerical coincidence 15 + 10 − 1 = 24 reflects:

    N(N+1)/2 + N(N−1)/2 − 1 = N² − 1

This is arithmetic (N²/2 + N²/2 − 1 = N² − 1), not a representation decomposition.
The 15 and 10-bar live in the symmetric/antisymmetric parts of 5⊗5, while the 24
lives in the traceless part of 5⊗5-bar. Different tensor spaces.

## 7. Summary: When can diquarks fill the 15?

**For J=0, color-3-bar diquarks: NEVER.**

The result is a theorem: J=0 forces L=S, making L+S even, which always gives
antisymmetric flavor (10-bar). This holds for both ordinary fermions (Pauli) and
SUSY composites (ε-tensor).

### Available loopholes

| # | Loophole | Status | Physical cost |
|---|----------|--------|---------------|
| 1 | Nonzero L (J ≥ 1) | Closed for J=0 | Must abandon scalar diquarks |
| 2 | Color-6 channel | Open | Exotic; no standard baryons |
| 3 | Broken SU(5) flavor | Open | Ad hoc; not an irrep |
| 4 | SUSY reinterpretation | Partial | Need color-6 or different composite |
| 5 | Relativistic corrections | Closed | Spin-statistics theorem is exact |
| 6 | Extended composites | Reduces to #1 | Higher L means higher mass/spin |

### Most natural resolution

The composites that fill the 15 should **not be identified as diquarks**. If the
theory's scalars arise as Seiberg dual mesons (living in the adjoint 24 of
SU(5)_flavor), or through some other dual description where the representation
constraint arises differently, the Pauli/ε-tensor conflict is entirely avoided.
The 15 can appear as part of a decomposition in a different context (e.g., under a
subgroup, or in a theory where the "composites" are fundamental fields of a dual
description that are not subject to the same antisymmetrization).
"""

with open("/home/codexssh/phys3/results/pauli_15_10bar.md", "w") as f:
    f.write(output)

print("\n" + "=" * 70)
print("Results written to results/pauli_15_10bar.md")
print("=" * 70)
