# Pauli Exclusion and the 15 vs 10-bar of SU(5) Flavor

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
