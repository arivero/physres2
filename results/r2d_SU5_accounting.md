# SU(5) Representation Accounting — Results

## CRITICAL CORRECTION: Charge Census

The assembly document had an **incorrect charge census**. Corrected:

| |Q| | 0 | 1/3 | 2/3 | 1 | 4/3 | 2 |
|-----|---|-----|-----|---|-----|---|
| WRONG (assembly) | 16 | 2×4 | 2×10 | 2×4 | 2×1 | — |
| CORRECT | **12** | **2×6** | **2×9** | **2×2** | **2×3** | **2×1** |

Total still 54 states. The |Q| = 2 sector was missing entirely.

### Full charge table by representation

| Q | −2 | −4/3 | −1 | −2/3 | −1/3 | 0 | +1/3 | +2/3 | +1 | +4/3 | +2 |
|---|-----|------|-----|------|------|---|------|------|-----|------|-----|
| 24 | 0 | 3 | 1 | 0 | 3 | 10 | 3 | 0 | 1 | 3 | 0 |
| 15 | 0 | 0 | 0 | 6 | 3 | 1 | 0 | 3 | 1 | 0 | 1 |
| 15̄ | 1 | 0 | 1 | 3 | 0 | 1 | 3 | 6 | 0 | 0 | 0 |
| **Total** | **1** | **3** | **2** | **9** | **6** | **12** | **6** | **9** | **2** | **3** | **1** |

## Neutral State Decomposition

**12 neutral states** (not 16):

| Source | Count | Identification |
|--------|-------|---------------|
| (8,1)₀ ⊂ 24 | 8 | SU(3)_C adjoint (gluon directions) |
| (1,3)₀ ⊂ 24, T₃=0 | 1 | Neutral SU(2)_L adjoint (W₃/Z-like) |
| (1,1)₀ ⊂ 24 | 1 | U(1)_Y generator |
| (1,3)₁ ⊂ 15, T₃=−1 | 1 | Color-singlet, electrically neutral scalar |
| (1,3̄)₋₁ ⊂ 15̄, T₃=+1 | 1 | Conjugate of above |
| **Total** | **12** | |

This resolves the 12-vs-16 discrepancy: it was an error in the census, not missing physics.
r² + s² − 1 = 9 + 4 − 1 = 12 = 4N matches exactly.

## Diophantine Equation Mapping

| Equation | Value | SU(5) states | Match? |
|----------|-------|-------------|--------|
| rs = 2N = 6 | 6 | Q=+1/3: 3 from (3̄,2)₅/₆⊂24 + 3 from (3̄,2)₋₁/₆⊂15̄ | **Exact** |
| r(r+1)/2 = 2N = 6 | 6 | Q=−2/3: 6 from (6,1)₋₂/₃⊂15 | **Partial** (total Q=−2/3 is 9, not 6) |
| r²+s²−1 = 4N = 12 | 12 | Q=0: all 12 neutral states | **Exact** |

Note: Eq 2 matches the **symmetric color-sextet** (6,1)₋₂/₃ exactly (dim 6 of SU(3) = r(r+1)/2 for r=3). The extra 3 states at Q=−2/3 come from (3̄,2)₋₁/₆⊂15̄ at T₃=−1/2, which are not symmetric down-down composites but rather mixed composites from the conjugate representation.

## Charge Operator

Q = Y + T₃ with (α, β) = (1, 1) confirmed correct with standard SU(5) hypercharge normalization.

## Anomaly Cancellation

| Rep | A(R) |
|-----|------|
| 24 (adjoint, real) | 0 |
| 15 (symmetric) | N+4 = 9 |
| 15̄ | −9 |
| **Total** | **0 ✓** |

24 ⊕ 15 ⊕ 15̄ is **anomaly-free**. Guaranteed by: real rep (24) + vector-like pair (15 ⊕ 15̄).

### Dynkin indices
- T(5) = 1/2, T(10) = 3/2, T(15) = 7/2, T(24) = 5
- T(24 ⊕ 15 ⊕ 15̄) = 5 + 7/2 + 7/2 = 12

## ±4/3 States: Correction

The assembly claimed "±4/3 states live exclusively in (1,3)₁ ⊂ 15" — this is **WRONG**.

Actual Q = ±4/3 states:
- Q = −4/3: 3 states from **(3,2)₋₅/₆ ⊂ 24** at T₃ = −1/2
- Q = +4/3: 3 states from **(3̄,2)₊₅/₆ ⊂ 24** at T₃ = +1/2

The ±4/3 states are **color triplets in the adjoint 24**, not in the 15.
The (1,3)₁ ⊂ 15 instead contains the Q = +2 state (at T₃ = +1).
