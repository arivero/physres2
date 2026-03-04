# SO(32) Decomposition with U(1) Charges — Results

## Embedding Chain

**SO(32) → SU(16) × U(1)_A → SU(5) × SU(3) × U(1)_A × U(1)_B**

Fundamental: 32 = (5,3) ⊕ (5̄,3̄) ⊕ (1,1) ⊕ (1,1)
Dim check: 15 + 15 + 1 + 1 = 32 ✓

The 16 of SU(16) decomposes as (5,3) + (1,1) under SU(5) × SU(3), with U(1)_B from the traceless generator in SU(16) Cartan orthogonal to SU(5) × SU(3).

## Adjoint Decomposition: ∧²(32) = 496

| SU(5)×SU(3) rep | Dim | Q_A | Q_B | Origin |
|---|---|---|---|---|
| (24, 8) | 192 | 0 | 0 | R₁ ⊗ R₂ |
| (24, 1) | 24 | 0 | 0 | R₁ ⊗ R₂ |
| (1, 8) | 8 | 0 | 0 | R₁ ⊗ R₂ |
| (1, 1) | 1 | 0 | 0 | R₁ ⊗ R₂ |
| (1, 1) | 1 | 0 | 0 | R₃ ⊗ R₄ |
| (10, 6) | 60 | +2 | +2 | ∧²(R₁) |
| (10̄, 6̄) | 60 | −2 | −2 | ∧²(R₂) |
| (15, 3̄) | 45 | +2 | +2 | ∧²(R₁) |
| (15̄, 3) | 45 | −2 | −2 | ∧²(R₂) |
| (5, 3)₊ | 15 | +2 | −14 | R₁ ⊗ R₃ |
| (5, 3)₋ | 15 | 0 | +16 | R₁ ⊗ R₄ |
| (5̄, 3̄)₊ | 15 | −2 | +14 | R₂ ⊗ R₃ |
| (5̄, 3̄)₋ | 15 | 0 | −16 | R₂ ⊗ R₄ |
| **Total** | **496** | | | |

Dim check: 192+24+8+1+1+60+60+45+45+15+15+15+15 = 496 ✓

## Derivation: ∧²(R₁ ⊕ R₂ ⊕ R₃ ⊕ R₄)

∧²(5,3) = [∧²(5) ⊗ S²(3)] ⊕ [S²(5) ⊗ ∧²(3)] = (10,6) ⊕ (15,3̄)
  Dim: 60 + 45 = 105 = C(15,2) ✓

(5,3)⊗(5̄,3̄) = (5⊗5̄)⊗(3⊗3̄) = (24⊕1)⊗(8⊕1) = (24,8)⊕(24,1)⊕(1,8)⊕(1,1)
  Dim: 192+24+8+1 = 225 ✓

## Adjoint Sector (Q_A = Q_B = 0)

(24,8) + (24,1) + (1,8) + 2×(1,1) = 192 + 24 + 8 + 2 = **226 states**

Of these:
- 34 = dim[SU(5) × SU(3) × U(1)² adjoint]
- 192 = (24,8) are gauge bosons in coset SU(15)/[SU(5)×SU(3)×U(1)_B]

Note: (24,1) + (1,8) + (1,1) + (24,8) = 225 = dim[SU(15) adjoint] + 1

## Electric Charges (from SU(5) → SU(3)_c × SU(2)_L × U(1)_Y)

Charges come from SU(5) factor only. SU(3)_horiz provides multiplicity.

Most exotic charges:
- **|Q| = 4/3**: from (3,2)₋₅/₆ ⊂ 24 and (3̄,2)₊₅/₆ ⊂ 24, multiplied by SU(3)_horiz
- **|Q| = 2**: from (1,3)₊₁ ⊂ 15 and (1,3̄)₋₁ ⊂ 15̄, multiplied by SU(3)_horiz

## State Count Outside Adjoint Sector

496 − 226 = **270 "matter" states** carrying nonzero U(1) charges:
- (10,6) ⊕ (10̄,6̄) = 120 at |Q_A| = 2
- (15,3̄) ⊕ (15̄,3) = 90 at |Q_A| = 2
- 2×(5,3) ⊕ 2×(5̄,3̄) = 60 at mixed Q_A

These 270 = 496 − 226 states would need to be projected out or made massive in a realistic model.
The previous initial count of "291 extra states" was incorrect — the actual split is 226 adjoint + 270 charged.

Actually: if we identify the SM gauge group within SU(5) × SU(3), the SM adjoint is only 12-dimensional (8 gluons + W⁺W⁻W³ + B). So the "extra" count relative to SM is 496 − 12 = 484, of which 214 are in the extended adjoint sector and 270 carry U(1) charges.
