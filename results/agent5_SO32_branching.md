# Agent 5: SO(32) Adjoint Branching Rules — Results

## Embedding Chain

SO(32) ⊃ SU(16) × U(1)_A ⊃ SU(15) × U(1)_A × U(1)_B ⊃ SU(5) × SU(3) × U(1)_A × U(1)_B

Key embedding: **15** of SU(15) → **(5, 3)** of SU(5) × SU(3) (tensor product embedding)

## Step-by-Step Decomposition

### SO(32) → SU(16) × U(1)_A

496 → (120)₊₂ + (255)₀ + (1)₀ + (120̄)₋₂

Dimension check: 120 + 255 + 1 + 120 = 496 ✓

### SU(16) → SU(15) × U(1)_B

- 255 → 224₀ + 15₊₁₆ + 15̄₋₁₆ + 1₀
- 120 → 105₊₂ + 15₋₁₄
- 120̄ → 105̄₋₂ + 15̄₊₁₄

### SU(15) → SU(5) × SU(3) (tensor product)

Using ∧²(V⊗W) = (∧²V ⊗ S²W) + (S²V ⊗ ∧²W):

- 224 → (24, 8) + (24, 1) + (1, 8)
- 105 → (10, 6) + (15, 3̄)
- 105̄ → (10̄, 6̄) + (15̄, 3)
- 15 → (5, 3)
- 15̄ → (5̄, 3̄)

## Complete Decomposition: 496 under SU(5) × SU(3) × U(1)²

| SU(5) | SU(3) | dim | (q_A, q_B) | Origin |
|--------|--------|-----|------------|--------|
| **24** | **8** | 192 | (0, 0) | adj of SU(15) |
| **24** | **1** | 24 | (0, 0) | adj of SU(15) |
| **1** | **8** | 8 | (0, 0) | adj of SU(15) |
| **5** | **3** | 15 | (0, +16) | adj of SU(16) |
| **5̄** | **3̄** | 15 | (0, −16) | adj of SU(16) |
| **1** | **1** | 1 | (0, 0) | adj of SU(16) |
| **1** | **1** | 1 | (0, 0) | U(1)_A generator |
| **10** | **6** | 60 | (+2, +2) | ∧²(16) |
| **15** | **3̄** | 45 | (+2, +2) | ∧²(16) |
| **5** | **3** | 15 | (+2, −14) | ∧²(16) |
| **10̄** | **6̄** | 60 | (−2, −2) | ∧²(16̄) |
| **15̄** | **3** | 45 | (−2, −2) | ∧²(16̄) |
| **5̄** | **3̄** | 15 | (−2, +14) | ∧²(16̄) |

**Total: 192 + 24 + 8 + 15 + 15 + 1 + 1 + 60 + 45 + 15 + 60 + 45 + 15 = 496 ✓**

## Self-Check

Two independent computational routes (via SU(15) intermediate, and direct SU(16) → SU(5) × SU(3)) produce **identical decomposition tables**. All dimension checks pass at every level.

## Key Identities Used

1. ∧²(V ⊕ W) = ∧²V ⊕ (V ⊗ W) ⊕ ∧²W
2. ∧²(V ⊗ W) = (∧²V ⊗ S²W) ⊕ (S²V ⊗ ∧²W)
3. V ⊗ V* = adj ⊕ 1 for SU(n)
4. SO(2n) ⊃ SU(n) × U(1): **2n** → **n**₊₁ + **n̄**₋₁

## Notable Content

- The **(15, 3̄)** and **(15̄, 3)** pieces (45 + 45 = 90 states) carry SU(3) color triplet/antitriplet structure
- The **(10, 6)** and **(10̄, 6̄)** pieces (60 + 60 = 120 states) carry SU(3) sextet structure
- The adjoint sector (24,8) + (24,1) + (1,8) + 2×(1,1) has 226 states at zero U(1) charges
