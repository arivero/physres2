# SU(6) Branching Rules

Branching rules for selected SU(6) representations under two subgroup chains.
All dimension checks are included. U(1) charges use the conventions stated in each section.

## Part 1: SU(6) → SU(5) × U(1)

**Convention.** Split index 6 from {1,…,5}. The U(1) generator is

    T = diag(1, 1, 1, 1, 1, −5) / √30

so each fundamental index carries charge +1 (indices 1–5) or −5 (index 6).
For a multi-index state the charge is the sum of individual charges.
Conjugate representations have all charges negated.

SU(5) antisymmetric tensor representations appearing below:

| k | ∧^k(5) | SU(5) label | dim |
|---|--------|-------------|-----|
| 0 | ∧^0(5) | 1 | 1 |
| 1 | ∧^1(5) | 5 | 5 |
| 2 | ∧^2(5) | 10 | 10 |
| 3 | ∧^3(5) | 10̄ | 10 |
| 4 | ∧^4(5) | 5̄ | 5 |
| 5 | ∧^5(5) | 1 | 1 |

### SU(6) rep **6** (∧¹(6) — fundamental, dim 6)  —  dimension check: ✓

| SU(5) rep | dim | U(1) charge |
|-----------|-----|-------------|
| 1 | 1 | -5 |
| 5 | 5 | +1 |
| **Total** | **6** | |

Contains SU(5) singlet(s) at U(1) charge(s): [-5].

### SU(6) rep **6̄** (∧⁵(6) — anti-fundamental, dim 6)  —  dimension check: ✓

| SU(5) rep | dim | U(1) charge |
|-----------|-----|-------------|
| 5̄ | 5 | -1 |
| 1 | 1 | +5 |
| **Total** | **6** | |

Contains SU(5) singlet(s) at U(1) charge(s): [5].

### SU(6) rep **15** (∧²(6) — antisymmetric 2-tensor, dim 15)  —  dimension check: ✓

| SU(5) rep | dim | U(1) charge |
|-----------|-----|-------------|
| 5 | 5 | -4 |
| 10 | 10 | +2 |
| **Total** | **15** | |

Contains no SU(5) singlets.

### SU(6) rep **20** (∧³(6) — antisymmetric 3-tensor, dim 20)  —  dimension check: ✓

| SU(5) rep | dim | U(1) charge |
|-----------|-----|-------------|
| 10 | 10 | -3 |
| 10̄ | 10 | +3 |
| **Total** | **20** | |

Contains no SU(5) singlets.

### SU(6) rep **21** (Sym²(6) — symmetric 2-tensor, dim 21)  —  dimension check: ✓

| SU(5) rep | dim | U(1) charge |
|-----------|-----|-------------|
| 15 | 15 | +2 |
| 5 | 5 | -4 |
| 1 | 1 | -10 |
| **Total** | **21** | |

Contains SU(5) singlet(s) at U(1) charge(s): [-10].

### SU(6) rep **35** (adjoint (6 ⊗ 6̄ − 1), dim 35)  —  dimension check: ✓

| SU(5) rep | dim | U(1) charge |
|-----------|-----|-------------|
| 24 | 24 | +0 |
| 1 | 1 | +0 |
| 5 | 5 | +6 |
| 5̄ | 5 | -6 |
| **Total** | **35** | |

Contains SU(5) singlet(s) at U(1) charge(s): [0].

### SU(6) rep **56** (Sym³(6) — symmetric 3-tensor, dim 56)  —  dimension check: ✓

| SU(5) rep | dim | U(1) charge |
|-----------|-----|-------------|
| 35(SU5) | 35 | +3 |
| 15(SU5) | 15 | -3 |
| 5 | 5 | -9 |
| 1 | 1 | -15 |
| **Total** | **56** | |

Contains SU(5) singlet(s) at U(1) charge(s): [-15].

## Part 2: SU(6) → SU(3)_A × SU(3)_B × U(1)

**Convention.** Split indices {1,2,3} → SU(3)_A, {4,5,6} → SU(3)_B.
The U(1) generator is

    T = diag(1, 1, 1, −1, −1, −1) / √6

so A-type indices carry charge +1 and B-type indices carry charge −1.

SU(3) representations that appear:

| label | description | dim |
|-------|-------------|-----|
| 1 | singlet | 1 |
| 3 | fundamental | 3 |
| 3̄ | anti-fundamental | 3 |
| 6 | symmetric 2-tensor | 6 |
| 8 | adjoint | 8 |
| 10 | symmetric 3-tensor | 10 |

### SU(6) rep **6** (∧¹(6) — fundamental, dim 6)  —  dimension check: ✓

| SU(3)_A rep | SU(3)_B rep | dim | U(1) charge |
|-------------|-------------|-----|-------------|
| 1 | 3 | 3 | -1 |
| 3 | 1 | 3 | +1 |
| | **Total** | **6** | |

Contains no SU(3)_A × SU(3)_B singlets.

### SU(6) rep **6̄** (∧⁵(6) — anti-fundamental, dim 6)  —  dimension check: ✓

| SU(3)_A rep | SU(3)_B rep | dim | U(1) charge |
|-------------|-------------|-----|-------------|
| 1 | 3 | 3 | +1 |
| 3 | 1 | 3 | -1 |
| | **Total** | **6** | |

Contains no SU(3)_A × SU(3)_B singlets.

### SU(6) rep **15** (∧²(6) — antisymmetric 2-tensor, dim 15)  —  dimension check: ✓

| SU(3)_A rep | SU(3)_B rep | dim | U(1) charge |
|-------------|-------------|-----|-------------|
| 1 | 3̄ | 3 | -2 |
| 3 | 3 | 9 | +0 |
| 3̄ | 1 | 3 | +2 |
| | **Total** | **15** | |

Contains no SU(3)_A × SU(3)_B singlets.

### SU(6) rep **20** (∧³(6) — antisymmetric 3-tensor, dim 20)  —  dimension check: ✓

| SU(3)_A rep | SU(3)_B rep | dim | U(1) charge |
|-------------|-------------|-----|-------------|
| 1 | 1 | 1 | -3 |
| 3 | 3̄ | 9 | -1 |
| 3̄ | 3 | 9 | +1 |
| 1 | 1 | 1 | +3 |
| | **Total** | **20** | |

Contains SU(3)_A × SU(3)_B singlet(s) at U(1) charge(s): [-3, 3].

### SU(6) rep **21** (Sym²(6) — symmetric 2-tensor, dim 21)  —  dimension check: ✓

| SU(3)_A rep | SU(3)_B rep | dim | U(1) charge |
|-------------|-------------|-----|-------------|
| 6 | 1 | 6 | +2 |
| 3 | 3 | 9 | +0 |
| 1 | 6 | 6 | -2 |
| | **Total** | **21** | |

Contains no SU(3)_A × SU(3)_B singlets.

### SU(6) rep **35** (adjoint (6 ⊗ 6̄ − 1), dim 35)  —  dimension check: ✓

| SU(3)_A rep | SU(3)_B rep | dim | U(1) charge |
|-------------|-------------|-----|-------------|
| 8 | 1 | 8 | +0 |
| 1 | 8 | 8 | +0 |
| 1 | 1 | 1 | +0 |
| 3 | 3̄ | 9 | +2 |
| 3̄ | 3 | 9 | -2 |
| | **Total** | **35** | |

Contains SU(3)_A × SU(3)_B singlet(s) at U(1) charge(s): [0].

### SU(6) rep **56** (Sym³(6) — symmetric 3-tensor, dim 56)  —  dimension check: ✓

| SU(3)_A rep | SU(3)_B rep | dim | U(1) charge |
|-------------|-------------|-----|-------------|
| 10 | 1 | 10 | +3 |
| 6 | 3 | 18 | +1 |
| 3 | 6 | 18 | -1 |
| 1 | 10 | 10 | -3 |
| | **Total** | **56** | |

Contains no SU(3)_A × SU(3)_B singlets.

## Part 3: Dimension Check Summary

| SU(6) rep | dim | Part 1 sum | Part 1 | Part 2 sum | Part 2 |
|-----------|-----|------------|--------|------------|--------|
| 6 | 6 | 6 | OK | 6 | OK |
| 6̄ | 6 | 6 | OK | 6 | OK |
| 15 | 15 | 15 | OK | 15 | OK |
| 20 | 20 | 20 | OK | 20 | OK |
| 21 | 21 | 21 | OK | 21 | OK |
| 35 | 35 | 35 | OK | 35 | OK |
| 56 | 56 | 56 | OK | 56 | OK |

## Part 4: Notable Features

### SU(6) → SU(5) × U(1)

| SU(6) rep | SU(5) reps present | Contains singlet? | U(1) charge spectrum |
|-----------|-------------------|-------------------|----------------------|
| 6 | 1, 5 | yes, q = [-5] | [-5, 1] |
| 6̄ | 5̄, 1 | yes, q = [5] | [-1, 5] |
| 15 | 5, 10 | no | [-4, 2] |
| 20 | 10, 10̄ | no | [-3, 3] |
| 21 | 15, 5, 1 | yes, q = [-10] | [-10, -4, 2] |
| 35 | 24, 1, 5, 5̄ | yes, q = [0] | [-6, 0, 6] |
| 56 | 35(SU5), 15(SU5), 5, 1 | yes, q = [-15] | [-15, -9, -3, 3] |

### SU(6) → SU(3)_A × SU(3)_B × U(1)

| SU(6) rep | (A,B) reps present | Contains (1,1) singlet? | U(1) charge spectrum |
|-----------|-------------------|------------------------|----------------------|
| 6 | (1,3), (3,1) | no | [-1, 1] |
| 6̄ | (1,3), (3,1) | no | [-1, 1] |
| 15 | (1,3̄), (3,3), (3̄,1) | no | [-2, 0, 2] |
| 20 | (1,1), (3,3̄), (3̄,3), (1,1) | yes, q = [-3, 3] | [-3, -1, 1, 3] |
| 21 | (6,1), (3,3), (1,6) | no | [-2, 0, 2] |
| 35 | (8,1), (1,8), (1,1), (3,3̄), (3̄,3) | yes, q = [0] | [-2, 0, 2] |
| 56 | (10,1), (6,3), (3,6), (1,10) | no | [-3, -1, 1, 3] |

### Physical interpretation notes

**6 = (5)(+1) + (1)(−5):** The fundamental splits into the SU(5) quintet plus one singlet with large negative U(1) charge.

**15 = (10)(+2) + (5)(−4):** The antisymmetric 2-tensor of SU(6) branches to the SU(5) decuplet at charge +2 and the quintet at charge −4. No SU(5) singlet. This is the representation relevant for antisymmetric diquark operators.

**20 = (10)(+3) + (10̄)(−3) :** The antisymmetric 3-tensor is self-conjugate as an SU(6) branching set: the SU(5) 10 and its conjugate 10̄ appear with opposite charges. No SU(5) singlet.

**21 = (15)(+2) + (5)(−4) + (1)(−10):** The symmetric 2-tensor contains one SU(5) singlet at charge −10 (high magnitude). The 15 of SU(5) (symmetric 2-tensor of SU(5)) appears at charge +2.

**35 = (24)(0) + (1)(0) + (5)(+6) + (5̄)(−6):** The adjoint decomposes into the SU(5) adjoint (24) and one singlet, both at zero U(1) charge, plus the quintet/antiquintet pair at ±6. The zero-charge singlet is the U(1) generator itself.

**56 = (35(SU5))(+3) + (15(SU5))(−3) + (5)(−9) + (1)(−15):** The symmetric cubic contains one SU(5) singlet at charge −15. The SU(5) representations 35 and 15 here are the symmetric 3- and 2-tensors of SU(5) (not the same as ∧²(5)=10 or ∧³(5)=10̄).

**SU(3)×SU(3) split of 20:** The antisymmetric 3-tensor decomposes as (1,1)+( 3,3̄)+(3̄,3)+(1,1) at charges ±3 and ±1 — wait, let us be precise.

- (1, 1), dim=1, q=-3
- (3, 3̄), dim=9, q=-1
- (3̄, 3), dim=9, q=+1
- (1, 1), dim=1, q=+3

The two (1,1) singlets in the 20 under SU(3)×SU(3) are the totally antisymmetric tensors ε_{abc} and ε_{def} acting as scalars on each SU(3) factor — the volume forms.

**SU(3)×SU(3) split of 35 (adjoint):** The decomposition (8,1)+(1,8)+(1,1)+(3,3̄)+(3̄,3) shows the two SU(3) sub-algebras as (8,1) and (1,8), the U(1) as (1,1), and the off-diagonal coset directions as (3,3̄)+(3̄,3). This is the standard pattern for a symmetric SU(3)×SU(3)×U(1) embedding.

---

*All branching rules derived by explicit index splitting. Dimension checks: sum of component dimensions equals SU(6) rep dimension in all cases.*