# Agent 6: SU(5) Charge Assignment — Results

## Branching Rules: SU(5) → SU(3) × SU(2) × U(1)

| SU(5) | Decomposition | Dim check |
|-------|---------------|-----------|
| **5** | (3,1)₋₁/₃ ⊕ (1,2)₁/₂ | 3+2=5 ✓ |
| **5̄** | (3̄,1)₁/₃ ⊕ (1,2)₋₁/₂ | 3+2=5 ✓ |
| **10** | (3̄,1)₋₂/₃ ⊕ (3,2)₁/₆ ⊕ (1,1)₁ | 3+6+1=10 ✓ |
| **10̄** | (3,1)₂/₃ ⊕ (3̄,2)₋₁/₆ ⊕ (1,1)₋₁ | 3+6+1=10 ✓ |
| **15** | (6,1)₋₂/₃ ⊕ (3,2)₁/₆ ⊕ (1,3)₁ | 6+6+3=15 ✓ |
| **15̄** | (6̄,1)₂/₃ ⊕ (3̄,2)₋₁/₆ ⊕ (1,3)₋₁ | 6+6+3=15 ✓ |
| **24** | (8,1)₀ ⊕ (3,2)₋₅/₆ ⊕ (3̄,2)₅/₆ ⊕ (1,3)₀ ⊕ (1,1)₀ | 8+6+6+3+1=24 ✓ |

## Tensor Product Cross-Checks

- **24** = **5** ⊗ **5̄** − **1**: Verified ✓
- **10** = ∧²(**5**): Verified ✓
- **15** = Sym²(**5**): Verified ✓

## Valid Charge Assignments Q = α Y₁ + β T₃

Constraint: charges in 24 ⊕ 15 ⊕ 15̄ must be exclusively {0, ±1/3, ±2/3, ±1, ±4/3}.

4 nontrivial solutions (related by conjugation):

| α | β | Produces |
|---|---|---------|
| +1 | +1/3 | {0, ±1/3, ±2/3, ±1, ±4/3} |
| +1 | −1/3 | same spectrum |
| −1 | +1/3 | same spectrum |
| −1 | −1/3 | same spectrum |

## Charge Census for 24 ⊕ 15 ⊕ 15̄ (54 states)

| Q | Multiplicity | Sources |
|---|-------------|---------|
| 0 | 16 | (8,1)₀ from 24 + (1,3)₀ from 24 + (1,1)₀ from 24 + (3,2) from 15 + (3̄,2) from 15̄ |
| ±1/3 | 4 each | (1,3)₀ from 24 + (3,2) from 15 / (3̄,2) from 15̄ |
| ±2/3 | 10 each | (3̄,2) from 24 + (6̄,1) from 15̄ + (1,3) from 15 / conjugates |
| ±1 | 4 each | (3,2) from 24 + (1,3) from 15 / conjugates |
| **±4/3** | **1 each** | **(1,3)₁ ⊂ 15 and (1,3)₋₁ ⊂ 15̄** |

Total: 16 + 8 + 20 + 8 + 2 = 54 ✓

## Key Finding: The ±4/3 States

The Q = ±4/3 states live **exclusively** in:
- Q = +4/3: (1, 3)₁ ⊂ **15**, at T₃ = +1
- Q = −4/3: (1, 3)₋₁ ⊂ **15̄**, at T₃ = −1

These are **SU(3) color singlets** and **SU(2) weak isospin triplets**.
The **24** adjoint contributes NO states with |Q| = 4/3.
There is exactly **1 state each** with Q = +4/3 and Q = −4/3.
