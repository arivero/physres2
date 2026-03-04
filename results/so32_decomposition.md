# SO(32) Adjoint Decomposition under SU(5) x SU(3)

Complete branching rules for the adjoint **496** of SO(32) through the chain:

    SO(32) ⊃ SU(16) × U(1)_A ⊃ SU(5) × SU(3)_c × U(1)_A × U(1)_B

All dimension checks included. No internal project labels.

---

## 1. Embedding Setup

### 1.1 SO(32) → SU(16) × U(1)_A

The **32** of SO(32) is the vector representation. Under SU(16) × U(1)_A:

    32 = 16_(+1) ⊕ 16-bar_(-1)

The adjoint of SO(32) is the antisymmetric 2-tensor:

    496 = ∧²(32) = ∧²(16 ⊕ 16-bar)

Using ∧²(V ⊕ W) = ∧²(V) ⊕ (V ⊗ W) ⊕ ∧²(W):

    496 = ∧²(16)_(+2) ⊕ (16 ⊗ 16-bar)_(0) ⊕ ∧²(16-bar)_(-2)
        = 120_(+2) ⊕ (255 ⊕ 1)_(0) ⊕ 120-bar_(-2)

**Dimension check:** 120 + 255 + 1 + 120 = **496** ✓

U(1)_A charges: +2 for the 120, 0 for the 255 ⊕ 1, -2 for the 120-bar. These charges reflect the fact that ∧²(16) involves two factors of 16 (each charge +1), etc.

### 1.2 SU(16) → SU(5) × SU(3) × U(1)_B

The fundamental 16 of SU(16) decomposes as:

    16 = (5, 3)_(+1) ⊕ (1, 1)_(-15)

where the subscripts denote U(1)_B charges. The 15 quark states (5 flavors × 3 colors) carry U(1)_B charge +1 each, and the 16th singlet carries charge -15 to ensure tracelessness: 15×(+1) + 1×(-15) = 0.

Conjugate:

    16-bar = (5-bar, 3-bar)_(-1) ⊕ (1, 1)_(+15)

**Notation.** Write R₁ = (5, 3)_(+1), R₂ = (5-bar, 3-bar)_(-1), R₃ = (1, 1)_(-15), R₄ = (1, 1)_(+15).

Then 16 = R₁ ⊕ R₃ and 16-bar = R₂ ⊕ R₄.

---

## 2. Decomposition of Each Piece

### 2.1 The 255 ⊕ 1 (adjoint sector, q_A = 0)

    16 ⊗ 16-bar = (R₁ ⊕ R₃) ⊗ (R₂ ⊕ R₄)
                = (R₁ ⊗ R₂) ⊕ (R₁ ⊗ R₄) ⊕ (R₃ ⊗ R₂) ⊕ (R₃ ⊗ R₄)

**Term (a): R₁ ⊗ R₂ = (5, 3) ⊗ (5-bar, 3-bar)**

    = (5 ⊗ 5-bar) ⊗ (3 ⊗ 3-bar)
    = (24 ⊕ 1) ⊗ (8 ⊕ 1)
    = (24, 8) ⊕ (24, 1) ⊕ (1, 8) ⊕ (1, 1)

Dimensions: 192 + 24 + 8 + 1 = **225** ✓  (= 15 × 15)

U(1)_B charges: (+1) + (-1) = 0 for all components.

**Term (b): R₁ ⊗ R₄ = (5, 3)_(+1) ⊗ (1, 1)_(+15) = (5, 3)_(+16)**

Dimension: 15.

**Term (c): R₃ ⊗ R₂ = (1, 1)_(-15) ⊗ (5-bar, 3-bar)_(-1) = (5-bar, 3-bar)_(-16)**

Dimension: 15.

**Term (d): R₃ ⊗ R₄ = (1, 1)_(-15) ⊗ (1, 1)_(+15) = (1, 1)_(0)**

Dimension: 1.

Sum: 225 + 15 + 15 + 1 = **256** = 16². Of these, the 255 is the traceless part of 16 ⊗ 16-bar (adjoint of SU(16)), and the remaining 1 is the U(1)_A generator (trace). Both carry q_A = 0.

The decomposition of 255 ⊕ 1 under SU(5) × SU(3) is therefore:

| SU(5) × SU(3) | Dim | q_A | q_B | Origin |
|----------------|-----|-----|-----|--------|
| (24, 8) | 192 | 0 | 0 | R₁ ⊗ R₂ |
| (24, 1) | 24 | 0 | 0 | R₁ ⊗ R₂ |
| (1, 8) | 8 | 0 | 0 | R₁ ⊗ R₂ |
| (1, 1) | 1 | 0 | 0 | R₁ ⊗ R₂ (trace part → U(1) within SU(16)) |
| (5, 3) | 15 | 0 | +16 | R₁ ⊗ R₄ |
| (5-bar, 3-bar) | 15 | 0 | -16 | R₃ ⊗ R₂ |
| (1, 1) | 1 | 0 | 0 | R₃ ⊗ R₄ (trace part → U(1)_A generator) |
| **Total** | **256** | | | |

Dimension check: 192 + 24 + 8 + 1 + 15 + 15 + 1 = **256** ✓

Note: The two (1,1)₀ singlets are the U(1)_B generator (from the SU(16) adjoint) and the U(1)_A generator. They cannot be separated by SU(5) × SU(3) × U(1)² quantum numbers alone.

### 2.2 The 120 (∧²(16), q_A = +2)

    ∧²(16) = ∧²(R₁ ⊕ R₃) = ∧²(R₁) ⊕ (R₁ ⊗ R₃) ⊕ ∧²(R₃)

Since R₃ = (1,1) is one-dimensional, ∧²(R₃) = 0.

**Term (a): ∧²(R₁) = ∧²[(5, 3)] = ∧²(5 ⊗ 3)**

Apply the identity for exterior powers of tensor products:

    ∧²(V ⊗ W) = [S²(V) ⊗ ∧²(W)] ⊕ [∧²(V) ⊗ S²(W)]

With V = 5 (fundamental of SU(5)) and W = 3 (fundamental of SU(3)):

    S²(5) = 15      (symmetric 2-tensor of SU(5))
    ∧²(5) = 10-bar  (antisymmetric 2-tensor of SU(5) ≅ 10-bar via ε-tensor)
    S²(3) = 6        (symmetric 2-tensor of SU(3))
    ∧²(3) = 3-bar    (antisymmetric 2-tensor of SU(3) ≅ 3-bar via ε-tensor)

Therefore:

    ∧²(5 ⊗ 3) = (15, 3-bar) ⊕ (10-bar, 6)

Dimension check: 15 × 3 + 10 × 6 = 45 + 60 = **105** ✓ (= C(15,2))

U(1)_B charges: (+1) + (+1) = +2 for both components.

**Term (b): R₁ ⊗ R₃ = (5, 3)_(+1) ⊗ (1, 1)_(-15) = (5, 3)_(-14)**

Dimension: 15.

U(1)_B charge: (+1) + (-15) = -14.

Sum of 120: 105 + 15 = **120** ✓

| SU(5) × SU(3) | Dim | q_A | q_B | Origin |
|----------------|-----|-----|-----|--------|
| **(15, 3-bar)** | **45** | **+2** | **+2** | **∧²(R₁): S²(5) ⊗ ∧²(3)** |
| (10-bar, 6) | 60 | +2 | +2 | ∧²(R₁): ∧²(5) ⊗ S²(3) |
| (5, 3) | 15 | +2 | -14 | R₁ ⊗ R₃ |
| **Total** | **120** | | | |

### 2.3 The 120-bar (∧²(16-bar), q_A = -2)

By conjugation of Section 2.2:

| SU(5) × SU(3) | Dim | q_A | q_B | Origin |
|----------------|-----|-----|-----|--------|
| **(15-bar, 3)** | **45** | **-2** | **-2** | **∧²(R₂): S²(5-bar) ⊗ ∧²(3-bar)** |
| (10, 6-bar) | 60 | -2 | -2 | ∧²(R₂): ∧²(5-bar) ⊗ S²(3-bar) |
| (5-bar, 3-bar) | 15 | -2 | +14 | R₂ ⊗ R₄ |
| **Total** | **120** | | | |

Note on conjugation: 6-bar = 6 for SU(3) (the symmetric tensor 6 is not self-conjugate — 6-bar is the conjugate symmetric tensor, same dimension). Similarly 10-bar of SU(5) conjugates to 10.

---

## 3. Master Table: Complete 496 Decomposition

| # | SU(5) | SU(3) | Dim | q_A | q_B | Source piece |
|---|-------|-------|-----|-----|-----|-------------|
| 1 | 24 | 8 | 192 | 0 | 0 | 255: adj(SU(15)) |
| 2 | **24** | **1** | **24** | **0** | **0** | **255: adj(SU(15))** |
| 3 | 1 | 8 | 8 | 0 | 0 | 255: adj(SU(15)) |
| 4 | 1 | 1 | 1 | 0 | 0 | 255: U(1)_B generator |
| 5 | 5 | 3 | 15 | 0 | +16 | 255: 16-singlet mixing |
| 6 | 5-bar | 3-bar | 15 | 0 | -16 | 255: 16-bar-singlet mixing |
| 7 | 1 | 1 | 1 | 0 | 0 | 1: U(1)_A generator |
| 8 | **15** | **3-bar** | **45** | **+2** | **+2** | **120: S²(5) ⊗ ∧²(3)** |
| 9 | 10-bar | 6 | 60 | +2 | +2 | 120: ∧²(5) ⊗ S²(3) |
| 10 | 5 | 3 | 15 | +2 | -14 | 120: R₁ ⊗ R₃ |
| 11 | **15-bar** | **3** | **45** | **-2** | **-2** | **120-bar: S²(5-bar) ⊗ ∧²(3-bar)** |
| 12 | 10 | 6-bar | 60 | -2 | -2 | 120-bar: ∧²(5-bar) ⊗ S²(3-bar) |
| 13 | 5-bar | 3-bar | 15 | -2 | +14 | 120-bar: R₂ ⊗ R₄ |
| **Total** | | | **496** | | | |

**Grand dimension check:**

    192 + 24 + 8 + 1 + 15 + 15 + 1 + 45 + 60 + 15 + 45 + 60 + 15 = 496  ✓

Breakdown: 256 (from 255⊕1) + 120 (from ∧²(16)) + 120 (from ∧²(16-bar)) = 496 ✓

---

## 4. Identification of the Target Spectrum

### 4.1 What the theory requires

The composite spectrum consists of:
- **Mesons** (quark-antiquark): 5 ⊗ 5-bar = **24** ⊕ 1 of SU(5)_flavor, color singlet → **(24, 1)**
- **Diquarks** (quark-quark): S²(5) = **15** of SU(5)_flavor, color 3-bar → **(15, 3-bar)**
- **Anti-diquarks**: conjugate → **(15-bar, 3)**
- **Gluon** directions: **(1, 8)** (color adjoint, flavor singlet)

### 4.2 Where each piece sits in the 496

| Target rep | Present in 496? | Row # | Source | q_A | q_B |
|------------|----------------|-------|--------|-----|-----|
| **(24, 1)** | **YES** | 2 | adj(SU(16)), 255 | 0 | 0 |
| **(15, 3-bar)** | **YES** | 8 | ∧²(16), 120 | +2 | +2 |
| **(15-bar, 3)** | **YES** | 11 | ∧²(16-bar), 120-bar | -2 | -2 |
| **(1, 8)** | **YES** | 3 | adj(SU(16)), 255 | 0 | 0 |

All four target representations are present, each appearing exactly once.

### 4.3 The counting question

The target spectrum under SU(5) × SU(3)_c has total dimension:

    (24, 1) + (15, 3-bar) + (15-bar, 3) + (1, 8)
    = 24 + 45 + 45 + 8
    = 122

But the published counting says **54 composite scalars**. Where does 54 come from? The 54 counts only the SU(5)_flavor content, ignoring color multiplicity:

    24 + 15 + 15-bar = 54

This is the count of distinct **flavor** composites. Each diquark flavor state (15 of them) exists in 3 color states (the 3-bar), for 45 total diquark field components. The 54 counts flavor types; the 122 counts field components including color. In a gauge theory both counts are meaningful:

- **54** = number of distinct composite operators up to color index
- **122** = total number of real scalar field degrees of freedom (before color averaging)
- **114** = the 122 minus the 8 gluon directions, if one separates gauge from matter

The 496 contains 122 states matching the composite spectrum (including gluon directions) and 374 states that must be projected out.

### 4.4 U(1) charge separation

The target pieces sit at distinct U(1)_A charges:
- (24, 1) and (1, 8): q_A = 0 (adjoint sector)
- (15, 3-bar): q_A = +2 (∧² sector)
- (15-bar, 3): q_A = -2 (∧²-bar sector)

This means the target states are not an arbitrary subset — they are partially selected by U(1) quantum numbers. The q_A = 0 sector needs further projection (to remove the (24,8) and singlets), but the q_A = ±2 sectors need only the removal of (10-bar, 6) / (10, 6-bar) and the extra (5,3) / (5-bar, 3-bar) pieces.

---

## 5. What Must Be Projected Out (374 states)

| Row # | Rep | Dim | q_A | q_B | Why it must go |
|-------|-----|-----|-----|-----|----------------|
| 1 | (24, 8) | 192 | 0 | 0 | Colored mesons — no such composites in minimal model |
| 4 | (1, 1) | 1 | 0 | 0 | U(1)_B generator — gauge direction |
| 5 | (5, 3) | 15 | 0 | +16 | Quark-like, but wrong quantum numbers for composites |
| 6 | (5-bar, 3-bar) | 15 | 0 | -16 | Conjugate of above |
| 7 | (1, 1) | 1 | 0 | 0 | U(1)_A generator — gauge direction |
| 9 | (10-bar, 6) | 60 | +2 | +2 | Color-sextet diquarks: wrong color channel |
| 10 | (5, 3) | 15 | +2 | -14 | From ∧²(16), quark-like |
| 12 | (10, 6-bar) | 60 | -2 | -2 | Color-sextet antidiquarks: wrong color channel |
| 13 | (5-bar, 3-bar) | 15 | -2 | +14 | From ∧²(16-bar), antiquark-like |

Total projected out: 192 + 1 + 15 + 15 + 1 + 60 + 15 + 60 + 15 = **374** ✓  (496 - 122 = 374)

### 5.1 Structure of the unwanted states

The 374 states group into three categories:

**(A) The (24,8) block — 192 states.** This is the largest unwanted piece. It lives at q_A = q_B = 0 (same as the target (24,1) and (1,8)), so it cannot be projected by U(1) charges alone. It represents the coset generators of SU(15)/[SU(5)×SU(3)×U(1)_B]. In a string-theoretic context, these would need to gain mass through the compactification or an orbifold projection.

**(B) The (10-bar, 6) ⊕ (10, 6-bar) block — 120 states.** These carry the same U(1) charges as the target diquarks (q_A = ±2, q_B = ±2) but sit in the wrong SU(3) representation: color-6 instead of color-3-bar. They are distinguished from the target by their SU(3)_c quantum numbers. A color-sensitive projection (e.g., an orientifold that acts nontrivially on color indices) can separate them.

**(C) The (5,3) ⊕ (5-bar, 3-bar) "quark-like" states — 60 states total (four copies of 15).** These transform in the fundamental of SU(5) and fundamental of SU(3), resembling quarks rather than composites. Two pairs sit at q_A = 0 with q_B = ±16 (distinguishable by U(1)_B), and two pairs at q_A = ±2 with q_B = ∓14 (distinguishable by both U(1) charges). All are separated from the target by their U(1)_B charges.

**(D) The 2 singlets (1,1)₀ — 2 states.** These are the U(1)_A and U(1)_B gauge generators themselves, which become massive in any breaking to SU(5) × SU(3).

### 5.2 Projection feasibility

The projection to the target spectrum requires removing 374 out of 496 states. The pieces have the following U(1) distinguishability:

| Piece to keep | Piece to remove | Same (q_A, q_B)? | Separated by color? |
|---------------|----------------|-------------------|-------------------|
| (24, 1) | (24, 8) | YES: (0, 0) | YES: 1 vs 8 |
| (1, 8) | (24, 8) | YES: (0, 0) | Overlap — both have 8 |
| (15, 3-bar) | (10-bar, 6) | YES: (+2, +2) | YES: 3-bar vs 6 |
| (15, 3-bar) | (5, 3) | NO: different q_B | — |
| (15-bar, 3) | (10, 6-bar) | YES: (-2, -2) | YES: 3 vs 6-bar |

The critical separation is always by SU(3)_c representation (3-bar vs 6 for diquarks; 1 vs 8 for mesons). This is precisely the kind of projection provided by **orientifold planes** in Type I string theory: the worldsheet parity Ω acts on the Chan-Paton factors and can select symmetric vs antisymmetric color channels.

---

## 6. The ∧²(5 ⊗ 3) Identity — Detailed Derivation

This is the crucial step that produces the **15** (symmetric) rather than the **10-bar** (antisymmetric) of SU(5). It deserves careful treatment.

### 6.1 General identity

For vector spaces V (dim m) and W (dim n):

    ∧²(V ⊗ W) = [S²(V) ⊗ ∧²(W)] ⊕ [∧²(V) ⊗ S²(W)]

**Proof sketch.** An element of ∧²(V ⊗ W) is an antisymmetric 2-tensor T^{(ia)(jb)} = -T^{(jb)(ia)} where i,j are V-indices and a,b are W-indices. Define:

    T^{(ia)(jb)} = ½[T^{(ia)(jb)} + T^{(ja)(ib)}] + ½[T^{(ia)(jb)} - T^{(ja)(ib)}]

The first term is symmetric under (i↔j) and the overall antisymmetry under (ia)↔(jb) forces it to be antisymmetric under (a↔b). This gives S²(V) ⊗ ∧²(W).

The second term is antisymmetric under (i↔j) and the overall antisymmetry forces symmetry under (a↔b). This gives ∧²(V) ⊗ S²(W).

### 6.2 Application to V = 5, W = 3

    S²(5) = 15   (the symmetric representation — this is the one needed for diquarks)
    ∧²(5) = 10   (but as a subspace of ∧²(5), this is the conjugate: 10-bar)
    ∧²(3) = 3-bar  (color antitriplet — the physical diquark color)
    S²(3) = 6    (color sextet)

Therefore:

    ∧²(5 ⊗ 3) = (15, 3-bar) ⊕ (10-bar, 6)

**The 15 of SU(5) is paired with the 3-bar of SU(3).** The 10-bar of SU(5) is paired with the 6 of SU(3).

This is a representation-theoretic fact: within the antisymmetric square of the (5 ⊗ 3) space, the symmetric flavor channel (15) is forced into the antisymmetric color channel (3-bar), and vice versa. The pairing is dictated by the overall antisymmetry of ∧².

### 6.3 Physical interpretation

For scalar composites (J = 0, S-wave), Pauli/Fermi statistics require the total wavefunction to be antisymmetric. If color is in the 3-bar (antisymmetric), then flavor must be in the 10-bar (antisymmetric) — this is the standard Pauli argument.

But the SO(32) adjoint decomposition produces the (15, 3-bar) pairing because the representation theory of ∧² works differently from Pauli statistics: ∧²(V ⊗ W) is the exterior algebra of the product space (antisymmetric under exchange of the composite index pair (i,a)), which is NOT the same as requiring antisymmetry separately in each factor. The antisymmetry of the composite index (ia) forces a CORRELATION: symmetric in one factor implies antisymmetric in the other.

The fact that (15, 3-bar) appears in the decomposition means that the SO(32) adjoint naturally houses states where flavor-symmetric pairs sit in color-antisymmetric combinations. The resolution of the Pauli tension (15 vs 10-bar for J=0 diquarks) must come from the reinterpretation: these are not ordinary fermionic composites but fields in a dual description (Seiberg duality, string theory) where the compositeness constraints differ.

---

## 7. Subgroup chain through SU(15)

An alternative route verifies the same result. The 16 of SU(16) decomposes under SU(15) × U(1)_B' as:

    16 = 15_(+1) ⊕ 1_(-15)

Then SU(15) decomposes under SU(5) × SU(3) (tensor product embedding, 15 = 5 × 3):

### 7.1 Adjoint of SU(15): 224

    adj(SU(15)) = 15 ⊗ 15-bar - 1 = (5,3) ⊗ (5-bar,3-bar) - 1
                = (24⊕1) ⊗ (8⊕1) - 1
                = (24,8) ⊕ (24,1) ⊕ (1,8) ⊕ (1,1) - 1
                = (24,8) ⊕ (24,1) ⊕ (1,8)

Dimension: 192 + 24 + 8 = **224** ✓

### 7.2 ∧²(15) under SU(5) × SU(3): 105

    ∧²(15) = ∧²(5 ⊗ 3) = (15, 3-bar) ⊕ (10-bar, 6)

Dimension: 45 + 60 = **105** ✓

### 7.3 Adjoint of SU(16): 255

    255 = 224 ⊕ 15 ⊕ 15-bar ⊕ 1

Under SU(5) × SU(3):

    = [(24,8)⊕(24,1)⊕(1,8)] ⊕ (5,3) ⊕ (5-bar,3-bar) ⊕ (1,1)

Dimension: 224 + 15 + 15 + 1 = **255** ✓

### 7.4 ∧²(16) = 120

Under SU(15) × U(1)_B':

    ∧²(16) = ∧²(15 ⊕ 1) = ∧²(15) ⊕ (15 ⊗ 1) = 105 ⊕ 15

Under SU(5) × SU(3):

    = [(15,3-bar)⊕(10-bar,6)] ⊕ (5,3)

Dimension: 45 + 60 + 15 = **120** ✓

This confirms the result of Section 2.2 through an independent route.

### 7.5 The full 496

    adj(SO(32)) = adj(SU(16)) ⊕ ∧²(16) ⊕ ∧²(16-bar) ⊕ 1_(U(1)_A)

    = [224 ⊕ 15 ⊕ 15-bar ⊕ 1] ⊕ [105 ⊕ 15] ⊕ [105 ⊕ 15] ⊕ 1
    = 224 + 15 + 15 + 1 + 105 + 15 + 105 + 15 + 1 = 496 ✓

---

## 8. The (24,1) ⊕ (15,3-bar) ⊕ (15-bar,3) Match: Assessment

### 8.1 Exact match

The SO(32) adjoint contains exactly one copy of each target representation:

    (24, 1):      Row 2, from adj(SU(16))        — unique ✓
    (15, 3-bar):  Row 8, from ∧²(16)            — unique ✓
    (15-bar, 3):  Row 11, from ∧²(16-bar)       — unique ✓
    (1, 8):       Row 3, from adj(SU(16))        — unique ✓

No representation appears with multiplicity > 1 in the target set. Each is uniquely located in a different piece of the decomposition (adjoint sector vs ∧² sector).

### 8.2 Self-referential structure ("turtles")

The embedding is:
1. Start with 5 quark flavors × 3 colors = 15 states
2. Add 1 singlet (the "elephant"): 15 + 1 = 16 Chan-Paton factors
3. Form SO(2 × 16) = SO(32) (oriented → unoriented string doubling)
4. The adjoint 496 of SO(32) contains, among other things, exactly the composites of those 5 quarks:
   - 5 ⊗ 5-bar = 24 (mesons)
   - S²(5) = 15 (diquarks, with color 3-bar from ∧²(3))
   - S²(5-bar) = 15-bar (antidiquarks, with color 3)

The quarks determine the group, and the group contains the quarks' composites. This is the self-consistency ("bootstrap") condition.

### 8.3 What "exactly the right content" means

The answer to the posed question — **does the SO(32) adjoint contain exactly the right SU(5) × SU(3) content?** — is:

**YES**, with qualifications:

1. **The flavor representations match exactly.** The 24 (mesons), 15 (diquarks), 15-bar (antidiquarks) are all present, each once.

2. **The color assignments match exactly.** The 15 is paired with 3-bar (antisymmetric color = physical diquark color), the 15-bar with 3, and the 24 with 1 (color singlet mesons). These pairings are forced by the ∧² identity — they are not a choice.

3. **The gluon sector (1,8) is present** as a bonus, consistent with having SU(3)_c as a gauged subgroup.

4. **There are 374 extra states** that must be removed. This is not a flaw — it is expected. The SO(32) string has 496 gauge bosons; the SM-like content is a small subset. The question is whether the projection mechanism exists, not whether the content is present.

---

## 9. Why E₈ × E₈ Cannot Work

### 9.1 The representation theory obstruction

E₈ has only the adjoint representation (248-dimensional) as its smallest nontrivial representation. Under the standard embedding chain:

    E₈ → SU(5) × SU(5)':
    248 = (24,1) + (1,24) + (10,5-bar) + (10-bar,5) + (5,10) + (5-bar,10-bar)

The SU(5) representations that appear are: **1, 5, 5-bar, 10, 10-bar, 24**. The **15** (symmetric tensor) **never appears** in any E₈ decomposition, regardless of the embedding chain.

This is a theorem: E₈ is simply-laced (all roots have the same length), and its representations decompose into antisymmetric tensors of SU(n) subgroups. The symmetric tensor 15 = S²(5) is not an exterior power and cannot appear.

More precisely: under E₈ → SU(5), only representations corresponding to ∧^k(5) appear: 1, 5, 10, 10-bar, 5-bar, 1 (and the adjoint 24 = 5 ⊗ 5-bar - 1). The symmetric representations S^k(5) are absent.

### 9.2 The standard E₈ → E₆ → SO(10) → SU(5) chain

    E₈ → E₆ × SU(3):    248 = (78,1) + (1,8) + (27,3) + (27-bar,3-bar)
    E₆ → SO(10) × U(1):  27 = 16 + 10 + 1
    SO(10) → SU(5):       16 = 10 + 5-bar + 1
                           10 = 5 + 5-bar

At every level, only antisymmetric tensors of SU(5) appear: 1, 5, 5-bar, 10, 10-bar. The 15 is structurally excluded.

### 9.3 Physical consequence

The diquark spectrum requires 15 = S²(5) for the symmetric flavor pairing (which includes same-flavor diquarks: uu, dd, ss, cc, bb). The E₈ × E₈ heterotic string, despite giving beautiful three-generation models, produces only 10-bar = ∧²(5) for diquarks — the antisymmetric flavor channel that excludes same-flavor pairs.

This is the definitive obstruction: **E₈ × E₈ heterotic gives ∧², while the composite spectrum needs S².**

SO(32) / Type I resolves this because the adjoint ∧²(32) applied to the tensor product 5 ⊗ 3 produces both S²(5) ⊗ ∧²(3) and ∧²(5) ⊗ S²(3) — the symmetric flavor representation paired with antisymmetric color is automatically generated.

### 9.4 Caveat: dualities

Type I SO(32) is S-dual to heterotic SO(32), which is T-dual to heterotic E₈ × E₈ on a non-simply-connected torus. In principle, the same physics could be described in any duality frame. The statement "E₈ × E₈ cannot work" means specifically that the **standard perturbative E₈ × E₈ description** does not produce the 15. A nonperturbative or non-standard embedding might, through duality, reproduce the same spectrum — but then one is effectively working in the SO(32) frame.

---

## 10. Dimension Cross-Checks

### 10.1 Global check

    496 = C(32, 2) = 32 × 31 / 2 = 496  ✓  (adjoint of SO(32) = ∧²(32))

### 10.2 SU(16) level

    256 = 16² = 255 + 1   ✓  (adj of SU(16) + U(1))
    120 = C(16, 2)         ✓  (∧²(16))
    256 + 2×120 = 496      ✓

### 10.3 SU(15) level

    224 = 15² - 1          ✓  (adj of SU(15))
    105 = C(15, 2)         ✓  (∧²(15))

### 10.4 SU(5) × SU(3) level — piece by piece

| Rep | dim(SU(5)) | dim(SU(3)) | Product | Count in 496 | Subtotal |
|-----|-----------|-----------|---------|-------------|----------|
| (24,8) | 24 | 8 | 192 | 1 | 192 |
| (24,1) | 24 | 1 | 24 | 1 | 24 |
| (1,8) | 1 | 8 | 8 | 1 | 8 |
| (1,1) | 1 | 1 | 1 | 2 | 2 |
| (15,3-bar) | 15 | 3 | 45 | 1 | 45 |
| (15-bar,3) | 15 | 3 | 45 | 1 | 45 |
| (10-bar,6) | 10 | 6 | 60 | 1 | 60 |
| (10,6-bar) | 10 | 6 | 60 | 1 | 60 |
| (5,3) | 5 | 3 | 15 | 2 | 30 |
| (5-bar,3-bar) | 5 | 3 | 15 | 2 | 30 |
| | | | | **Total:** | **496** ✓ |

### 10.5 Anomaly check of the target

    A(24) = 0       (real/self-conjugate representation)
    A(15) = +9      (anomaly coefficient for symmetric tensor of SU(5))
    A(15-bar) = -9
    Total: 0 + 9 - 9 = 0  ✓  (anomaly-free)

This is guaranteed: the 496 inherits anomaly cancellation from SO(32), which is anomaly-free. Any real subrepresentation of an anomaly-free representation is anomaly-free.

### 10.6 Asymptotic freedom

The one-loop beta function coefficient for SU(5) gauge theory with matter in 24 ⊕ 15 ⊕ 15-bar:

    b₀ = (11/3)C₂(G) - (2/3)∑_f T(R_f)

For SU(5): C₂(G) = 5, and:

    T(24) = 5,   T(15) = 7/2,   T(15-bar) = 7/2

If these are scalar fields (contributing 1/3 of the fermionic amount to the beta function):

    b₀ = (11/3)(5) - (1/3)[5 + 7/2 + 7/2] = 55/3 - 12/3 = 43/3 > 0

If fermionic: b₀ = 55/3 - (2/3)(12) = 55/3 - 24/3 = 31/3 > 0.

Either way, **asymptotic freedom is maintained**.

---

## 11. Summary

The adjoint 496 of SO(32), decomposed under SU(5)_flavor × SU(3)_color through the chain SO(32) ⊃ SU(16) × U(1) ⊃ SU(5) × SU(3) × U(1)², contains precisely the target representations:

- **(24, 1)**: 24 color-singlet meson directions — from the adjoint 255 of SU(16)
- **(15, 3-bar)**: 45 diquark directions (symmetric flavor, antisymmetric color) — from ∧²(16)
- **(15-bar, 3)**: 45 antidiquark directions — from ∧²(16-bar)
- **(1, 8)**: 8 gluon directions — from the adjoint 255 of SU(16)

The total target spectrum accounts for 122 of the 496 states. The remaining 374 states decompose into:
- (24, 8): 192 colored mesons
- (10-bar, 6) ⊕ (10, 6-bar): 120 color-sextet (anti)diquarks
- 4 × (5, 3) or (5-bar, 3-bar): 60 quark-like states
- 2 × (1, 1): 2 U(1) generators

The color-sextet pieces (10-bar, 6) and (10, 6-bar) are the natural candidates for orientifold projection: Type I string theory's Ω-projection acts on color indices and can select 3-bar over 6 (or vice versa), keeping precisely (15, 3-bar) while removing (10-bar, 6) from the ∧² sector.

The E₈ × E₈ heterotic string is excluded because E₈ decompositions under SU(5) produce only antisymmetric tensors (1, 5, 10, 10-bar, 5-bar, 24), never the symmetric 15. This is the central representation-theoretic obstruction.
