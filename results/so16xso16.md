# SO(16) x SO(16): Adjoint Decomposition and the Absence of (15, 3bar)

Complete analysis of whether the symmetric representation (15, 3bar) of SU(5) x SU(3) can appear in representations of SO(16) x SO(16), compared with SO(32) and E8 x E8.

---

## 1. Embedding SU(5) x SU(3) into a Single SO(16)

### 1.1 Why (5,3) cannot fit in the vector 16

The vector 16 of SO(16) is a **real** representation. The tensor product (5,3) of SU(5) x SU(3) has dimension 15, which almost fits (15 + 1 = 16). But (5,3) is **complex**: the 5 of SU(5) and 3 of SU(3) are both complex representations. A real representation of SO(16) must contain complex conjugate pairs, requiring (5,3) + (5bar,3bar) = 30 real dimensions. Since 30 > 16, the tensor product embedding is impossible.

**Dimensional obstruction:**
- SO(32): dim(vector) = 32 >= 30 + 2. Can fit (5,3) + (5bar,3bar) + 2 singlets. CHECK.
- SO(16): dim(vector) = 16 < 30. Cannot fit (5,3) + (5bar,3bar). FAILS.

This is the root cause of everything that follows.

### 1.2 The only available embedding: via SU(8)

The maximal SU(n) subgroup of SO(16) is SU(8), with 16 = 8 + 8bar.

Under SU(5) x SU(3) x U(1) in SU(8), the fundamental 8 decomposes as:

    8 = (5,1)_{+3} + (1,3)_{-5}

(U(1) charge chosen for tracelessness: 5(+3) + 3(-5) = 0.)

Therefore:

    16 = (5,1)_{+3} + (1,3)_{-5} + (5bar,1)_{-3} + (1,3bar)_{+5}

Dimension check: 5 + 3 + 5 + 3 = **16**. CHECK.

The 5 and 3 live in **separate direct summands** of the 8, not as a tensor product.

---

## 2. Task 1: Adjoint of a Single SO(16)

The adjoint is 120 = wedge^2(16). Label the four pieces of 16:

    A = (5,1),   B = (1,3),   Abar = (5bar,1),   Bbar = (1,3bar)

Then wedge^2(A + B + Abar + Bbar) = sum of wedge^2 and cross terms:

| Source | SU(5) x SU(3) | Dim | U(1) charge |
|--------|---------------|-----|-------------|
| wedge^2(A) | (10bar, 1) | 10 | +6 |
| wedge^2(B) | (1, 3bar) | 3 | -10 |
| wedge^2(Abar) | (10, 1) | 10 | -6 |
| wedge^2(Bbar) | (1, 3) | 3 | +10 |
| A tensor B | (5, 3) | 15 | -2 |
| A tensor Abar | (24, 1) + (1, 1) | 25 | 0 |
| A tensor Bbar | (5, 3bar) | 15 | +8 |
| B tensor Abar | (5bar, 3) | 15 | -8 |
| B tensor Bbar | (1, 8) + (1, 1) | 9 | 0 |
| Abar tensor Bbar | (5bar, 3bar) | 15 | -2 |
| **Total** | | **120** | |

Dimension check: 10 + 3 + 10 + 3 + 15 + 25 + 15 + 15 + 9 + 15 = **120** = C(16,2). CHECK.

### Result: (15, 3bar) is ABSENT from wedge^2(16) of SO(16).

The SU(5) representations appearing are: **10bar, 10, 24, 5, 5bar, 1**. Only antisymmetric tensors and their products. The symmetric tensor 15 = Sym^2(5) does not appear.

**Why**: The cross term A tensor B = (5,1) tensor (1,3) = (5,3), which is a rank-1 tensor product, not Sym^2(5) tensor anything. The identity wedge^2(V tensor W) = Sym^2(V) tensor wedge^2(W) + wedge^2(V) tensor Sym^2(W) requires V and W to be **factors of the same vector space**. Here 5 and 3 are separate summands of 8, not factors of a single 15-dimensional space.

---

## 3. Task 2: Split Embedding (SU(5) in first SO(16), SU(3) in second)

### 3.1 First factor: 120_1 under SU(5)

Embedding: SU(5) in SU(8) in SO(16)_1, with 8 = 5 + 1 + 1 + 1.

    16_1 = 5 + 5bar + 6 x (1)   [under SU(5)]

wedge^2(16_1) under SU(5):

| Source | SU(5) rep | Dim |
|--------|-----------|-----|
| wedge^2(5) | 10bar | 10 |
| wedge^2(5bar) | 10 | 10 |
| wedge^2(6 singlets) | 15 singlets | 15 |
| 5 tensor 5bar | 24 + 1 | 25 |
| 5 tensor 6 singlets | 6 copies of 5 | 30 |
| 5bar tensor 6 singlets | 6 copies of 5bar | 30 |
| **Total** | | **120** |

SU(5) reps appearing: **10bar, 10, 24, 1, 5, 5bar**. No 15.

### 3.2 Second factor: 120_2 under SU(3)

Embedding: SU(3) in SU(8) in SO(16)_2, with 8 = 3 + 1 + 1 + 1 + 1 + 1.

    16_2 = 3 + 3bar + 10 x (1)   [under SU(3)]

wedge^2(16_2) under SU(3):

| Source | SU(3) rep | Dim |
|--------|-----------|-----|
| wedge^2(3) | 3bar | 3 |
| wedge^2(3bar) | 3 | 3 |
| wedge^2(10 singlets) | 45 singlets | 45 |
| 3 tensor 3bar | 8 + 1 | 9 |
| 3 tensor 10 singlets | 10 copies of 3 | 30 |
| 3bar tensor 10 singlets | 10 copies of 3bar | 30 |
| **Total** | | **120** |

SU(3) reps appearing: **3bar, 3, 8, 1**. No 6.

### 3.3 The adjoint of SO(16) x SO(16)

    adj = (120_1, 1) + (1, 120_2)

Under SU(5) x SU(3):
- **(120_1, 1)** gives (R, 1) for each SU(5) rep R in 120_1: (10bar, 1), (10, 1), (24, 1), (5, 1), (5bar, 1), (1, 1)
- **(1, 120_2)** gives (1, R') for each SU(3) rep R' in 120_2: (1, 3bar), (1, 3), (1, 8), (1, 1)

**All representations are of the form (R, 1) or (1, R').**

This is a trivial consequence of the product structure: adj(G1 x G2) = adj(G1) + adj(G2), and each summand transforms nontrivially under only one factor. When SU(5) in G1 and SU(3) in G2, no mixed representation (R_SU5, R'_SU3) with both nontrivial can appear in the adjoint.

### Result: (15, 3bar) is ABSENT. In fact, no representation (R, R') with both R and R' nontrivial appears at all.

---

## 4. Task 3: Comparison Table

| Group | Contains (15, 3bar) in adjoint? | Mechanism | Key obstruction |
|-------|:-----:|-----------|-----------------|
| **SO(32)** | **YES** | 32 = (5,3) + (5bar,3bar) + 2 singlets. wedge^2(5 tensor 3) = (15, 3bar) + (10bar, 6) | None. The identity wedge^2(V tensor W) = Sym^2(V) tensor wedge^2(W) + ... generates 15 from 5 tensor 3. |
| **SO(16) x SO(16)** | **NO** | 5 and 3 live in separate factors (or separate summands within one factor). No single vector space has 5 tensor 3 structure. | **Dimensional**: (5,3) + (5bar,3bar) needs 30 real dimensions; the vector 16 has only 16. **Structural**: adjoint of product group is sum of adjoints, never mixed. |
| **E8 x E8** | **NO** | E8 reps decompose into only antisymmetric tensors of SU(5): 1, 5, 10, 10bar, 5bar, 24. Simply-laced roots exclude symmetric tensors. | **Algebraic**: E8 is simply-laced; all reps are exterior powers of fundamental representations of SU(n) subgroups. Sym^2(5) = 15 is structurally excluded. |

### The gradient from SO(32) to E8 x E8

The three cases form a hierarchy of obstructions:

1. **SO(32)**: flavor and color are **entangled** in the vector (5 tensor 3 = 15 states per slot). The exterior square identity does the rest.
2. **SO(16) x SO(16)**: flavor and color are **separated** (either in different factors, or in different summands of one factor). The tensor product structure 5 tensor 3 never appears in any vector representation, so the exterior square identity never fires.
3. **E8 x E8**: even if one could somehow embed SU(5) x SU(3), the representation theory of E8 forbids symmetric tensors entirely. This is a deeper algebraic obstruction beyond the dimensional one.

---

## 5. Task 4: Beyond the Adjoint -- All SO(16) x SO(16) Representations

### 5.1 Which SO(16) representations contain 15 of SU(5)?

For the embedding SU(5) in SU(8) in SO(16):

**Representations built from exterior algebra** (including spinors):
- Vector 16: contains 5, 5bar, 1. **No 15.**
- Adjoint 120 = wedge^2(16): contains 10bar, 10, 24, 5, 5bar, 1. **No 15.**
- wedge^3(16) = 560: only wedge^k(5) for various k. **No 15.**
- All wedge^k(16): only produce wedge^p(5) = 1, 5, 10bar, 10, 5bar. **No 15.**
- Spinor 128_s (even forms of SU(8)): decomposes as

        128_s = sum over even k of wedge^k(8)
              = 1 + wedge^2(8) + wedge^4(8) + wedge^6(8) + wedge^8(8)

    Under SU(5), each wedge^k(8) = wedge^k(5 + 1 + 1 + 1) contains only wedge^p(5) for p <= k. **No 15.**

- Spinor 128_c (odd forms): same analysis. **No 15.**

**The symmetric tensor representation**:

    Sym^2(16) = 136 = traceless 135 + trace 1

There are two cases depending on the embedding.

**Case A: SU(5) only (split embedding).** With 16 = 5 + 5bar + 6 singlets:

    Sym^2(16) contains Sym^2(5) = 15 of SU(5).

The 135 contains **(15)** as an SU(5) irrep. This is the smallest SO(16) representation containing the symmetric 15 of SU(5).

**Case B: SU(5) x SU(3) together via SU(8).** With 16 = 8 + 8bar, 8 = (5,1) + (1,3):

    Sym^2(16) = Sym^2(8) + 8 tensor 8bar + Sym^2(8bar)

Under SU(5) x SU(3):

| Source | Rep | Dim |
|--------|-----|-----|
| Sym^2(8) | **(15, 1)** | 15 |
| Sym^2(8) | (5, 3) | 15 |
| Sym^2(8) | (1, 6) | 6 |
| Sym^2(8bar) | (15bar, 1) | 15 |
| Sym^2(8bar) | (5bar, 3bar) | 15 |
| Sym^2(8bar) | (1, 6bar) | 6 |
| 8 tensor 8bar | (24, 1) + (1, 1) | 25 |
| 8 tensor 8bar | (5, 3bar) | 15 |
| 8 tensor 8bar | (5bar, 3) | 15 |
| 8 tensor 8bar | (1, 8) + (1, 1) | 9 |
| **Total** | | **136** |

The 15 of SU(5) appears as **(15, 1)** -- paired with a color **singlet**, not with 3bar. The trace subtraction (to get 135 from 136) removes one (1,1) from the 8 tensor 8bar sector and does not affect the (15, 1). So **the 135 contains (15, 1) but NOT (15, 3bar)**.

**The 1344 representation.** The product 16 tensor 120 (dim 1920) contains (15, 3bar) via the chain:

    16 has (5,1);  120 has (5, 3bar)  [from (5,1) tensor (1,3bar)]
    (5,1) tensor (5, 3bar) = (5 tensor 5, 3bar) = (15 + 10bar, 3bar)

The product 16 tensor 120 decomposes as 16 + 560 + 1344 in SO(16) irreps. The 560 = wedge^3(16) contains only exterior powers of 5, so (15, 3bar) is NOT in 560. Therefore **(15, 3bar) lives in the 1344** (the hook-tableau [2,1,0,...,0] irrep of SO(16)).

The 1344 is a mixed-symmetry tensor, not an exterior power, not a spinor, and not a symmetric tensor. It does not appear in any known string massless spectrum.

### 5.2 Does (15, 3bar) appear in any (R1, R2)?

For the **split embedding** (SU(5) in one factor, SU(3) in the other), irreps are (R1, R2). Under SU(5) x SU(3):

    (R1, R2)|_{SU(5) x SU(3)} = [R1|_{SU(5)}] x [R2|_{SU(3)}]

So (15, 3bar) appears in (R1, R2) iff:
- R1|_{SU(5)} contains 15, AND
- R2|_{SU(3)} contains 3bar.

From the analysis above:
- R1 must be the **135** (symmetric traceless 2-tensor) or higher symmetric power. Not an exterior power, not a spinor, not the adjoint.
- R2 can be the adjoint 120, spinors, or many other reps (3bar is common, being wedge^2(3)).

**Answer: (15, 3bar) CAN appear in (135, R2) of SO(16) x SO(16)** for any R2 containing 3bar of SU(3). But note that this is a **factored** product: (15, 1) x (1, 3bar) = (15, 3bar). The flavor and color are NOT entangled -- they arise from independent representations of separate gauge factors. This is qualitatively different from SO(32), where (15, 3bar) arises as a single irreducible piece of wedge^2(5 tensor 3), with flavor and color indices entangled within the same vector space.

For SU(5) x SU(3) **both in the same SO(16)**, the smallest irrep containing (15, 3bar) is the **1344**, a massive string state.

### 5.3 The bifundamental (16_1, 16_2)

    (16_1, 16_2) under SU(5) x SU(3) = (5 + 5bar + 6 singlets) x (3 + 3bar + 10 singlets)

This contains: (5,3), (5,3bar), (5bar,3), (5bar,3bar), (5,1), (5bar,1), (1,3), (1,3bar), (1,1).

**No (15, 3bar).** The 15 of SU(5) requires Sym^2(5), which is quadratic in the 5. A single copy of the fundamental 16_1 provides only linear pieces (5 and 5bar), never Sym^2(5).

### 5.4 Tensor products of bifundamentals

What about (16_1, 16_2) tensor (16_1, 16_2) = (16_1 tensor 16_1, 16_2 tensor 16_2)?

Under SO(16)_1: 16_1 tensor 16_1 = Sym^2(16_1) + wedge^2(16_1) = 136 + 120 = 256.
Under SO(16)_2: 16_2 tensor 16_2 = Sym^2(16_2) + wedge^2(16_2) = 136 + 120 = 256.

The Sym^2(16_1) piece contains 15 of SU(5), and the wedge^2(16_2) piece contains 3bar of SU(3).

So the tensor square **(16_1, 16_2)^{tensor 2}** does contain (15, 3bar) -- but this is a 256 x 256 = 65536 dimensional reducible representation. Its irreducible decomposition under SO(16) x SO(16) includes (135, 120), and this irrep contains (15, 3bar).

**But this is a massive string state**, not part of any known low-energy spectrum.

### 5.5 Physical spectrum of the SO(16) x SO(16) string

The non-supersymmetric SO(16) x SO(16) heterotic string (Alvarez-Gaume, Ginsparg, Moore, Vafa, 1986) has massless spectrum:

- **Bosonic**: graviton + B-field + dilaton + gauge bosons in (120, 1) + (1, 120)
- **Fermionic**: spinors (128_s, 1) + (1, 128_s) (or 128_c, depending on GSO)
- **Tachyonic**: the bifundamental (16, 16) is tachyonic in the non-compact theory (but can be stabilized by compactification)

The **135** does not appear in the massless spectrum. It arises only at massive string levels (mass of order the string scale). Therefore, (15, 3bar) is absent from the physically relevant low-energy spectrum of SO(16) x SO(16).

---

## 6. The Structural Reason: Entanglement of Flavor and Color

The central mechanism that produces (15, 3bar) in SO(32) is the identity:

    wedge^2(V tensor W) = [Sym^2(V) tensor wedge^2(W)] + [wedge^2(V) tensor Sym^2(W)]

This requires V and W to be **factors of the same vector space**, i.e., the fundamental representation must contain the tensor product V tensor W as a single irreducible piece.

In SO(32): 32 contains (5,3) = 5 tensor 3 as one block. The wedge^2 identity applies with V = 5, W = 3, and Sym^2(5) = 15 emerges paired with wedge^2(3) = 3bar.

In SO(16) x SO(16): No single 16-dimensional vector can contain 5 tensor 3 (which needs 30 real dimensions for the complex pair). The 5 and 3 can only appear as **separate summands** (in the same factor, via SU(8)) or in **separate factors**. In either case:
- Separate summands: the cross-term is A tensor B = (5,1) tensor (1,3) = (5,3), a rank-1 product, not Sym^2(5) tensor wedge^2(3).
- Separate factors: the adjoint is a direct sum, never mixed.

The identity wedge^2(V tensor W) fundamentally requires the flavor-color **entanglement** that only a sufficiently large single vector representation provides.

### Minimum group size

For the tensor product embedding (5,3) + (5bar,3bar) + singlets to fit in the vector representation of SO(2n), we need:

    2n >= 2 x dim(5,3) + 2 = 2 x 15 + 2 = 32

(The factor of 2 is for the complex conjugate pair; the +2 is for at least one real singlet pair.) Therefore **SO(32) is the minimum orthogonal group** that can accommodate the embedding, and it does so with exactly 0 dimensions to spare (modulo the 2 singlets needed for the SU(16) structure).

---

## 7. Summary

| | SO(32) | SO(16) x SO(16) (same factor) | SO(16) x SO(16) (split) | E8 x E8 |
|---|:---:|:---:|:---:|:---:|
| **Vector contains (5,3)?** | YES: 32 = (5,3)+... | NO: 16 too small | N/A | N/A |
| **(15,3bar) in adjoint?** | YES | NO | NO | NO |
| **(15,3bar) in spinors?** | N/A | NO | NO | NO |
| **(15,1) in Sym^2?** | N/A | YES (135) | YES (135) | NO |
| **(15,3bar) in any irrep?** | YES (adjoint 120) | YES (1344) | YES (factored: 135 x R2) | NO (never) |
| **Smallest rep with (15,3bar)** | 120 (adjoint) | 1344 (hook tableau) | (135, 120) = dim 16200 | Does not exist |
| **Obstruction level** | None | Dimensional (30 > 16) | Structural (product) + Dimensional | Algebraic (simply-laced) |
| **Physical relevance** | Massless gauge bosons | Massive strings only | Massive strings only | Excluded entirely |

### The hierarchy of obstructions

1. **SO(32)**: No obstruction. The wedge^2(V tensor W) identity generates (15, 3bar) automatically in the adjoint 496. This is the massless gauge boson spectrum.

2. **SO(16) x SO(16)**: The (15, 3bar) exists in principle, but in qualitatively different ways depending on the embedding:
   - **Split embedding** (SU(5) in one factor, SU(3) in other): appears in (135, R2) as a factored product (15,1) x (1,3bar), with no flavor-color entanglement. Requires the symmetric tensor 135, not available from exterior algebra or spinors.
   - **Same-factor embedding** (SU(5) x SU(3) both in one SO(16)): the 135 gives only (15, 1). Getting (15, 3bar) requires the 1344 (hook-tableau irrep), far larger.
   - In neither case does (15, 3bar) appear in the massless string spectrum. It requires massive string states at the string scale.

3. **E8 x E8**: The (15, 3bar) is **algebraically forbidden**. E8 is simply-laced: all representations decompose into only antisymmetric tensors of SU(n) subgroups. Sym^2(5) = 15 never appears in any E8 representation, at any mass level, in any embedding.

### The discriminant

The appearance of (15, 3bar) in the **massless** spectrum is unique to SO(32). This provides a sharp representation-theoretic criterion:

> The symmetric flavor representation 15 = Sym^2(5) of SU(5), paired with the antisymmetric color representation 3bar = wedge^2(3) of SU(3), appears in the massless gauge boson spectrum **only** for the Type I / SO(32) string. Neither SO(16) x SO(16) nor E8 x E8 can produce it at the massless level.

---

## 8. Dimension Cross-Checks

### 8.1 Single SO(16) factor

    C(16, 2) = 120   [adjoint = wedge^2(16)]
    10 + 3 + 10 + 3 + 15 + 25 + 15 + 15 + 9 + 15 = 120  CHECK

### 8.2 Split embedding, first factor

    10 + 10 + 15 + 25 + 30 + 30 = 120  CHECK

### 8.3 Split embedding, second factor

    3 + 3 + 45 + 9 + 30 + 30 = 120  CHECK

### 8.4 Sym^2(16)

    16 x 17 / 2 = 136  CHECK

Split embedding (SU(5) only):

    15 + 15 + 21 + 25 + 30 + 30 = 136  CHECK

Same-factor embedding (SU(5) x SU(3) via SU(8)):

    Sym^2(8) = 36:   15 + 15 + 6 = 36  CHECK
    Sym^2(8bar) = 36: 15 + 15 + 6 = 36  CHECK
    8 x 8bar = 64:    25 + 15 + 15 + 9 = 64  CHECK
    Total: 36 + 36 + 64 = 136  CHECK

### 8.5 Bifundamental

    16 x 16 = 256  CHECK
    15 + 15 + 50 + 15 + 15 + 50 + 18 + 18 + 60 = 256  CHECK

### 8.6 Spinor 128_s

    sum over even k of C(8,k) = 1 + 28 + 70 + 28 + 1 = 128  CHECK
    Each C(8,k) decomposes correctly under SU(5):
      C(8,0) = 1, C(8,2) = 28, C(8,4) = 70, C(8,6) = 28, C(8,8) = 1
    Total under SU(5): all wedge^p(5) with multiplicities from C(3,q), sum = 128  CHECK
