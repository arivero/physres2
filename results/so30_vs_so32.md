# SO(30) vs SO(32): When Does the Symmetric 15 Appear?

Comparison of adjoint decompositions under SU(5) x SU(3) for different embeddings. The question: does the symmetric representation 15 = Sym^2(5) appear in wedge^2, and does this depend on the ambient SO group or only on the embedding structure?

---

## 1. SO(30) Adjoint Decomposition (Task 1)

### Setup

The adjoint of SO(30) is wedge^2(30), with dimension C(30,2) = 435.

Embedding: SO(30) contains SU(15) x U(1), with the vector decomposing as

    30 = 15_{+1} + 15-bar_{-1}

where 15 of SU(15) = (5,3) of SU(5) x SU(3) via the tensor product embedding.

The adjoint decomposes as

    adj(SO(30)) = adj(SU(15)) + 1_{U(1)} + wedge^2(15) + wedge^2(15-bar)
               = 224 + 1 + 105 + 105 = 435

Dimension check: (15^2 - 1) + 1 + C(15,2) + C(15,2) = 224 + 1 + 105 + 105 = 435 = C(30,2).

### Piece 1: wedge^2((5,3))

Applying the fundamental identity

    wedge^2(V tensor W) = [Sym^2(V) tensor wedge^2(W)] + [wedge^2(V) tensor Sym^2(W)]

with V = 5, W = 3:

    Sym^2(5) = 15,  wedge^2(3) = 3-bar,  wedge^2(5) = 10-bar,  Sym^2(3) = 6

Therefore:

    wedge^2((5,3)) = (15, 3-bar) + (10-bar, 6)

Dimensions: 15 x 3 + 10 x 6 = 45 + 60 = 105 = C(15,2).

### Piece 2: (5,3) tensor (5-bar,3-bar) = adj(SU(15)) + 1

    (5 tensor 5-bar) x (3 tensor 3-bar) = (24 + 1) x (8 + 1)
                                         = (24,8) + (24,1) + (1,8) + (1,1)

Dimensions: 192 + 24 + 8 + 1 = 225 = 15^2. The adj(SU(15)) is the 224-dimensional traceless part; the (1,1) is the U(1) within SU(15) x U(1) inside SU(16).

### Piece 3: wedge^2((5-bar,3-bar))

By conjugation of Piece 1:

    wedge^2((5-bar,3-bar)) = (15-bar, 3) + (10, 6-bar)

Dimensions: 45 + 60 = 105.

### Complete SO(30) decomposition

| SU(5) x SU(3) | Dim | Source |
|----------------|-----|--------|
| (24, 8) | 192 | adj(SU(15)) |
| (24, 1) | 24 | adj(SU(15)) |
| (1, 8) | 8 | adj(SU(15)) |
| (1, 1) | 1 | U(1) generator |
| **(15, 3-bar)** | **45** | **wedge^2(15)** |
| (10-bar, 6) | 60 | wedge^2(15) |
| **(15-bar, 3)** | **45** | **wedge^2(15-bar)** |
| (10, 6-bar) | 60 | wedge^2(15-bar) |
| **Total** | **435** | |

**Result: YES, the (15, 3-bar) appears in SO(30), identically to SO(32).**

The difference from SO(32) is purely in the adjoint sector: SO(32) has a 256-dimensional piece (255 + 1) from 16 x 16-bar, including extra (5,3) and (5-bar,3-bar) states mixed with the singlets. SO(30) has a 225-dimensional piece from 15 x 15-bar. The wedge^2 sector is identical.

---

## 2. SO(10) Adjoint Decomposition (Task 2)

### Setup

adj(SO(10)) = wedge^2(10), dim = C(10,2) = 45.

Embedding: 10 = (5,1) + (5-bar,1) under SU(5) x {trivial}. No color tensor product.

    wedge^2(5 + 5-bar) = wedge^2(5) + (5 x 5-bar) + wedge^2(5-bar)

### Piece by piece

    wedge^2(5) = 10-bar     (dim 10)
    5 x 5-bar = 24 + 1      (dim 25)
    wedge^2(5-bar) = 10      (dim 10)
    Total: 10 + 25 + 10 = 45

### Complete SO(10) decomposition

| SU(5) | Dim | Source |
|-------|-----|--------|
| 24 | 24 | adj(SU(5)) |
| 1 | 1 | U(1) |
| **10-bar** | **10** | **wedge^2(5)** |
| **10** | **10** | **wedge^2(5-bar)** |
| **Total** | **45** | |

**Result: NO 15 appears. Only the antisymmetric 10-bar = wedge^2(5).**

Without the tensor product structure 5 tensor 3, the identity that produces Sym^2(5) = 15 is never activated. Plain wedge^2(5) gives only the antisymmetric representation.

---

## 3. Why the Tensor Product Structure Matters (Task 3)

The key identity is

    wedge^2(V tensor W) = [Sym^2(V) tensor wedge^2(W)] + [wedge^2(V) tensor Sym^2(W)]

and its companion

    Sym^2(V tensor W) = [Sym^2(V) tensor Sym^2(W)] + [wedge^2(V) tensor wedge^2(W)]

### Proof (index-level)

An element of (V tensor W)^{tensor 2} carries indices T^{ia,jb} where i,j are V-indices and a,b are W-indices. Decompose under the exchange (ia) <-> (jb):

    T^{ia,jb} = 1/2 [T^{ia,jb} + T^{ib,ja}] + 1/2 [T^{ia,jb} - T^{ib,ja}]
                         (i<->j symmetric)            (i<->j antisymmetric)

For the ANTISYMMETRIC part of the full composite index (ia):

    T^{ia,jb} = -T^{jb,ia}  (overall antisymmetry)

Combine with the (i<->j) decomposition:
- The (i<->j)-symmetric part, having T^{ia,jb} + T^{ja,ib} symmetric in V-indices, must be antisymmetric in W-indices (a<->b) to maintain overall antisymmetry. This gives **Sym^2(V) tensor wedge^2(W)**.
- The (i<->j)-antisymmetric part, having T^{ia,jb} - T^{ja,ib} antisymmetric in V-indices, must be symmetric in W-indices to maintain overall antisymmetry. This gives **wedge^2(V) tensor Sym^2(W)**.

### Why this matters

When V = 5 (SU(5) flavor) and W = 3 (SU(3) color):

    wedge^2(5 tensor 3) = [Sym^2(5) tensor wedge^2(3)] + [wedge^2(5) tensor Sym^2(3)]
                        = (15, 3-bar)                   + (10-bar, 6)

The overall antisymmetry is satisfied in BOTH terms, but distributed differently:

| Term | Flavor symmetry | Color symmetry | Overall |
|------|----------------|----------------|---------|
| (15, 3-bar) | Symmetric (Sym^2) | Antisymmetric (wedge^2) | Antisymmetric |
| (10-bar, 6) | Antisymmetric (wedge^2) | Symmetric (Sym^2) | Antisymmetric |

Both terms are individually consistent with overall antisymmetry. The tensor product V tensor W creates a CORRELATION between the symmetry types of V and W: symmetric in one forces antisymmetric in the other.

### Contrast with direct sum

When the fundamental is V + V-bar (direct sum, no tensor product):

    wedge^2(V + V-bar) = wedge^2(V) + V tensor V-bar + wedge^2(V-bar)

Here wedge^2(V) = wedge^2(5) = 10-bar. There is no mechanism to produce Sym^2(V) = 15 because V-indices are never "mixed" with other indices through a tensor product. The antisymmetry acts on V-indices alone, forcing wedge^2(V).

**The tensor product 5 tensor 3 is the structural prerequisite for the 15 to appear in wedge^2.**

---

## 4. Comparison Table (Task 4)

### SO(16) decomposition

With 16 = (5,1) + (1,3) + (5-bar,1) + (1,3-bar), adj(SO(16)) = wedge^2(16), dim = C(16,2) = 120.

    wedge^2(A+B+C+D) = wedge^2(A) + wedge^2(B) + wedge^2(C) + wedge^2(D)
                      + A tensor B + A tensor C + A tensor D
                      + B tensor C + B tensor D + C tensor D

where A=(5,1), B=(1,3), C=(5-bar,1), D=(1,3-bar).

| SU(5) x SU(3) | Dim | Source |
|----------------|-----|--------|
| (10-bar, 1) | 10 | wedge^2(A) |
| (1, 3-bar) | 3 | wedge^2(B) |
| (10, 1) | 10 | wedge^2(C) |
| (1, 3) | 3 | wedge^2(D) |
| (5, 3) | 15 | A tensor B |
| (24, 1) + (1, 1) | 25 | A tensor C |
| (5, 3-bar) | 15 | A tensor D |
| (5-bar, 3) | 15 | B tensor C |
| (1, 8) + (1, 1) | 9 | B tensor D |
| (5-bar, 3-bar) | 15 | C tensor D |
| **Total** | **120** | |

**NO (15, 3-bar) appears.** The cross terms A tensor B = (5,1) tensor (1,3) = (5,3) are fundamentals of SU(5) times fundamentals of SU(3), not symmetric tensors. Without the tensor product embedding (where 5 and 3 share the same index space), Sym^2(5) is never produced.

### Master comparison table

| Embedding | Group | Fund. decomp. | dim(adj) | wedge^2 has (15,R)? | wedge^2 has (10-bar,R')? |
|-----------|-------|---------------|----------|---------------------|--------------------------|
| (5,3)+(5-bar,3-bar) | SO(30) | 30 | 435 | **YES**: (15,3-bar) | YES: (10-bar,6) |
| (5,3)+(5-bar,3-bar)+(1,1)^2 | SO(32) | 32 | 496 | **YES**: (15,3-bar) | YES: (10-bar,6) |
| (5,1)+(5-bar,1) | SO(10) | 10 | 45 | **NO** | YES: (10-bar,1) |
| (5,1)+(1,3)+(5-bar,1)+(1,3-bar) | SO(16) | 16 | 120 | **NO** | YES: (10-bar,1) |
| (5,3)+(5-bar,3-bar)+(1,1) | SO(31) | 31 | 465 | **YES**: (15,3-bar) | YES: (10-bar,6) |

### The answer

**The appearance of the symmetric 15 is a consequence of the tensor product embedding 5 tensor 3, independent of the ambient group.** Any SO(n) whose vector representation contains a (5,3) summand under SU(5) x SU(3) will have (15, 3-bar) in its adjoint. The ambient group (SO(30), SO(31), SO(32), ...) changes only the "adjoint sector" (from n x n-bar) and possibly extra cross terms with singlets, but the wedge^2 sector always contains the same (15, 3-bar) + (10-bar, 6) decomposition.

The tensor product is the mechanism. The group is the container.

---

## 5. Clean Projection: Can (15, 3-bar) Appear Without (10-bar, 6)? (Task 5)

### Answer: NO.

This is impossible within any simple Lie group of rank at most 16.

### Argument

**Step 1: Source analysis.** The (15, 3-bar) representation can only arise through the identity

    wedge^2(V tensor W) = [Sym^2(V) tensor wedge^2(W)] + [wedge^2(V) tensor Sym^2(W)]

applied with V = 5 of SU(5) and W = 3 of SU(3). Both terms are nonzero whenever dim(V) >= 2 and dim(W) >= 2. This identity is an exact decomposition; the two pieces span all of wedge^2(V tensor W). Therefore **(15, 3-bar) and (10-bar, 6) always appear together** from this source.

**Step 2: Uniqueness of source.** The 15 = Sym^2(5) of SU(5) has Dynkin label (2,0,0,0). It arises ONLY from the symmetric square of the fundamental representation. Checked explicitly:

| Product | Result | Contains 15? |
|---------|--------|-------------|
| 5 tensor 5 | 15 + 10-bar | YES (Sym^2 part) |
| 5 tensor 5-bar | 24 + 1 | No |
| 5 tensor 10 | 40 + 10-bar | No |
| 5 tensor 15 | 35 + 40 | No |
| 5 tensor 24 | 70 + 45 + 5 | No |
| 10 tensor 10 | 50 + 45 + 5-bar | No |
| 24 tensor 24 | 200 + 126 + 126-bar + 75 + 24 + 24 + 1 | No |

No SU(5) tensor product other than Sym^2(5) produces the 15.

**Step 3: Exhaustion of simple Lie groups.**

*Exceptional groups:*
- G2 (rank 2), F4 (rank 4): Cannot contain SU(5) x SU(3) (rank 7). Not relevant.
- E6 (rank 6, adj = 78): Under SU(5), the adjoint decomposes into {1, 5, 5-bar, 10, 10-bar, 24} only. No 15.
- E7 (rank 7, adj = 133): Same restriction via E7 -> E6 -> SO(10) -> SU(5). No 15.
- E8 (rank 8, adj = 248): Simply-laced; only exterior powers of SU(n) appear. No 15.

*SU(n):*
- adj(SU(n)) = n tensor n-bar - 1. Under the tensor product embedding SU(15) -> SU(5) x SU(3):
  adj(SU(15)) = (24,8) + (24,1) + (1,8). No (15, 3-bar).
- The adjoint of SU(n) never contains wedge^2 or Sym^2 of the fundamental.

*Sp(2n):*
- adj(Sp(2n)) contains Sym^2(n), not wedge^2(n).
- Sym^2(V tensor W) = [Sym^2(V) tensor Sym^2(W)] + [wedge^2(V) tensor wedge^2(W)]
- For V = 5, W = 3: Sym^2((5,3)) = **(15, 6)** + (10-bar, 3-bar)
- This gives (15, 6), NOT (15, 3-bar). The color pairing is wrong.
- Interestingly, Sp(30) gives (10-bar, 3-bar) which is the Pauli-correct diquark, while SO(30) gives (15, 3-bar) which is the Pauli-violating but sBootstrap-needed one.

*SO(2n):*
- adj(SO(2n)) contains wedge^2(n) and wedge^2(n-bar).
- wedge^2((5,3)) = (15, 3-bar) + (10-bar, 6). Both ALWAYS appear together.

### SO vs Sp comparison

The companion identities create a remarkable SO/Sp duality in the color pairings:

| Group family | adj contains | Identity | Flavor-color pairings |
|-------------|-------------|----------|----------------------|
| SO(2n) | wedge^2(n) | wedge^2(V tensor W) | **(15, 3-bar)** + (10-bar, 6) |
| Sp(2n) | Sym^2(n) | Sym^2(V tensor W) | (15, 6) + **(10-bar, 3-bar)** |

SO and Sp swap which SU(5) representation is paired with which SU(3) representation. The bolded entries are:
- SO: symmetric flavor + antisymmetric color = (15, 3-bar), the sBootstrap diquark
- Sp: antisymmetric flavor + antisymmetric color = (10-bar, 3-bar), the Pauli-correct diquark

The orientifold projection in Type I string theory selects SO(32), not Sp(32), and therefore selects the (15, 3-bar) pairing.

### Conclusion

No simple Lie group of any rank produces (15, 3-bar) in its adjoint without simultaneously producing (10-bar, 6). The two representations are born together from the wedge^2(V tensor W) identity and cannot be separated by any group-theoretic branching rule. Their separation requires a projection mechanism external to the Lie algebra structure (orientifold, orbifold, or other stringy mechanism).

---

## 6. Dimension Cross-Checks

### SO(30)
    C(30,2) = 435
    224 + 1 + 105 + 105 = 435
    (192 + 24 + 8 + 1) + (45 + 60) + (45 + 60) = 435

### SO(32)
    C(32,2) = 496
    255 + 1 + 120 + 120 = 496
    (192 + 24 + 8 + 1 + 15 + 15 + 1) + (45 + 60 + 15) + (45 + 60 + 15) = 496

### SO(10)
    C(10,2) = 45
    24 + 1 + 10 + 10 = 45

### SO(16)
    C(16,2) = 120
    10 + 3 + 10 + 3 + 15 + 25 + 15 + 15 + 9 + 15 = 120

### SO(31)
    C(31,2) = 465
    (192 + 24 + 8 + 1) + (45 + 60) + (45 + 60) + 15 + 15 = 465

### Sp(30) (for comparison)
    15 * 31 = 465
    (192 + 24 + 8) + 1 + (90 + 30) + (90 + 30) = 465
    adj(SU(15)) + U(1) + Sym^2(15) + Sym^2(15-bar)
    224 + 1 + 120 + 120 = 465

### wedge^2 vs Sym^2 identity check
    wedge^2(15) = C(15,2) = 105 = 45 + 60
    Sym^2(15) = C(16,2) = 120 = 90 + 30
    Total: 105 + 120 = 225 = 15^2

All dimension checks pass.
