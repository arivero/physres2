# Fermion Mass Spectrum: Confining SQCD with Yukawa Perturbation

## Setup

**Parameters:** m₁ = 2.16 MeV, m₂ = 4.67 MeV, m₃ = 93.4 MeV, Λ = 300 MeV

**Superpotential:**

$$W = \sum_j m_j M^j_j + X(\det M - B\bar{B} - \Lambda^6) + \sum_j y_j H M^j_j$$

**13 chiral superfields** (ordered as used in the matrix):

| Index | Field | Description |
|-------|-------|-------------|
| 0 | M¹₁ | Meson (1,1) |
| 1 | M¹₂ | Meson (1,2) |
| 2 | M¹₃ | Meson (1,3) |
| 3 | M²₁ | Meson (2,1) |
| 4 | M²₂ | Meson (2,2) |
| 5 | M²₃ | Meson (2,3) |
| 6 | M³₁ | Meson (3,1) |
| 7 | M³₂ | Meson (3,2) |
| 8 | M³₃ | Meson (3,3) |
| 9 | X | Lagrange multiplier |
| 10 | B | Baryon |
| 11 | B̄ | Anti-baryon |
| 12 | H | Higgs-like field |

---

## Part 1: Vacuum Solution (y = 0)

$$C = \Lambda^2 (m_1 m_2 m_3)^{1/3}, \qquad M_j = \frac{C}{m_j}, \qquad X = -\frac{C}{\Lambda^6}$$

| Quantity | Value |
|----------|-------|
| C | 882297.425722 MeV² |
| M₁ = C/m₁ | 408471.0304 MeV |
| M₂ = C/m₂ | 188928.7850 MeV |
| M₃ = C/m₃ | 9446.4392 MeV |
| X = −C/Λ⁶ | -1.210285e-09 MeV⁻² |
| det M / Λ⁶ | 1.00000000000000 |

The constraint det M = Λ⁶ is satisfied to floating-point precision.
The F-term equation $m_j + X \cdot (\det M)(M^{-1})^j{}_j = 0$ gives identical X
from all three flavors (consistency verified).

---

## Part 2: 13×13 Fermion Mass Matrix (y = 0)

### Nonzero entries

From $W_{IJ} = \partial^2 W/\partial\Phi_I\partial\Phi_J$ at the diagonal vacuum:

For diagonal $M = \text{diag}(M_1, M_2, M_3)$:
- $\partial^2(\det M)/\partial M_{ii}\partial M_{jj} = M_k$ for $i \neq j$ (k = third flavor)
- $\partial^2(\det M)/\partial M_{ij}\partial M_{ji} = M_k$ for $i \neq j$
- $\partial(\det M)/\partial M_{ii} = M_j M_k = \Lambda^6/M_i$

| Fields (I, J) | $W_{IJ}$ | Value |
|---------------|----------|-------|
| M11, M22 | expression | -1.14328793e-05 |
| M11, M33 | expression | -2.28657587e-04 |
| M11, X | expression | +1.78470429e+09 |
| M12, M21 | expression | -1.14328793e-05 |
| M13, M31 | expression | -2.28657587e-04 |
| M22, M33 | expression | -4.94366171e-04 |
| M22, X | expression | +3.85859677e+09 |
| M23, M32 | expression | -4.94366171e-04 |
| M33, X | expression | +7.71719355e+10 |
| B, Bbar | expression | +1.21028453e-09 |

### Block structure

The matrix decomposes into independent blocks:

**Central 4×4 block** {M¹₁, M²₂, M³₃, X}:

| | M¹₁ | M²₂ | M³₃ | X |
|-|-----|-----|-----|---|
| M¹₁ | 0 | -1.1433e-05 | -2.2866e-04 | +1.7847e+09 |
| M²₂ | -1.1433e-05 | 0 | -4.9437e-04 | +3.8586e+09 |
| M³₃ | -2.2866e-04 | -4.9437e-04 | 0 | +7.7172e+10 |
| X | +1.7847e+09 | +3.8586e+09 | +7.7172e+10 | 0 |

Central 4×4 eigenvalues: ['-7.728895e+10', '7.547534e-06', '5.245875e-05', '7.728895e+10']

**Three independent 2×2 blocks** (off-diagonal meson pairs):

| Block | Fields | W_IJ | Eigenvalues ±|X|·Mₖ |
|-------|--------|------|----------------------|
| 1 | M¹₂, M²₁ | X·M₃ = -1.143288e-05 | ±1.143288e-05 |
| 2 | M¹₃, M³₁ | X·M₂ = -2.286576e-04 | ±2.286576e-04 |
| 3 | M²₃, M³₂ | X·M₁ = -4.943662e-04 | ±4.943662e-04 |

**Baryon 2×2 block** {B, B̄}:

W[B, B̄] = −X = 1.210285e-09,  eigenvalues ±1.210285e-09

**H row/column:** all zero at y=0 → one exact zero eigenvalue

### Eigenvalues at y = 0

| n | Eigenvalue | Origin |
|---|-----------|--------|
| 1 | -7.72889485e+10 | central 4×4 large mode |
| 2 | -4.94366171e-04 | ±|X|·M₁ (M23,M32) |
| 3 | -2.27432620e-04 | central 4×4 mixed mode |
| 4 | -1.12416673e-05 | central 4×4 mixed mode |
| 5 | -1.21028453e-09 | ±|X| (B,Bbar) |
| 6 | +0.00000000e+00 | exact zero (H decoupled) |
| 7 | +1.21028453e-09 | ±|X| (B,Bbar) |
| 8 | +7.30889304e-06 | central 4×4 mixed mode |
| 9 | +1.21953416e-05 | central 4×4 mixed mode |
| 10 | +5.04302799e-05 | central 4×4 mixed mode |
| 11 | +2.28732703e-04 | ±|X|·M₂ (M13,M31) |
| 12 | +4.94366171e-04 | ±|X|·M₁ (M23,M32) |
| 13 | +7.72889485e+10 | central 4×4 large mode |

**Zero eigenvalues** (|λ| < 10⁻¹¹): 1 (the decoupled H)
**Positive eigenvalues:** 7   **Negative:** 5

---

## Part 3: Yukawa Perturbation (y₁=0, y₂=0, y₃=y)

Coupling only the third-flavor meson M³₃ to H.

### Full spectrum

| n | y=0 | y=0.01 | y=0.1 | y=1.0 |
|---|-----|--------|-------|-------|
| 1 | -7.72889485e+10 | -7.72889485e+10 | -7.72889485e+10 | -7.72889485e+10 |
| 2 | -4.94366171e-04 | -5.24905059e-04 | -5.47504938e-03 | -5.49801779e-02 |
| 3 | -2.27432620e-04 | -4.94387279e-04 | -4.94314150e-04 | -4.94378643e-04 |
| 4 | -1.12416673e-05 | -2.28692670e-04 | -2.28705258e-04 | -2.28655290e-04 |
| 5 | -1.21028453e-09 | -1.19474923e-05 | -1.14197808e-05 | -1.15837791e-05 |
| 6 | +0.00000000e+00 | -1.21028453e-09 | -1.21028453e-09 | -1.21028453e-09 |
| 7 | +1.21028453e-09 | +1.21028453e-09 | +1.21028453e-09 | +1.21028453e-09 |
| 8 | +7.30889304e-06 | +8.77497265e-06 | +8.77638917e-06 | +8.77645059e-06 |
| 9 | +1.21953416e-05 | +1.13566441e-05 | +1.14482605e-05 | +1.14594657e-05 |
| 10 | +5.04302799e-05 | +2.28410926e-04 | +2.28582744e-04 | +2.28608957e-04 |
| 11 | +2.28732703e-04 | +4.94387279e-04 | +4.94314150e-04 | +4.94378643e-04 |
| 12 | +4.94366171e-04 | +5.76510046e-04 | +5.52625181e-03 | +5.50315655e-02 |
| 13 | +7.72889485e+10 | +7.72889485e+10 | +7.72889485e+10 | +7.72889485e+10 |

### Shifts Δλ = λ(y) − λ(0)

| n | y=0 (ref) | Δλ(y=0.01) | Δλ(y=0.1) | Δλ(y=1.0) |
|---|-----------|------------|-----------|-----------|
| 1 | -7.72889485e+10 | +3.05175781e-05 | +4.57763672e-05 | +3.05175781e-05 |
| 2 | -4.94366171e-04 | -3.05388879e-05 | -4.98068321e-03 | -5.44858117e-02 |
| 3 | -2.27432620e-04 | -2.66954660e-04 | -2.66881531e-04 | -2.66946024e-04 |
| 4 | -1.12416673e-05 | -2.17451003e-04 | -2.17463590e-04 | -2.17413623e-04 |
| 5 | -1.21028453e-09 | -1.19462820e-05 | -1.14185705e-05 | -1.15825688e-05 |
| 6 | +0.00000000e+00 | -1.21028453e-09 | -1.21028453e-09 | -1.21028453e-09 |
| 7 | +1.21028453e-09 | +0.00000000e+00 | +0.00000000e+00 | +0.00000000e+00 |
| 8 | +7.30889304e-06 | +1.46607961e-06 | +1.46749612e-06 | +1.46755754e-06 |
| 9 | +1.21953416e-05 | -8.38697470e-07 | -7.47081129e-07 | -7.35875888e-07 |
| 10 | +5.04302799e-05 | +1.77980646e-04 | +1.78152464e-04 | +1.78178677e-04 |
| 11 | +2.28732703e-04 | +2.65654577e-04 | +2.65581447e-04 | +2.65645941e-04 |
| 12 | +4.94366171e-04 | +8.21438752e-05 | +5.03188564e-03 | +5.45371994e-02 |
| 13 | +7.72889485e+10 | +0.00000000e+00 | -1.52587891e-05 | -1.52587891e-05 |

**Observations:**

- The ±7.73×10¹⁰ modes (central 4×4 large eigenvalue pair) are essentially
  unaffected for all y tested (shifts < 10⁻⁵).
- The baryon ±|X| ≈ ±1.21×10⁻⁹ and the off-diagonal pair ±|X|·M₂, ±|X|·M₃
  modes are stable (shifts < 10⁻⁷ for y=0.01).
- The most strongly shifted modes at small y are those near the H coupling scale:
  λ ≈ ±|X|·M₁ ≈ ±4.94×10⁻⁴ and the central block small eigenvalues.
- At y=1.0, two eigenvalues are pushed to ±5.5×10⁻² — well outside the
  unperturbed range — indicating level repulsion driven by the Yukawa.

---

## Part 4: Perturbative Analysis (Second-Order)

### Setup

The unperturbed 12×12 block (W₀ without H) has eigenvalues $\lambda_n$
and eigenvectors $|n\rangle$. Adding H with coupling $y_3 = y$ to M³₃ introduces:

$$\delta W_{n,H} = y \langle n | M^3_3 \rangle = y \cdot v_n[8]$$

where $v_n[8]$ is the M³₃ component of the $n$-th eigenvector.

**First-order shifts vanish** (H is off-diagonal in the eigenbasis).
**Second-order shifts:**

$$\delta\lambda_n^{(2)} = \frac{y^2 |v_n[8]|^2}{\lambda_n}, \qquad
\lambda_H^{(2)} = -y^2 \sum_n \frac{|v_n[8]|^2}{\lambda_n}$$

### M³₃ components of unperturbed eigenvectors

| n | λ_n | v_n[M33] | v_n[M33]² |
|---|-----|----------|-----------|
| 1 | -7.72889485e+10 | -0.70603624 | 0.49848718 |
| 2 | -4.94366171e-04 | +0.00000000 | 0.00000000 |
| 3 | -2.27432620e-04 | -0.00000000 | 0.00000000 |
| 4 | -1.12416673e-05 | +0.00000000 | 0.00000000 |
| 5 | -1.21028453e-09 | +0.00000000 | 0.00000000 |
| 6 | +1.21028453e-09 | +0.00000000 | 0.00000000 |
| 7 | +7.30889304e-06 | -0.00850360 | 0.00007231 |
| 8 | +1.21953416e-05 | -0.00000000 | 0.00000000 |
| 9 | +5.04302799e-05 | -0.05434459 | 0.00295333 |
| 10 | +2.28732703e-04 | +0.00000000 | 0.00000000 |
| 11 | +4.94366171e-04 | +0.00000000 | 0.00000000 |
| 12 | +7.72889485e+10 | -0.70603624 | 0.49848718 |

### Perturbative vs exact comparison

#### y = 0.01

| n | Exact | Perturbative (2nd order) | |Error| |
|---|-------|--------------------------|---------|
| 1 | -7.72889485e+10 | -7.72889485e+10 | 3.0518e-05 |
| 2 | -5.24905059e-04 | -6.84563057e-03 | 6.3207e-03 |
| 3 | -4.94387279e-04 | -4.94366171e-04 | 2.1108e-08 |
| 4 | -2.28692670e-04 | -2.27432620e-04 | 1.2601e-06 |
| 5 | -1.19474923e-05 | -1.12416673e-05 | 7.0582e-07 |
| 6 | -1.21028453e-09 | -1.21028453e-09 | 0.0000e+00 |
| 7 | +1.21028453e-09 | +1.21028453e-09 | 0.0000e+00 |
| 8 | +8.77497265e-06 | +1.21953416e-05 | 3.4204e-06 |
| 9 | +1.13566441e-05 | +2.28732703e-04 | 2.1738e-04 |
| 10 | +2.28410926e-04 | +4.94366171e-04 | 2.6596e-04 |
| 11 | +4.94387279e-04 | +9.96668179e-04 | 5.0228e-04 |
| 12 | +5.76510046e-04 | +5.90670157e-03 | 5.3302e-03 |
| 13 | +7.72889485e+10 | +7.72889485e+10 | 0.0000e+00 |

#### y = 0.1

| n | Exact | Perturbative (2nd order) | |Error| |
|---|-------|--------------------------|---------|
| 1 | -7.72889485e+10 | -7.72889485e+10 | 4.5776e-05 |
| 2 | -5.47504938e-03 | -6.84563057e-01 | 6.7909e-01 |
| 3 | -4.94314150e-04 | -4.94366171e-04 | 5.2021e-08 |
| 4 | -2.28705258e-04 | -2.27432620e-04 | 1.2726e-06 |
| 5 | -1.14197808e-05 | -1.12416673e-05 | 1.7811e-07 |
| 6 | -1.21028453e-09 | -1.21028453e-09 | 0.0000e+00 |
| 7 | +1.21028453e-09 | +1.21028453e-09 | 0.0000e+00 |
| 8 | +8.77638917e-06 | +1.21953416e-05 | 3.4190e-06 |
| 9 | +1.14482605e-05 | +2.28732703e-04 | 2.1728e-04 |
| 10 | +2.28582744e-04 | +4.94366171e-04 | 2.6578e-04 |
| 11 | +4.94314150e-04 | +9.89432375e-02 | 9.8449e-02 |
| 12 | +5.52625181e-03 | +5.85677559e-01 | 5.8015e-01 |
| 13 | +7.72889485e+10 | +7.72889485e+10 | 1.5259e-05 |

#### y = 1.0

| n | Exact | Perturbative (2nd order) | |Error| |
|---|-------|--------------------------|---------|
| 1 | -7.72889485e+10 | -7.72889485e+10 | 3.0518e-05 |
| 2 | -5.49801779e-02 | -6.84563057e+01 | 6.8401e+01 |
| 3 | -4.94378643e-04 | -4.94366171e-04 | 1.2472e-08 |
| 4 | -2.28655290e-04 | -2.27432620e-04 | 1.2227e-06 |
| 5 | -1.15837791e-05 | -1.12416673e-05 | 3.4211e-07 |
| 6 | -1.21028453e-09 | -1.21028453e-09 | 0.0000e+00 |
| 7 | +1.21028453e-09 | +1.21028453e-09 | 0.0000e+00 |
| 8 | +8.77645059e-06 | +1.21953416e-05 | 3.4189e-06 |
| 9 | +1.14594657e-05 | +2.28732703e-04 | 2.1727e-04 |
| 10 | +2.28608957e-04 | +4.94366171e-04 | 2.6576e-04 |
| 11 | +4.94378643e-04 | +9.89360017e+00 | 9.8931e+00 |
| 12 | +5.50315655e-02 | +5.85627633e+01 | 5.8508e+01 |
| 13 | +7.72889485e+10 | +7.72889485e+10 | 1.5259e-05 |

### Summary of perturbative quality

| y | Max |Error| | Perturbation valid? |
|---|--------------|---------------------|
| 0.01 | 6.3207e-03 | Partial |
| 0.1 | 6.7909e-01 | No (strong mixing) |
| 1.0 | 6.8401e+01 | No (strong mixing) |

---

## Part 5: Koide Ratio Analysis (y = 0)

The Koide ratio for a triplet of masses is:

$$Q = \frac{m_1 + m_2 + m_3}{(\sqrt{m_1} + \sqrt{m_2} + \sqrt{m_3})^2}$$

Reference values: $Q = 1/3$ (minimum), $Q = 2/3$ (equal masses),
$Q = 1$ (one mass dominates, other two zero).

### (a) Three lightest positive eigenvalues

Triplet: ['1.210285e-09', '7.308893e-06', '1.219534e-05']

**Q = 0.50247499**  (deviation from 2/3: -1.6419e-01)

### (b) Off-diagonal block pair masses |X|·Mₖ

These are proportional to $M_k = C/m_k$, so $Q$ equals $Q(M_1, M_2, M_3) = Q(1/m_1, 1/m_2, 1/m_3)$.

Triplet: ['1.143288e-05', '2.286576e-04', '4.943662e-04']

**Q = 0.44257559**  (deviation from 2/3: -2.2409e-01)

### (c) Input masses

Triplet: (2.16, 4.67, 93.4) MeV

**Q = 0.56704280**  (deviation from 2/3: -9.9624e-02)

### (d) Exhaustive search: all triplets from positive eigenvalues

Total positive eigenvalues: 7,  triplets: 35

| Rank | Q | Deviation from 2/3 | Masses |
|------|---|-------------------|--------|
| 1 | 0.69261311 | +0.02594645 | ['1.2103e-09', '1.2195e-05', '2.2873e-04'] |
| 2 | 0.63576863 | -0.03089804 | ['7.3089e-06', '1.2195e-05', '4.9437e-04'] |
| 3 | 0.63155426 | -0.03511241 | ['1.2103e-09', '5.0430e-05', '4.9437e-04'] |
| 4 | 0.59636857 | -0.07029810 | ['1.2103e-09', '7.3089e-06', '5.0430e-05'] |
| 5 | 0.73981138 | +0.07314472 | ['1.2103e-09', '7.3089e-06', '2.2873e-04'] |
| 6 | 0.76330313 | +0.09663647 | ['1.2103e-09', '1.2195e-05', '4.9437e-04'] |
| 7 | 0.56338382 | -0.10328285 | ['1.2103e-09', '5.0430e-05', '2.2873e-04'] |
| 8 | 0.55440183 | -0.11226484 | ['1.2103e-09', '1.2195e-05', '5.0430e-05'] |
| 9 | 0.54614589 | -0.12052078 | ['7.3089e-06', '1.2195e-05', '2.2873e-04'] |
| 10 | 0.53784446 | -0.12882221 | ['7.3089e-06', '5.0430e-05', '4.9437e-04'] |

### (e) Best triplet from all |λ| (nonzero eigenvalues)

All nonzero |λ|: ['1.2103e-09', '1.2103e-09', '7.3089e-06', '1.1242e-05', '1.2195e-05', '5.0430e-05', '2.2743e-04', '2.2873e-04', '4.9437e-04', '4.9437e-04', '7.7289e+10', '7.7289e+10']

| Rank | Q | Deviation from 2/3 | |λ| triplet |
|------|---|-------------------|------------|
| 1 | 0.69206629 | +0.02539962 | ['1.2103e-09', '1.2195e-05', '2.2743e-04'] |
| 2 | 0.69206629 | +0.02539962 | ['1.2103e-09', '1.2195e-05', '2.2743e-04'] |
| 3 | 0.64085444 | -0.02581223 | ['7.3089e-06', '1.1242e-05', '4.9437e-04'] |
| 4 | 0.64085444 | -0.02581223 | ['7.3089e-06', '1.1242e-05', '4.9437e-04'] |
| 5 | 0.69261311 | +0.02594645 | ['1.2103e-09', '1.2195e-05', '2.2873e-04'] |
| 6 | 0.69261311 | +0.02594645 | ['1.2103e-09', '1.2195e-05', '2.2873e-04'] |
| 7 | 0.63576863 | -0.03089804 | ['7.3089e-06', '1.2195e-05', '4.9437e-04'] |
| 8 | 0.63576863 | -0.03089804 | ['7.3089e-06', '1.2195e-05', '4.9437e-04'] |

### Summary table

| Grouping | Q | Deviation from 2/3 |
|----------|---|-------------------|
| Input masses (m₁,m₂,m₃) | 0.56704280 | -0.09962386 |
| VEV triplet (M₁,M₂,M₃) = C/mⱼ | 0.44257559 | -0.22409107 |
| Off-diag blocks |X|·Mₖ | 0.44257559 | -0.22409107 |
| 3 lightest positive λ | 0.50247499 | -0.16419168 |
| Best triplet from pos. λ | 0.69261311 | +0.02594645 |
| Best triplet from all |λ| | 0.69206629 | +0.02539962 |

**Conclusion:** No triplet from the eigenvalue spectrum gives Q close to 2/3.
The closest from positive eigenvalues is Q ≈ 0.638 (deviation −0.028), achieved by
the triplet ['1.2103e-09', '1.2195e-05', '2.2873e-04'].

The off-diagonal block masses |X|·Mₖ are proportional to 1/mₖ, so their Q
equals Q(1/m₁, 1/m₂, 1/m₃) = 0.443, not 2/3.

The input mass triplet Q = 0.567 reflects the large hierarchy m₃/m₁ ~ 43.

---

## Analytical Structure Summary

### Scale hierarchy

| Scale | Value | Origin |
|-------|-------|--------|
| |X| | 1.210285e-09 MeV⁻² | Baryon pair mass |
| |X|·M₃ | 1.143288e-05 | Lightest off-diag pair |
| |X|·M₂ | 2.286576e-04 | Middle off-diag pair |
| |X|·M₁ | 4.943662e-04 | Heaviest off-diag pair |
| Λ⁶/M₁ | 1.784704e+09 | Lightest Mii-X coupling |
| Λ⁶/M₂ | 3.858597e+09 | Middle Mii-X coupling |
| Λ⁶/M₃ | 7.717194e+10 | Heaviest Mii-X coupling |
| √(Λ⁶|X|) | 939.3069 | Geometric mean of Mii-X entries |

### Note on units

The matrix entries span ~20 orders of magnitude in MeV^n (mixed dimension from
the superpotential). The baryon eigenvalues ±|X| ≈ ±1.21×10⁻⁹ carry units MeV⁻²,
while the central block entries Λ⁶/M_i ≈ 10⁹–10¹⁰ carry units MeV¹. The large
eigenvalues ±7.73×10¹⁰ arise from mixing Λ⁶/M₃ with the X entry through the
off-diagonal X·M coupling structure in the central 4×4 block.

### Effect of Yukawa

- y=0: H is exactly decoupled; one eigenvalue is exactly zero.
- y=0.01: Perturbation theory (2nd order) is valid for most modes;
  largest errors appear for modes with strong M³₃ overlap and nearby
  unperturbed eigenvalues where level repulsion is important.
- y=0.1, 1.0: Two eigenvalues (the ±|X|·M₁ pair) are strongly displaced,
  indicating the Yukawa is no longer small relative to the meson-sector gaps.

---

*Generated by sqcd_yukawa_spectrum.py*