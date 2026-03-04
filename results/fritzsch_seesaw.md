# Fritzsch Texture, CKM Mixing, and Seesaw Meson VEVs

## Part (a): Fritzsch Relations — Characteristic Polynomial Verification

The Fritzsch nearest-neighbor Hermitian mass matrix is:

```
M = [[0,  A,  0],
     [A,  0,  B],
     [0,  B,  C]]
```

The characteristic polynomial is:

    lambda^3 - C lambda^2 - (A^2 + B^2) lambda + A^2 C = 0

By Vieta's formulas, if the eigenvalues are (+m_1, -m_2, +m_3) with m_i > 0:

| Vieta relation | LHS | RHS |
|---|---|---|
| sum lambda_i | m_1 - m_2 + m_3 | C |
| sum lambda_i lambda_j | m_1 m_2 + m_2 m_3 - m_1 m_3 | A^2 + B^2 |
| prod lambda_i | -m_1 m_2 m_3 | -A^2 C |

The **exact** Fritzsch relations are:

    C = m_3 - m_2 + m_1
    A^2 = m_1 m_2 m_3 / C = m_1 m_2 m_3 / (m_3 - m_2 + m_1)
    B^2 = m_1 m_2 + m_2 m_3 - m_1 m_3 - A^2

The **approximate** (leading order in mass ratios) relations are:

    |A|^2 ~ m_1 m_2,   |B|^2 ~ m_2 m_3,   C ~ m_3 - m_2

**Numerical verification (exact relations):**

| Parameter | Up sector | Down sector |
|---|---|---|
| C | 171487.16 MeV | 4091.27 MeV |
| A^2 (exact) | 2774.44 | 445.64 |
| A^2 (approx m_1 m_2) | 2754.00 | 436.18 |
| B^2 (exact) | 2.199 x 10^8 | 3.709 x 10^5 |
| B^2 (approx m_2 m_3) | 2.203 x 10^8 | 3.904 x 10^5 |
| Eigenvalue max error | 2.3 x 10^{-13} | 1.4 x 10^{-14} |

The exact relations reproduce input masses to machine precision (< 10^{-13} relative error). The approximate relations have ~1% (up) and ~3% (down) errors in B due to the non-negligible m_1/m_2 ratios.


## Part (b): CKM Angles

### GST/Oakes derivation

For a 2x2 Fritzsch block M = [[0, a],[a, D]], the diagonalizing rotation has:

    sin theta = sqrt(m_light / (m_light + m_heavy)) ~ sqrt(m_light / m_heavy)

The CKM (1,2) element combines up and down sector rotations:

    V_us ~ phi_d - phi_u ~ sqrt(m_d/m_s) - sqrt(m_u/m_c)

In the limit m_u -> 0: **sin theta_C ~ sqrt(m_d/m_s)** (GST/Oakes relation).

### Numerical results

| Method | sin theta_12 | theta_12 | sin theta_23 | theta_23 | sin theta_13 | theta_13 |
|---|---|---|---|---|---|---|
| GST/Oakes | 0.2236 | 12.92 deg | 0.1495 | 8.60 deg | 0.0334 | 1.92 deg |
| Fritzsch approx | 0.1824 | 10.51 deg | 0.0636 | 3.64 deg | 0.0299 | 1.71 deg |
| Full diag (exact) | 0.1780 | 10.25 deg | 0.0589 | 3.38 deg | 0.0031 | 0.18 deg |
| PDG 2024 | 0.2256 | 13.04 deg | 0.0415 | 2.38 deg | 0.0035 | 0.20 deg |

**Key observations:**

1. **GST/Oakes gives the best theta_12**: sqrt(m_d/m_s) = 0.2236 vs PDG 0.2256 (0.9% deviation, 0.12 deg).

2. **The Fritzsch texture overshoots theta_23** by ~42-50% (predicts 3.4-3.6 deg vs PDG 2.38 deg). This is a well-known failure mode of the original Fritzsch texture.

3. **The full diagonalization undershoots theta_12** relative to the approximate formula. This happens because the 3x3 exact diagonalization is not well approximated by sequential 2x2 rotations when mass ratios are only O(10-100).

4. **theta_13 from full diag** (0.18 deg) matches PDG (0.20 deg) reasonably, but the Fritzsch approximate formula for theta_13 is wildly off (1.7 deg vs 0.20 deg) because the product approximation breaks down.

5. The full CKM matrix from exact Fritzsch diagonalization is:

```
|V_CKM| = [[0.9840  0.1780  0.0031]
            [0.1775  0.9823  0.0589]
            [0.0136  0.0574  0.9983]]
```


## Part (c): Seesaw Map — Matrix Inversion

**The seesaw is a matrix inversion, not just eigenvalue inversion.**

For N_f = N_c = 3 SQCD with quantum-modified constraint det(Phi) = Lambda^6 and mass deformation W_mass = Tr(m_q Phi), the F-term condition gives:

    Phi = alpha * m_q^{-T}

where alpha = Lambda^2 (det m_q)^{1/3} is fixed by the constraint.

If m_q = U diag(m_k) V^dag, then:

    Phi = alpha * (V diag(1/m_k) U^dag)^T = alpha * U* diag(1/m_k) V^T

For real symmetric Fritzsch (U = V):

    Phi = U diag(alpha/m_k) U^T

**The eigenvectors of Phi are identical to those of m_q. The off-diagonal structure in flavor space is preserved. Only the eigenvalues are inverted (with overall rescaling).**

This is simply the statement that the matrix inverse shares eigenvectors with the original matrix, with inverted eigenvalues.


## Part (d): Numerical Meson VEV Matrix (Down Sector)

Scale: alpha = Lambda^2 (m_d m_s m_b)^{1/3} = 1.099 x 10^7 MeV^2

Meson eigenvalues:
- alpha/m_d = 2,354,357 MeV (lightest quark -> largest meson VEV)
- alpha/m_s = 117,718 MeV
- alpha/m_b = 2,630 MeV (heaviest quark -> smallest meson VEV)

Meson VEV matrix in flavor basis (MeV):

```
         d              s              b
d    2,247,799      471,206       -70,306
s      471,206      219,569       -31,954
b      -70,306      -31,954         7,338
```

**Off-diagonal ratios:**

| Ratio | Value | Compare to |
|---|---|---|
| Phi_ds / Phi_dd | 0.210 | sqrt(m_d/m_s) = 0.224 |
| Phi_sb / Phi_ss | -0.146 | sqrt(m_s/m_b) = 0.149 |
| Phi_db / Phi_dd | -0.031 | sqrt(m_d/m_b) = 0.033 |

The off-diagonal meson VEV ratios are close to (but not identical with) the Cabibbo-like mass ratios. This is expected: the off-diagonal structure of Phi inherits the mixing pattern of the Fritzsch texture.


## Part (e): Determinant Consistency

**Algebraic proof:**

    det(Phi) = det(U)^2 * prod(alpha/m_k)
             = 1 * alpha^3 / (m_d m_s m_b)
             = [Lambda^6 m_d m_s m_b] / (m_d m_s m_b)
             = Lambda^6

This is independent of the mixing matrix U. The Seiberg constraint det(Phi) = Lambda^6 is automatically satisfied whether or not Phi has off-diagonal entries.

**Numerical verification:**

    det(Phi) = 7.290000 x 10^{14} MeV^6
    Lambda^6 = 7.290000 x 10^{14} MeV^6
    Relative error: 8.6 x 10^{-16}

The off-diagonal meson VEVs from the Fritzsch texture are fully consistent with the quantum-modified constraint of SQCD.


## Part (f): Cabibbo Angle Comparison

### theta_12 (Cabibbo angle)

| Method | sin theta_C | theta_C | Deviation from PDG |
|---|---|---|---|
| GST/Oakes: sqrt(m_d/m_s) | 0.2236 | 12.92 deg | -0.12 deg |
| Fritzsch: \|sqrt(m_d/m_s) - sqrt(m_u/m_c)\| | 0.1824 | 10.51 deg | -2.53 deg |
| Full diag | 0.1780 | 10.25 deg | -2.79 deg |
| PDG | 0.2257 | 13.04 deg | --- |

### Heavy angles

| Method | sin theta_23 | theta_23 | sin theta_13 | theta_13 |
|---|---|---|---|---|
| Fritzsch approx | 0.0636 | 3.64 deg | 0.0299 | 1.71 deg |
| Full diag | 0.0589 | 3.38 deg | 0.0031 | 0.18 deg |
| PDG | 0.0415 | 2.38 deg | 0.0035 | 0.20 deg |

### Individual sector mixing angles

| Sector | 1-2 rotation | 2-3 rotation |
|---|---|---|
| Down | sqrt(m_d/m_s) = 0.2236 | sqrt(m_s/m_b) = 0.1495 |
| Up | sqrt(m_u/m_c) = 0.0412 | sqrt(m_c/m_t) = 0.0859 |

### Summary

1. **GST/Oakes is the winner for theta_12**: sqrt(m_d/m_s) = 0.224 matches |V_us|_{PDG} = 0.226 to 0.9%.

2. **The Fritzsch texture overshoots theta_23** by ~42-50%. This is the historical failure of the original (1978) Fritzsch ansatz and is why more general textures were later explored.

3. **theta_13 from full diagonalization** (0.18 deg) agrees with PDG (0.20 deg) to 11%, but this is somewhat accidental in the Fritzsch texture. The approximate formula fails badly for theta_13.

4. The GST/Oakes relation sin theta_C = sqrt(m_d/m_s) encodes the core physics: the Cabibbo angle is dominated by the down-sector mass hierarchy, with the up-sector contribution (sqrt(m_u/m_c) = 0.041) being a small correction that actually *worsens* the fit when included.
