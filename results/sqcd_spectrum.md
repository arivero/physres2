# SQCD Seiberg Effective Theory: Mass Spectrum

## Setup

SU(3)_c SQCD with quark mass matrix $m = \text{diag}(m_u, m_d, m_s, m_c, m_b)$.

Inputs (MS-bar at 2 GeV):
$m_u = 2.16$ MeV, $m_d = 4.67$ MeV, $m_s = 93.4$ MeV, $m_c = 1270$ MeV, $m_b = 4180$ MeV, $\Lambda = 300$ MeV.

Koide convention throughout: $Q = (m_1+m_2+m_3)/(\sqrt{m_1}+\sqrt{m_2}+\sqrt{m_3})^2$, so $Q(e,\mu,\tau) = 2/3$.

---

## 1. Seiberg Seesaw: Meson VEVs ($N_f = N_c = 3$)

For the confining phase with quantum-modified constraint $\det M - B\bar{B} = \Lambda^6$:

$$W = \text{Tr}(m\,M) + X(\det M - \Lambda^6)$$

The SUSY vacuum at $B = \bar{B} = 0$ gives:

$$m_j M_j = C \quad \text{(constant for all $j$)}$$

where $C = \Lambda^2 (\prod_k m_k)^{1/3}$ and $M_j = C/m_j$.

This is the **Seiberg seesaw**: $M_j \propto 1/m_j$.

### Triple (u,d,s)

| Parameter | Value |
|-----------|-------|
| $C$ | $882297$ MeV$^2$ |
| $X$ | $-1.210 \times 10^{-9}$ MeV$^{-4}$ |
| $M_u$ | $408471$ MeV |
| $M_d$ | $188929$ MeV |
| $M_s$ | $9446$ MeV |
| $Q(m_u,m_d,m_s)$ | 0.5670 |
| $Q(1/m_u,1/m_d,1/m_s) = Q(M)$ | 0.4426 |

### Triple (d,s,b) -- the dual Koide triple

| Parameter | Value |
|-----------|-------|
| $C$ | $10{,}994{,}847$ MeV$^2$ |
| $M_d$ | $2{,}354{,}357$ MeV |
| $M_s$ | $117{,}718$ MeV |
| $M_b$ | $2630$ MeV |
| $Q(m_d,m_s,m_b)$ | 0.7314 |
| $Q(1/m_d,1/m_s,1/m_b) = Q(M)$ | **0.6652** |
| $|Q - 2/3|$ | $0.00146$ |
| Deviation | **0.22%** from 2/3 |

This confirms the dual Koide observation. The Seiberg seesaw maps it to a standard Koide condition on the meson VEVs.

### Triple (s,c,b)

| Parameter | Value |
|-----------|-------|
| $C$ | $71{,}233{,}585$ MeV$^2$ |
| $M_s$ | $762{,}672$ MeV |
| $M_c$ | $56{,}089$ MeV |
| $M_b$ | $17{,}042$ MeV |
| $Q(M)$ | 0.5430 |

### All triples -- Koide table

| Triple | $Q(m)$ | $Q(1/m)$ | $|Q(1/m)-2/3|$ | % dev |
|--------|--------|----------|-----------------|-------|
| $(u,d,s)$ | 0.5670 | 0.4426 | 0.2241 | 33.6% |
| $(u,d,c)$ | 0.8281 | 0.4942 | 0.1725 | 25.9% |
| $(u,d,b)$ | 0.8980 | 0.5046 | 0.1621 | 24.3% |
| $(u,s,c)$ | 0.6242 | 0.7197 | 0.0530 | 8.0% |
| $(u,s,b)$ | 0.7444 | 0.7417 | 0.0750 | 11.3% |
| $(u,c,b)$ | 0.5265 | 0.8853 | 0.2187 | 32.8% |
| $(d,s,c)$ | 0.6073 | 0.6389 | 0.0278 | 4.2% |
| **(d,s,b)** | **0.7314** | **0.6652** | **0.0015** | **0.22%** |
| $(d,c,b)$ | 0.5197 | 0.8394 | 0.1728 | 25.9% |
| $(s,c,b)$ | 0.4585 | 0.5430 | 0.1237 | 18.6% |

The triple $(d,s,b)$ is the only one with $Q(1/m)$ within 1% of 2/3. It is the **unique near-Koide triple** in the inverse-mass sector.

### Seesaw meson VEVs for all triples

| Triple | $M_1$ (MeV) | $M_2$ (MeV) | $M_3$ (MeV) | $Q(M)$ |
|--------|-------------|-------------|-------------|--------|
| $(u,d,s)$ | 408471 | 188929 | 9446 | 0.4426 |
| $(u,d,c)$ | 974945 | 450938 | 1658 | 0.4942 |
| $(u,d,b)$ | 1450233 | 670772 | 749 | 0.5046 |
| $(u,s,c)$ | 2646408 | 61202 | 4501 | 0.7197 |
| $(u,s,b)$ | 3936538 | 91038 | 2034 | 0.7417 |
| $(u,c,b)$ | 9395791 | 15980 | 4855 | 0.8853 |
| $(d,s,c)$ | 1582759 | 79138 | 5820 | 0.6389 |
| $(d,s,b)$ | 2354357 | 117718 | 2630 | **0.6652** |
| $(d,c,b)$ | 5619417 | 20664 | 6278 | 0.8394 |
| $(s,c,b)$ | 762672 | 56089 | 17042 | 0.5430 |

---

## 2. Fermion Mass Matrix

The 10 chiral fields at the vacuum are: 9 meson components $M^i_j$ ($3\times 3$) + Lagrange multiplier $X$.

The fermion mass matrix $\mathcal{W}_{IJ} = \partial^2 W/\partial \Phi_I \partial \Phi_J$ at the vacuum decomposes into blocks:

### Off-diagonal 2x2 blocks

For each pair $(M^i_j, M^j_i)$ with $i \neq j$, the third index $k$ gives:

$$\mathcal{W} = \begin{pmatrix} 0 & -X M_k \\ -X M_k & 0 \end{pmatrix}$$

Eigenvalues: $\pm |X| M_k$.

| Pair | $k$ | Mass $= |X| M_k$ | Value (MeV) |
|------|-----|-------------------|-------------|
| $(M_{ud}, M_{du})$ | $s$ | $|X| M_s$ | $1.143 \times 10^{-5}$ |
| $(M_{us}, M_{su})$ | $d$ | $|X| M_d$ | $2.287 \times 10^{-4}$ |
| $(M_{ds}, M_{sd})$ | $u$ | $|X| M_u$ | $4.944 \times 10^{-4}$ |

These masses are $\sim 10^{-5}$--$10^{-4}$ MeV because $|X| = 1.21 \times 10^{-9}$ MeV$^{-4}$ is extremely small.

### Diagonal + X block (4x4)

The entries are:
- $\mathcal{W}[M_{ii}, M_{jj}] = X M_k$ (where $\{i,j,k\}$ distinct)
- $\mathcal{W}[M_{ii}, X] = \prod_{j \neq i} M_j$ (cofactor)
- $\mathcal{W}[X, X] = 0$

$$\mathcal{W}_{4\times4} = \begin{pmatrix}
0 & X M_s & X M_d & M_d M_s \\
X M_s & 0 & X M_u & M_u M_s \\
X M_d & X M_u & 0 & M_u M_d \\
M_d M_s & M_u M_s & M_u M_d & 0
\end{pmatrix}$$

Eigenvalues:

| Eigenvalue | Value (MeV) |
|------------|-------------|
| $\lambda_1$ | $+7.548 \times 10^{-6}$ |
| $\lambda_2$ | $+5.246 \times 10^{-5}$ |
| $\lambda_3$ | $+7.729 \times 10^{10}$ |
| $\lambda_4$ | $-7.729 \times 10^{10}$ |

The heavy mass comes from the cofactor entries (the $M_i M_j$ terms dominate):
$$m_{\text{heavy}} \approx \sqrt{\sum_{i<j} (M_i M_j)^2} = 7.729 \times 10^{10} \text{ MeV}$$

### Complete spectrum

| Mode | Mass (MeV) | Multiplicity | Origin |
|------|-----------|-------------|--------|
| Light 1 | $7.55 \times 10^{-6}$ | 1 | 4x4 block |
| Light 2 | $1.14 \times 10^{-5}$ | 2 | $(M_{ud},M_{du})$ pair |
| Light 3 | $5.25 \times 10^{-5}$ | 1 | 4x4 block |
| Light 4 | $2.29 \times 10^{-4}$ | 2 | $(M_{us},M_{su})$ pair |
| Light 5 | $4.94 \times 10^{-4}$ | 2 | $(M_{ds},M_{sd})$ pair |
| Heavy | $7.73 \times 10^{10}$ | 2 | 4x4 block ($\pm$) |

Total: 10 mass eigenvalues from 10 chiral fields.

---

## 3. Scalar Mass-Squared Matrix

At the SUSY vacuum where all $F_I = \partial W/\partial \Phi_I = 0$, the scalar mass-squared matrix for the complex scalars is:

$$(M^2_{\text{scalar}})_{IJ} = \sum_K \mathcal{W}^*_{IK} \mathcal{W}_{KJ} = (\mathcal{W}^\dagger \mathcal{W})_{IJ}$$

Since $\mathcal{W}$ is real and symmetric: $M^2_{\text{scalar}} = \mathcal{W}^2$.

The scalar mass eigenvalues are $|m_{\text{fermion},n}|$ (the absolute values of the fermion mass eigenvalues). The spectrum is **exactly degenerate** between scalars and fermions, as guaranteed by unbroken SUSY.

---

## 4. Supertrace

$$\text{STr}[M^2] = \sum_n (-1)^{2J_n}(2J_n+1)\, m_n^2 = 2\sum_n m^2_{\text{scalar},n} - 2\sum_n m^2_{\text{fermion},n} = 0$$

This is **identically zero** -- not a numerical result but a consequence of unbroken supersymmetry at the confining vacuum. All higher supertraces $\text{STr}[M^{2k}]$ also vanish.

---

## 5. ISS Regime ($N_f = 4$, $N_c = 3$)

Taking $(u,d,s,c)$ as the active flavors, with $N_c^{\text{mag}} = N_f - N_c = 1$.

Parameters: $h = \sqrt{3}$, $\mu^2 = 8464$ MeV$^2$, $\phi_0 = 92.0$ MeV.

### Tree-level spectrum

**Sector A (pseudo-modulus):**

| State | Mass |
|-------|------|
| Goldstino | 0 |
| Heavy fermion | $\sqrt{2}\, h\, \phi_0 = 225.3$ MeV |
| Pseudo-modulus scalar | $m^2 = 0$ (tree-level flat) |

**Sector B (link fields), per broken flavor $a = d, s, c$:**

| Flavor | $m_q$ (MeV) | $F_a$ (MeV$^2$) | $m^2_{B+}$ (MeV$^2$) | $m^2_{B-}$ (MeV$^2$) | $m_F$ (MeV) |
|--------|-------------|-----------------|----------------------|----------------------|------------|
| $d$ | 4.67 | 14655 | 50772 | 4.35 | 159.3 |
| $s$ | 93.4 | 14567 | 50618 | 158 | 159.3 |
| $c$ | 1270 | 13390 | 48581 | 2196 | 159.3 |

The fermion mass is flavor-universal: $m_F = h\,\phi_0 = 159.3$ MeV.
The scalar masses are split: $m^2_{B\pm} = h^2 \phi_0^2 \pm h\, F_a$.

**Sector C (broken mesons $\chi_a$):**
- Tree-level: **flat** (massless)
- One-loop CW: $m_{\text{CW}} \approx 24$--$27$ MeV for all flavors (nearly degenerate)

### ISS Supertrace

$$\text{STr}[M^2]_{\text{tree}} = (m^2_{B+} + m^2_{B-}) - 2 m_F^2 = 0 \quad \text{(per flavor, exact)}$$

Despite broken SUSY ($F_a \neq 0$), the tree-level supertrace vanishes because the theory has canonical Kahler potential and F-term breaking only.

### Koide transmission

The ISS CW spectrum does **not** transmit Koide: the CW masses are nearly degenerate ($Q_{\text{CW}} \to 1/3$, the degenerate limit) because all quark masses are small compared to $h\mu^2$.

The Seiberg seesaw (SUSY vacuum, not ISS vacuum) **does** transmit Koide via $M_j \propto 1/m_j$ and scale invariance of $Q$.

---

## 6. $N_f = 5$ classification

| $N_f$ | Phase | Mechanism |
|-------|-------|-----------|
| 3 ($= N_c$) | Confining | Quantum-modified moduli space, $\det M = \Lambda^6$ |
| 4 | Free magnetic ($N_c < N_f < 3N_c/2$) | ISS metastable SUSY breaking |
| 5 | Conformal window ($3N_c/2 < N_f < 3N_c$) | Magnetic dual with $N_c^{\text{mag}} = 2$ |

The Seiberg seesaw $M_j = C/m_j$ applies only for $N_f = N_c = 3$. For 5 flavors, the sBootstrap approach uses overlapping $N_c = 3$ blocks: $(u,d,s)$, $(d,s,b)$, $(s,c,b)$, etc.

---

## Key result: Dual Koide uniqueness

Among all $\binom{5}{3} = 10$ flavor triples, the inverse-mass Koide $Q(1/m_i, 1/m_j, 1/m_k)$ is uniquely closest to $2/3$ for the **down-type triple $(d,s,b)$**:

$$Q(1/m_d, 1/m_s, 1/m_b) = 0.6652 \quad (0.22\% \text{ from } 2/3)$$

The next-closest is $(d,s,c)$ at 4.2% deviation. The dual Koide picks out the down-type quarks with a 19:1 margin.

Under the Seiberg seesaw, this becomes a standard Koide condition on the meson VEVs:

$$Q(M_d, M_s, M_b) = 0.6652$$

where $M_j = \Lambda^2 (m_d m_s m_b)^{1/3}/m_j$.

---

## Computation files

- `sqcd_spectrum_final.py` -- main computation (correct Koide convention)
- `sqcd_spectrum_v2.py`, `sqcd_spectrum_calc.py` -- earlier versions (wrong Koide convention in `koide_Q`; the standard Koide formula is $Q = \sum m / (\sum \sqrt{m})^2$, not $Q = (\sum \sqrt{m})^2 / (3\sum m)$)
