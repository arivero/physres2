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

where $C^3 = \Lambda^6 \prod_k m_k$ and $M_j = C/m_j$. The cubic has three roots
$C_k = \omega^k C_0$ with $\omega = e^{2\pi i/3}$ and $C_0 = \Lambda^2(\prod_k m_k)^{1/3}$:
these are the three supersymmetric vacua. The tables below use the real root $C_0$.
All three vacua leave to infinity in the chiral limit $m_u = m_d \to 0$ (Section 7).

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

## 7. Chiral limit $m_u = m_d = 0$: no finite vacuum, Affleck–Dine–Seiberg runaway (issue #2)

**Status.** Known physics (Affleck–Dine–Seiberg 1984; Seiberg 1994;
Intriligator–Seiberg, hep-th/9509066, §§4.1–4.2). It corrects the earlier
reading of this limit as a dynamical O'Raifeartaigh mechanism with a
Goldstino at a finite vacuum. The soft-mass stabilization in 7.5 is a
standard construction evaluated here; it rests on the stated extra input.
Script: `calculations/sqcd_chiral_limit_runaway.py` (Python 3, sympy + numpy).
Manuscript: `sbootstrap_v4d.tex`, subsection "The chiral limit of SQCD and
the Affleck–Dine–Seiberg runaway" (`subsec:sqcd_runaway`).

### 7.1 Constrained rank-one F-term system

$$W = m_s M_{33} + X(\det M - B\tilde B - \Lambda_3^6),\qquad m_s \neq 0 .$$

F-terms:

- $W_{M_{ij}} = m_s\,\delta_{i3}\delta_{j3} + X\,\mathrm{cof}(M)_{ij}$ (nine equations)
- $W_B = -X\tilde B$, $W_{\tilde B} = -XB$
- $W_X = \det M - B\tilde B - \Lambda_3^6$

No finite common zero. $W_{M_{33}} = 0$ requires $X \neq 0$. Then
$W_B = W_{\tilde B} = 0$ give $B = \tilde B = 0$, and $W_X = 0$ gives
$\det M = \Lambda_3^6 \neq 0$. The cofactor matrix
$\Lambda_3^6 (M^{-1})^T$ has rank 3 and cannot equal the rank-1 matrix
$-m_s E_{33}/X$.

Ideal-membership certificate (verified symbolically by the script):

$$m_s\Lambda_3^6 = \Lambda_3^6 W_{M_{33}} - \mathrm{cof}(M)_{33}\Big(\sum_j M_{1j}W_{M_{1j}} + B\,W_B - X\,W_X\Big).$$

It uses $\sum_j M_{1j}\,\mathrm{cof}(M)_{1j} = \det M$ (row-1 Laplace expansion;
row 1 carries no mass term).

### 7.2 The runaway trajectory: all F-terms

$$M = \mathrm{diag}(t,\,t,\,\Lambda_3^6/t^2),\qquad B = \tilde B = 0,\qquad X = -m_s/t^2 .$$

| F-term | value on the trajectory |
|---|---|
| $W_{M_{11}}$, $W_{M_{22}}$ | $-m_s\Lambda_3^6/t^3$ |
| $W_{M_{33}}$ | $m_s + X t^2 = 0$ |
| $W_{M_{ij}}$, $i \neq j$ | 0 (cofactors of a diagonal matrix) |
| $W_X$ | $t\cdot t\cdot\Lambda_3^6/t^2 - \Lambda_3^6 = 0$ |
| $W_B$, $W_{\tilde B}$ | 0 |

### 7.3 Large-field potential in weakly Higgsed coordinates

For $|t| \gg \Lambda_3^2$ the light squarks have D-flat VEVs
$Q = \tilde Q^T = \sqrt{t}\,\mathrm{diag}(1,1)$ embedded in colour. Two
flavours Higgs $SU(3)$ completely, at weak coupling. The induced light-meson
Kähler potential is $K = 2\,\mathrm{Tr}\sqrt{M^\dagger M}$, so
$K_{t\bar t} = 1/(2|t|)$ and $g^{11} = g^{22} = 2|t|$.

Constrained description. The physical F-term is the projection of
$dW_{\rm tree} = m_s\,dM_{33}$ onto the tangent space of $\det M = \Lambda_3^6$,
$V = |a|^2 - |\langle a,n\rangle|^2/|n|^2$ with $n = d\det M$. With
$|n_{11}| = |n_{22}| = |\Lambda_3|^6/|t|$ and $|n_{33}| = |t|^2$:

$$V = \frac{4|m_s|^2|\Lambda_3|^{12}}{|t|^5}\left[1 + \frac{4|\Lambda_3|^{12}}{g^{33}|t|^5}\right]^{-1} \to \frac{4|m_s|^2|\Lambda_3|^{12}}{|t|^5}.$$

The limit holds for any $g^{33}$ that falls more slowly than $|t|^{-5}$;
the script checks $g^{33} = 1$ and $g^{33} = 2|M_{33}|$. A canonical meson
metric on the light entries would give $2|m_s|^2|\Lambda_3|^{12}/|t|^6$:
meson-coordinate derivatives alone give the wrong power.

Two-flavour ADS description. Integrating out $M_{33}$ and $X$
holomorphically on $B = 0$, $M_{3a} = M_{a3} = 0$:

$$W_{\rm eff} = m_s M_{33}\big|_{\det M = \Lambda_3^6} = \frac{\Lambda_2^7}{\det M_{\rm light}},\qquad \Lambda_2^7 = m_s\Lambda_3^6 .$$

In squark coordinates each of the four nonzero squarks has
$|\partial W/\partial Q| = |\Lambda_2|^7/v^5$, so
$V = 4|\Lambda_2|^{14}/v^{10}$ with $v^2 = |t|$. This equals the constrained
result. On the D-flat manifold $V$ is homogeneous of degree $-10$ in $v$, so
it has no stationary point in the weakly coupled region. $V > 0$ at every
finite point and $V \to 0$ along the runaway. Near $|t| \sim \Lambda_3^2$
the Kähler metric is not calculable; no local minimum there has a derivation.

### 7.4 Massive control: three branches, equal-light-mass limit

With all three masses nonzero, $X\,\mathrm{cof}(M) = -\mathrm{diag}(m_i)$
gives $X\det M\,(M^{-1})^T = -\mathrm{diag}(m_i)$, so $M$ is diagonal with
$m_i M_{ii} = C \equiv -X\Lambda_3^6$ and $C^3 = \Lambda_3^6 m_u m_d m_s$.
For $m_u = m_d = m$:

$$M_{uu} = M_{dd} = t_k = \omega^k\Lambda_3^2 (m_s/m)^{1/3},\qquad M_{ss} = \Lambda_3^6/t_k^2 = \omega^k\Lambda_3^2(m/m_s)^{2/3},\qquad X = -m_s/t_k^2 .$$

Every branch lies exactly on the runaway trajectory of 7.2 at $t = t_k$ and is
the two-flavour ADS vacuum $t_k^3 = \Lambda_2^7/m$. As $m \to 0$ at fixed
$m_s$, the light meson VEVs $M_{uu} = M_{dd}$ diverge as $m^{-1/3}$ on all
three branches; $M_{ss}$ and $X$ vanish as $m^{2/3}$.

Numbers ($\Lambda_3 = 300$ MeV, $m_s = 93.4$ MeV; identical moduli on the
three branches, phases $\omega^k$; F-term residuals $< 10^{-14}$):

| $m$ (MeV) | $\lvert M_{uu}\rvert = \lvert M_{dd}\rvert$ (MeV²) | $\lvert M_{ss}\rvert$ (MeV²) | $\lvert X\rvert$ (MeV⁻⁴) |
|---|---|---|---|
| 4.67 | $2.443\times10^5$ | $1.221\times10^4$ | $1.565\times10^{-9}$ |
| 1 | $4.083\times10^5$ | $4.372\times10^3$ | $5.601\times10^{-10}$ |
| $10^{-1}$ | $8.797\times10^5$ | $9.419\times10^2$ | $1.207\times10^{-10}$ |
| $10^{-2}$ | $1.895\times10^6$ | $2.029\times10^2$ | $2.600\times10^{-11}$ |
| $10^{-3}$ | $4.083\times10^6$ | $4.372\times10^1$ | $5.601\times10^{-12}$ |
| $10^{-4}$ | $8.797\times10^6$ | $9.419$ | $1.207\times10^{-12}$ |
| $10^{-6}$ | $4.083\times10^7$ | $0.4372$ | $5.601\times10^{-14}$ |

Fitted exponents: $M_{uu} \sim m^{-0.3333}$, $M_{ss} \sim m^{0.6667}$,
$X \sim m^{0.6667}$.

### 7.5 What has a derivation; a stabilized deformation

- $m_u = m_d = 0$, $m_s \neq 0$: no vacuum at finite fields; runaway.
- $m_u, m_d, m_s \neq 0$: three finite supersymmetric vacua ($F = 0$).
- Neither case contains a Goldstino. The manuscript statements placing a
  Goldstino or $\langle F\rangle \neq 0$ at a broken vacuum in this chiral
  limit are now marked as sBootstrap hypotheses.

A finite vacuum requires extra input. Controlled example: soft squark masses
$\tilde m^2(|Q|^2 + |\tilde Q|^2)$ for the two light flavours. On the D-flat
direction $V = 4|\Lambda_2|^{14}/v^{10} + 4\tilde m^2 v^2$, stationary at

$$v_0^{12} = 5|m_s|^2|\Lambda_3|^{12}/\tilde m^2 .$$

Weak coupling needs $v_0 \gg \Lambda_3$, i.e. $\tilde m \ll m_s$
($v_0/\Lambda_3 = 2.47$ at $\tilde m = m_s/100$). The script builds the full
two-flavour theory (12 complex squarks, $SU(3)$ D-terms with $g = 1$,
8 gauginos) at $\tilde m = 10^{-2} m_s$ and finds:

| quantity | result |
|---|---|
| stationarity | $\max\lvert\partial V\rvert = 2\times10^{-16}$, D-flat |
| scalar zero modes | 12 = 8 eaten + 3 pions ($SU(2)_L\times SU(2)_R\to SU(2)_V$) + 1 $R$-axion |
| scalar triplet | $12\tilde m^2/5$ |
| radial singlet | $12\tilde m^2$ |
| heavy scalars | $\mathcal O(g^2v_0^2)$ (vector multiplets) |
| negative modes | none |
| light fermions | $\tilde m/\sqrt5$ (triplet), $\sqrt5\,\tilde m$ (singlet) |
| massless fermions | none (explicit breaking: no Goldstino) |

The Hessian is cross-checked against finite differences. In this vacuum
$m_\pi^2 - m_\psi^2 = -\tilde m^2/5$, opposite in sign to the manuscript's
$m_\pi^2 - m_\mu^2 = f_\pi^2$. For $\tilde m \sim f_\pi \approx m_s$,
$v_0/\Lambda_3 = 1.15$: the minimum sits at strong coupling and is not
controlled. No stabilized finite vacuum of the chiral-limit theory is
claimed.

---

## Computation files

- `sqcd_spectrum_final.py` -- main computation (correct Koide convention)
- `sqcd_chiral_limit_runaway.py` -- Section 7: rank-one F-term certificate, runaway potential, ADS matching, massive branches, soft-mass stabilization
- `sqcd_spectrum_v2.py`, `sqcd_spectrum_calc.py` -- earlier versions (wrong Koide convention in `koide_Q`; the standard Koide formula is $Q = \sum m / (\sum \sqrt{m})^2$, not $Q = (\sum \sqrt{m})^2 / (3\sum m)$)
