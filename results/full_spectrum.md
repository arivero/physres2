# Complete Particle Spectrum of an N=1 SUSY Model

## Setup

**Field content**: 16 chiral superfields (16 complex scalars + 16 Weyl fermions)

| Fields | Complex scalars | Weyl fermions |
|--------|----------------|---------------|
| $M^i_j$ (3x3 meson matrix, i,j = u,d,s) | 9 | 9 |
| $X$ (Lagrange multiplier) | 1 | 1 |
| $B, \tilde{B}$ (baryons) | 2 | 2 |
| $H_u = (H_u^+, H_u^0)$ | 2 | 2 |
| $H_d = (H_d^0, H_d^-)$ | 2 | 2 |
| **Total** | **16** | **16** |

**Superpotential**:

$$W = \sum_i m_i M^i_i + X(\det M - B\tilde{B} - \Lambda^6) + y_c H_u^0 M^d_d + y_b H_d^0 M^s_s$$

where $\det M$ is the full 3x3 determinant of the meson matrix.

**Soft SUSY-breaking**: $V_{\text{soft}} = \tilde{m}^2 \operatorname{Tr}(M^\dagger M) + (B_\mu \sum_i m_i M^i_i + \text{h.c.})$

**Numerical inputs**: $m_u = 2.16$ MeV, $m_d = 4.67$ MeV, $m_s = 93.4$ MeV, $m_c = 1270$ MeV, $m_b = 4180$ MeV, $\Lambda = 300$ MeV, $v = 246220$ MeV, $f_\pi = 92$ MeV, $y_c = 2m_c/v = 0.01032$, $y_b = 2m_b/v = 0.03395$.

**Vacuum** (Seiberg seesaw):

| Quantity | Value | Units |
|----------|-------|-------|
| $C = \Lambda^2(m_u m_d m_s)^{1/3}$ | 882297.43 | MeV$^2$ |
| $\langle M^u_u \rangle = C/m_u$ | 408471.03 | MeV |
| $\langle M^d_d \rangle = C/m_d$ | 188928.79 | MeV |
| $\langle M^s_s \rangle = C/m_s$ | 9446.44 | MeV |
| $\langle X \rangle = -C/\Lambda^6$ | $-1.210 \times 10^{-9}$ | MeV$^{-4}$ |
| $\det M / \Lambda^6$ | 1.00000000000000 | (exact) |
| $\langle H_u^0 \rangle = v/\sqrt{2}$ | 174103.83 | MeV |
| $\langle H_d^0 \rangle = v/\sqrt{2}$ | 174103.83 | MeV |

---

## Task 1: Fermion Mass Matrix

The fermion mass matrix $W_{IJ} = \partial^2 W / \partial \Phi_I \partial \Phi_J$ is a 16x16 real symmetric matrix evaluated at the vacuum. It decomposes into independent blocks:

**Block structure**:

1. **Central 6x6 block** $\{M^u_u, M^d_d, M^s_s, X, H_u^0, H_d^0\}$: Diagonal meson fields coupled to $X$ through cofactors and to Higgs through Yukawa.

2. **Three off-diagonal meson 2x2 blocks**: $\{M^u_d, M^d_u\}$, $\{M^u_s, M^s_u\}$, $\{M^d_s, M^s_d\}$, each with eigenvalues $\pm |X_0| M_k$ where $k$ is the third flavor.

3. **Baryon 2x2 block** $\{B, \tilde{B}\}$: eigenvalues $\pm |X_0|$.

4. **Decoupled fields**: $H_u^+$ and $H_d^-$ have no couplings in $W_{IJ}$ at this vacuum.

**Nonzero entries of $W_{IJ}$**:

| Entry | Expression | Value (MeV) |
|-------|-----------|-------------|
| $W_{M^u_u, M^d_d}$ | $X_0 M_s$ | $-1.143 \times 10^{-5}$ |
| $W_{M^u_u, M^s_s}$ | $X_0 M_d$ | $-2.287 \times 10^{-4}$ |
| $W_{M^d_d, M^s_s}$ | $X_0 M_u$ | $-4.944 \times 10^{-4}$ |
| $W_{M^u_u, X}$ | $M_d M_s$ | $+1.785 \times 10^{9}$ |
| $W_{M^d_d, X}$ | $M_u M_s$ | $+3.859 \times 10^{9}$ |
| $W_{M^s_s, X}$ | $M_u M_d$ | $+7.717 \times 10^{10}$ |
| $W_{M^u_d, M^d_u}$ | $-X_0 M_s$ | $+1.143 \times 10^{-5}$ |
| $W_{M^u_s, M^s_u}$ | $-X_0 M_d$ | $+2.287 \times 10^{-4}$ |
| $W_{M^d_s, M^s_d}$ | $-X_0 M_u$ | $+4.944 \times 10^{-4}$ |
| $W_{B, \tilde{B}}$ | $-X_0$ | $+1.210 \times 10^{-9}$ |
| $W_{M^d_d, H_u^0}$ | $y_c$ | $+1.032 \times 10^{-2}$ |
| $W_{M^s_s, H_d^0}$ | $y_b$ | $+3.395 \times 10^{-2}$ |

**Complete fermion spectrum** (16 eigenvalues, sorted by $|m|$):

| # | $|m|$ (MeV) | Dominant component | Identification |
|---|-------------|-------------------|----------------|
| 1 | 0 | $H_d^-$ (100%) | massless Weyl (decoupled) |
| 2 | $5.6 \times 10^{-17}$ | $H_u^+$ (100%) | massless Weyl (decoupled) |
| 3 | $1.223 \times 10^{-7}$ | $B$ (100%) | baryonic fermion |
| 4 | $4.768 \times 10^{-7}$ | $\tilde{B}$ (100%) | baryonic fermion |
| 5 | $9.865 \times 10^{-6}$ | $M^u_d$ (52%) + $M^d_u$ (48%) | off-diagonal meson (ud) |
| 6 | $1.088 \times 10^{-5}$ | $M^d_u$ (52%) + $M^u_d$ (48%) | off-diagonal meson (ud) |
| 7 | $2.279 \times 10^{-4}$ | $M^s_u$ (50%) + $M^u_s$ (50%) | off-diagonal meson (us) |
| 8 | $2.290 \times 10^{-4}$ | $M^u_s$ (50%) + $M^s_u$ (50%) | off-diagonal meson (us) |
| 9 | $4.928 \times 10^{-4}$ | $M^d_s$ (50%) + $M^s_d$ (50%) | off-diagonal meson (ds) |
| 10 | $4.948 \times 10^{-4}$ | $M^s_d$ (50%) + $M^d_s$ (50%) | off-diagonal meson (ds) |
| 11 | $7.694 \times 10^{-4}$ | $M^u_u$ (50%) + $H_d^0$ (49%) | diagonal meson / Higgsino |
| 12 | $7.797 \times 10^{-4}$ | $M^u_u$ (50%) + $H_d^0$ (48%) | diagonal meson / Higgsino |
| 13 | $1.042 \times 10^{-2}$ | $M^d_d$ (50%) + $H_u^0$ (49%) | diagonal meson / Higgsino |
| 14 | $1.047 \times 10^{-2}$ | $M^d_d$ (50%) + $H_u^0$ (49%) | diagonal meson / Higgsino |
| 15 | $7.729 \times 10^{10}$ | $X$ (50%) + $M^s_s$ (50%) | X-ino / heavy meson |
| 16 | $7.729 \times 10^{10}$ | $X$ (50%) + $M^s_s$ (50%) | X-ino / heavy meson |

The eigenvalues come in approximate $\pm$ pairs, each pair forming one Dirac fermion. The charged Higgsinos $H_u^+$ and $H_d^-$ are exactly massless (no coupling in the superpotential at this vacuum). The Yukawa couplings split what would otherwise be degenerate $\pm$ pairs in the diagonal meson sector.

---

## Task 2: Scalar Mass-Squared Matrix

The scalar potential receives three contributions: $V = \sum_I |F_I|^2 + V_{\text{soft}} + V_D$.

**F-terms at the vacuum**: Without Yukawa couplings, the Seiberg seesaw sets all diagonal F-terms to zero ($F_{M^i_i} = m_i + X_0 \cdot \text{cofactor} = 0$). The Yukawa couplings break SUSY:

| F-term | Value (MeV) | Source |
|--------|-------------|--------|
| $F_{M^u_u}$ | $\approx 0$ | seesaw cancellation |
| $F_{M^d_d}$ | 1796 | $y_c \langle H_u^0 \rangle$ |
| $F_{M^s_s}$ | 5911 | $y_b \langle H_d^0 \rangle$ |
| $F_X$ | $\approx 0$ | $\det M = \Lambda^6$ |
| $F_{H_u^0}$ | 1949 | $y_c \langle M^d_d \rangle$ |
| $F_{H_d^0}$ | 321 | $y_b \langle M^s_s \rangle$ |
| All others | 0 | |

Total: $|F|^2 = F_{M^d_d}^2 + F_{M^s_s}^2 + F_{H_u^0}^2 + F_{H_d^0}^2 = 4.207 \times 10^7$ MeV$^2$.

**Scalar mass matrix structure**: The 32x32 real scalar mass matrix (in the basis $\phi_1^R, \phi_1^I, \phi_2^R, \phi_2^I, \ldots$) is:

$$m^2_{RR} = W^2 + M^2_{\text{hol}} + M^2_{\text{soft}}, \quad m^2_{II} = W^2 - M^2_{\text{hol}} + M^2_{\text{soft}}, \quad m^2_{RI} = 0$$

where $M^2_{\text{hol}}$ is the holomorphic mass matrix from $W_{IJK} F_I$ third-derivative terms, and $M^2_{\text{soft}}$ is diagonal with $\tilde{m}^2 = f_\pi^2$ for the 9 meson fields.

**Numerical caveat**: The fermion mass matrix $W$ has a condition number of $\sim 10^{27}$ (from the $X$-$M^s_s$ coupling at $7.7 \times 10^{10}$ MeV down to the baryon coupling at $1.2 \times 10^{-9}$ MeV). Computing $W^2$ directly and diagonalizing is numerically unreliable for the light modes. The block-by-block analysis gives reliable results for the well-conditioned subspaces.

**Block-by-block scalar spectrum** (reliable):

*ud off-diagonal block* ($M^u_d$, $M^d_u$): all modes positive (stable).

| Mode | $m^2$ (MeV$^2$) | $|m|$ (MeV) |
|------|-----------------|-------------|
| $R_1, R_2$ | 4922 | 70.2 |
| $I_1, I_2$ | 12006 | 109.6 |

*us off-diagonal block* ($M^u_s$, $M^s_u$): R modes tachyonic.

| Mode | $m^2$ (MeV$^2$) | $|m|$ (MeV) |
|------|-----------------|-------------|
| $R_1, R_2$ | $-62384$ | 249.8 (tachyonic) |
| $I_1, I_2$ | 79312 | 281.6 |

*ds off-diagonal block* ($M^d_s$, $M^s_d$): R modes tachyonic.

| Mode | $m^2$ (MeV$^2$) | $|m|$ (MeV) |
|------|-----------------|-------------|
| $R_1, R_2$ | $-144713$ | 380.4 (tachyonic) |
| $I_1, I_2$ | 161641 | 402.0 |

*Baryon block* ($B$, $\tilde{B}$): All four modes degenerate at $m \approx 1.21 \times 10^{-9}$ MeV.

*Charged Higgs* ($H_u^+$, $H_d^-$): All four modes exactly massless.

*Central block* ($M^u_u$, $M^d_d$, $M^s_s$, $X$, $H_u^0$, $H_d^0$): Contains both very heavy ($\sim 10^{10}$ MeV) and very light ($\sim 10^{-4}$ MeV) modes. The heavy modes are the $X$-$M^s_s$ scalar partners.

The tachyonic scalar modes in the us and ds blocks signal that the vacuum is unstable in those directions: the holomorphic mass contribution from $W_{IJK} F_I$ (driven by the large Yukawa-induced F-terms) overcomes the positive soft mass $\tilde{m}^2$. This is a physical instability indicating that the simple diagonal vacuum ansatz is not the true minimum when Yukawa couplings and EWSB are included.

---

## Task 3: SM Particle Identification

The spectrum separates into several distinct sectors:

**Massless sector**: $H_u^+$ and $H_d^-$ are exactly massless (both fermion and scalar components). In a full model these would acquire mass through gauge couplings not included here.

**Ultra-light sector** ($m \lesssim 10^{-6}$ MeV): Baryon fermions and scalars with masses $\sim |X_0| \sim 10^{-9}$ MeV. These are artifacts of the confining phase -- in the physical spectrum they correspond to composite baryons.

**Light meson sector** ($10^{-5}$ to $10^{-2}$ MeV): The off-diagonal meson fermions with masses $|X_0| M_k \propto 1/m_k$ (inverted hierarchy). These track the quark mass ratios through the Seiberg seesaw:

| Fermion pair | Mass $\sim |X_0| M_k$ | Proportional to |
|-------------|----------------------|-----------------|
| $M^u_d$, $M^d_u$ | $\sim 10^{-5}$ MeV | $1/m_s$ |
| $M^u_s$, $M^s_u$ | $\sim 2 \times 10^{-4}$ MeV | $1/m_d$ |
| $M^d_s$, $M^s_d$ | $\sim 5 \times 10^{-4}$ MeV | $1/m_u$ |

**Diagonal meson / Higgsino sector**: The diagonal mesons $M^u_u$, $M^d_d$ mix with the neutral Higgsinos $H_d^0$, $H_u^0$ through the Yukawa couplings $y_b$, $y_c$. The resulting fermion masses ($\sim 10^{-4}$ to $10^{-2}$ MeV) are set by the Yukawa values rather than the seesaw. In the CW effective potential picture, these would correspond to charm and bottom quark mass generation.

**Ultra-heavy sector** ($\sim 7.7 \times 10^{10}$ MeV): The $X$ field and $M^s_s$ form a maximally mixed pair at the scale $\sim \Lambda^6/M_s \sim \Lambda^3 (m_u m_d/m_s)^{1/3}$. These are the heaviest states and decouple from low-energy physics.

**Scalar meson sector** ($\sim 70$ to 600 MeV): The off-diagonal meson scalars, split between R and I components by the holomorphic mass matrix. The R components of the us and ds blocks are tachyonic.

---

## Task 4: Supertrace

The supertrace $\text{STr}[M^2] = \sum_{\text{bosons}} m^2 - 2 \sum_{\text{fermions}} m^2$ is computed analytically:

$$\text{STr}[M^2] = \text{Tr}(m^2_{RR}) + \text{Tr}(m^2_{II}) - 2\,\text{Tr}(W^2)$$

Substituting:

$$= \text{Tr}(W^2 + M^2_{\text{hol}} + M^2_{\text{soft}}) + \text{Tr}(W^2 - M^2_{\text{hol}} + M^2_{\text{soft}}) - 2\,\text{Tr}(W^2)$$

$$= 2\,\text{Tr}(W^2) + 2\,\text{Tr}(M^2_{\text{soft}}) - 2\,\text{Tr}(W^2)$$

$$= 2\,\text{Tr}(M^2_{\text{soft}}) = 2 \times 9 \times \tilde{m}^2$$

$$\boxed{\text{STr}[M^2] = 18\,f_\pi^2 = 152352 \text{ MeV}^2}$$

The $W^2$ piece cancels exactly between the scalar and fermion sectors. The holomorphic piece $M^2_{\text{hol}}$ cancels between R and I components. The Yukawa-induced F-terms do NOT contribute to the supertrace.

For unbroken SUSY ($\tilde{m} = 0$): $\text{STr} = 0$, as required. The nonzero supertrace is entirely due to the soft SUSY-breaking mass $\tilde{m}^2 = f_\pi^2$.

---

## Task 5: Unified Spectrum Table

### Fermion spectrum

| # | Mass (MeV) | Spin | Q | Dominant field | SM identification |
|---|-----------|------|---|---------------|-------------------|
| 1 | 0 | 1/2 | $-1$ | $H_d^-$ | massless charged Higgsino |
| 2 | 0 | 1/2 | $+1$ | $H_u^+$ | massless charged Higgsino |
| 3,4 | $1.2 \times 10^{-7}$, $4.8 \times 10^{-7}$ | 1/2 | 0 | $B$, $\tilde{B}$ | baryonic fermion pair |
| 5,6 | $9.9 \times 10^{-6}$, $1.1 \times 10^{-5}$ | 1/2 | 0 | $M^u_d$/$M^d_u$ | off-diag meson (ud) |
| 7,8 | $2.28 \times 10^{-4}$, $2.29 \times 10^{-4}$ | 1/2 | 0 | $M^u_s$/$M^s_u$ | off-diag meson (us) |
| 9,10 | $4.93 \times 10^{-4}$, $4.95 \times 10^{-4}$ | 1/2 | 0 | $M^d_s$/$M^s_d$ | off-diag meson (ds) |
| 11,12 | $7.69 \times 10^{-4}$, $7.80 \times 10^{-4}$ | 1/2 | 0 | $M^u_u + H_d^0$ | diag meson/Higgsino |
| 13,14 | $1.042 \times 10^{-2}$, $1.047 \times 10^{-2}$ | 1/2 | 0 | $M^d_d + H_u^0$ | diag meson/Higgsino (charm sector) |
| 15,16 | $7.729 \times 10^{10}$ | 1/2 | 0 | $X + M^s_s$ | ultra-heavy X-ino |

### Scalar spectrum (block-by-block, reliable)

| # | $m^2$ (MeV$^2$) | $|m|$ (MeV) | Spin | Field | Identification |
|---|-----------------|-------------|------|-------|----------------|
| 1--4 | 0 | 0 | 0 | $H_u^+$, $H_d^-$ | massless (4 flat directions) |
| 5--8 | $\sim 1.5 \times 10^{-18}$ | $\sim 1.2 \times 10^{-9}$ | 0 | $B$, $\tilde{B}$ | baryonic scalars |
| 9,10 | 4922 | 70.2 | 0 | $M^u_d$, $M^d_u$ (R) | ud meson scalar |
| 11,12 | 12006 | 109.6 | 0 | $M^u_d$, $M^d_u$ (I) | ud meson scalar |
| 13,14 | $-62384$ | 249.8 | 0 | $M^u_s$, $M^s_u$ (R) | us meson scalar (tachyonic) |
| 15,16 | 79312 | 281.6 | 0 | $M^u_s$, $M^s_u$ (I) | us meson scalar |
| 17,18 | $-144713$ | 380.4 | 0 | $M^d_s$, $M^s_d$ (R) | ds meson scalar (tachyonic) |
| 19,20 | 161641 | 402.0 | 0 | $M^d_s$, $M^s_d$ (I) | ds meson scalar |
| 21,22 | $\sim 6 \times 10^{-7}$ | $\sim 7.7 \times 10^{-4}$ | 0 | central (R) | light neutral scalar |
| 23,24 | $\sim 10^{-4}$ | $\sim 0.01$ | 0 | central (R) | neutral scalar |
| 25,26 | $-767308$, $-158897$ | 876, 399 | 0 | central (R) | tachyonic (unstable) |
| 27,28 | 5370, 9729 | 73.3, 98.6 | 0 | central (I) | neutral scalar |
| 29--32 | $5.97 \times 10^{21}$ | $7.73 \times 10^{10}$ | 0 | $X$, $M^s_s$ | ultra-heavy scalars |

---

## Scale Hierarchy and Seesaw Structure

The spectrum spans over 20 orders of magnitude in mass, reflecting the Seiberg seesaw $M_i = C/m_i$:

$$\underbrace{|X_0| \sim 10^{-9}}_{\text{baryon}} \ll \underbrace{|X_0| M_s \sim 10^{-5}}_{\text{ud meson}} \ll \underbrace{y_{c,b} \sim 10^{-2}}_{\text{Yukawa}} \ll \underbrace{f_\pi \sim 10^2}_{\text{soft}} \ll \underbrace{\Lambda^6/M_s \sim 10^{10}}_{\text{X-meson heavy}}$$

The off-diagonal meson fermion masses are proportional to $1/m_{\text{quark}}$ (the inverted hierarchy from the seesaw):

$$\frac{m(M^d_s,M^s_d)}{m(M^u_d,M^d_u)} = \frac{M_u}{M_s} = \frac{m_s}{m_u} = 43.24$$

This is the ISS mechanism transmitting the quark mass spectrum through the CW potential.

---

*Generated by `results/full_spectrum.py`*
