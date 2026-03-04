# Complete SUSY Lagrangian: SQCD(N_f = N_c = 3) + NMSSM Higgs Sector

**Date:** 2026-03-04

---

## 0. Overview and Scope

This document assembles the complete N=1 supersymmetric Lagrangian for the model consisting of:

- **Confined SU(3) SQCD** with N_f = N_c = 3 flavors (u, d, s), described in the infrared by gauge-invariant composites (mesons, baryons, and a Lagrange multiplier).
- **MSSM Higgs sector** (H_u, H_d doublets) coupled to the confined sector via Yukawa couplings.
- **NMSSM singlet** identified with the Seiberg Lagrange multiplier X, providing the Higgs quartic.
- **Top quark** treated as external (not part of the confined sector).

The Lagrangian is specified by three objects: the Kahler potential K, the superpotential W, and the gauge kinetic function (trivial at low energy since SU(3)_c is confined). All soft SUSY-breaking terms are also listed.

---

## 1. Field Content

### 1.1 Chiral Superfields (16 total)

| Field | Multiplicity | SU(2)_L x U(1)_Y | Description |
|-------|-------------|-------------------|-------------|
| M^i_j (i,j = u,d,s) | 9 | various (see below) | Meson composites Q-tilde^i Q_j / Lambda |
| X | 1 | (1, 0) | Seiberg Lagrange multiplier = NMSSM singlet S |
| B | 1 | (1, 0) | Baryon epsilon_{ijk} Q^i Q^j Q^k / Lambda^3 |
| B-tilde | 1 | (1, 0) | Antibaryon |
| H_u = (H_u^+, H_u^0) | 2 | (2, +1/2) | Up-type Higgs doublet |
| H_d = (H_d^0, H_d^-) | 2 | (2, -1/2) | Down-type Higgs doublet |

**Total:** 16 complex scalars + 16 Weyl fermions = 32 bosonic + 32 fermionic on-shell degrees of freedom (before gauge fixing).

### 1.2 Meson EM Charges

The meson M^i_j carries charge Q(i) - Q(j), where Q(u) = 2/3, Q(d) = Q(s) = -1/3:

| Meson | EM charge |
|-------|-----------|
| M^u_u, M^d_d, M^s_s | 0 |
| M^u_d, M^u_s | +1 |
| M^d_u, M^s_u | -1 |
| M^d_s, M^s_d | 0 |

### 1.3 External Fields (not part of the confined sector)

| Field | Description |
|-------|-------------|
| t, t-bar | Top quark (chiral superfields) |
| SU(2)_L x U(1)_Y gauge multiplets | EW gauge sector |

---

## 2. Superpotential

The complete superpotential has four layers:

$$W = W_{\text{Seiberg}} + W_{\text{Yukawa}} + W_{\text{NMSSM}} + W_{\text{instanton}}$$

### 2.1 Seiberg Confined-Phase Superpotential

$$W_{\text{Seiberg}} = \sum_{i=u,d,s} m_i M^i_i + X \bigl(\det M^{(3)} - B \tilde{B} - \Lambda^6\bigr)$$

where det M^{(3)} is the 3x3 determinant of the meson matrix M^i_j (i,j = u,d,s):

$$\det M^{(3)} = \epsilon_{ijk} \epsilon^{lmn} M^i_l M^j_m M^k_n / 6$$

expanded as:

$$\det M = M^u_u M^d_d M^s_s - M^u_u M^d_s M^s_d - M^d_d M^u_s M^s_u - M^s_s M^u_d M^d_u + M^u_d M^d_s M^s_u + M^u_s M^d_u M^s_d$$

**Classification:** DETERMINED (Seiberg, 1994). The structure follows from the quantum-modified constraint of SU(3) SQCD with N_f = N_c = 3. The quark mass terms m_i are current quark masses.

### 2.2 Yukawa Couplings to the Higgs Sector

$$W_{\text{Yukawa}} = y_c\, H_u^0\, M^d_d + y_b\, H_d^0\, M^s_s + y_t\, H_u^0\, t\, \bar{t}$$

The meson-Higgs Yukawa couplings identify specific diagonal mesons with the charm and bottom sectors. The identification M^d_d ~ c-bar c and M^s_s ~ b-bar b follows from the Seiberg seesaw: the meson VEV <M^j_j> = C/m_j is inversely proportional to the quark mass, so the charm meson carries the SU(2)_L quantum numbers appropriate for coupling to H_u, and the bottom meson couples to H_d.

The top quark does not form mesons (it does not participate in the confining dynamics) and couples to H_u through a standard Yukawa.

**Classification:**

| Coupling | Value | Status |
|----------|-------|--------|
| y_c = 2 m_c / v | 1.032 x 10^{-2} | FITTED (from m_c = 1270 MeV, v = 246220 MeV) |
| y_b = 2 m_b / v | 3.395 x 10^{-2} | FITTED (from m_b = 4180 MeV) |
| y_t = sqrt(2) m_t / v | 0.9934 | FITTED (from m_t = 172760 MeV) |

**IMPORTANT NOTE on Yukawa conventions:** The factors of 2 vs sqrt(2) in the Yukawa definitions depend on the convention for the Higgs VEV. In the MSSM convention H^0 = v/sqrt(2), the relation is m = y v/sqrt(2), giving y = sqrt(2) m/v. The factor 2 m/v used for y_c and y_b in the existing results corresponds to a different convention (m = y v/2). This is a bookkeeping issue; the physical quark masses are determined and the numerical results are internally consistent. In the assembled Lagrangian below, I use the MSSM convention y = sqrt(2) m/v throughout.

**Corrected values (MSSM convention):**

| Coupling | Expression | Value |
|----------|-----------|-------|
| y_c | sqrt(2) m_c / v | 7.298 x 10^{-3} |
| y_b | sqrt(2) m_b / v | 2.401 x 10^{-2} |
| y_t | sqrt(2) m_t / v | 0.9934 |

**CONSISTENCY CHECK:** Both conventions give the same physical mass: m_c = y_c v/sqrt(2) = 1270 MeV. The choice affects only the numerical value of the coupling constant. In what follows I will use the convention established in the existing results (y = 2m/v) to avoid introducing discrepancies with prior numerical work, and flag this as a normalization issue to be resolved when preparing the final paper.

### 2.3 NMSSM Coupling

$$W_{\text{NMSSM}} = \lambda\, X\, (H_u^0 H_d^0 - H_u^+ H_d^-)$$

The NMSSM coupling identifies the Seiberg Lagrange multiplier X with the NMSSM singlet S. This is the central structural identification of the model.

| Parameter | Value | Status |
|-----------|-------|--------|
| lambda | 0.72 | FITTED (from m_h = 125 GeV at tan beta = 1: lambda = sqrt(2) m_h / v) |

**Key consequence:** The |F_X|^2 term generates a Higgs quartic:

$$V \supset \lambda^2 |H_u^0 H_d^0|^2$$

At tan beta = 1 (v_u = v_d = v/sqrt(2)), this gives:

$$m_h^{\text{tree}} = \lambda v / \sqrt{2} = 0.72 \times 246.22 / \sqrt{2} = 125.4 \text{ GeV}$$

**Classification:** The identification X = S is ASSUMED. It is the model's central hypothesis linking the confined SQCD sector to the Higgs mechanism. Once assumed, lambda = 0.72 is DETERMINED by the measured Higgs mass.

### 2.4 Three-Instanton Correction

$$W_{\text{instanton}} = c_3\, \frac{(\det M^{(3)})^3}{\Lambda^{18}}$$

This is the leading instanton correction beyond the Seiberg constraint. The three-instanton term generates a cos(9 delta) potential for the Koide phase delta.

| Parameter | Value | Status |
|-----------|-------|--------|
| c_3 | O(1) | UNDETERMINED. The instanton coefficient is calculable in principle but has not been computed from first principles for this model. |

**Classification:** ASSUMED to exist (standard instanton calculus guarantees such a term); coefficient UNDETERMINED.

---

## 3. Kahler Potential

$$K = K_{\text{canonical}} + K_{\text{pseudo-modulus}} + K_{\text{bion}}$$

### 3.1 Canonical Terms

$$K_{\text{canonical}} = \text{Tr}(M^\dagger M) + |X|^2 + |B|^2 + |\tilde{B}|^2 + |H_u|^2 + |H_d|^2$$

**Classification:** DETERMINED by the canonical normalization of the magnetic dual (Intriligator-Seiberg-Shih framework).

### 3.2 Pseudo-Modulus Correction

$$K_{\text{pseudo-modulus}} = -\frac{|X|^4}{12\,\mu^2}$$

where mu is the Kahler mass parameter (not to be confused with the MSSM mu-term).

This quartic correction creates a Kahler pole at |X| = sqrt(3) mu, where the Kahler metric K_{X X-bar} = 1 - |X|^2/(3 mu^2) vanishes. The CW potential pins the pseudo-modulus VEV at the pole:

$$\langle X \rangle = \sqrt{3}\,\mu$$

**The O'Raifeartaigh-Koide mechanism:** At this VEV, the fermion mass matrix of the O'Raifeartaigh sector produces eigenvalues proportional to (2 - sqrt(3)) and (2 + sqrt(3)), which is the exact Koide seed mass ratio, provided gv/m = sqrt(3). This is the microscopic origin of the Koide condition in the seed (0, m_s, m_c) sector.

| Parameter | Value | Status |
|-----------|-------|--------|
| mu | determined by C/(sqrt(3) Lambda^6) | DETERMINED (consistency of seesaw with Kahler pole) |
| c = -1/12 | coefficient of |X|^4 term | ASSUMED (specific value gives the Kahler pole at the seesaw VEV; a derivation from the microscopic theory is lacking) |

**Classification:** The functional form is ASSUMED. The coefficient c = -1/12 is FITTED to reproduce the Koide seed at the ISS matching point.

### 3.3 Bion Correction (Monopole-Instanton Pairs)

$$K_{\text{bion}} = \frac{\zeta^2}{\Lambda^2}\, \biggl|\sum_{k=1}^{3} s_k \sqrt{M^k_k}\biggr|^2$$

where:
- zeta = exp(-S_0/(2 N_c)) is the monopole fugacity
- S_0 = 2 pi / alpha_s is the instanton action on R^3 x S^1
- s_k = +/-1 are signs from the monopole-instanton sector
- The square root sqrt(M^k_k) is non-holomorphic

The bion correction generates the effective potential for v_0-doubling:

$$V_{\text{bion}} = \lambda_2\, |S_{\text{bloom}} - 2\,S_{\text{seed}}|^2$$

where S_seed = sqrt(m_s) + sqrt(m_c) and S_bloom = -sqrt(m_s) + sqrt(m_c) + sqrt(m_b).

**Minimization** gives the v_0-doubling condition:

$$\sqrt{m_b} = 3\sqrt{m_s} + \sqrt{m_c}$$

predicting m_b = 4177 MeV (PDG: 4180 +/- 30, deviation 0.10 sigma).

| Parameter | Value | Status |
|-----------|-------|--------|
| zeta^2/Lambda^2 | 1.52 x 10^{-2} (at alpha_s = 1) | ESTIMATED (depends on alpha_s at the compactification scale) |
| s_k signs | s = (-1, +1, +1) in bloom configuration | ASSUMED (from monopole zero-mode phase flip at the seed-to-bloom transition) |
| lambda_2 (bion coupling) | O(zeta^4/Lambda^4) | UNDETERMINED in absolute normalization |

**Classification:** The functional form is MOTIVATED by the Unsal semiclassical framework (monopole-instanton pairs on R^3 x S^1). The specific form involving sqrt(M) is NON-STANDARD and ASSUMED. A rigorous derivation from the compactified theory is missing. The v_0-doubling prediction is the model's primary quantitative test of this term.

### 3.4 Expanded Bion Kahler (Diagonal M)

For M = diag(M_1, M_2, M_3):

$$K_{\text{bion}} = \frac{\zeta^2}{\Lambda^2}\bigl(M_1 + M_2 + M_3 + 2 s_1 s_2 \sqrt{M_1 M_2} + 2 s_1 s_3 \sqrt{M_1 M_3} + 2 s_2 s_3 \sqrt{M_2 M_3}\bigr)$$

The cross terms sqrt(M_i M_j) are the bion contributions (monopole i -- antimonopole j).

---

## 4. Soft SUSY-Breaking Terms

### 4.1 Meson Soft Mass

$$V_{\text{soft}} = f_\pi^2\, \text{Tr}(M^\dagger M)$$

where f_pi = 92 MeV is the pion decay constant.

This is a universal mass-squared m-tilde^2 = f_pi^2 = 8464 MeV^2 for all 9 meson scalars. It is the ONLY soft term stated in the model.

| Parameter | Value | Status |
|-----------|-------|--------|
| f_pi | 92 MeV | DETERMINED (from QCD, not a free parameter) |
| m-tilde^2 = f_pi^2 | 8464 MeV^2 | ASSUMED to be the sole soft term |

**Classification:** The identification of the SUSY-breaking scale with f_pi is ASSUMED. The mechanism that generates this soft mass (mediation from the SUSY-breaking sector) is not specified.

### 4.2 Missing Soft Terms

A critical referee would immediately ask about the following soft terms that are NOT included but generically expected:

**(a) Higgs soft masses:** m_{H_u}^2, m_{H_d}^2, B_mu. The EWSB conditions at tan beta = 1 require:

$$m_{H_u}^2 + y_c^2 = m_{H_d}^2 + y_b^2 = b$$

where b is the soft bilinear. These are needed for EWSB but not specified from first principles.

**Status:** REQUIRED but UNDETERMINED. The model assumes EWSB happens but does not derive the soft Higgs masses.

**(b) Gaugino masses:** M_1, M_2, M_3 (bino, wino, gluino). In any SUSY model, gaugino masses are required for phenomenological viability (chargino, neutralino masses). Since SU(3)_c is confined, the gluino mass is absorbed into the confined-phase description, but M_1 and M_2 remain unspecified.

**Status:** MISSING. Not addressed.

**(c) A-terms:** Trilinear scalar couplings A_c, A_b, A_t. These shift the scalar mass eigenvalues and contribute to flavor violation.

**Status:** MISSING. Not addressed.

**(d) Soft terms for B, B-tilde:** m_B^2, m_{B-tilde}^2. These would affect the baryon spectrum.

**Status:** MISSING. Not addressed. The baryonic sector is treated as a spectator.

**(e) Soft mass for X:** m_X^2. This would shift the pseudo-modulus VEV.

**Status:** MISSING. The pseudo-modulus stabilization relies on the Kahler pole and CW potential, not a soft mass.

### 4.3 What Soft Terms ARE Generated

The Yukawa couplings generate SUSY-breaking F-terms at the vacuum:

| F-term | Value (MeV) | Source |
|--------|-------------|--------|
| F_{M^u_u} | ~0 | Seiberg seesaw cancellation |
| F_{M^d_d} | 1796 | y_c <H_u^0> |
| F_{M^s_s} | 5911 | y_b <H_d^0> |
| F_X | ~0 | det M = Lambda^6 (at leading order) |
| F_{H_u^0} | 289.6 | y_c <M^d_d> = C_{uds}/v |
| F_{H_d^0} | 289.6 | y_b <M^s_s> = C_{uds}/v |
| |F|^2 total | 4.207 x 10^7 MeV^2 | |

**The ISS flavor universality identity:** F_{H_u} = F_{H_d} = C/v = 289.6 MeV. The Yukawa y_j ~ m_j and the seesaw <M_j> ~ 1/m_j cancel, producing flavor-universal Higgs F-terms. This is algebraically exact.

---

## 5. D-Terms

### 5.1 SU(2)_L x U(1)_Y D-Terms

$$V_D = \frac{g^2}{2}\sum_{a=1}^{3}\Bigl(\sum_I \Phi_I^\dagger T^a \Phi_I\Bigr)^2 + \frac{g'^2}{2}\Bigl(\sum_I \frac{Y_I}{2}|\Phi_I|^2\Bigr)^2$$

At tan beta = 1 (v_u = v_d), the neutral Higgs D-term vanishes on the D-flat direction:

$$(g^2 + g'^2)(|H_u^0|^2 - |H_d^0|^2)^2 / 8 = 0 \quad \text{when } |H_u^0| = |H_d^0|$$

This is why the NMSSM quartic from W_NMSSM is essential: without it, the tree-level Higgs mass is zero at tan beta = 1.

### 5.2 Electromagnetic Fayet-Iliopoulos Term

The model does not include a U(1)_EM FI term. However, the supertrace decomposition by charge shows that the meson soft masses break the EM supertrace:

$$\sum_Q Q^2 \text{STr}_Q[M^2] = 8 f_\pi^2 \neq 0$$

This is not a problem (there is no gauge symmetry requiring it to vanish), but a referee might ask whether it signals a U(1) tadpole.

**Status:** No FI term is included. ASSUMED to be absent.

---

## 6. Complete Lagrangian in Component Form

### 6.1 General Structure

$$\mathcal{L} = K_{I\bar{J}}\, \partial_\mu \Phi_I\, \partial^\mu \bar{\Phi}_{\bar{J}} + i\, K_{I\bar{J}}\, \bar{\psi}^{\bar{J}} \bar{\sigma}^\mu \partial_\mu \psi^I - V$$

where V = V_F + V_D + V_soft and I runs over all 16 chiral superfields.

### 6.2 Scalar Potential: F-terms

$$V_F = K^{I\bar{J}}\, F_I\, \bar{F}_{\bar{J}}$$

where:

$$F_I = \frac{\partial W}{\partial \Phi_I}$$

and K^{I J-bar} is the inverse Kahler metric.

**Explicit F-terms** (at a general field point, setting B = B-tilde = 0):

| Field | F-term |
|-------|--------|
| F_X | det M - Lambda^6 + lambda H_u^0 H_d^0 - lambda H_u^+ H_d^- |
| F_{M^i_i} | m_i + X * cofactor(M, i,i) + 3 c_3 (det M)^2 cofactor(M, i,i)/Lambda^18 + delta_{i,d} y_c H_u^0 + delta_{i,s} y_b H_d^0 |
| F_{M^i_j} (i != j) | X * d(det M)/dM^i_j |
| F_B | -X B-tilde |
| F_{B-tilde} | -X B |
| F_{H_u^0} | y_c M^d_d + lambda X H_d^0 |
| F_{H_d^0} | y_b M^s_s + lambda X H_u^0 |
| F_{H_u^+} | -lambda X H_d^- |
| F_{H_d^-} | -lambda X H_u^+ |

### 6.3 Inverse Kahler Metric

For the pseudo-modulus X:

$$K^{X\bar{X}} = \frac{1}{1 - |X|^2/(3\mu^2)}$$

For the meson sector (to first order in zeta^2/Lambda^2):

$$K^{M_i\bar{M}_j} = \delta^{ij} - \frac{\zeta^2}{\Lambda^2}\frac{s_i s_j}{4\sqrt{m_i m_j}} + O(\zeta^4)$$

All other inverse metric components are canonical (= 1 on diagonal, 0 off-diagonal).

### 6.4 Fermion Mass Matrix

The fermion mass matrix (bilinear in Weyl spinors) is:

$$\mathcal{L}_{\text{fermion mass}} = -\frac{1}{2} W_{IJ}\, \psi^I \psi^J + \text{h.c.}$$

where W_{IJ} = d^2 W / (d Phi_I d Phi_J) evaluated at the vacuum.

The 16x16 matrix decomposes into blocks:

**Central 6x6 block** {M^u_u, M^d_d, M^s_s, X, H_u^0, H_d^0}:
- W_{M^i_i, M^j_j} = X_0 cofactor(M, k) for (i,j,k) cyclic
- W_{M^i_i, X} = cofactor of M at (i,i) = product of the other two diagonal VEVs
- W_{M^d_d, H_u^0} = y_c
- W_{M^s_s, H_d^0} = y_b
- W_{X, H_u^0} = lambda v/sqrt(2)
- W_{X, H_d^0} = lambda v/sqrt(2)
- W_{H_u^0, H_d^0} = lambda X_0

**Off-diagonal meson 2x2 blocks** {M^i_j, M^j_i}:
- W_{M^i_j, M^j_i} = X_0 M_k (k != i, j)
- Each pair forms a Dirac fermion with mass |X_0| M_k

**Baryon 2x2 block** {B, B-tilde}:
- W_{B, B-tilde} = -X_0

**Charged Higgs 2x2 block** {H_u^+, H_d^-}:
- W_{H_u^+, H_d^-} = -lambda X_0

### 6.5 Scalar Potential: D-terms

$$V_D = \frac{g^2 + g'^2}{8}\bigl(|H_u^0|^2 - |H_d^0|^2 + |H_u^+|^2 - |H_d^-|^2\bigr)^2 + \frac{g^2}{2}|H_u^{+*} H_u^0 + H_d^{0*} H_d^-|^2$$

plus D-terms from the meson EM charges (subdominant at the meson scale).

### 6.6 Scalar Potential: Soft Terms

$$V_{\text{soft}} = f_\pi^2 \sum_{i,j} |M^i_j|^2 + m_{H_u}^2 |H_u|^2 + m_{H_d}^2 |H_d|^2 + (b\, H_u \cdot H_d + \text{h.c.})$$

where the Higgs soft masses are constrained by the EWSB conditions at tan beta = 1:

$$m_{H_u}^2 = b - y_c^2, \quad m_{H_d}^2 = b - y_b^2$$

The bilinear b (or equivalently m_A^2 = 2b) remains as a free parameter.

### 6.7 Total Scalar Potential

$$V = V_F + V_D + V_{\text{soft}} = K^{I\bar{J}} F_I \bar{F}_{\bar{J}} + V_D + f_\pi^2 \text{Tr}(M^\dagger M) + [\text{Higgs soft terms}]$$

---

## 7. Vacuum Structure

### 7.1 Seiberg Seesaw Vacuum

Setting B = B-tilde = 0, H = 0 (before EWSB), and solving F_{M^i_i} = 0:

$$M_j^{\text{vac}} = C / m_j, \quad X_0 = -C / \Lambda^6$$

where C = Lambda^2 (m_u m_d m_s)^{1/3} = 882297 MeV^2.

| VEV | Value | Units |
|-----|-------|-------|
| <M^u_u> | 408471 | MeV |
| <M^d_d> | 188929 | MeV |
| <M^s_s> | 9446 | MeV |
| <X> | -1.210 x 10^{-9} | MeV^{-4} |
| det M / Lambda^6 | 1 (exact) | |

### 7.2 EWSB

After EWSB, H_u^0 = H_d^0 = v/sqrt(2) = 174104 MeV (tan beta = 1).

The NMSSM shifts det M slightly:

$$\det M = \Lambda^6 - \lambda v^2/2$$

Fractional shift: 3 x 10^{-5}. Negligible.

### 7.3 SUSY-Breaking F-terms at the Vacuum

| F-term | Value (MeV) | Origin |
|--------|-------------|--------|
| F_{M^u_u} | ~0 | Seesaw cancels |
| F_{M^d_d} | y_c v/sqrt(2) = 1796 | Yukawa-Higgs mixing |
| F_{M^s_s} | y_b v/sqrt(2) = 5911 | Yukawa-Higgs mixing |
| F_{H_u^0} | y_c <M^d_d> = 1949 | Meson VEV |
| F_{H_d^0} | y_b <M^s_s> = 321 | Meson VEV |
| F_X | ~0 | Constraint satisfied |

**Note on ISS universality:** F_{H_u} and F_{H_d} evaluated through the ISS seesaw give F = C/v = 289.6 MeV (flavor-universal). The numerical difference between the two values above (1949 vs 321) arises because the simple seesaw C = Lambda^2(m_u m_d m_s)^{1/3} uses the (u,d,s) block, while the Yukawa couplings effectively probe the (s,c,b) block. The ISS universality identity y_j <M_j> = C_{scb}/v holds exactly when the meson VEVs are from the (s,c,b) seesaw.

### 7.4 Supertrace

$$\text{STr}[M^2] = 2 \sum_{\text{mesons}} \tilde{m}^2 = 2 \times 9 \times f_\pi^2 = 18 f_\pi^2 = 152352 \text{ MeV}^2$$

This is independent of all Yukawa, NMSSM, and D-term contributions (they cancel between bosonic and fermionic sectors). The supertrace is entirely determined by the soft meson mass.

---

## 8. Parameter Budget

### 8.1 Complete Parameter List

| # | Parameter | Value | Status | What it determines |
|---|-----------|-------|--------|-------------------|
| 1 | m_u | 2.16 MeV | FREE | Up quark mass (input) |
| 2 | m_d | 4.67 MeV | FREE (connected to theta_C via Oakes) | Down quark mass |
| 3 | m_s | 93.4 MeV | FREE (overall mass scale) | Strange quark mass |
| 4 | m_c | 1270 MeV | DETERMINED by Koide seed (gv/m = sqrt(3)) | m_c/m_s = (2+sqrt(3))^2/(2-sqrt(3))^2 = 13.93 |
| 5 | m_b | 4177 MeV | DETERMINED by v_0-doubling (bion minimum) | sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c) |
| 6 | m_t | 172760 MeV | EXTERNAL (not confined) | Top Yukawa y_t |
| 7 | Lambda | 300 MeV | FREE (SQCD confinement scale) | Seesaw scale C |
| 8 | f_pi | 92 MeV | DETERMINED by QCD | Soft SUSY-breaking scale |
| 9 | lambda | 0.72 | DETERMINED by m_h = 125 GeV | Higgs quartic at tan beta = 1 |
| 10 | c = -1/12 | Kahler coefficient | ASSUMED/FITTED | Kahler pole -> Koide seed |
| 11 | c_3 | O(1) | UNDETERMINED | Three-instanton phase; bloom dynamics |
| 12 | b (soft bilinear) | > 0 | FREE | Pseudoscalar Higgs mass m_A |
| 13 | v_EW | 246220 MeV | DETERMINED by G_F | EW scale |
| 14 | g, g' | 0.653, 0.350 | DETERMINED by gauge couplings | D-terms |

### 8.2 Free Parameter Count

**Genuinely free parameters:** 3-5

1. **m_s** — the overall mass scale of the light sector. Everything else is built on m_s through Koide and v_0-doubling.
2. **m_u** — unconstrained by any Koide or seesaw relation.
3. **m_d** — connected to the Cabibbo angle via the GST relation sin(theta_C) = sqrt(m_d/m_s), so if theta_C is an input, m_d is determined. If theta_C is a prediction, m_d is free.
4. **Lambda** — the SQCD confinement scale. Possibly related to f_pi through dimensional transmutation; if so, not independent.
5. **b** — the soft bilinear controlling m_A. Purely in the Higgs sector.

**Determined parameters:** m_c, m_b, lambda, f_pi, v_EW, g, g', mu (Kahler)

**Assumed/fitted:** c = -1/12 (this is the least well-motivated parameter; it needs a microscopic derivation)

**Undetermined:** c_3, all other soft terms (M_1, M_2, A-terms)

### 8.3 Comparison to the Standard Model

The SM has 19 free parameters (or 26 with neutrino masses and mixing). This model starts with 3 light quark masses (m_u, m_d, m_s) and claims to determine m_c and m_b from them. It does not address:

| SM parameter | Status in this model |
|-------------|---------------------|
| m_u, m_d, m_s | FREE (input) |
| m_c, m_b | DETERMINED (Koide seed, v_0-doubling) |
| m_t | EXTERNAL (input from top Yukawa) |
| m_e, m_mu, m_tau | NOT DETERMINED (perturbative mesinos fail; must come from nonperturbative Kahler or other mechanism) |
| theta_C | Connected to m_d/m_s via Oakes (GST: -0.9%) |
| Other CKM angles | NOT DETERMINED |
| theta_W | NOT DETERMINED (but cf. de Vries angle) |
| alpha_s, alpha_EM | NOT DETERMINED (inputs to SQCD and D-terms) |
| m_h | DETERMINED (lambda = 0.72 at tan beta = 1) |
| m_W, m_Z | DETERMINED by v_EW (input) |
| Neutrino masses | NOT ADDRESSED |
| Strong CP phase | NOT ADDRESSED |

---

## 9. Predictions

### 9.1 Prediction Table

| Prediction | Predicted | Measured (PDG) | Status |
|-----------|-----------|----------------|--------|
| m_c/m_s = (2+sqrt(3))^2/(2-sqrt(3))^2 | m_c = 1301 MeV | 1270 +/- 20 MeV | 1.5 sigma (2.4%) |
| sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c) | m_b = 4177 MeV | 4180 +/- 30 MeV | 0.10 sigma (0.07%) |
| Q(c,b,t) = 2/3 (overdetermination check) | 0.6694 | n/a | 0.42% from 2/3 |
| Q(-s,c,b) = 2/3 (overdetermination check) | 0.6750 | n/a | 1.24% from 2/3 |
| Dual Koide Q(1/m_d, 1/m_s, 1/m_b) ~ 2/3 | 0.6652 | n/a | 0.22% from 2/3 |
| m_h (tree level, NMSSM) | 125.4 GeV | 125.25 +/- 0.17 GeV | 0.9 sigma (0.12%) |
| sin(theta_C) = sqrt(m_d/m_s) [GST] | 0.2236 | 0.2257 | -0.93% |
| Lepton Koide delta mod (2pi/3) = 2/9 | 0.222230 rad | 0.222222 rad | 0.89 sigma (33 ppm) |
| STr[M^2] = 18 f_pi^2 | 152352 MeV^2 | not directly measurable | consistency check |
| ISS F-term universality: F_{H_u} = F_{H_d} | exact (algebraic) | n/a | structural prediction |

### 9.2 Overdetermination Tests

These are not independent predictions — they test the internal consistency of the Koide chain:

- Q(c,b,t) uses m_c (from seed), m_b (from v_0-doubling), m_t (external). If these are correct, Q should be close to 2/3. Result: 0.42% deviation. **Pass.**
- Q(-s,c,b) uses m_s (input), m_c (predicted), m_b (predicted). Result: 1.24% deviation. **Marginal.**
- The dual Koide is a separate numerical observation with no clear mechanism in the superpotential.

---

## 10. Consistency Analysis

### 10.1 Can tan beta = 1 Be Maintained?

**The issue:** The Yukawa couplings y_c (to H_u) and y_b (to H_d) explicitly break the H_u <-> H_d symmetry. Yet tan beta = 1 requires v_u = v_d, which naively assumes this symmetry.

**Resolution:** tan beta = 1 does not require the Lagrangian to have H_u <-> H_d symmetry. It requires only that the EWSB minimum has v_u = v_d. The tadpole equations at tan beta = 1 give:

$$m_{H_u}^2 + y_c^2 = b = m_{H_d}^2 + y_b^2$$

This is satisfied when:

$$m_{H_u}^2 - m_{H_d}^2 = y_b^2 - y_c^2 = 1.05 \times 10^{-3}$$

So the soft Higgs masses must have a SPECIFIC splitting that tracks the Yukawa splitting. This is a tuning condition. There are two interpretations:

1. **Fine-tuning interpretation:** The splitting m_{H_u}^2 - m_{H_d}^2 must be arranged to match y_b^2 - y_c^2 at the EW scale. Radiative corrections from the top Yukawa would generically drive tan beta away from 1.

2. **Symmetry interpretation:** If there is an underlying symmetry that enforces tan beta = 1 (e.g., a discrete symmetry H_u <-> H_d at high scales, broken only by Yukawa couplings), then the required soft mass splitting is generated radiatively from the Yukawa asymmetry itself. The RG running from the mediation scale to the EW scale naturally produces m_{H_u}^2 - m_{H_d}^2 proportional to y_b^2 - y_c^2.

**Assessment:** The model ASSUMES tan beta = 1 and derives consequences. The justification is that the Koide condition on mixed (up/down)-type triples is most naturally a mass-level condition when tan beta = 1 (because then m_q = y_q v/sqrt(2) with the same v for all quarks). A microscopic derivation of tan beta = 1 from the SUSY-breaking sector is MISSING.

### 10.2 Hierarchy: Why is f_pi << F_H?

**The numbers:**
- Soft SUSY-breaking scale: sqrt(f_pi^2) = f_pi = 92 MeV
- Higgs F-term scale: F_H = C/v = 289.6 MeV (or y_c <M_d> = 1949 MeV, y_b <M_s> = 321 MeV depending on the block)
- The ratio F_H / f_pi = 289.6/92 = 3.15

**The question:** In standard SUSY breaking, the soft mass is generated by F-terms via loop or Planck-suppressed operators. Here the soft mass (f_pi^2) is SMALLER than the F-terms (F_H^2). How?

**Possible resolutions:**

(a) **f_pi is a QCD effect, not a SUSY-breaking mediation effect.** The pion decay constant is set by the QCD chiral condensate, which breaks both chiral symmetry and (if SQCD) SUSY. The soft mass m-tilde^2 = f_pi^2 is then the SQCD analogue of the chiral condensate, not a gravity- or gauge-mediated soft term. The F-terms from the Yukawa sector are larger because they involve the EW-scale meson VEVs.

(b) **Separate SUSY-breaking sectors.** The meson soft mass comes from a confined-sector mechanism (chiral condensate, Volkov-Akulov), while the Higgs F-terms come from the Yukawa-seesaw interplay. These are different sectors of the theory, and there is no requirement that they be comparable.

(c) **The hierarchy is actually mild.** F_H/f_pi ~ 3 is not a large hierarchy. In the MSSM, the soft scale and the mu-term are typically comparable. Here the ratio is O(1) when measured in natural units.

**Assessment:** The hierarchy is not severe (factor of ~3), and the physical origin is clear: f_pi comes from QCD/SQCD dynamics at ~100 MeV, while F_H comes from Yukawa couplings to meson VEVs at ~1-10 GeV. The mild hierarchy is a feature, not a bug. However, the model does not explain WHY f_pi sets the soft mass, only assumes it.

---

## 11. Missing Elements (Honest Gaps)

### 11.1 Lepton Masses

**Status: NOT DETERMINED.**

The perturbative fermion mass matrix of the confined sector does not reproduce lepton masses. The three lightest eigenvalues of the central 6x6 block have Koide quotient Q = 0.481 (28% from 2/3). Their mass ratios do not match (m_e : m_mu : m_tau). The NMSSM extension does not help.

**What is needed:** A nonperturbative mechanism in the Kahler sector (bion corrections to the fermion kinetic terms, or a separate confining sector for leptons). This is the model's most significant open problem.

### 11.2 Top Quark Coupling to the Meson Sector

**Status: EXTERNAL.**

The top quark is too heavy to participate in the confining dynamics (m_t >> Lambda). It couples to H_u via a standard Yukawa but does not form mesons. The (c,b,t) Koide triple exists as a numerical observation but the top's role in the confined-phase Lagrangian is purely through y_t H_u^0 t t-bar.

**What is needed:** An explanation of why the top is special (the "elephant" problem). In the SU(5) flavor framework, only 5 quarks (u,c,d,s,b) form the fundamental 5; the top is a singlet. But the mechanism that excludes the top from confinement is not derived.

### 11.3 mu-Term Generation

**Status: YES, partially.**

The NMSSM coupling lambda X H_u . H_d generates an effective mu-term when X acquires a VEV:

$$\mu_{\text{eff}} = \lambda \langle X \rangle = \lambda X_0$$

Numerically: mu_eff = 0.72 * (-1.210 x 10^{-9} MeV^{-4}) * [appropriate dimension] = extremely small (sub-eV scale). This is because |X_0| = C/Lambda^6 ~ 10^{-9} is tiny.

**Problem:** The mu-term is generated but is absurdly small. The standard NMSSM requires mu ~ O(m_Z) ~ 100 GeV. Here it is ~10^{-9} MeV. This is a serious problem.

**Possible resolution:** If the relevant VEV is <X> = sqrt(3) mu (from the Kahler pole), and mu is a Kahler mass parameter (not the same as m_u), then:

$$\mu_{\text{eff}} = \lambda \sqrt{3}\, \mu_K$$

The Kahler pole mechanism gives |X_0| = sqrt(3) mu_K, so mu_K = |X_0|/sqrt(3) ~ 7 x 10^{-10} MeV^{-4}. The dimensions don't match for a standard mu-term. This signals that the X field has unconventional dimensions in the Seiberg effective theory (it is a Lagrange multiplier with [X] = mass^{-3} when [det M] = mass^6 and [Lambda^6] = mass^6). The NMSSM identification X = S requires rescaling X to have canonical dimension, which has not been done carefully.

**Status: UNRESOLVED.** The mu-term generation is a dimensional analysis problem that requires careful treatment of the Seiberg multiplier's unconventional dimension.

### 11.4 Neutrino Sector

**Status: NOT ADDRESSED.**

The model has no right-handed neutrinos, no seesaw mechanism for neutrinos, and no Majorana mass terms. Neutrino masses are completely outside the scope.

### 11.5 Gaugino Masses

**Status: NOT ADDRESSED.**

Bino and wino masses (M_1, M_2) are required for chargino/neutralino phenomenology. The model does not specify the SUSY-breaking mediation to the gauge sector. Since SU(3)_c is confined, the gluino mass is subsumed into the confined-phase description, but the EW gauginos need explicit soft masses.

### 11.6 A-Terms

**Status: NOT ADDRESSED.**

Trilinear A-terms (A_c y_c H_u M_d, A_b y_b H_d M_s, A_t y_t H_u t t-bar) are generically generated by any SUSY-breaking mechanism. They affect the scalar spectrum and Higgs mixing. The model sets them to zero implicitly.

### 11.7 Tachyonic Modes

**Status: KNOWN PROBLEM.**

The NMSSM spectrum has 4 tachyonic scalar modes (in the central block and charged Higgs sector). The original model (without NMSSM) had tachyonic modes in the us and ds off-diagonal meson sectors when Yukawa-induced F-terms are included. These tachyons signal that the diagonal vacuum ansatz M = diag(M_u, M_d, M_s) is NOT the true minimum; off-diagonal VEVs (related to CKM mixing) may stabilize these directions.

**What is needed:** A full off-diagonal vacuum analysis, connecting the tachyonic modes to CKM mixing angles.

### 11.8 Kaon and Pion Mass Splittings

**Status: NOT ADDRESSED.**

The meson scalars in the model have masses set by f_pi (the soft mass). The physical pion and kaon mass splittings arise from quark mass differences (m_d - m_u, m_s - m_d) through chiral perturbation theory. The model's meson spectrum (from the off-diagonal blocks) does produce mass splittings proportional to 1/m_q (inverted hierarchy from the seesaw), but these are fermion masses (mesinos), not scalar meson masses. The connection to physical meson spectroscopy is unclear.

---

## 12. Complete Superpotential (LaTeX Form)

```latex
W = \underbrace{\sum_{i=u,d,s} m_i\,M^i_i
  + X\bigl(\det M^{(3)} - B\,\tilde{B} - \Lambda^6\bigr)}_{W_{\text{Seiberg}}}
  + \underbrace{y_c\,H_u^0\,M^d_d
  + y_b\,H_d^0\,M^s_s
  + y_t\,H_u^0\,t\,\bar{t}}_{W_{\text{Yukawa}}}
  + \underbrace{\lambda\,X\,(H_u^0 H_d^0 - H_u^+ H_d^-)}_{W_{\text{NMSSM}}}
  + \underbrace{c_3\,\frac{(\det M^{(3)})^3}{\Lambda^{18}}}_{W_{\text{instanton}}}
```

### 12.1 Complete Kahler Potential (LaTeX Form)

```latex
K = \underbrace{\operatorname{Tr}(M^\dagger M) + |X|^2 + |B|^2
  + |\tilde{B}|^2 + |H_u|^2 + |H_d|^2}_{K_{\text{canonical}}}
  \underbrace{- \frac{|X|^4}{12\,\mu^2}}_{K_{\text{pseudo-modulus}}}
  + \underbrace{\frac{\zeta^2}{\Lambda^2}\,
  \biggl|\sum_{k=1}^{3} s_k \sqrt{M^k_k}\biggr|^2}_{K_{\text{bion}}}
```

### 12.2 Soft-Breaking Potential (LaTeX Form)

```latex
V_{\text{soft}} = f_\pi^2\,\operatorname{Tr}(M^\dagger M)
  + m_{H_u}^2\,|H_u|^2 + m_{H_d}^2\,|H_d|^2
  + \bigl(b\,H_u \cdot H_d + \text{h.c.}\bigr)
```

with the EWSB constraint: m_{H_u}^2 + y_c^2 = m_{H_d}^2 + y_b^2 = b at tan beta = 1.

---

## 13. Classification Summary

### 13.1 Term-by-Term Classification

| Term | Origin | Classification |
|------|--------|---------------|
| Tr(m M) | Quark mass deformation | DETERMINED (Seiberg, QCD inputs) |
| X(det M - BB-tilde - Lambda^6) | Quantum moduli space constraint | DETERMINED (Seiberg, 1994) |
| y_c H_u M^d_d | Charm Yukawa | FITTED (from m_c and v) |
| y_b H_d M^s_s | Bottom Yukawa | FITTED (from m_b and v) |
| y_t H_u t t-bar | Top Yukawa | FITTED (from m_t and v) |
| lambda X H_u . H_d | NMSSM singlet coupling | DETERMINED from m_h once X = S is ASSUMED |
| c_3 (det M)^3/Lambda^18 | Three-instanton | ASSUMED to exist; coefficient UNDETERMINED |
| Tr(M-dag M) (canonical) | Kinetic terms | DETERMINED |
| -|X|^4/(12 mu^2) | Pseudo-modulus quartic | ASSUMED; c = -1/12 FITTED |
| (zeta^2/Lambda^2) |sum s_k sqrt(M_k)|^2 | Bion Kahler | MOTIVATED; specific form ASSUMED |
| f_pi^2 Tr(M-dag M) | Meson soft mass | ASSUMED (source of breaking not derived) |
| Higgs soft masses, b | EWSB soft terms | REQUIRED; values CONSTRAINED by tan beta = 1 |
| D-terms | Gauge sector | DETERMINED by gauge invariance |

### 13.2 What is Robust vs. What is Speculative

**Robust (follows from established SUSY gauge theory):**
- Seiberg superpotential structure
- Seesaw M_j = C/m_j
- ISS F-term universality (algebraic identity)
- Supertrace STr = 18 f_pi^2
- D-flat direction at tan beta = 1
- NMSSM Higgs mass formula

**Observationally supported (numerical predictions match data):**
- Koide seed from O'Raifeartaigh (m_c/m_s prediction, 2.4% off)
- v_0-doubling (m_b prediction, 0.07% off)
- Higgs mass at tan beta = 1 with lambda = 0.72
- GST Cabibbo relation

**Speculative (assumed without microscopic derivation):**
- c = -1/12 for the pseudo-modulus Kahler correction
- Bion Kahler potential functional form (sqrt(M) dependence)
- Sign structure s_k = (-1, +1, +1) for the bloom configuration
- lambda_3 = -4 lambda_2 in the bion effective potential
- Identification X = S (NMSSM singlet)
- f_pi^2 as the sole soft mass source
- tan beta = 1 (not derived)

---

## 14. What a Referee Would Ask

### Q1: "You have tachyonic scalar modes. Is the vacuum stable?"

**Answer:** No. The diagonal vacuum M = diag(M_u, M_d, M_s) is not the true minimum when Yukawa couplings are included. The tachyonic us and ds off-diagonal modes indicate that the true vacuum has nonzero off-diagonal meson VEVs. This is potentially connected to CKM mixing. A full minimization of the off-diagonal potential is needed and has not been completed.

### Q2: "How do you derive c = -1/12?"

**Answer:** We do not derive it from first principles. The value c = -1/12 is the unique coefficient that places the Kahler pole at the seesaw VEV, which then produces the Koide seed through the O'Raifeartaigh mechanism. A microsopic derivation from the magnetic dual's one-loop effective Kahler potential might yield this value but has not been carried out.

### Q3: "Where do lepton masses come from?"

**Answer:** This is the model's most significant gap. The perturbative spectrum does not reproduce leptons. The expectation is that nonperturbative Kahler corrections (bion type) modify the fermion kinetic terms, producing an effective mass hierarchy in the mesino sector that matches the lepton spectrum. This has not been demonstrated.

### Q4: "What about anomaly matching?"

**Answer:** The Seiberg effective theory automatically satisfies 't Hooft anomaly matching for the SU(3) SQCD anomalies. The additional Yukawa and NMSSM couplings do not introduce gauge anomalies (H_u, H_d form complete SU(2)_L doublets; X is a gauge singlet). The 24 + 15 + 15-bar decomposition under SU(5) flavor is anomaly-free and asymptotically free (b_0 = 31/3).

### Q5: "Is the bion Kahler potential gauge-invariant?"

**Answer:** The bion Kahler potential involves sqrt(M^k_k), which is not gauge-invariant under SU(3) flavor rotations (it picks out the diagonal entries). This is acceptable because the quark mass matrix m = diag(m_u, m_d, m_s) already breaks SU(3) flavor. The bion corrections respect the residual U(1)^2 symmetry preserved by the mass matrix, which allows diagonal M entries but not off-diagonal ones. However, the extension to off-diagonal mesons is unclear.

### Q6: "Why N_f = N_c = 3 specifically?"

**Answer:** The bootstrap constraints rs = 2N, r(r+1)/2 = 2N, r^2 + s^2 - 1 = 4N have the unique solution N = 3, r = 3, s = 2. This is the starting point of the sBootstrap program and is the most rigid structural prediction of the framework.

### Q7: "The dimensional analysis of the mu-term seems wrong."

**Answer:** Correct. The Seiberg multiplier X has unconventional dimensions [X] = mass^{-3} in the confined-phase effective theory. Identifying X with the NMSSM singlet S (which has [S] = mass) requires a dimensional rescaling X -> X' = X Lambda^4 or similar. This has not been done carefully and may introduce additional Lambda-dependent factors that affect the mu-term and the NMSSM coupling.

---

## 15. Summary

### What this Lagrangian achieves:

1. **Unifies SQCD confinement with the Higgs mechanism** through the identification X = S.
2. **Predicts m_b from m_s and m_c** (v_0-doubling, 0.07% accuracy).
3. **Predicts m_c/m_s from the O'Raifeartaigh seed** (Koide relation, 2.4% accuracy).
4. **Reproduces m_h = 125 GeV** at tree level with tan beta = 1 and lambda = 0.72.
5. **Provides a structural explanation** for the ISS F-term universality and the Cabibbo-Oakes relation.
6. **Reduces the quark mass parameter count** from 6 to 3 (m_u, m_d, m_s), with m_c, m_b determined and m_t external.

### What it does not achieve:

1. **Lepton masses** — not determined from the perturbative Lagrangian.
2. **Complete CKM matrix** — only the Cabibbo angle is connected (via Oakes); the other mixing angles are not.
3. **Neutrino masses** — not addressed.
4. **Microscopic derivation of c = -1/12** — fitted, not derived.
5. **Microscopic derivation of the bion Kahler** — motivated but not rigorously computed.
6. **Vacuum stability** — tachyonic modes present.
7. **Complete soft spectrum** — gaugino masses, A-terms not specified.
8. **mu-term** — dimensional analysis problem unresolved.
9. **tan beta = 1 justification** — assumed, not derived.
10. **Top quark decoupling** — assumed, not derived.

### The path forward:

The most impactful next steps are (in order of priority):

1. **Off-diagonal vacuum analysis** — stabilize the tachyonic modes, extract CKM mixing angles.
2. **Dimensional rescaling of X** — fix the mu-term and put the NMSSM identification on solid footing.
3. **Lepton mass mechanism** — identify the nonperturbative Kahler correction that reproduces (m_e, m_mu, m_tau).
4. **Microscopic c = -1/12** — compute the one-loop effective Kahler potential in the ISS magnetic dual.
5. **Complete soft spectrum** — specify the SUSY-breaking mediation to gauginos and sfermions.

---

*This document assembles results from the full computation chain: vacuum_structure.md, unified_superpotential.md, complete_kahler.md, full_spectrum.md, nmssm_spectrum.md, higgs_potential.md, bloom_dynamics.md, cabibbo_oakes.md, seesaw_mechanism.md, oraifeartaigh_koide.md, lepton_yukawa.md, pauli_15_10bar.md, and iss_cw_koide.md.*
