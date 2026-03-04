# Anomalous Dimension and Wave Function Renormalization for Mesons in N_f = N_c SQCD

## Setup

N = 1 SQCD with SU(3) gauge group and N_f = N_c = 3 flavors. Confinement scale Lambda = 300 MeV. The composite mesons M^i_j = Q^i Q-bar_j appear as elementary fields in the confined description via the quantum-modified moduli space constraint det M - B B-bar = Lambda^6.

---

## Part (a): Anomalous Dimension at Confinement

### The would-be fixed point is unphysical

The NSVZ exact beta function gives, at a would-be fixed point (beta = 0):

  gamma* = 1 - 3 N_c / N_f = 1 - 3(3)/3 = -2

This is negative and large, confirming that N_f = N_c = 3 is far below the conformal window (N_f = 4.5 to 9 for SU(3)). There is no perturbative IR fixed point.

### R-charge and scaling dimension

If one naively applies the a-maximization / unitarity bound analysis as though the theory were conformal:

  R(Q) = 1 - N_c/N_f = 1 - 1 = 0
  R(M = Q Q-bar) = 2 R(Q) = 0
  Delta(M) = (3/2) R(M) = 0

This violates the unitarity bound Delta >= 1 for scalar operators, signaling that the theory confines rather than flowing to a fixed point. The R-charge assignment R(Q) = 0 corresponds to M being a Lagrange multiplier (the field X in the confined description), not a propagating degree of freedom.

### Confining description: M as an elementary field

In the Seiberg confined description, M is a fundamental chiral superfield with canonical dimension 1. But in the UV, M corresponds to the composite Q Q-bar with UV dimension 2. The "effective" nonperturbative anomalous dimension bridging UV to IR is:

  gamma_eff = Delta_IR - Delta_UV_free = 1 - 2 = -1

This is NOT a perturbative anomalous dimension. It represents the full nonperturbative transmutation from a composite of UV dimension 2 to an elementary field of IR dimension 1.

### Wave function renormalization Z_M

The wave function renormalization relates UV and IR normalizations:

  M_IR = Z_M^{1/2} * M_UV

With the parametrization Z_M = (Lambda / M_UV)^{2 gamma}, where M_UV is the UV cutoff and gamma is the anomalous dimension:

For gamma = 0.5 and M_UV = M_Planck = 2.4 x 10^18 GeV:

  Lambda / M_Planck = 1.25 x 10^{-19}
  Z_M = (1.25 x 10^{-19})^{1.0} = 1.25 x 10^{-19}
  1/Z_M = 8.0 x 10^{18}

Table of Z_M for various gamma:

| gamma | Z_M | 1/Z_M | 1/sqrt(Z_M) |
|-------|-----|-------|-------------|
| 0.1 | 1.66 x 10^{-4} | 6.03 x 10^3 | 77.7 |
| 0.2 | 2.75 x 10^{-8} | 3.64 x 10^7 | 6.03 x 10^3 |
| 0.3 | 4.55 x 10^{-12} | 2.20 x 10^{11} | 4.69 x 10^5 |
| 0.5 | 1.25 x 10^{-19} | 8.00 x 10^{18} | 2.83 x 10^9 |
| 0.7 | 3.43 x 10^{-27} | 2.91 x 10^{26} | 1.71 x 10^{13} |
| 1.0 | 1.56 x 10^{-38} | 6.40 x 10^{37} | 8.00 x 10^{18} |

The enormous range of Z_M (spanning 38 orders of magnitude for gamma from 0 to 1) shows that even a modest anomalous dimension produces a huge hierarchy between holomorphic and physical masses.

### Interpretation

The question "what is the anomalous dimension of Q Q-bar at the confinement scale for N_f = N_c?" is not well-posed in the standard perturbative sense. The theory is not near a fixed point, so there is no single number gamma that characterizes the RG flow. What exists is:

1. A UV region (mu >> Lambda) where M = Q Q-bar / mu has dimension 1 and gamma is small and perturbatively calculable.
2. A confinement transition at mu ~ Lambda where the description changes to confined mesons.
3. An IR region (mu << Lambda) where M is an elementary field with canonical dimension 1.

The full wave function renormalization Z_M integrates all of these effects. It is not captured by any single value of gamma.

---

## Part (b): Mesino Mass from Anomalous Dimension

### Tree-level mesino masses

At the Seiberg seesaw vacuum with (u,d,s) flavors:

  C = Lambda^2 (m_u m_d m_s)^{1/3} = 882,297 MeV^2
  X_0 = -C / Lambda^6 = -1.210 x 10^{-9} MeV^{-4}

The off-diagonal mesino masses (Dirac pairs):

| Mesino | Mass formula | Tree-level mass |
|--------|-------------|-----------------|
| psi_{ud} | \|X_0\| M_s | 11.4 eV |
| psi_{us} | \|X_0\| M_d | 228.7 eV |
| psi_{ds} | \|X_0\| M_u | 494.4 eV |

### Physical mass with wave function renormalization

The physical mesino mass involves two factors of Z_M^{-1/2} (one per leg of the bilinear):

  m_phys = m_tree / Z_M

With Z_M = (Lambda / M_UV)^{2 gamma}:

| gamma | m(psi_{ud}) phys | m(psi_{us}) phys |
|-------|-----------------|-----------------|
| 0.1 | 0.069 MeV | 1.38 MeV |
| 0.12 | 0.51 MeV | 10.2 MeV |
| 0.15 | 5.3 MeV | 106 MeV |
| 0.2 | 416 MeV | 8,326 MeV |
| 0.23 | 4,700 MeV | 94,000 MeV |
| 0.26 | 51,000 MeV | 1.0 x 10^6 MeV |
| 0.3 | 2.5 x 10^6 MeV | 5.0 x 10^7 MeV |

### Phenomenological viability threshold

For the charged mesinos (psi_{ud} and psi_{us}) to be heavier than the LEP limit (approximately 100 GeV):

  gamma > 0.263 (for psi_{ud})
  gamma > 0.229 (for psi_{us})

So gamma > 0.26 suffices to make ALL charged mesinos heavier than 100 GeV.

### Lepton identification attempt

If mesinos are to be identified with leptons (a speculative possibility examined in mesino_phenomenology.md):

  gamma = 0.123 would give m(psi_{ud}) = m_e = 0.511 MeV
  gamma = 0.150 would give m(psi_{us}) = m_mu = 105.66 MeV

However, a UNIVERSAL gamma cannot simultaneously reproduce both masses. The tree-level mass ratio is:

  m_tree(psi_{us}) / m_tree(psi_{ud}) = M_d / M_s = m_s / m_d = 20.0

But the target ratio is:

  m_mu / m_e = 206.8

A universal Z_M preserves mass ratios, so the predicted ratio remains 20.0 regardless of gamma. The discrepancy factor is 206.8 / 20.0 = 10.3.

To fix this, one would need FLAVOR-DEPENDENT anomalous dimensions, meaning Z_M depends on the specific meson component. This is possible in principle (the anomalous dimension of Q^i Q-bar_j depends on the quark masses m_i, m_j through their couplings to the strong dynamics), but there is no controlled calculation of this effect.

---

## Part (c): Literature on Anomalous Dimensions for N_f = N_c SQCD

### 1. Lattice SQCD results

Direct lattice simulations of N = 1 SQCD with N_f = N_c = 3 remain computationally challenging. The main difficulty is preserving supersymmetry on the lattice; fine-tuning of scalar and gaugino masses is required to reach the SUSY continuum limit.

Existing lattice work relevant to SQCD composite operators:

- **Costa et al. (2019)** [arXiv:1812.06770, Phys. Rev. D 99, 074512]: Computed renormalization and mixing of composite operators (quark and squark bilinears) in lattice SQCD with SU(N_c) gauge group. Found significant mixing between quark bilinears and gluino bilinears, complicating extraction of anomalous dimensions. Calculations performed with massive matter fields to regulate IR singularities in squark bilinear renormalization.

- **Costa et al. (2018)** [arXiv:1812.00066]: One-loop lattice study of composite bilinear operators in SQCD, finding perturbative renormalization factors and mixing coefficients.

- **Bergner et al. (2023)** [arXiv:2302.10514]: Studied the spectrum of QCD with one flavor (N_f = 1), finding connections to SUSY dynamics. While not N_f = N_c = 3 directly, this provides the closest lattice data on how anomalous dimensions behave in confined SQCD-like theories.

No lattice computation has yet extracted the nonperturbative anomalous dimension of Q Q-bar at confinement for N_f = N_c = 3.

### 2. Analytic estimates

For the confining phase below the conformal window, the anomalous dimension is not accessible by perturbation theory at the confinement scale. Known analytic results:

- **NSVZ formula**: Gives gamma* = -2 at the would-be fixed point for N_f/N_c = 1, which is unphysical (signals confinement, not conformality).

- **Seiberg's analysis**: For N_f = N_c, the R-charge of Q is R(Q) = 0, giving Delta(M) = 0. This violates unitarity bounds and indicates M becomes a Lagrange multiplier. The physical interpretation is that the description changes from composite operators to elementary confined fields.

- **Interpolation from conformal window**: Banks-Zaks analysis gives gamma* approximately 0 at N_f = 3 N_c (upper edge) and gamma* = 1 at N_f = 3 N_c / 2 (lower edge, where the unitarity bound is saturated). Extrapolating below the window is unreliable.

- **Instanton contributions**: For N_f = N_c, the one-instanton amplitude generates the det M term in the superpotential (Affleck-Dine-Seiberg). Multi-instanton contributions to the Kahler potential are not protected by holomorphy and are incalculable in general.

### 3. a-theorem and conformal bootstrap bounds

- **a-theorem (Komargodski-Schwimmer 2011)**: Provides a_UV > a_IR, constraining the number of degrees of freedom. For N_f = N_c = 3 SU(3) SQCD, a_UV (quarks + gauginos) > a_IR (mesons + baryons + X). This is satisfied by the Seiberg confined spectrum but does not directly constrain gamma.

- **Conformal bootstrap** [arXiv:1603.01995]: Upper bound on mass anomalous dimension for SU(N) gauge theories with many flavors. For SU(12) at the conformal window edge, gamma_m* < 1.29 assuming no SU(N)-breaking relevant operators. These bounds apply to the conformal window, not to the confined phase below it.

- **Nakayama-Ohtsuki (2021)** [JHEP 10 (2021) 114]: Searched for gauge theories with conformal bootstrap, finding kinks and bounds on scalar operator dimensions in theories with SU(N) global symmetry. At large N, bounds approach free fermion bilinear dimensions Delta = D - 1 = 3.

None of these methods directly constrain gamma for N_f = N_c in the confined phase.

### 4. What is known

The honest answer: **very little is known rigorously about the anomalous dimension of Q Q-bar at the confinement scale for N_f = N_c SQCD.** The key facts are:

- The holomorphic information (superpotential, R-charges, moduli space) is exactly known.
- The Kahler potential is NOT protected and is essentially unknown beyond leading order.
- The anomalous dimension is a Kahler potential quantity, so it is not determined by holomorphy.
- Lattice SQCD has not yet reached the precision needed to extract this quantity.
- The a-theorem and conformal bootstrap provide inequalities, not equalities, and apply primarily to the conformal window.

The only constraints are:
1. The meson kinetic term must be positive (Z_M > 0).
2. The a-theorem is satisfied (checked for the Seiberg spectrum).
3. For N_f slightly above 3N_c/2, gamma* approaches 1 from below (unitarity bound).
4. Below the conformal window, gamma ceases to have a fixed-point interpretation and becomes a running quantity.

---

## Part (d): Alternative Mechanisms for Mesino Mass

### (1) Soft SUSY-breaking from spurion F_S = f_pi^2

If SUSY is broken in the QCD sector with F-term F_S = f_pi^2 = (92 MeV)^2 = 8464 MeV^2:

**Scalar soft mass**:
  m^2_soft ~ |F_S|^2 / Lambda^2 = 796 MeV^2
  m_soft ~ 28 MeV

**Mesino soft mass** (from A-term or direct superpotential coupling):
  m_mesino ~ F_S / Lambda = 28.2 MeV

This is 5 orders of magnitude above the tree-level holomorphic masses (eV scale), but still well below the electron mass of 0.511 MeV for the lightest charged mesino. It does not resolve the phenomenological problem.

### (2) Gauge-mediated contributions through the confined sector

**Electromagnetic contribution** (mesinos carry EM charge +/-1 or 0):
  delta_m ~ (alpha_em / 4pi) * Lambda = 0.17 MeV = 170 keV

**QCD contribution** (if mesons carry color, but in the confined phase they do not):
  delta_m ~ (alpha_s / 4pi) * Lambda = 7.2 MeV

**Gauge mediation with F_S as messenger**:
  m ~ (alpha_em / 4pi) * F_S / Lambda = 0.016 MeV = 16 keV

Gauge-mediated contributions give masses in the keV-to-MeV range. Insufficient to reach the electroweak scale.

### (3) Anomaly-mediated contributions (AMSB)

The AMSB contribution to the mesino mass is:

  m_AMSB = gamma_M * m_{3/2}

where gamma_M is the meson anomalous dimension and m_{3/2} is the gravitino mass.

For m_{3/2} = 30 TeV (typical gravity-mediation scale):

| gamma_M | m_AMSB |
|---------|--------|
| 0.01 | 300 MeV |
| 0.1 | 3 GeV |
| 0.5 | 15 GeV |
| 1.0 | 30 GeV |
| 2.0 | 60 GeV |

The gaugino AMSB mass before confinement:
  m_gaugino = (b_0 alpha_s / 4pi) m_{3/2} = (6 * 0.3 / 4pi) * 30 TeV = 4.3 GeV

AMSB gives the largest contributions among the alternatives, reaching tens of GeV for O(1) anomalous dimensions. However, AMSB requires a specific UV completion (sequestered sectors, conformal compensator) and the gravitino mass is a free parameter.

For the scalar meson AMSB mass-squared:
  m^2 ~ (gamma_M * alpha_s / 8pi) * m_{3/2}^2

For gamma_M = 1, m_{3/2} = 30 TeV: m_scalar ~ 3.3 GeV. This is small compared to the TeV scale.

### (4) Higher-dimensional Kahler operators

The Kahler potential receives corrections from operators of the form:

  K = |M|^2 (1 + c_1 |M|^2 / Lambda^2 + c_2 |M|^4 / Lambda^4 + ...)

At the seesaw vacuum, the meson VEVs are:

  M_u = 408,471 MeV >> Lambda = 300 MeV
  M_d = 188,929 MeV >> Lambda
  M_s = 9,446 MeV >> Lambda

The expansion parameter |M|^2 / Lambda^2 ranges from 10^3 to 10^6, far exceeding unity. The perturbative Kahler expansion is DIVERGENT at the seesaw vacuum.

This has two implications:
1. The canonical Kahler potential K = |M|^2 is an extremely poor approximation at the seesaw vacuum.
2. The physical meson masses and mesino masses receive non-perturbative corrections of unknown sign and magnitude from the full Kahler potential.

The breakdown of the Kahler expansion is itself evidence that large corrections to the mesino masses are EXPECTED. The question is whether these corrections are calculable.

---

## Summary Table

| Mechanism | Mesino mass estimate | Comment |
|-----------|---------------------|---------|
| Tree-level (holomorphic) | 11 - 494 eV | Exact (superpotential protected) |
| Soft SUSY-breaking (f_pi) | ~ 28 MeV | Requires specific breaking pattern |
| Gauge-mediated (EM) | ~ 16 keV | Small (loop-suppressed) |
| Gauge-mediated (QCD) | ~ 7 MeV | Mesons are color singlets |
| AMSB (gamma = 1, m_{3/2} = 30 TeV) | ~ 30 GeV | Depends on gravitino mass |
| Higher-dim Kahler | Non-perturbative | Expansion breaks down at seesaw vacuum |
| Anom. dim. gamma = 0.23 (Planck cutoff) | ~ 100 GeV | Viable for collider bounds |
| Anom. dim. gamma = 0.12 - 0.15 | ~ 0.5 - 106 MeV | Lepton mass range |
| Anom. dim. gamma = 0.26 | > 100 GeV (all charged) | Minimum for LEP exclusion |

## Key Conclusions

1. **The anomalous dimension gamma of Q Q-bar at confinement for N_f = N_c = 3 is not known.** It is not determined by holomorphy, not accessible perturbatively, and has not been computed on the lattice. The NSVZ fixed-point value (gamma* = -2) is unphysical. The R-charge analysis gives Delta(M) = 0, signaling a change of description rather than a finite anomalous dimension.

2. **Even a modest gamma ~ 0.12 - 0.26 produces enormous wave function renormalization** when running from the Planck scale to Lambda. This is because log(M_Planck / Lambda) = 43.6 in natural units, so Z_M = exp(-2 gamma * 43.6) which is exponentially sensitive to gamma.

3. **A universal gamma cannot reproduce lepton masses.** The tree-level mass ratio m(psi_{us})/m(psi_{ud}) = m_s/m_d = 20.0 is preserved by universal Z_M, but the target m_mu/m_e = 206.8 requires a factor of 10.3 from flavor-dependent anomalous dimensions.

4. **The Kahler expansion breaks down at the seesaw vacuum** because M_i >> Lambda for the light flavors. This is both a problem (loss of calculational control) and a feature (large corrections are expected and could bridge the gap between eV and MeV/GeV scales).

5. **Among the alternative mechanisms, AMSB is the most promising** for generating large mesino masses, giving m ~ gamma * m_{3/2} which can reach the tens-of-GeV range. However, it requires a specific UV completion with a gravitino mass as a free parameter. The soft SUSY-breaking spurion gives masses ~ f_pi^2 / Lambda ~ 28 MeV, which is parametrically determined but numerically insufficient.

6. **The honest assessment**: getting mesino masses to the phenomenologically required range (> 100 GeV for charged states, or ~ 0.5 - 106 MeV if identified with leptons) requires either (a) a nonperturbative Kahler potential computation that no existing technique can provide, or (b) a specific SUSY-breaking mechanism (AMSB or gauge mediation) with tuned parameters. The problem is open.
