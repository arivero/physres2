# AMSB Soft Spectrum for Confined SU(3) SQCD with N_f = N_c = 3

## Setup

SU(3) SQCD with N_f = N_c = 3 flavors (u, d, s) in the confining phase. Below the confinement scale Lambda = 300 MeV, the effective theory has:

- 9 complex mesons M^i_j with Seiberg seesaw VEVs M_j = C/m_j
- Baryons B, B-tilde
- Lagrange multiplier X (no kinetic term)
- NMSSM Higgs sector (H_u, H_d, S)

Numerical inputs: m_u = 2.16 MeV, m_d = 4.67 MeV, m_s = 93.4 MeV, Lambda = 300 MeV, f_pi = 92 MeV.

Seiberg seesaw parameters:

| Quantity | Value |
|----------|-------|
| C = Lambda^2 (m_u m_d m_s)^{1/3} | 882,297 MeV^2 |
| M_u = C/m_u | 408,471 MeV |
| M_d = C/m_d | 188,929 MeV |
| M_s = C/m_s | 9,446 MeV |
| X_0 = -C/Lambda^6 | -1.210 x 10^{-9} MeV^{-4} |

---

## Part (a): AMSB Scalar Soft Masses

### Anomalous dimension in the UV electric theory

The composite meson M^i_j = Q^i Q-bar_j has anomalous dimension gamma_M = 2 gamma_Q at one loop, where gamma_Q is the quark anomalous dimension.

For SU(N_c) with fundamental matter:

    gamma_Q = -C_2(fund) g^2 / (8 pi^2) = -C_2(fund) alpha_s / (2 pi)

    C_2(fund) = (N_c^2 - 1) / (2 N_c) = 4/3

    gamma_M = 2 gamma_Q = -C_2(fund) alpha_s / pi

One-loop values:

| alpha_s | gamma_Q | gamma_M |
|---------|---------|---------|
| 0.3 | -0.0637 | -0.1273 |
| 0.5 | -0.1061 | -0.2122 |
| 1.0 | -0.2122 | -0.4244 |

Note: The NSVZ exact beta function gives a would-be fixed point gamma* = 1 - 3N_c/N_f = -2 for N_f = N_c = 3. This is below the unitarity bound (gamma > -1 for scalars), confirming that N_f = N_c lies outside the conformal window. The theory confines rather than flowing to a fixed point. Consequently, gamma_M is not a fixed-point value but a running quantity evaluated at mu = Lambda. Its value at strong coupling (alpha_s ~ 1) is the relevant number.

### Beta function coefficient

    b_0 = 3 N_c - N_f = 6

    beta_g = -b_0 g^3 / (16 pi^2)

### Running of the anomalous dimension

    dot{gamma}_M = d gamma_M / d ln mu
                 = (d gamma_M / d g) * beta_g
                 = [-C_2 g / (4 pi^2)] * [-b_0 g^3 / (16 pi^2)]
                 = C_2 b_0 g^4 / (64 pi^4)
                 = C_2 b_0 alpha_s^2 / (4 pi^2)     [per quark]
                 = 2 * C_2 b_0 alpha_s^2 / (4 pi^2)  [for composite M]
                 = C_2 b_0 alpha_s^2 / (2 pi^2)

Derivation detail:

    d gamma_Q / d g = -C_2 * 2g / (8 pi^2) = -C_2 g / (4 pi^2)

    dot{gamma}_Q = [-C_2 g / (4 pi^2)] * [-b_0 g^3 / (16 pi^2)]
                 = C_2 b_0 g^4 / (64 pi^4)
                 = C_2 b_0 alpha_s^2 / (4 pi^2)

    dot{gamma}_M = 2 dot{gamma}_Q = C_2 b_0 alpha_s^2 / (2 pi^2)
                 = (4/3)(6) alpha_s^2 / (2 pi^2) = 8 alpha_s^2 / (2 pi^2)

### AMSB scalar soft mass-squared

The AMSB formula for scalar soft masses:

    m^2_phi(AMSB) = -(1/4) dot{gamma}_phi m^2_{3/2}

For the composite meson:

    m^2_M(AMSB) = -(1/4) dot{gamma}_M m^2_{3/2}
                = -C_2 b_0 alpha_s^2 / (8 pi^2) m^2_{3/2}
                = -0.1013 alpha_s^2 m^2_{3/2}

Important: dot{gamma}_M > 0 (since the anomalous dimension gamma_M is negative and becomes MORE negative as g increases, so d gamma/d ln mu > 0 in the asymptotically free direction). Therefore **m^2_AMSB < 0** (tachyonic). This is the well-known tachyonic slepton problem of AMSB, now occurring for the meson scalars.

### Numerical results

| alpha_s | dot{gamma}_M | m^2_AMSB (m_{3/2}=30 TeV) | m^2_AMSB (m_{3/2}=50 TeV) | m^2_AMSB (m_{3/2}=100 TeV) |
|---------|-------------|---------------------------|---------------------------|----------------------------|
| 0.3 | 0.0365 | -(2865 GeV)^2 | -(4775 GeV)^2 | -(9549 GeV)^2 |
| 0.5 | 0.1013 | -(4775 GeV)^2 | -(7958 GeV)^2 | -(15915 GeV)^2 |
| 0.8 | 0.2594 | -(7639 GeV)^2 | -(12732 GeV)^2 | -(25465 GeV)^2 |
| 1.0 | 0.4053 | -(9549 GeV)^2 | -(15915 GeV)^2 | -(31831 GeV)^2 |

All entries are **tachyonic** (negative m^2). The magnitude ranges from 3 TeV to 32 TeV depending on alpha_s and m_{3/2}.

---

## Part (b): AMSB Contributions to Mesino Masses

### Gaugino masses

In the confined phase, the magnetic gauge group is SU(N_f - N_c) = SU(0), i.e., trivial. There are **no gauginos** in the IR effective theory. The AMSB gaugino mass formula M_lambda = (b_0 g^2 / (16 pi^2)) m_{3/2} does not apply.

In the UV electric theory (above Lambda), the SU(3) gaugino acquires an AMSB mass:

    M_gaugino(UV) = (b_0 alpha_s / (4 pi)) m_{3/2}

| alpha_s | M_g (30 TeV) | M_g (50 TeV) | M_g (100 TeV) |
|---------|-------------|-------------|--------------|
| 0.3 | 4.3 TeV | 7.2 TeV | 14.3 TeV |
| 0.5 | 7.2 TeV | 11.9 TeV | 23.9 TeV |
| 1.0 | 14.3 TeV | 23.9 TeV | 47.7 TeV |

After confinement, the gaugino is absorbed into the composite states (mesinos). The UV gaugino AMSB mass is transmitted to the confined spectrum through matching at the confinement scale.

### Fermion masses: no direct AMSB contribution

**Critical point.** The AMSB mechanism generates SOFT SUSY-BREAKING terms:

1. Scalar mass-squared: m^2_i = -(1/4) dot{gamma}_i m^2_{3/2}
2. B-term: B_{ij} = -(gamma_i + gamma_j) m_{3/2}
3. A-term: A_{ijk} = -(gamma_i + gamma_j + gamma_k) m_{3/2}
4. Gaugino mass: M_a = (beta_{g_a}/g_a) m_{3/2}

**None of these generate fermion mass terms.** The fermion masses are determined by the holomorphic superpotential W_{IJ} = d^2 W / d Phi_I d Phi_J, which is protected by holomorphy and does not receive soft SUSY-breaking corrections.

The claim "delta m_mesino = gamma_M * m_{3/2}" conflates the B-term (which is a scalar bilinear parameter) with a fermion mass. The B-term enters the scalar potential as:

    V supset (B_{ij} mu_{ij}) phi_i phi_j + h.c.

It shifts the SCALAR mass-squared matrix, not the fermion mass matrix.

### B-term values

    B = -gamma_M * m_{3/2}

| alpha_s | B (30 TeV) | B (50 TeV) | B (100 TeV) |
|---------|-----------|-----------|------------|
| 0.3 | 3.8 TeV | 6.4 TeV | 12.7 TeV |
| 0.5 | 6.4 TeV | 10.6 TeV | 21.2 TeV |
| 1.0 | 12.7 TeV | 21.2 TeV | 42.4 TeV |

The B-term enters the scalar meson bilinear as B * mu_eff, where mu_eff = X_0 M_k ~ 10^{-5} MeV (the holomorphic mass parameter). The product B * mu_eff ~ (10 TeV) * (10^{-5} MeV) ~ 10^{-1} MeV^2, which is negligible compared to m^2_AMSB ~ (10 TeV)^2.

### A-term values

    A = -(gamma_X + 3 gamma_M) m_{3/2}

Taking gamma_X = 0 (X has no gauge coupling in the confined theory):

| alpha_s | A (30 TeV) | A (50 TeV) | A (100 TeV) |
|---------|-----------|-----------|------------|
| 0.3 | 11.5 TeV | 19.1 TeV | 38.2 TeV |
| 0.5 | 19.1 TeV | 31.8 TeV | 63.7 TeV |
| 1.0 | 38.2 TeV | 63.7 TeV | 127.3 TeV |

### Indirect effect: vacuum destabilization

While AMSB does not directly generate mesino masses, the tachyonic scalar soft masses m^2_AMSB ~ -(10 TeV)^2 completely overwhelm the tree-level soft mass f_pi^2 = (92 MeV)^2 by a factor of ~ 10^{10}. The Seiberg seesaw vacuum is **catastrophically destabilized**.

The scalar fields roll to new VEVs determined by the balance of:
- Tachyonic AMSB masses (negative m^2)
- Quartic terms from |F|^2
- Kahler potential corrections
- Possible D-term contributions

At the new vacuum, ALL masses (including mesino fermion masses) are reorganized. The mesino masses at the new vacuum are not calculable without the full scalar potential.

---

## Part (c): Summary Spectrum Table

Reference values: alpha_s = 1.0 (strong coupling at confinement), m_{3/2} = 50 TeV.

### AMSB soft terms

| Soft term | Formula | Value |
|-----------|---------|-------|
| Scalar m^2_AMSB | -(1/4) dot{gamma}_M m^2_{3/2} | -(15,915 GeV)^2 (tachyonic) |
| B-term | -gamma_M m_{3/2} | +21,221 GeV |
| A-term (det M) | -3 gamma_M m_{3/2} | +63,662 GeV |
| UV gaugino mass | (b_0 alpha_s / 4pi) m_{3/2} | 23,873 GeV |
| Gaugino mass (IR) | 0 | 0 (no gauge group) |

### Full spectrum comparison

| Field | Tree-level mass | AMSB soft term | Total status |
|-------|----------------|----------------|--------------|
| Scalar mesons (off-diag) | 92 MeV (from f_pi^2) | m^2 = -(15,915 GeV)^2 | **Tachyonic**, destabilizes vacuum |
| Scalar mesons (diag) | 92 MeV (from f_pi^2) | m^2 = -(15,915 GeV)^2 | **Tachyonic**, destabilizes vacuum |
| Mesino psi_{ud} (Q = +/-1) | 11.4 eV | No direct AMSB | Remains at 11.4 eV (but vacuum invalid) |
| Mesino psi_{us} (Q = +/-1) | 228.7 eV | No direct AMSB | Remains at 228.7 eV (but vacuum invalid) |
| Mesino psi_{ds} (Q = 0) | 494.4 eV | No direct AMSB | Remains at 494.4 eV (but vacuum invalid) |
| Diag mesinos (2 modes) | 7.5 keV, 52.5 keV | No direct AMSB | Same (but vacuum invalid) |
| Scalar baryons | ~300 MeV | m^2 = -(15,915 GeV)^2 | **Tachyonic** |
| Baryon fermions | ~1.2 x 10^{-9} MeV | No direct AMSB | Unchanged |
| X (Lagrange mult.) | No kinetic term | N/A | Auxiliary field |
| UV gaugino (above Lambda) | 0 | 23,873 GeV | Confined into composites |

### Scale hierarchy

    f_pi        = 92 MeV          (tree-level soft mass)
    M_FCNC      = 745 MeV         (K-Kbar bound, g=1)
    LEP bound   = 100 GeV         (charged fermion)
    |m_AMSB|    = 15,915 GeV      (scalar soft mass magnitude)
    B-term      = 21,221 GeV
    UV M_gaugino = 23,873 GeV
    m_{3/2}     = 50,000 GeV
    A-term      = 63,662 GeV

The AMSB contributions are 5 orders of magnitude above the tree-level soft mass and 4 orders of magnitude above the FCNC bound.

---

## Part (d): Phenomenological Implications

### (d1) Are the mesinos pushed above the LEP bound (100 GeV)?

**No, not by AMSB directly.** AMSB does not generate fermion masses. The mesino masses are holomorphic quantities (W_{IJ}) that are not shifted by soft SUSY-breaking. At the Seiberg seesaw vacuum (if it survives), the mesino masses remain at their tree-level values of 11-494 eV.

However, the seesaw vacuum does NOT survive: the tachyonic AMSB scalar masses (of order 16 TeV) overwhelm the tree-level scalar potential (of order 92 MeV) by 10^{10}. The meson fields roll to a completely different vacuum.

At the new vacuum, the mesino masses are determined by:
- The new VEVs (set by balancing AMSB tachyonic masses against other contributions)
- The superpotential evaluated at the new VEVs
- Possible new mass terms generated by the vacuum rearrangement

Without additional stabilization (D-terms, Kahler corrections), the new vacuum is M = 0 (all meson VEVs vanish), where the mesino mass matrix degenerates and all mesino masses become **exactly zero**. This is worse than the seesaw vacuum.

### (d2) Are the scalar mesons pushed above FCNC bounds (745 MeV)?

**Not in the useful sense.** The AMSB scalar m^2 has magnitude (15,915 GeV)^2, which is far above (745 MeV)^2. But the sign is **wrong**: m^2_AMSB < 0. Tachyonic scalars do not decouple from low-energy physics; they condense, destabilizing the vacuum.

To satisfy FCNC bounds, the scalar mesons need **positive** m^2 >> (745 MeV)^2. This requires additional positive contributions:

1. **D-terms from gauged SU(5) flavor**: If SU(5)_F is gauged, the D-term contribution m^2_D ~ g_F^2 D / (16 pi^2) can be positive and large. For g_F ~ 1 and D ~ Lambda^2: m^2_D ~ (Lambda / 4pi)^2 ~ (24 MeV)^2, too small. One needs D ~ m^2_{3/2} to compete with AMSB.

2. **Deflected AMSB**: Adding messenger fields that couple to the meson sector can generate positive threshold corrections that overcome the tachyonic AMSB masses (Pomarol & Rattazzi, 1999; Rattazzi & Sarid, 2000).

3. **Kahler potential corrections**: Higher-dimensional Kahler operators can generate positive scalar masses of order Lambda^2. At the seesaw vacuum where M >> Lambda, these corrections can be enormous but are incalculable.

### (d3) What value of m_{3/2} is needed for LEP-safe mesinos?

Since AMSB does not directly generate fermion masses, the question has no direct answer in the AMSB framework. The fermion mass problem must be addressed through one of:

**Scenario A: Vacuum rearrangement.** If the vacuum is stabilized (by D-terms or Kahler corrections) at VEVs of order m_{3/2}, then mesino masses are parametrically:

    m_mesino ~ (coupling) * m_{3/2}

For the coupling set by the det M superpotential, this gives:

    m_mesino ~ alpha_s/(4pi) * m_{3/2}   (loop factor)

The LEP bound then requires:

| alpha_s | m_{3/2} for m_mesino > 100 GeV |
|---------|-------------------------------|
| 0.3 | > 4.2 TeV |
| 1.0 | > 1.3 TeV |

**Scenario B: Wave function renormalization.** The anomalous dimension gamma produces Z_M = (Lambda/M_UV)^{2 gamma}, giving m_phys = m_hol / Z_M. For gamma = 0.26 with M_UV = M_Planck, all charged mesinos are above 100 GeV (see mesino_anomdim.md). This mechanism is independent of m_{3/2}.

**Scenario C: UV gaugino mass transmission.** The AMSB gaugino mass M_g ~ (b_0 alpha_s / 4pi) m_{3/2} is generated above Lambda. After confinement, this mass is distributed among the composite fermions. If the mesino receives a fraction of the gaugino mass:

    m_mesino ~ M_g ~ (b_0 alpha_s / 4pi) m_{3/2}

This gives the same parametric estimate as Scenario A. For LEP safety with alpha_s = 1:

    m_{3/2} > 100 GeV * (4 pi) / (6 * 1.0) = 209 GeV

This is a very mild requirement. The issue is not the magnitude but whether the gaugino mass actually transmits to the mesinos through the confinement transition. There is no controlled calculation of this matching.

---

## Error Analysis of the Problem Statement

The problem statement suggests "delta m_mesino ~ -(1/3) * 50 TeV ~ -17 TeV" from AMSB. This claim is **incorrect** for the following reasons:

1. **The formula delta m_mesino = gamma_M * m_{3/2} is not an AMSB fermion mass.** The AMSB B-term B = -gamma_M * m_{3/2} is a SCALAR bilinear parameter. It enters the scalar mass-squared matrix as B * mu, where mu is the holomorphic mass. It does NOT generate a fermion mass term.

2. **The value gamma_M ~ -1/3 used in the problem statement.** At alpha_s = 1, the one-loop result is gamma_M = -C_2 alpha_s / pi = -(4/3)/pi = -0.424, not -1/3. At alpha_s = pi/3 ~ 1.05, gamma_M = -4/9, which is closer to -1/3 at alpha_s ~ 0.75.

3. **Even if a fermion mass were generated,** it would be proportional to the holomorphic mass mu = X_0 M_k times the anomalous dimension, not to m_{3/2} alone. The shift would be:

       delta mu / mu = -(gamma_i + gamma_j) * (m_{3/2} / mu)

   which diverges because mu ~ 10^{-5} MeV << m_{3/2} ~ 10^7 MeV. This signals a breakdown of perturbation theory in the AMSB expansion, not a physical 17 TeV mass.

4. **The correct large effect is the scalar m^2_AMSB**, which is of order -(16 TeV)^2 and completely destabilizes the seesaw vacuum. The physical consequence is vacuum rearrangement, not a direct fermion mass.

---

## Key Conclusions

1. **The AMSB scalar mass-squared for mesons is TACHYONIC**, with magnitude (3-32 TeV)^2 depending on alpha_s and m_{3/2}. This is the standard tachyonic slepton problem of AMSB applied to the meson sector.

2. **AMSB does not directly generate mesino fermion masses.** The fermion spectrum is determined by the holomorphic superpotential, not by soft terms. The mesino masses at the seesaw vacuum remain at their tree-level values (11-494 eV).

3. **The seesaw vacuum is catastrophically destabilized** by AMSB-scale tachyonic scalar masses. The scalar potential at the seesaw vacuum is dominated by the AMSB contribution (10^{10} times larger than the tree-level soft mass).

4. **The AMSB scale exceeds the FCNC bound by 4 orders of magnitude** in magnitude. If the sign problem is resolved (by D-terms, deflection, or Kahler corrections), AMSB at m_{3/2} ~ 50 TeV easily satisfies FCNC constraints for the scalar mesons.

5. **For LEP-safe mesinos,** one needs either (a) vacuum rearrangement that gives mesino masses at the AMSB scale (requiring m_{3/2} > 1-4 TeV), or (b) wave function renormalization with gamma > 0.26 (independent of m_{3/2}), or (c) gaugino mass transmission through confinement (requiring m_{3/2} > 200 GeV but with uncontrolled matching).

6. **The tachyonic AMSB scalar masses may be phenomenologically useful**: they destabilize the seesaw vacuum, forcing a vacuum rearrangement that could simultaneously resolve the mesino mass problem (by generating large fermion masses at the new vacuum) and the FCNC problem (by decoupling the scalar mesons). But computing the new vacuum requires knowledge of the full Kahler potential, which is incalculable in the strongly coupled regime.

---

*Computation performed March 2026. Uses one-loop perturbative formulae for anomalous dimensions, which are unreliable at alpha_s ~ 1. The NSVZ exact beta function confirms that N_f = N_c = 3 is outside the conformal window (gamma* = -2 violates unitarity), so the perturbative estimates should be treated as order-of-magnitude guides, not precision calculations.*
