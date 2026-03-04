# FCNC Constraints on Off-Diagonal Meson Scalars in Confined SQCD

## Model Setup

N=1 SUSY with confined SU(3) SQCD (N_f = N_c = 3). The low-energy spectrum
contains 9 complex meson scalars M^i_j and their fermionic superpartners
(mesinos). The meson scalars have universal soft mass m_tilde^2 = f_pi^2 = (92 MeV)^2.

### Diagonal meson VEVs (Seiberg seesaw)

| Field | VEV (MeV) |
|-------|-----------|
| M^u_u | 408471 |
| M^d_d | 188929 |
| M^s_s | 9446 |
| X | -1.2103e-09 MeV^{-4} |

### Mesino masses at the seesaw vacuum

| Pair | Mass (MeV) | Expression |
|------|-----------|------------|
| d-s  | 4.9437e-04 | |X| M_u |
| u-s  | 2.2866e-04 | |X| M_d |
| u-d  | 1.1433e-05 | |X| M_s |

The mesino masses are extremely small (O(10^-4) MeV), far below the kaon mass.

---

## (a) Tree-Level K-Kbar Mixing

The off-diagonal meson M^d_s mediates Delta S = 2 transitions at tree level.
The Wilson coefficient for the (LL)(LL) operator is:

C_1 = g^2 / (2 m_S^2)

where g is the effective d-s-scalar coupling and m_S is the scalar mass.

### Results for m_S = f_pi = 92 MeV

| Coupling | C_1 (GeV^-2) | UTfit bound (GeV^-2) | Ratio | Status |
|----------|-------------|---------------------|-------|--------|
| g = 1 | 5.9074e-11 | 9.0e-13 | 6.56e+01 | EXCLUDED |
| g = NDA (7.3) | 3.1095e-09 | 9.0e-13 | 3.46e+03 | EXCLUDED |
| g = V_ds (0.225) | 2.9906e-12 | 9.0e-13 | 3.32e+00 | EXCLUDED |

The K-Kbar mass difference contribution Dm_K = C_1 * (8/3) f_K^2 m_K B_K exceeds
the experimental value by a factor of 4e+14 for g = 1.

---

## (b) Mesino Box Diagrams

Box diagrams with mesino propagators give:

C_1^(box) = g^4 / (16 pi^2 m_mesino^2)

With m_mesino(d-s) = 4.9437e-04 MeV, these contributions are
even more severely excluded than tree-level exchange:

- g = 1: C_1^(box) = 2.5911e+04 MeV^-2, ratio to bound = 2.88e+10
- g = NDA: C_1^(box) = 7.1793e+07 MeV^-2, ratio to bound = 7.98e+13

At the M = 0 vacuum, mesinos are massless (all F-term couplings vanish),
but the effective coupling to quarks also vanishes.

---

## (c) Wilson Coefficient Comparison

### Effective 4-quark coupling at the kaon scale

For off-diagonal mesons with m_S^2 = f_pi^2 = 8464 MeV^2:

- C_1(g=1) = 1/(2 * 92.0^2) = 5.907372e-05 MeV^-2 = 5.907372e-11 GeV^-2
- UTfit bound: |Re C_1| < 9.0e-13 GeV^-2
- **Ratio: 7e+01**

Maximum coupling allowed at m_S = 92.0 MeV: g < 1.23e-01

---

## (d) D-Dbar and B-Bbar Constraints

| System | M_FCNC (g=1) | M_FCNC (g=NDA) | Model m_S | Tension factor |
|--------|-------------|---------------|-----------|---------------|
| K-Kbar (Delta S=2) | 0.7 GeV | 5 GeV | 0.092 GeV | 8x (g=1) |
| D-Dbar (Delta C=2) | 1.1 GeV | 8 GeV | 0.092 GeV | 12x (g=1) |
| B_d-Bbar (Delta B=2) | 0.1 GeV | 1 GeV | 0.092 GeV | 2x (g=1) |
| B_s-Bbar (Delta B=2) | 0.0 GeV | 0 GeV | 0.092 GeV | 0x (g=1) |

K-Kbar is the most stringent constraint by far.

---

## (e) Required M_FCNC Scale

The critical question: is m_tilde^2 = f_pi^2 compatible with FCNC bounds?

**No.** The K-Kbar constraint requires:

- For g = 1: M_FCNC > 0.7 GeV (f_pi = 0.092 GeV)
- For g = NDA: M_FCNC > 5 GeV

The off-diagonal meson scalars at m_S = f_pi are **8x too light** (g = 1)
or **59x too light** (g = NDA) to satisfy K-Kbar mixing constraints.

Required soft mass squared: m_tilde^2 > 5.56e+05 MeV^2 (g = 1), much larger than f_pi^2 = 8464 MeV^2.

---

## (f) Mass at the Stabilized Vacuum

At the M = 0 vacuum (true global minimum of the scalar potential):
- Off-diagonal scalar mass^2 = 2 m_tilde^2 = 2 f_pi^2 = 16928 MeV^2
- Scalar mass = sqrt(2) f_pi = 130.1 MeV

At the seesaw vacuum (metastable):
- SUSY contribution: |X M_u|^2 = 2.4440e-07 MeV^2 (negligible)
- Soft contribution: f_pi^2 = 8464 MeV^2 (dominant)
- Total: m_scalar ~ 92.0 MeV

In **neither** vacuum do the off-diagonal scalars acquire masses large enough
to satisfy FCNC bounds. The tachyonic stabilization at M = 0 gives
m_scalar ~ 130 MeV, which is only marginally larger than f_pi.

---

## Possible Resolutions

1. **Non-universal soft masses**: If SUSY-breaking generates m_tilde^2(offdiag) >> f_pi^2
   while keeping m_tilde^2(diag) ~ f_pi^2, the off-diagonal scalars can be pushed above
   M_FCNC. Required enhancement: ~66x (for g = 1).

2. **Off-diagonal mesons decouple at Lambda**: If these states are confined above the
   QCD scale and not present in the low-energy spectrum, M_FCNC ~ Lambda = 300.0 MeV,
   which is marginal for g = 1 (requires further suppression).

3. **Flavour alignment / GIM mechanism**: If the meson-quark coupling matrix is aligned
   with the quark mass matrix, FCNC couplings are CKM-suppressed. With g_eff ~ V_ds ~ 0.225,
   M_FCNC drops to 0.2 GeV, still requiring scalars well above f_pi.

4. **Gauging flavour SU(3)**: D-term contributions from a gauged flavour symmetry could
   generate large masses for off-diagonal scalars proportional to the D-term VEV.

5. **RG running**: Above the confinement scale, the Kahler potential receives large anomalous
   dimension corrections that could enhance off-diagonal scalar masses.

---

## Key Numbers Summary

| Quantity | Value |
|----------|-------|
| f_pi | 92.0 MeV |
| m_tilde^2 = f_pi^2 | 8464 MeV^2 |
| m_scalar (M = 0 vacuum) | 130.1 MeV |
| m_scalar (seesaw vacuum) | 92.0 MeV |
| m_mesino(d-s) at seesaw | 4.9437e-04 MeV |
| M_FCNC (K-Kbar, g=1) | 745 MeV = 0.7 GeV |
| M_FCNC (K-Kbar, g=NDA) | 5408 MeV = 5 GeV |
| M_FCNC (D-Dbar, g=1) | 1078 MeV = 1.1 GeV |
| M_FCNC (Bd-Bbar, g=1) | 147 MeV = 0.15 GeV |
| Tension factor (K, g=1) | 8x |
| Tension factor (K, g=NDA) | 59x |
| Max coupling at m_S = f_pi | 1.23e-01 |

---

*Generated by fcnc_constraints.py*
