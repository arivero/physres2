# Two-Higgs-Doublet Scalar Potential at tan(beta) = 1

## Setup

We consider an N=1 SUSY two-Higgs-doublet model with Higgs superpotential terms

W_H = mu_H (H_u . H_d) + y_c H_u M_c + y_b H_d M_b

where H_u, H_d are the MSSM Higgs doublets and M_c, M_b are singlet meson fields
with VEVs determined by an ISS seesaw: <M_j> = C/m_j. The soft SUSY-breaking terms
are

V_soft = m_Hu^2 |H_u|^2 + m_Hd^2 |H_d|^2 + (B_mu H_u . H_d + h.c.)

and the D-term quartic is

V_D = (g^2 + g'^2)/8 (|H_u|^2 - |H_d|^2)^2

We set mu_H = 0 (pure meson-mediated SUSY breaking) and work at tan(beta) = v_u/v_d = 1 throughout.

**Numerical inputs:**
- v = 246.22 GeV, m_Z = 91.1876 GeV
- g = 0.653, g' = 0.350
- m_c = 1270 MeV, m_b = 4180 MeV, m_t = 172.76 GeV
- C = 882297 MeV^2 (SQCD seesaw scale)
- y_c = 2 m_c/v = 0.01032, y_b = 2 m_b/v = 0.03395


## 1. The full scalar potential

The F-terms derived from W_H (with mu_H = 0) give:

|F_{M_c}|^2 = y_c^2 |H_u^0|^2,    |F_{M_b}|^2 = y_b^2 |H_d^0|^2

These contribute to the potential as effective mass terms, not quartic couplings.
Defining effective squared masses m_1^2 = m_Hd^2 + y_b^2 and m_2^2 = m_Hu^2 + y_c^2,
the full neutral scalar potential is

V = m_2^2 |H_u^0|^2 + m_1^2 |H_d^0|^2 - b(H_u^0 H_d^0 + h.c.) + (g^2+g'^2)/8 (|H_u^0|^2 - |H_d^0|^2)^2

where b > 0 is the soft bilinear (related to B_mu by a sign convention).

The crucial point: the meson Yukawa couplings produce mass-dimension-two terms (y_j^2 |H|^2),
not quartic couplings. The ONLY quartic is the D-term. This is the standard MSSM structure,
with the Yukawa F-terms absorbed into redefined soft masses.


## 2. EWSB conditions at tan(beta) = 1

Parameterizing H_u^0 = v_u/sqrt(2), H_d^0 = v_d/sqrt(2), the tadpole equations are

dV/dv_u = 0:   m_2^2 v_u - b v_d + (g_Z^2/4) v_u (v_u^2 - v_d^2) = 0
dV/dv_d = 0:   m_1^2 v_d - b v_u - (g_Z^2/4) v_d (v_u^2 - v_d^2) = 0

At tan(beta) = 1 (v_u = v_d = v/sqrt(2)), the D-term vanishes identically, and both
equations reduce to

m_1^2 = m_2^2 = b

This gives two constraints on the soft parameters:

1. **Soft mass splitting tracks Yukawa splitting:**
   m_Hu^2 - m_Hd^2 = y_b^2 - y_c^2 = 1.05 x 10^{-3}

2. **Each soft mass is determined by b:**
   m_Hu^2 = b - y_c^2,   m_Hd^2 = b - y_b^2

The bilinear b remains as a free parameter that controls the pseudoscalar mass.


## 3. CP-even Higgs mass matrix

In the (h_d, h_u) basis (real fluctuations around the VEVs), the CP-even mass matrix is

M^2 = (m_A^2 + m_Z^2)/2 * [[1, -1], [-1, 1]]

where m_A^2 = m_1^2 + m_2^2 = 2b is the pseudoscalar mass squared.

The eigenvalues are

- m_h^2 = 0  (the D-flat direction h = (h_u + h_d)/sqrt(2))
- m_H^2 = m_A^2 + m_Z^2  (the D-term direction H = (h_u - h_d)/sqrt(2))

**The lightest CP-even Higgs is massless at tree level.**

This is the famous MSSM result at tan(beta) = 1: the D-term quartic, proportional to
(|H_u|^2 - |H_d|^2)^2, vanishes along the D-flat direction h_u = h_d. Since the meson
F-terms contribute only mass terms (not quartics), they cannot lift this flat direction.


## 4. Complete tree-level spectrum

| State | Mass^2 formula | Notes |
|-------|---------------|-------|
| h (CP-even) | 0 | Massless: D-flat direction |
| H (CP-even) | m_A^2 + m_Z^2 | D-term massive |
| A (CP-odd) | m_A^2 = 2b | Physical pseudoscalar |
| G^0 | 0 | Goldstone, eaten by Z |
| H^pm (charged) | m_A^2 + m_W^2 | m_W = g v/2 = 80.39 GeV |
| G^pm | 0 | Goldstones, eaten by W |

Numerical spectrum for representative values of m_A:

| m_A (GeV) | m_H (GeV) | m_H^pm (GeV) | b (GeV^2) |
|-----------|-----------|--------------|-----------|
| 50 | 104.0 | 94.7 | 1250 |
| 100 | 135.3 | 128.3 | 5000 |
| 200 | 219.8 | 215.5 | 20000 |
| 300 | 313.6 | 310.6 | 45000 |
| 500 | 508.3 | 506.4 | 125000 |
| 1000 | 1004.2 | 1003.2 | 500000 |


## 5. Soft masses from F-term SUSY breaking

The meson F-terms at the EWSB minimum give

F_{M_c} = y_c v_u / sqrt(2) = y_c v/2 = m_c = 1270 MeV
F_{M_b} = y_b v_d / sqrt(2) = y_b v/2 = m_b = 4180 MeV

**The F-term of each meson singlet equals the corresponding quark mass.**

The Higgs F-terms involve the meson VEVs:

F_{H_u} = y_c <M_c> = (2 m_c/v)(C/m_c) = 2C/v = 7.167 MeV
F_{H_d} = y_b <M_b> = (2 m_b/v)(C/m_b) = 2C/v = 7.167 MeV

The ISS seesaw produces a remarkable cancellation: the Yukawa y_j proportional to m_j and the
meson VEV proportional to 1/m_j combine to give a flavor-universal result. The Higgs F-terms
do not depend on which quark flavor mediates them.

The SUSY-breaking scale is

F = C/v = 3.583 MeV = 0.003583 GeV

The effective soft mass-squared contributions from these F-terms are

y_c^2 = 1.064 x 10^{-4} (dimensionless, from |F_{M_c}|^2 / v^2)
y_b^2 = 1.153 x 10^{-3} (dimensionless, from |F_{M_b}|^2 / v^2)


## 6. The tree-level Higgs mass and the MSSM bound

The MSSM tree-level upper bound on the lightest Higgs mass is

m_h <= m_Z |cos(2 beta)|

At tan(beta) = 1:

cos(2 beta) = cos(pi/2) = 0

so the bound gives m_h = 0 identically. This is exactly what our explicit diagonalization
confirms: the lightest CP-even eigenvalue is zero.

**Physical interpretation:** At tan(beta) = 1, the two Higgs doublets contribute equally to
EWSB. The D-term potential (g^2+g'^2)/8 (|H_u|^2 - |H_d|^2)^2 only provides a quartic
along the relative-size direction |H_u|^2 - |H_d|^2. The overall-size direction
|H_u|^2 + |H_d|^2 has no quartic at all -- it is a D-flat direction. Without a quartic
restoring force, the curvature of V along this direction is zero: hence m_h = 0.

The meson Yukawa F-terms (y_c^2 |H_u|^2, y_b^2 |H_d|^2) only add mass terms that are
already absorbed into the EWSB conditions. They do not generate a quartic coupling that
could lift the flat direction.

To achieve m_h = 125 GeV, two options exist:

**(a) Radiative correction from stop loops (pure MSSM):**
At tan(beta) = 1, the leading-log one-loop correction is

Delta m_h^2 = (3 y_t^4 v^2) / (32 pi^2) * ln(m_stop^2 / m_t^2)

with y_t = sqrt(2) m_t/v = 0.992. To achieve m_h = 125.25 GeV requires
ln(m_stop^2/m_t^2) = 28.1, giving m_stop = 2 x 10^8 GeV. This is unphysical.

The suppression factor sin^4(beta) = 1/4 at tan(beta) = 1 makes the top-loop correction
four times less effective than at large tan(beta), where sin(beta) -> 1.

**(b) NMSSM singlet coupling (lambda S H_u . H_d):**
This generates an extra quartic lambda^2 |H_u . H_d|^2, contributing

Delta m_h^2 = lambda^2 v^2 sin^2(2 beta) / 2 = lambda^2 v^2 / 2

at sin(2 beta) = 1. For m_h = 125.25 GeV:

lambda = sqrt(2) m_h / v = 0.719

This is a perturbative coupling. The NMSSM naturally accommodates tan(beta) = 1 with
a moderate singlet coupling.


## 7. ISS seesaw identity: flavor universality

The most structurally significant result is the flavor universality of the Higgs F-terms.
The ISS seesaw gives meson VEVs <M_j> = C/m_j, and the Yukawa couplings are y_j = 2m_j/v.
Their product is

y_j <M_j> = (2 m_j / v)(C / m_j) = 2C/v

independent of the quark flavor j. Numerically:

- y_c <M_c> = 0.01032 x 694.72 = 7.1667 MeV
- y_b <M_b> = 0.03395 x 211.08 = 7.1667 MeV
- 2C/v = 2 x 882297 / 246220 = 7.1667 MeV

This identity holds to machine precision and is algebraically exact.

Consequences:
1. F_{H_u} = F_{H_d} = 2C/v regardless of which quark mediates.
2. The Higgs sector is flavor-blind despite coupling to different quark flavors.
3. The soft bilinear b is the only free parameter controlling the heavy Higgs spectrum.
4. FCNC constraints are automatically satisfied by the flavor-universal structure.


## Output files

- Python script: `results/higgs_potential.py`
- This writeup: `results/higgs_potential.md`
