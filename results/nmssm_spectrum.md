# NMSSM Extension of Seiberg Effective Theory Spectrum

## Setup

**Original superpotential:**

$$W = \sum_i m_i M^i_i + X(\det M - B\tilde{B} - \Lambda^6) + y_c H_u^0 M^d_d + y_b H_d^0 M^s_s$$

**NMSSM extension:** Add $\lambda X (H_u^0 H_d^0 - H_u^+ H_d^-)$ with $\lambda = 0.72$.

**16 chiral superfields:** $M^i_j$ (9), $X$ (1), $B, \tilde{B}$ (2), $H_u^+, H_u^0, H_d^0, H_d^-$ (4)

**Parameters:**

| Parameter | Value |
|-----------|-------|
| $m_u$ | 2.16 MeV |
| $m_d$ | 4.67 MeV |
| $m_s$ | 93.4 MeV |
| $\Lambda$ | 300.0 MeV |
| $v$ | 246220.0 MeV |
| $y_c = 2m_c/v$ | 1.03159776e-02 |
| $y_b = 2m_b/v$ | 3.39533750e-02 |
| $\lambda$ | 0.72 |
| $f_\pi$ | 92.0 MeV |
| $\tilde{m}^2 = f_\pi^2$ | 8464.0 MeV$^2$ |

---

## Task 1: Modified F-terms

The NMSSM coupling $\lambda X H_u \cdot H_d$ modifies the F-terms:

| F-term | Expression | Value (MeV) |
|--------|-----------|-------------|
| $F_X$ | $\det M - \Lambda^6 + \lambda H_u^0 H_d^0$ | -3.814697e-06 |
| $F_{H_u^0}$ | $y_c M^d_d + \lambda X H_d^0$ | 1.948966e+03 |
| $F_{H_d^0}$ | $y_b M^s_s + \lambda X H_u^0$ | 3.207351e+02 |
| $F_{H_u^+}$ | $-\lambda X H_d^-$ | 0 (at vacuum) |
| $F_{H_d^-}$ | $-\lambda X H_u^+$ | 0 (at vacuum) |
| $F_{M^u_u}$ | $m_u + X \cdot M_d M_s$ | 4.884981e-15 |
| $F_{M^d_d}$ | $m_d + X \cdot M_u M_s + y_c v/\sqrt{2}$ | 1.796051e+03 |
| $F_{M^s_s}$ | $m_s + X \cdot M_u M_d + y_b v/\sqrt{2}$ | 5.911413e+03 |

**Key change:** $F_{H_u^+}$ and $F_{H_d^-}$ are now structurally nonzero
(proportional to $\lambda X$), though they vanish at the vacuum where
$\langle H_u^+ \rangle = \langle H_d^- \rangle = 0$.

---

## Task 2: Modified vacuum

The constraint $\partial W / \partial X = 0$ now reads:

$$\det M - B\tilde{B} - \Lambda^6 + \lambda H_u^0 H_d^0 = 0$$

At $\langle B \rangle = \langle \tilde{B} \rangle = 0$, this shifts $\det M$ by:

$$\det M = \Lambda^6 - \lambda v^2/2 = 7.290000e+14 - 2.182474e+10 = 7.289782e+14 \text{ MeV}^6$$

The shift is a fraction $2.993792e-05$ of $\Lambda^6$,
so the Seiberg vacuum survives with:

| Quantity | Original | NMSSM-modified | Fractional change |
|----------|----------|----------------|-------------------|
| $C$ | 882297.425722 | 882288.620918 | -9.979406e-06 |
| $M_u$ | 408471.0304 | 408466.9541 | -9.979406e-06 |
| $M_d$ | 188928.7850 | 188926.8996 | -9.979406e-06 |
| $M_s$ | 9446.4392 | 9446.3450 | -9.979406e-06 |
| $X$ | -1.21028453e-09 | -1.21030869e-09 | 1.995911e-05 |

The NMSSM shift to the vacuum is negligible ($\sim 10^{-5}$ fractional change).

---

## Task 3: Fermion mass matrix

The 16x16 matrix $W_{IJ}$ now includes NMSSM cross-terms:

| New entry | Expression | Value (MeV) |
|-----------|-----------|-------------|
| $W_{X, H_u^0}$ | $\lambda \langle H_d^0 \rangle = \lambda v/\sqrt{2}$ | +125354.7588 |
| $W_{X, H_d^0}$ | $\lambda \langle H_u^0 \rangle = \lambda v/\sqrt{2}$ | +125354.7588 |
| $W_{H_u^0, H_d^0}$ | $\lambda X$ | -8.71422257e-10 |
| $W_{H_u^+, H_d^-}$ | $-\lambda X$ | +8.71422257e-10 |

### Block-by-block eigenvalues

| Block | Fields | Eigenvalues (MeV) |
|-------|--------|-------------------|
| central | Muu, Mdd, Mss, X, Hu0, Hd0 | 7.695421e-04, 7.795732e-04, 1.041764e-02, 1.046708e-02, 7.728741e+10, 7.728741e+10 |
| ud_offdiag | Mud, Mdu | 1.143299e-05, 1.143299e-05 |
| us_offdiag | Mus, Msu | 2.286599e-04, 2.286599e-04 |
| ds_offdiag | Mds, Msd | 4.943711e-04, 4.943711e-04 |
| baryon | B, Bt | 1.210309e-09, 1.210309e-09 |
| charged_Higgs | Hu+, Hd- | 8.714223e-10, 8.714223e-10 |

### Comparison: NMSSM vs original fermion eigenvalues

| Block | Mode | Original |m| | NMSSM |m| | Shift |
|-------|------|-----------|----------|-------|
| central | 1 | 7.694586e-04 | 7.695421e-04 | +8.343151e-08 |
| central | 2 | 7.796563e-04 | 7.795732e-04 | -8.318573e-08 |
| central | 3 | 1.041754e-02 | 1.041764e-02 | +1.019891e-07 |
| central | 4 | 1.046720e-02 | 1.046708e-02 | -1.235797e-07 |
| central | 5 | 7.728741e+10 | 7.728741e+10 | +2.032776e-01 |
| central | 6 | 7.728741e+10 | 7.728741e+10 | +2.032776e-01 |
| ud_offdiag | 1 | 1.143299e-05 | 1.143299e-05 | +0.000000e+00 |
| ud_offdiag | 2 | 1.143299e-05 | 1.143299e-05 | +0.000000e+00 |
| us_offdiag | 1 | 2.286599e-04 | 2.286599e-04 | +0.000000e+00 |
| us_offdiag | 2 | 2.286599e-04 | 2.286599e-04 | +0.000000e+00 |
| ds_offdiag | 1 | 4.943711e-04 | 4.943711e-04 | +0.000000e+00 |
| ds_offdiag | 2 | 4.943711e-04 | 4.943711e-04 | +0.000000e+00 |
| baryon | 1 | 1.210309e-09 | 1.210309e-09 | +0.000000e+00 |
| baryon | 2 | 1.210309e-09 | 1.210309e-09 | +0.000000e+00 |
| charged_Higgs | 1 | 0.000000e+00 | 8.714223e-10 | +8.714223e-10 (NEW) |
| charged_Higgs | 2 | 0.000000e+00 | 8.714223e-10 | +8.714223e-10 (NEW) |

---

## Task 4: Charged Higgsino mass

In the original model, $H_u^+$ and $H_d^-$ were **exactly massless**
(decoupled from the fermion mass matrix).

The NMSSM coupling gives them a mass through $W_{H_u^+, H_d^-} = -\lambda X$:

$$m_{\text{charged Higgsino}} = |\lambda X| = |\lambda \cdot C / \Lambda^6| = 8.714223e-10 \text{ MeV} = 8.714223e-13 \text{ GeV}$$

This is extremely small (8.71e-10 MeV) because $|X| = |C/\Lambda^6| \sim 10^{-9}$ MeV$^{-4}$.

---

## Task 5: Scalar mass-squared matrix

The 32x32 scalar mass-squared matrix is constructed block by block:

| Block | Fields | Scalar m$^2$ eigenvalues (MeV$^2$) | Tachyons? |
|-------|--------|-------------------------------------|-----------|
| central | Muu, Mdd, Mss, X, Hu0, Hd0 | -6.0779e+03, -2.0282e-06, +2.2494e-06, +1.0761e-04, +1.0833e-04, +9.0761e+03, +6.7552e+04, +1.4118e+06, +5.9733e+21, +5.9733e+21, +5.9733e+21, +5.9733e+21 | Yes (2) |
| ud_offdiag | Mud, Mdu | +8.4640e+03, +8.4640e+03, +8.4640e+03, +8.4640e+03 | No |
| us_offdiag | Mus, Msu | +8.4633e+03, +8.4633e+03, +8.4647e+03, +8.4647e+03 | No |
| ds_offdiag | Mds, Msd | +8.4624e+03, +8.4624e+03, +8.4656e+03, +8.4656e+03 | No |
| baryon | B, Bt | +1.4648e-18, +1.4648e-18, +1.4648e-18, +1.4648e-18 | No |
| charged_Higgs | Hu+, Hd- | -2.7466e-06, -2.7466e-06, +2.7466e-06, +2.7466e-06 | Yes (2) |

---

## Task 6: Tachyonic modes

| Model | Total tachyonic scalar modes |
|-------|----------------------------|
| Original (no NMSSM) | 0 |
| NMSSM ($\lambda = 0.72$) | 4 |

---

## Task 7: Higgs mass

The NMSSM coupling $\lambda X H_u^0 H_d^0$ generates a Higgs quartic through $|F_X|^2$:

$$V \supset |F_X|^2 \supset \lambda^2 |H_u^0|^2 |H_d^0|^2$$

For $\tan\beta = 1$ ($\langle H_u^0 \rangle = \langle H_d^0 \rangle = v/\sqrt{2}$),
the CP-even Higgs mass at tree level is:

$$m_h = \frac{\lambda v}{\sqrt{2}} = \frac{0.72 \times 246.22 \text{ GeV}}{\sqrt{2}} = 125.3548 \text{ GeV}$$

This matches the physical Higgs mass $m_h = 125.25$ GeV to better than 0.1%.

---

## Task 8: Supertrace

The supertrace is computed analytically:

$$\text{STr}[M^2] = \text{Tr}(m^2_{\text{scalar},R} + m^2_{\text{scalar},I}) - 2\text{Tr}(m^2_{\text{fermion}})$$

$$= 2\text{Tr}(W^2) + 2\text{Tr}(M^2_{\text{soft}}) - 2\text{Tr}(W^2) = 2\text{Tr}(M^2_{\text{soft}})$$

$$= 2 \times 9 \times \tilde{m}^2 = 18 f_\pi^2 = 152352 \text{ MeV}^2$$

**The NMSSM coupling does NOT change the supertrace.** The result $18 f_\pi^2$
is identical to the original model. This is because the supertrace depends only on
soft SUSY-breaking terms, and the NMSSM coupling is an F-term interaction whose
holomorphic contribution cancels between real and imaginary scalar components.

---

## Summary of key results

| Quantity | Original model | NMSSM extension |
|----------|---------------|-----------------|
| Charged Higgsino mass | 0 (decoupled) | 8.714223e-10 MeV |
| $W_{H_u^0, H_d^0}$ | 0 | $\lambda X = -8.714223e-10$ MeV |
| $W_{X, H_u^0}$ | 0 | $\lambda v/\sqrt{2} = 125354.7588$ MeV |
| $W_{X, H_d^0}$ | 0 | $\lambda v/\sqrt{2} = 125354.7588$ MeV |
| Higgs mass (tree level) | --- | $\lambda v/\sqrt{2} = 125.3548$ GeV |
| Tachyonic scalar modes | 0 | 4 |
| STr[$M^2$] | $18 f_\pi^2 = 152352$ MeV$^2$ | $18 f_\pi^2 = 152352$ MeV$^2$ (unchanged) |
| Vacuum shift (det M) | $\Lambda^6$ | $\Lambda^6 - \lambda v^2/2$ (shift $2.993792e-05$) |

The NMSSM coupling $\lambda X H_u \cdot H_d$ with $\lambda = 0.72$:

1. Gives mass to the previously massless charged Higgsinos ($m = |\lambda X| \sim 10^{-9}$ MeV)
2. Introduces X--Higgs cross-couplings of order $\lambda v/\sqrt{2} \sim 125355$ MeV into the central fermion block
3. Generates a tree-level Higgs mass $m_h = \lambda v/\sqrt{2} = 125.4$ GeV
4. Negligibly shifts the Seiberg seesaw vacuum ($\sim 10^{-5}$ fractional change)
5. Does not change the supertrace STr[$M^2$] = $18 f_\pi^2$

---

*Generated by nmssm_spectrum.py*