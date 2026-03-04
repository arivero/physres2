# ISS Coleman-Weinberg Potential and Koide Transmission

## Setup

SQCD with $N_c = 3$, $N_f = 4$. ISS window: $N_c < N_f < \frac{3}{2}N_c$, i.e., $3 < 4 < 4.5$.

Magnetic dual: $N_c^{\text{mag}} = N_f - N_c = 1$ magnetic color.

Fields:
- Meson singlet: $\Phi^i_j$ ($4 \times 4$ matrix)
- Magnetic quarks: $q^i$, $\tilde{q}_j$ (single magnetic color, 4 flavors)

Superpotential:
$$W = h \, q^i \Phi_i^j \tilde{q}_j - h\mu^2 \operatorname{Tr}\Phi + \operatorname{Tr}(m\,\Phi)$$

Parameters:
| Parameter | Value | Units |
|-----------|-------|-------|
| $h$ (magnetic Yukawa) | $\sqrt{3} = 1.7321$ | |
| $\mu$ (SUSY-breaking) | 92.0 | MeV |
| $\mu^2$ | 8464 | MeV$^2$ |
| $h\mu^2$ | 14660 | MeV$^2$ |
| $\Lambda$ (cutoff) | 300 | MeV |
| $m_u$ | 2.16 | MeV |
| $m_d$ | 4.67 | MeV |
| $m_s$ | 93.4 | MeV |
| $m_c$ | 1270 | MeV |

## 1. Tree-Level Spectrum

### Vacuum structure

The ISS metastable vacuum has the magnetic quark VEV aligned with the lightest flavor (up):
$$\langle q^1 \rangle = \langle \tilde{q}_1 \rangle = \phi_0 = \sqrt{\mu^2 - m_u/h} = 91.993 \text{ MeV}$$

The rank condition ($N_f - N_c^{\text{mag}} = 3$ directions have $F \neq 0$) gives SUSY breaking with F-terms:
$$F_a = h\mu^2 - m_a \quad (a = 2, 3, 4 \text{ for d, s, c})$$

### Field decomposition

Decompose into sectors at the ISS vacuum:

**Sector A (pseudo-modulus):** $X = \Phi^1_1$, $\rho = q^1$, $\tilde{\rho} = \tilde{q}_1$
- Goldstino: mass 0
- Heavy fermion pair: mass $\sqrt{2}\,h\,\phi_0 = 225.3$ MeV
- Pseudo-modulus $X$: tree-level FLAT direction

**Sector B (link fields):** For each broken flavor $a = d, s, c$:
$\psi_a = q^a$, $\tilde{\psi}_a = \tilde{q}_a$, $Y_a = \Phi^1_a$, $\tilde{Y}_a = \Phi^a_1$

| Flavor | $m^2_{B+}$ (MeV$^2$) | $m^2_{B-}$ (MeV$^2$) | $m_F$ (MeV) |
|--------|-----------------------|-----------------------|--------------|
| d | 50772 | 4.35 | 159.3 |
| s | 50618 | 158.0 | 159.3 |
| c | 48581 | 2196 | 159.3 |

Boson masses: $m^2_{B\pm} = h^2\phi_0^2 \pm h\,F_a$

Fermion mass: $m_F = h\,\phi_0 = 159.3$ MeV (flavor-universal)

Y fields: $m^2_Y = h^2\phi_0^2 = 25388$ MeV$^2$ (no SUSY-breaking splitting)

**Sector C (broken mesons):** $\chi_a = \Phi^a_a$ for $a = d, s, c$
- Tree-level MASSLESS (both scalar and fermion)
- Lifted by one-loop CW potential

### Supertrace check

$$\text{STr}[M^0] = 0, \quad \text{STr}[M^2] = 0 \quad \text{(verified for each flavor)}$$

$$\text{STr}[M^4] \neq 0 \quad \text{(drives the CW potential)}$$

| Flavor | STr$[M^4]$ (MeV$^4$) |
|--------|----------------------|
| u | $2.578 \times 10^9$ |
| d | $2.577 \times 10^9$ |
| s | $2.546 \times 10^9$ |
| c | $2.152 \times 10^9$ |

## 2. Coleman-Weinberg Potential

$$V_{\text{CW}} = \frac{1}{64\pi^2} \operatorname{STr}\left[M^4\left(\log\frac{M^2}{\Lambda^2} - \frac{3}{2}\right)\right]$$

### CW mass for broken mesons $\chi_a$

The CW potential lifts the $\chi_a$ flat directions. Computing $m^2_{\text{CW}} = d^2 V_{\text{CW}} / d\chi_a^2|_{\chi_a=0}$ numerically (5-point stencil):

| Flavor | $m_q$ (MeV) | $m^2_{\text{CW}}$ (MeV$^2$) | $m_{\text{CW}}$ (MeV) | $m_{\text{CW}}/m_q$ |
|--------|-------------|------------------------------|------------------------|----------------------|
| u | 2.16 | 745.3 | 27.300 | 12.64 |
| d | 4.67 | 744.9 | 27.294 | 5.844 |
| s | 93.4 | 733.5 | 27.082 | 0.290 |
| c | 1270 | 599.1 | 24.477 | 0.019 |

### Key observation: near-degeneracy

The CW masses are nearly DEGENERATE (within 10%) despite the UV quark masses spanning 3 orders of magnitude. This is because:

$$m^2_{\text{CW}} \propto F_a^2 = (h\mu^2 - m_a)^2 \approx (h\mu^2)^2\left(1 - \frac{2m_a}{h\mu^2}\right)$$

The fractional flavor dependence is $\Delta m^2 / m^2 \approx -2m_a/(h\mu^2)$:

| Flavor | $m_a/(h\mu^2)$ | Fractional correction |
|--------|-----------------|----------------------|
| u | 0.015% | negligible |
| d | 0.032% | negligible |
| s | 0.64% | small |
| c | 8.7% | moderate |

## 3. Koide Analysis

### UV quark Koide vs CW meson Koide

| Triple | $Q$ (UV quarks) | $Q$ (CW mesons) | $Q$ (Seesaw $M \sim 1/m$) |
|--------|----------------|-----------------|---------------------------|
| $(u,d,s)$ | 0.5670 | 0.3333 | 0.4426 |
| $(d,s,c)$ | 0.6073 | 0.3335 | 0.6389 |
| $(u,s,c)$ | 0.6242 | 0.3335 | 0.7197 |
| $(u,d,c)$ | 0.8281 | 0.3335 | 0.4942 |

The CW meson spectrum gives $Q \approx 1/3$ (degenerate limit), independent of the UV Koide ratio.

### Flavor-dependence of CW shift

The ISS CW maps quark masses through $f(m_a) = h\mu^2 - m_a$ (affine shift). This does NOT preserve the Koide condition because:

1. Koide $Q = 2/3$ is preserved by rescaling $m_a \to \alpha\,m_a$
2. Koide $Q = 2/3$ is NOT preserved by affine shifts $m_a \to m_0 - m_a$
3. The ISS CW mass depends on $m_a$ through an affine shift

Therefore: **the CW shift is flavor-dependent but pushes $Q$ toward 1/3, not 2/3**.

### Sensitivity to $\mu$

| $\mu$ (MeV) | $m_{\text{CW}}(u)$ | $m_{\text{CW}}(c)$ | $Q(u,d,s)$ | $Q(d,s,c)$ |
|-------------|---------------------|---------------------|-------------|-------------|
| 50 | 37.50 | 36.15 | 0.3333 | 0.3336 |
| 92 | 27.30 | 24.48 | 0.3333 | 0.3335 |
| 150 | 15.67 | 14.29 | 0.3333 | 0.3335 |
| 500 | 148.4 | 147.8 | 0.3333 | 0.3333 |
| 5000 | 1484 | 1484 | 0.3333 | 0.3333 |

For ALL values of $\mu$, the CW spectrum gives $Q \approx 1/3$. Increasing $\mu$ makes the spectrum more degenerate.

## 4. Connection to Seiberg Seesaw and Dual Koide

In the SUSY vacuum (not the ISS metastable vacuum), the Seiberg seesaw gives:
$$M_i^{\text{meson}} = \frac{\Lambda^2}{m_i^{\text{quark}}}$$

This maps the "dual Koide" condition to a standard Koide condition:
$$Q\left(\frac{1}{m_d}, \frac{1}{m_s}, \frac{1}{m_b}\right) = Q(M_d, M_s, M_b) = 0.6652$$

This is 0.22% from 2/3. The seesaw transmits dual Koide exactly (the Koide formula is scale-invariant).

However, this applies to the SUSY vacuum, not the ISS metastable vacuum where the CW potential operates.

## 5. Pseudo-Modulus and Koide Phase

The pseudo-modulus $X = \Phi^1_1$ is stabilized at $X = 0$ by the CW potential (Shih 2007).

The Koide phase $\delta$ parametrizes the mass eigenstates:
$$\sqrt{m_k} = \sqrt{M_0}\left(1 + \sqrt{2}\cos\left(\frac{2\pi k}{3} + \delta\right)\right)$$

At the CW vacuum ($X = 0$):
- SUSY is unbroken in sector A $\Rightarrow$ no CW potential for $X$ from sector A alone
- The CW potential for $X$ comes from sector B loops (mixed $X$-$\psi$ diagrams)
- The CW potential depends on $|X|$ only, not on the phase of $X$
- Therefore: **the Koide phase $\delta$ is a flat direction of the CW potential**

This confirms the finding in koide_phases.md: a $Z_9$-breaking term in the potential (beyond one-loop CW) is needed to select $\delta$.

## 6. Summary

### Main results

1. **Tree-level spectrum**: The ISS metastable vacuum has
   - 1 Goldstino (massless)
   - Heavy SUSY pairs at mass $\sim h\phi_0 \approx 160$ MeV
   - SUSY-breaking splitting $\pm h\,F_a$ in the link sector
   - Broken mesons $\chi_a$: tree-level FLAT (massless)

2. **CW masses are nearly degenerate**: $m_{\text{CW}} \approx 27$ MeV for all flavors (within 10%), because all quark masses are small compared to $h\mu^2 = 14660$ MeV$^2$.

3. **CW does NOT transmit Koide**: The CW spectrum gives $Q \approx 1/3$ (degenerate limit), regardless of the UV Koide ratio. This is because the CW mass depends on $m_a$ through an affine shift $f_a = h\mu^2 - m_a$, which does not preserve Koide.

4. **Seiberg seesaw DOES transmit dual Koide**: The SUSY vacuum seesaw $M_i \propto 1/m_i$ maps dual Koide to standard Koide exactly. But this is the wrong vacuum for ISS.

5. **Koide phase is CW-flat**: The one-loop CW potential depends on $|X|$ only, so $\delta$ is not selected. A $Z_9$-breaking mechanism beyond one-loop CW is needed.

### Implications for the sBootstrap program

The ISS CW potential is NOT the mechanism for Koide transmission from UV to IR. The relevant mechanisms are:

(a) **Bion Kahler potential** (complete_kahler.md): Non-perturbative monopole-instanton effects generate a Kahler potential $K_{\text{bion}} \propto |\sum_k s_k \sqrt{m_k}|^2$ that directly involves $\sqrt{m_k}$ and implements $v_0$-doubling.

(b) **D-term immunity** (unified_superpotential.md): Koide is a fermion mass condition. D-terms shift scalar masses but not fermion masses in $N=1$ SUSY, so the Koide condition is naturally protected.

(c) **Seiberg seesaw** (SUSY vacuum): The dual Koide $Q(1/m_d, 1/m_s, 1/m_b) = 0.665 \approx 2/3$ would become a standard Koide condition for seesaw meson masses.

The CW potential's role is limited to:
- Stabilizing the pseudo-modulus $X$ at the origin
- Lifting the $\chi_a$ flat directions (giving them mass $\sim 27$ MeV)
- NOT transmitting flavor structure from UV to IR

### Numerical verification

All computations verified with:
- 5-point stencil numerical differentiation (convergence tested)
- Supertrace checks: STr$[M^0] = 0$, STr$[M^2] = 0$ (exact)
- Multiple $\mu$ values confirming the degeneracy trend
- Dual Koide = Seesaw Koide identity (verified to machine precision)
