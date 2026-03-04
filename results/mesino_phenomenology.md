# Off-Diagonal Mesino Phenomenology at the Seiberg Seesaw Vacuum

## Setup

Confined SU(3) SQCD with $N_f = N_c = 3$, flavors $(u, d, s)$. The meson matrix $M^i_j = Q^i \bar{Q}_j$ is a $3 \times 3$ chiral superfield. At the Seiberg seesaw vacuum:

$$\langle M^i_j \rangle = \delta^i_j \frac{C}{m_i}, \qquad C = \Lambda^2 (m_u m_d m_s)^{1/3} = 882297 \text{ MeV}^2$$

$$X_0 = -\frac{C}{\Lambda^6} = -1.210 \times 10^{-9} \text{ MeV}^{-4}$$

The fermion mass matrix $W_{IJ} = \partial^2 W / \partial \Phi_I \partial \Phi_J$ decomposes into independent 2x2 blocks for the off-diagonal meson pairs.

---

## 1. Quantum Numbers

### EM charge assignment

The meson $M^i_j = Q^i \bar{Q}_j$ carries EM charge $Q_{\rm em}(M^i_j) = Q_i - Q_j$, where $Q(u) = +2/3$ and $Q(d) = Q(s) = -1/3$.

| Field | Composition | $Q_{\rm em}$ | $B$ | $S$ | Mass (MeV) |
|-------|-------------|:---:|:---:|:---:|:---:|
| $M^u_d$ | $u\bar{d}$ | $+1$ | 0 | 0 | $1.143 \times 10^{-5}$ |
| $M^d_u$ | $d\bar{u}$ | $-1$ | 0 | 0 | $1.143 \times 10^{-5}$ |
| $M^u_s$ | $u\bar{s}$ | $+1$ | 0 | $+1$ | $2.287 \times 10^{-4}$ |
| $M^s_u$ | $s\bar{u}$ | $-1$ | 0 | $-1$ | $2.287 \times 10^{-4}$ |
| $M^d_s$ | $d\bar{s}$ | $0$ | 0 | $+1$ | $4.944 \times 10^{-4}$ |
| $M^s_d$ | $s\bar{d}$ | $0$ | 0 | $-1$ | $4.944 \times 10^{-4}$ |

### Dirac structure

Each off-diagonal pair $(M^i_j, M^j_i)$ forms a **Dirac fermion** through the 2x2 mass matrix:

$$\mathcal{L}_{\rm mass} = |X_0| M_k \, (\psi_{M^i_j} \psi_{M^j_i} + \text{h.c.})$$

where $k$ is the third flavor not equal to $i$ or $j$. The Weyl fermion in $M^i_j$ (charge $Q_i - Q_j$) pairs with the Weyl fermion in $M^j_i$ (charge $Q_j - Q_i$) to form one Dirac fermion of mass $|X_0| M_k$ and charge $|Q_i - Q_j|$.

| Dirac fermion | Weyl pair | $|Q_{\rm em}|$ | Mass | Mass (eV) |
|-------|-----------|:---:|--------|:---:|
| $\psi_{ud}$ | $(M^u_d, M^d_u)$ | 1 | $|X_0| M_s$ | 11.4 |
| $\psi_{us}$ | $(M^u_s, M^s_u)$ | 1 | $|X_0| M_d$ | 228.7 |
| $\psi_{ds}$ | $(M^d_s, M^s_d)$ | 0 | $|X_0| M_u$ | 494.4 |

### Comparison with SM leptons

| Property | $\psi_{ud}$ | electron | $\psi_{us}$ | muon | $\psi_{ds}$ | $\nu_e$ |
|----------|:-----------:|:--------:|:-----------:|:----:|:-----------:|:------:|
| $Q_{\rm em}$ | $\pm 1$ | $-1$ | $\pm 1$ | $-1$ | 0 | 0 |
| $B$ | 0 | 0 | 0 | 0 | 0 | 0 |
| $L$ | 0 | 1 | 0 | 1 | 0 | 1 |
| $S$ | 0 | 0 | $\pm 1$ | 0 | $\pm 1$ | 0 |

**Matches**: EM charge and baryon number match for $\psi_{ud} \leftrightarrow e$ and $\psi_{us} \leftrightarrow \mu$.

**Failures**:
1. **Lepton number**: All mesinos have $L = 0$. If mesinos are leptons, $L$ must be an accidental/emergent symmetry of the low-energy theory.
2. **Strangeness**: $\psi_{us}$ carries $S = \pm 1$ while the muon has $S = 0$. Since strangeness is violated by the weak interaction (it is not a gauge quantum number), this is not necessarily fatal in a framework where flavor quantum numbers are emergent.

---

## 2. Electron g-2 Constraint

### Vacuum polarization from light charged fermions

A new charged fermion with mass $m_f \ll m_e$ and unit EM charge contributes to $a_e$ through vacuum polarization insertions at order $\alpha^2$:

$$\Delta a_e^{\rm VP} = \left(\frac{\alpha}{\pi}\right)^2 \left[\frac{1}{3} \ln\frac{m_e}{m_f} - \frac{25}{36} + \mathcal{O}\left(\frac{m_f}{m_e}\right)\right]$$

The contribution is logarithmic (not power-law) because the VP integral is dominated by momenta $q^2 \sim m_e^2$ where the new fermion loop is unsuppressed.

### u-d mesino ($m_f = 11.4$ eV, $Q = \pm 1$)

$$\frac{m_e}{m_f} = 4.47 \times 10^4, \qquad \ln\frac{m_e}{m_f} = 10.71$$

$$\Delta a_e = \left(\frac{\alpha}{\pi}\right)^2 \times 2.88 = 1.55 \times 10^{-5}$$

**Ratio to experimental precision**: $|\Delta a_e| / \delta a_e^{\rm exp} = 5.5 \times 10^7$

**Excluded by 7 orders of magnitude.**

### u-s mesino ($m_f = 229$ eV, $Q = \pm 1$)

$$\frac{m_e}{m_f} = 2.23 \times 10^3, \qquad \ln\frac{m_e}{m_f} = 7.71$$

$$\Delta a_e = \left(\frac{\alpha}{\pi}\right)^2 \times 1.88 = 1.01 \times 10^{-5}$$

**Ratio to experimental precision**: $|\Delta a_e| / \delta a_e^{\rm exp} = 3.6 \times 10^7$

**Excluded by 7 orders of magnitude.**

### d-s mesino ($m_f = 494$ eV, $Q = 0$)

Electrically neutral. No direct VP contribution to $a_e$. Not constrained by g-2.

### Running of $\alpha$

Beyond VP insertions in g-2 diagrams, a light charged fermion shifts the effective fine structure constant at the electron mass scale:

$$\Delta\alpha(m_e) = \frac{\alpha}{3\pi} \left[2\ln\frac{m_e}{m_f} - \frac{5}{3}\right]$$

For the u-d mesino: $\Delta\alpha / \alpha = 2.10$, meaning $\alpha$ would more than double its running contribution. This would destroy the agreement between the QED prediction of $a_e$ and the value of $\alpha$ measured from atomic recoil.

---

## 3. LEP and Collider Constraints

### Cross section at the Z pole

A new light charged Dirac fermion with unit EM charge would be pair-produced at LEP via $e^+e^- \to \gamma^*/Z \to \psi\bar{\psi}$. For $m_\psi \ll M_Z$:

$$\sigma(\psi\bar{\psi}) \sim \sigma(\mu^+\mu^-) \simeq 2 \text{ nb at } \sqrt{s} = M_Z$$

This would appear as an additional unit in the ratio $R_\ell = \Gamma_{\rm had}/\Gamma_\ell$, or equivalently as one additional light lepton generation. The measurement $N_\nu = 2.984 \pm 0.008$ (LEP+SLD) excludes this at $> 100\sigma$ significance.

The precise exclusion depends on the $SU(2)_L \times U(1)_Y$ quantum numbers of the mesino. If the mesino is an $SU(2)$ singlet with hypercharge $Y = -1$ (like $e_R$), it couples to the Z via $g_V = -1 + 2\sin^2\theta_W$, $g_A = 0$, giving a reduced but still easily observable cross section.

### Direct search limits

Any stable charged particle lighter than the electron would have been discovered in:
- Beam dump experiments (copious production, long-lived)
- Positronium spectroscopy (VP corrections)
- Lamb shift measurements in hydrogen
- Millicharged particle searches

**Verdict: A charged fermion at 11.4 eV or 229 eV with unit EM charge is experimentally impossible.**

---

## 4. Resolution Paths

### (a) Identification path: mesinos = leptons

| Test | Result |
|------|--------|
| EM charge match | $\checkmark$ for $\psi_{ud} \leftrightarrow e$, $\psi_{us} \leftrightarrow \mu$ |
| Baryon number | $\checkmark$ ($B = 0$ for both) |
| Strangeness | $\times$ ($\psi_{us}$ has $S \neq 0$; muon has $S = 0$) |
| Lepton number | Requires $L$ to be emergent |
| Mass scale | $\times$ (off by factors $10^4$--$10^5$) |
| Mass ratios | $\times$ ($m_s/m_d = 20$ vs $m_\mu/m_e = 207$) |

The mass scale failure is the most severe: off-diagonal mesino masses are set by $|X_0| M_k = |C/\Lambda^6| \times C/m_k$, which gives masses in the eV range, not the MeV range of leptons.

**To rescue this identification**, one needs the Kahler metric correction:

$$m_{\rm phys} = Z_M^{1/2} \times m_{\rm hol}$$

Required corrections:
- $Z_M^{1/2}(\psi_{ud} \to e) = m_e / (|X_0| M_s) = 4.5 \times 10^4$
- $Z_M^{1/2}(\psi_{us} \to \mu) = m_\mu / (|X_0| M_d) = 4.6 \times 10^5$

These values ($Z_M \sim 10^9$--$10^{11}$) are implausibly large, even for a strongly coupled theory.

### (b) Confinement path: mesinos are not asymptotic states

In the Seiberg effective theory for $N_f = N_c$, the mesons $M^i_j$ are the fundamental low-energy degrees of freedom **below** the confinement scale $\Lambda$. They are **not** further confined -- they **are** the confined spectrum. The mesinos (fermionic partners) are equally physical propagating particles.

This path does not resolve the problem. The mesinos cannot be dismissed as "confined away."

### (c) Non-perturbative mass modification

The most plausible resolution. The off-diagonal mesino masses $|X_0| M_k$ are computed from the **holomorphic** superpotential with a **canonical** Kahler potential. In a strongly coupled confining theory, corrections arise from:

1. **Non-canonical Kahler potential**: $K = Z_{ij}(M, M^\dagger) M^{i\dagger}_j M^j_i + \ldots$ where $Z_{ij}$ encodes the anomalous dimension of the composite operator $Q^i \bar{Q}_j$.

2. **Higher-order superpotential corrections**: Terms beyond the leading $\det M$ (instanton corrections at sub-leading order in $\Lambda$).

3. **SUSY-breaking contributions**: Soft masses, D-terms, and anomaly-mediated corrections modify the fermion spectrum.

4. **Off-diagonal VEV shifts**: If the true vacuum has $\langle M^i_j \rangle \neq 0$ for $i \neq j$ (CKM-like mixing), the mesino mass matrix is substantially reorganized.

**Key point**: The holomorphic superpotential is exact (protected by holomorphy), but the physical masses depend on the Kahler potential, which is not protected and receives corrections at every order in $\Lambda / M_{\rm UV}$.

---

## 5. Are Off-Diagonal Mesinos Physical?

### Status in the Seiberg effective theory

In the $N_f = N_c$ confined SQCD theory, the mesons $M^i_j$, baryons $B, \tilde{B}$, and the Lagrange multiplier $X$ are the **low-energy degrees of freedom**. They are physical propagating fields with well-defined masses. The off-diagonal mesinos are physical states that would appear in the spectrum.

### Experimental non-existence

The charged mesinos at 11.4 eV and 229 eV are excluded by:
- **Electron g-2**: Exceeds experimental precision by 7 orders of magnitude
- **LEP**: Would show up as additional light charged fermion generation
- **Vacuum polarization**: Would double the running of $\alpha$ below $m_e$
- **Positronium**: Spectroscopy excludes new light charged particles at this level

### Resolution

**The off-diagonal mesinos at the seesaw vacuum with canonical Kahler potential are physical states in the effective theory, but their masses as computed are not the physical masses.** The gap between the holomorphic mass ($\sim$ eV) and any viable physical mass ($\gtrsim m_e$ for charged states) must be filled by Kahler potential corrections.

Three scenarios remain open:

1. **The SQCD sector confines at a scale much higher than 300 MeV**, pushing all mesino masses above experimental thresholds. This changes the numerical values of $M_k$, $X_0$, and hence the mesino spectrum.

2. **The Kahler potential provides large wave function renormalization** that lifts the off-diagonal mesino masses to the lepton mass scale. This requires $Z_M \sim 10^9$--$10^{11}$, which is possible in principle for a strongly coupled theory but has no known calculational control.

3. **The off-diagonal meson directions are lifted by additional dynamics** (gauge interactions, anomalous U(1) masses, or additional superpotential terms from UV completion) that give them masses at or above the confinement scale, removing them from the low-energy spectrum entirely.

---

## Summary Table

| Question | Answer |
|----------|--------|
| Do quantum numbers match leptons? | Partially: $Q_{\rm em}$ and $B$ match; $L$ and $S$ fail |
| Can $\psi_{ud}$ be the electron? | Not at tree level (mass off by $10^4$, g-2 excluded) |
| Can $\psi_{us}$ be the muon? | Not at tree level (mass off by $10^5$, g-2 excluded, wrong $S$) |
| g-2 constraint on charged mesinos | Excluded by factor $\sim 10^7$ |
| Are mesinos physical in Seiberg theory? | Yes -- they are the confined degrees of freedom |
| Most viable resolution | Non-perturbative Kahler corrections or higher confinement scale |
| Is direct identification viable? | Only with large non-perturbative corrections to both masses and quantum numbers |

---

## Key Conclusion

The off-diagonal mesino spectrum at the Seiberg seesaw vacuum exhibits the correct **structure** for lepton identification (two charged Dirac fermions, one neutral Dirac fermion, correct charge assignments for the charged states) but fails quantitatively on masses by 4--5 orders of magnitude. The g-2 exclusion is not merely a precision test failure; it is a **catastrophic** exclusion at 7 orders of magnitude, meaning the charged mesinos at their holomorphic masses are absolutely ruled out as physical particles coupling to photons.

The resolution must come from Kahler potential physics, not from the superpotential. The holomorphic superpotential (which determines the mesino masses at tree level) is exact, but the physical fermion masses are $m_{\rm phys} = Z^{-1/2} W_{IJ}$, where $Z$ is the Kahler metric. In a strongly coupled confined theory, $Z$ can differ from unity by orders of magnitude, and its computation requires non-perturbative techniques (lattice, functional methods, or holographic duality) that go beyond the Seiberg effective action.

*Generated from mesino_phenomenology analysis, March 2026.*
