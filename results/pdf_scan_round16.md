# PDF Scan — Round 16
Date: 2026-03-04

Three PDFs scanned from /home/codexssh/phys3/sources/

---

## PDF 1: volkov_1973.pdf

**Title**: Higgs Effect for Goldstone Particles with Spin 1/2
**Authors**: D.V. Volkov and V.A. Soroka
**Journal**: JETP Letters 18 (1973), 529
**Institution**: Institute of Physics & Technology, Kharkov, USSR

### Summary

This is NOT the original Volkov-Akulov 1973 paper on nonlinear SUSY realization. It is the companion 1973 letter by Volkov and Soroka — the paper that first wrote down supergravity by gauging the Volkov-Akulov symmetry. The Volkov-Akulov nonlinear SUSY paper itself is cited as a preprint (ITF-73-51R) and a 1972 JETP Lett. letter.

The paper analyzes what happens when the VA symmetry is gauged. The symmetry group G is the direct product of the Poincare group with an internal symmetry group, extended by the nonlinear fermionic shifts:

```
psi -> psi' = psi + xi
x_mu -> x'_mu = x_mu - (a/2i)(xi+ sigma_mu psi - psi+ sigma_mu xi)
```

These are the VA transformations (eq. 1), reproduced from Volkov-Akulov (1973).

### Key Equations

**Eq. (1)**: The nonlinear SUSY transformation (VA transformations) — shifts in spinor space with x-dependent compensation.

**Eq. (4)**: The covariant differential forms (supervielbeins), including the gravitino field Phi(d) for spin 3/2:
```
omega_mu(d) = Dx_mu + (sigma/2i)[psi+ sigma_mu(Dpsi + 2f Phi(d)) - (Dpsi + 2f Phi(d))+ sigma_mu psi] + f W_mu(d)
omega^alpha_a(d) = (Dpsi + f Phi(d))^alpha_a
```

**Eq. (5)**: Curvature 2-forms for the Lorentz connection and internal gauge field.

**Eqs. (6)-(9)**: Invariant action combinations (in the language of exterior products):
- Eq. (6): Kinetic term for Goldstone (spin 1/2) particles
- Eq. (7): Kinetic term for spin 3/2 gauge field (gravitino)
- Eq. (8): Einstein gravity (spin 2)
- Eq. (9): Yang-Mills (spin 1)

**Eq. (10)**: New gravitino gauge form absorbing the goldstino:
```
f Phi'(d) = Dpsi + f Phi(d)
```

**Eq. (11)**: Invariant constraint omega^alpha(d) = 0 that eliminates the spin 3/2 field.

**Eq. (12)**: The key Higgs mechanism for the goldstino:
```
f W'(d) = omega_mu(d)
```
This redefinition of the metric/vierbein gauge form eliminates the Goldstone fermion — the goldstino is "eaten" by the graviton to give a massive gravitino. This is supergravity.

### Key Results

1. **Goldstino eaten by graviton**: The fermionic Goldstone particle (goldstino) is absorbed into the metric tensor gauge field when the VA symmetry is gauged, not into a spin 3/2 field. The graviton acquires a spin 3/2 superpartner with mass — this is the first supergravity.

2. **Spin 3/2 can be independently eliminated**: Imposing the invariant constraint omega^alpha(d) = 0 removes the gravitino from the spectrum, leaving only the graviton and Yang-Mills fields.

3. **Goldstone fermion retention requires gauge-group violation**: A massless goldstino can persist only if the gauge invariance is broken. This is the key obstruction to using VA goldstinos for weak interaction phenomenology.

4. **Cosmological term**: The invariant (6) acts as a cosmological constant term in the gravitational action.

### Connections to sBootstrap Program

- **Nonlinear SUSY / goldstino**: The VA model is the canonical example of SUSY broken with a massless Nambu-Goldstone fermion. In the sBootstrap context, this is the relevant framework if SM quarks are to be interpreted as goldstinos of a broken SUSY.
- **Higgs mechanism for fermions**: The mechanism here (goldstino eaten by a gauge boson of different quantum numbers) is structurally analogous to what would be needed if neutrinos or light quarks are to acquire mass from a "foreign" symmetry breaking.
- **No direct mass formula**: This paper does not contain Koide-type mass relations or any QCD/EW scale connection. It is purely about the kinematic structure of nonlinear SUSY gauging.
- **No hadronic SUSY**: No connection to Miyazawa or hadronic supersymmetry. The model is gravitational/supergravity context only.

---

## PDF 2: article_23569.pdf

**Contents**: This PDF is a journal proceedings page containing TWO separate short papers printed consecutively, plus the tail-end of a third paper's reference list. The actual content is:

### Paper 2a (pages 261-262 of the journal)

**Title**: Interaction of Goldstone Neutrino with Electromagnetic Field
**Authors**: V.P. Akulov and D.V. Volkov
**Institution**: Physico-Technical Institute, Ukrainian Academy of Sciences
**Published**: ZhETF Pis. Red. 17, No. 7, 367-369 (5 April 1973)
**Submitted**: 13 February 1973

### Summary of Paper 2a

Sequel to the original VA 1972 paper. The neutrino is modeled as a Goldstone particle of the expanded Poincare group (VA group). The paper derives the coupling of this Goldstone neutrino to the electromagnetic field and obtains observational bounds.

Starting from the VA action integral invariant under:
```
psi -> psi' = psi + zeta
x_mu -> x'_mu = x_mu - (a/2i)(zeta+ sigma_mu psi - psi+ sigma_mu zeta)
```
with single coupling constant a = l^4 (l has dimension of length).

### Key Equations (Paper 2a)

**Eq. (1)**: Invariant action in differential form language (same as original VA action but now including electromagnetic potential A_mu(x)):
```
S = integral [ D(A_mu dx^mu) wedge omega^nu wedge omega^rho D(A_a dx^a) wedge omega^nu' wedge omega^rho' / (epsilon_{mu nu rho lambda} omega^mu wedge omega^nu wedge omega^rho wedge omega^lambda) ] g_{nu nu'} g_{rho rho'}
```

**Eq. (2)**: Expanded action in conventional notation:
```
S = integral [ -(1/4) F_{mu nu} F^{mu nu} + a F_{mu nu} F^{mu rho} T^nu_rho + ... ] d^4x
```
where T^mu_nu = (i/2)(psi+ sigma^mu d_nu psi - d_nu psi+ sigma^mu psi) is the neutrino stress-energy tensor.

**Eq. (3)**: Differential cross section for the forbidden process gamma + gamma -> nu + nu-bar:
```
d sigma/dt = (a^2 / 8 pi 16) [ st + 3t^2 + 4 t^3/s + 2 t^4/s^2 ]
```

**Eq. (4)**: Total cross section:
```
sigma(s) = (a^2 / 80 pi) s^3
```
This grows as s^3 (energy cubed) — characteristic of a dimension-8 operator interaction.

**Eq. (5)**: Solar neutrino luminosity from Goldstone-neutrino pair emission:
```
Q_nu = (l^8 c^2 hbar / rho) * [Gamma(7) Gamma(6) zeta(7) zeta(6) 2^8 / 400 pi^5] * (T/hbar c)^13
```
evaluated at solar center conditions (rho ~ 140 g/cm^3, T = 1.3 keV).

**Eq. (6)**: Upper bound on the coupling length from solar energy loss:
```
l < 10^{-12} cm
```

**Eq. (7)**: Differential cross section for nu_mu p -> nu_mu p elastic scattering (from the VA Lagrangian for neutrino-fermion interaction):
```
d sigma(s,t)/dt = (a^2 / 16 pi) [(s-m^2)^2 + (s-2m^2)t + s/(s-m^2) t^2]
```

From CERN 1967 neutrino data on elastic scattering (nu_mu p -> nu_mu p), the bound tightens to:
```
l < 10^{-15} cm
```
which is close to the weak interaction length l_w = sqrt(G_F) ~ 0.66 x 10^{-16} cm.

### Key Results (Paper 2a)

1. **Process gamma + gamma -> nu + nu-bar is allowed** in the Goldstone-neutrino model (forbidden in standard V-A weak interaction theory).
2. **Cross section grows as s^3** (dimension-8 interaction), not s^1 (as in Fermi theory).
3. **Coupling constant bounded**: l < 10^{-15} cm from CERN neutrino data, approaching the weak scale.
4. **Universal interaction**: The same constant a = l^4 governs all neutrino interactions — neutrino self-coupling, neutrino-fermion, neutrino-photon.

### Paper 2b (page 263 of the journal)

**Title**: Influence of Low-Frequency Oscillations on Transport Processes Across a Magnetic Field
**Authors**: L.M. Kovrizhnykh
**Institution**: P.N. Lebedev Physics Institute, USSR Academy of Sciences
**Published**: ZhETF Pis. Red. 17, No. 7, 369-373 (5 April 1973)

This paper is plasma physics (anomalous transport in magnetized plasma due to drift waves) and is completely unrelated to SUSY or particle physics. It analyzes transport coefficients in a plasma cylinder with a potential wave of the form Phi = Phi_0 exp(i k_perp eta + i k_parallel zeta - i omega t). Not relevant to any sBootstrap topic. Documented here for completeness only.

### Connections to sBootstrap Program (Paper 2a)

- **Goldstone neutrino as massless fermion**: The VA model treats the neutrino as a massless Goldstone particle. This is the original motivation for identifying massless fermions with goldstinos of broken SUSY.
- **Universal coupling constant a = l^4**: The single parameter controlling all interactions has dimension [length]^4, related to the SUSY-breaking scale F by a = 1/F^2. The bound l < 10^{-15} cm gives F > (10^{-15} cm)^{-1/2} which is near the electroweak scale.
- **QCD/EW connection**: The bound l < 10^{-15} cm ~ l_w suggests that if this mechanism is operative, the SUSY-breaking scale is naturally the electroweak scale. This is a qualitative link between the Goldstone-fermion scale and the weak interaction scale.
- **No mass formula, no Koide connection, no hadronic SUSY**.
- **Dimension-8 operators**: The interaction terms F_{mu nu} F^{mu rho} T^nu_rho are dimension-8. These are the characteristic contact interactions of a nonlinear SUSY theory at low energy, suppressed by 1/F^2.

---

## PDF 3: Possible_universal_neutrino_interaction.pdf

**Title**: On the Possible Universal Neutrino Interaction
**Authors**: D.V. Volkov, V.P. Akulov
**Institution**: Physico-Technical Institute, Academy of Sciences of the Ukrainian SSR (Kharkov)
**Published**: Pis'ma ZhETF (JETP Lett.) 16, No. 11, pp. 621-624 (1972)
**Reprinted in**: Ukr. J. Phys. 2008, V. 53, Special Issue (with memorial text for Volkov)

This is the ORIGINAL Volkov-Akulov paper — the founding document of nonlinear supersymmetry and the goldstino. The PDF is a 2008 Ukrainian Journal of Physics reprint with a biographical memoir of Volkov appended.

### Summary

The paper proposes that the neutrino is the Goldstone particle of a spontaneously broken symmetry of the vacuum. The broken symmetry is an extension of the Poincare group by fermionic (anticommuting) shifts of both the spinor field and spacetime coordinates.

The free neutrino equation:
```
i sigma^mu d/dx^mu psi = 0   (Eq. 1)
```
is invariant under Poincare, chiral transformations, and constant spinor shifts psi -> psi + zeta (with x -> x). The authors replace this with the nonlinear transformation:

```
psi -> psi' = psi + zeta
psi+ -> psi+' = psi+ + zeta+
x_mu -> x'_mu = x_mu - (a/2i)(zeta+ sigma_mu psi - psi+ sigma_mu zeta)   (Eq. 3)
```

This is the VA group: 10 commuting (Poincare) + 4 anticommuting parameters. They note (footnote 2) this is a Lie group with anticommuting parameters, citing Berezin-Kats (1970) — one of the earliest uses of superspace language.

### Key Equations

**Eq. (1)**: Free massless neutrino (Weyl) equation.

**Eq. (2)**: Linear spinor shifts (abandoned in favor of eq. 3).

**Eq. (3)**: The VA nonlinear transformation — the defining equation of the Volkov-Akulov model. The constant a has dimension [length]^4.

**Eq. (4)**: The invariant 1-form (supervielbein):
```
omega_mu = dx_mu + (a/2i)(psi+ sigma_mu dpsi - dpsi+ sigma_mu psi)
```

**Eq. (5)**: Invariant action as exterior product of four omega forms:
```
S = (1/a) integral omega_0 wedge omega_1 wedge omega_2 wedge omega_3
```

**Eq. (6)**: Action in component form:
```
S = (1/sigma) integral |W| d^4x
```
where W_{mu nu} = delta_{mu nu} + sigma T_{mu nu} and sigma = a.

**Eq. (7)**: The stress-energy tensor entering the action:
```
T_{mu nu} = (1/2i)(psi+ sigma_mu d_nu psi - d_nu psi+ sigma_mu psi)
```

**Eq. (8)**: Full expanded action as series in T:
```
S = integral [ 1/sigma + T_{mu mu} + (a/2)(T_{mu mu} T_{nu nu} - T_{mu nu} T_{nu mu})
    + (a^2/3!) sum_p (-1)^p T_{mu mu} T_{nu nu} T_{rho rho}
    + (a^3/4!) sum_p (-1) T_{mu mu} T_{nu nu} T_{rho rho} T_{sigma sigma} ] d^4x
```
The term T_{mu mu} is the kinetic term; products of 2, 3, 4 tensors T give 4-, 6-, 8-fermion interactions respectively.

**Eq. (9)**: Neutrino-Dirac particle interaction action:
```
S = integral [ R_{mu mu} + a(R_{mu mu} T_{nu nu} - R_{mu nu} T_{nu mu})
    + (a^2/2) sum_p (-1)^p R_{mu mu} T_{nu nu} T_{rho rho}
    + (a^3/3!) sum_p (-1)^p R_{mu mu} T_{nu nu} T_{rho rho} T_{sigma sigma}
    + m phi-bar phi |W| ] d^4x
```
where R_{mu nu} = (1/2i)(phi-bar gamma_mu d_nu phi - d_nu phi-bar gamma_mu phi).

**Eq. (10)**: Dirac bilinear current tensor.

### Key Results

1. **First nonlinear SUSY realization**: The VA transformation (eq. 3) is the first explicit nonlinear realization of supersymmetry. The algebra closes only on-shell (for massless spinors) and generates a graded extension of the Poincare group.

2. **Single coupling constant**: All neutrino interactions (self-coupling, coupling to any fermion, coupling to any gauge field) are determined by the single dimensionful constant a = l^4. The interaction is universal in the sense that the coupling is geometrically fixed by the invariance requirement.

3. **Goldstone interpretation of zero-mass**: The masslessness of the neutrino is attributed to its being the Goldstone particle of the broken fermionic shift symmetry of the vacuum.

4. **Derivative interactions only**: The interaction Lagrangian contains only derivative couplings (through T_{mu nu}), so all contact interactions are suppressed by powers of (momentum/F) where F = a^{-1/2} is the SUSY-breaking scale.

5. **Weak/EM switching on via Higgs mechanism**: The paper notes that weak and electromagnetic interactions can be incorporated by introducing gauge fields for the approximate unitary symmetry group; the Goldstone-neutrino then vanishes (eaten by the spin 3/2 gauge boson) and a massive gravitino appears — anticipating the content of the companion paper (Volkov-Soroka 1973, PDF 1 above).

### Appended Biographical Text

The reprint includes a biographical memoir of Volkov (1925-1996), listing his main contributions:
- Theory of Regge poles (with Gribov, 1963)
- Goldstone neutrino / nonlinear SUSY (with Akulov, 1972-73)
- Phenomenological Lagrangian methods
- First supergravity (with Soroka, 1973)
- Twistor approach to superstrings and superbranes

### Connections to sBootstrap Program

- **This is the founding VA paper**: The nonlinear SUSY transformation (eq. 3) is the starting point for any treatment of goldstinos, including the scenario where SM quarks/leptons are goldstinos of spontaneously broken SUSY.

- **SUSY-breaking scale and electroweak scale**: The constraint from CERN neutrino data (in PDF 2) gives l < 10^{-15} cm, meaning F > ~ 100 GeV — directly at the electroweak scale. This is consistent with the sBootstrap identification of the SUSY-breaking VEV v_0 with the electroweak scale.

- **O'Raifeartaigh-Koide connection**: The O'Raifeartaigh model (used in the sBootstrap to generate the Koide seed) achieves F-term SUSY breaking, and its vacuum is described by a nonlinear SUSY in the low-energy limit — the goldstino from O'Raifeartaigh is exactly the VA goldstino in that effective theory. The coupling gv/m = sqrt(3) corresponds to a specific value of the VA coupling constant a.

- **Composite scalars / hadronic SUSY**: The VA model gives massless goldstinos but says nothing directly about hadronic SUSY or composite bosons. The connection to Miyazawa-type hadronic SUSY requires additional structure (e.g., identifying mesons as goldstino bilinears or as the scalar superpartners — which is the sBootstrap claim).

- **No mass formula**: The VA model is kinematic; no mass formula, no Koide relation, no QCD scale emerges from the original papers. The mass scale is put in by hand through the constant a.

- **Neutrino as goldstino**: The original motivation was to explain massless neutrinos, not quarks. Applying VA to quarks requires massive goldstinos (soft SUSY breaking) — which is precisely the sBootstrap setup where quark masses arise from Koide-constrained SUSY breaking.

---

## Cross-Paper Summary: Relation Among the Three Documents

All three papers are by the Kharkov group (Volkov, Akulov, Soroka) from 1972-1973 and form a coherent sequence:

1. **Volkov-Akulov 1972** (PDF 3): Original VA paper. Neutrino = Goldstone fermion. Nonlinear SUSY transformation. Universal interaction. Single coupling a = l^4.

2. **Akulov-Volkov 1973** (PDF 2, paper 2a): Application of VA to neutrino-EM coupling. Bounds on coupling constant; l_w = sqrt(G_F) emerges as natural scale.

3. **Volkov-Soroka 1973** (PDF 1): Gauging of VA symmetry. Goldstino eaten by graviton -> massive gravitino. First supergravity. Residual Yang-Mills fields.

The sequence thus goes: goldstino -> phenomenological bounds (EW scale) -> local SUSY (supergravity). The sBootstrap program enters at step 1, asking whether quarks (not neutrinos) play the role of VA goldstinos, with the SUSY-breaking scale at the electroweak scale and mass generation via Koide-structured O'Raifeartaigh breaking rather than pure VA kinematics.

### Items to Flag

- **No Koide formula** in any of the three papers.
- **No hadronic SUSY / Miyazawa** connection in any of the three papers.
- **No seesaw mechanism** in any of the three papers.
- **EW scale appears naturally**: The bound l < 10^{-15} cm in PDF 2 gives the VA coupling scale at ~ G_F^{1/2}, which is the Fermi scale. This is a genuine QCD/EW scale connection in the sense that the Goldstone-fermion interpretation of the neutrino, when constrained by neutrino scattering data, forces the SUSY-breaking scale to be at the electroweak scale.
- **Dimension-8 operators**: The leading VA interaction is dimension 8 (F_{mu nu} F^{mu rho} T^nu_rho), distinguishing VA neutrinos from ordinary neutrino interactions. This is a testable prediction.
- **Spin 3/2 field**: The gauged VA theory requires a spin 3/2 gauge field (gravitino). The invariant constraint omega^alpha(d) = 0 can eliminate it. If instead the gravitino is retained, it acquires mass ~ F/M_Pl, which is a soft mass formula involving the SUSY-breaking scale.
