# Complete SUSY Lagrangian: SQCD N_f = N_c = 3 with NMSSM Higgs

## Model Specification

The model is SQCD with N_f = N_c = 3 in the confining phase, with the Seiberg-dual
description of the meson moduli space, coupled to the MSSM Higgs sector via an NMSSM
singlet. The total field content is 16 chiral superfields.

**Chiral superfield content:**

| Field | Type | Dimension | Description |
|-------|------|-----------|-------------|
| M^i_j (i,j = u,d,s) | 9 meson superfields | [mass]^2 | Gauge-invariant composites Q̃^i Q_j / Λ |
| X | 1 singlet | [mass]^{-4} | Seiberg Lagrange multiplier / NMSSM singlet |
| B, B̃ | 2 baryons | [mass]^3 | ε_{ijk} Q^i Q^j Q^k / Λ^3 (and conjugate) |
| H_u = (H_u^+, H_u^0) | 2-component doublet | [mass]^1 | Up-type Higgs superfield |
| H_d = (H_d^0, H_d^-) | 2-component doublet | [mass]^1 | Down-type Higgs superfield |

Total: 9 + 1 + 2 + 4 = 16 chiral superfields.

---

## 1. The Complete Lagrangian

The N=1 SUSY Lagrangian in component notation is:

```
L = K_{IJ̄}(∂_μ Φ_I)(∂^μ Φ̄_J) + i K_{IJ̄} ψ̄_I σ̄^μ ∂_μ ψ_J
    - V_F - V_D - V_soft
    - (1/2) W_{IJ} ψ_I ψ_J + h.c.
    + L_gauge
```

where:
- `K_{IJ̄} = ∂^2 K / ∂Φ_I ∂Φ̄_J` is the Kähler metric (kinetic mixing)
- `V_F = K^{IJ̄} F_I F̄_J` with `F_I = -∂W/∂Φ_I` (F-term potential)
- `V_D` is the D-term potential (gauge sector)
- `V_soft` is the explicit SUSY breaking
- `W_{IJ} = ∂^2 W / ∂Φ_I ∂Φ_J` gives fermion masses

---

## 2. Kähler Potential

```
K = K_tree + K_X + K_bion
```

### 2a. Tree-level (canonical) Kähler

```
K_tree = Tr(M^† M) + |X|^2 + |B|^2 + |B̃|^2 + |H_u|^2 + |H_d|^2
```

**STATUS: DETERMINED** — canonical kinetic terms required by N=1 SUSY for
gauge-invariant composites in the Seiberg dual description. No free parameters.

In expanded form:
```
Tr(M^† M) = Σ_{i,j=u,d,s} |M^i_j|^2
|H_u|^2 = |H_u^+|^2 + |H_u^0|^2
|H_d|^2 = |H_d^0|^2 + |H_d^-|^2
```

### 2b. X pseudo-modulus Kähler correction

```
K_X = - |X|^4 / (12 μ^2)
```

where μ is the Kähler mass scale, determined by self-consistency of the vacuum:

```
μ = |C / (√3 Λ^6)|    with C = Λ^2 (m_u m_d m_s)^{1/3}
```

**STATUS: DETERMINED** by the requirement that the Kähler pole at |X| = √3 μ
coincides with the seesaw vacuum |X|_vac = C/Λ^6. This yields gv/m = √3 (the
O'Raifeartaigh–Koide condition) without free parameters.

Numerical value: μ = C/(√3 Λ^6) = 882297 / (1.732 × 7.29×10^{14}) = 6.99×10^{-10} MeV^{-4}.

**NOTE**: The task statement writes K_bion with c = -1/12 and p = 1.
That form, K_bion = c |Tr(M^†M)| / Λ^{2p-2}, is NOT the bion Kähler potential.
The -1/12 coefficient applies to the X quartic, K_X = -|X|^4/(12μ^2),
not to a power of Tr(M^†M). See Section 2c below.

### 2c. Bion (monopole-instanton) Kähler correction

```
K_bion = (ζ^2 / Λ^2) |Σ_{k=u,d,s} s_k √(M^k_k)|^2
```

where:
- `ζ = exp(-S_0 / (2 N_c)) = exp(-π/(3 α_s))` is the monopole fugacity
- `s_k = ±1` are signs from monopole-instanton zero-mode structure
- `S_0 = 2π/(α_s N_c)` is the monopole-instanton action

The sign assignments are:
- Seed configuration (0, m_s, m_c): s_u = 0 (or absent), s_d = +1, s_s = +1
- Bloom configuration (−m_s, m_c, m_b): s_s = −1, s_c = +1, s_b = +1

The sign s_k flips when the corresponding eigenvalue of M crosses zero during the
seed-to-bloom transition; the phase arises from the fermion zero-mode structure of
the third monopole-instanton activating at the boundary of moduli space.

**STATUS: DETERMINED** in functional form by the semiclassical bion analysis on
R^3 × S^1. The coefficient ζ^2/Λ^2 is determined by α_s at the confinement scale.
At α_s = 1: ζ^2/Λ^2 = exp(−2S_0/3) ≈ 1.52×10^{-2}.

**LIMITATION**: The matching between the compactified R^3×S^1 bion result and the
four-dimensional confining phase is not derived from first principles. The functional
form is established (bion_kahler.md); the precise coefficient in 4D requires
additional justification.

Explicitly for diagonal M with bloom signs (s_s = −1, s_c = +1, s_b = +1):
```
K_bion = (ζ^2/Λ^2) [ M_s + M_c + M_b
                      − 2√(M_s M_c) − 2√(M_s M_b) + 2√(M_c M_b) ]
```

The cross-terms √(M_i M_j) are the physical bion contributions (correlated
monopole–antimonopole pairs on adjacent roots of the extended SU(3) Dynkin diagram).

**Suppression factor at α_s = 0.3, N_f = 5 physical flavors**:
```
ζ^2/Λ^2 × [Π_f (m_f/Λ)]^2 = 1.33×10^{-5}
```
The light-quark mass suppression prevents c_bion from reaching -1/12 with physical
masses; the bion mechanism operates on the magnetic dual meson masses (all of order Λ)
where the suppression is removed.

### 2d. Kähler metric

The Kähler metric K_{IJ̄} = ∂^2 K / ∂Φ_I ∂Φ̄_J has the following nonzero components
in the diagonal meson sector (bloom configuration):

```
K_{M^k_k, M^k_k} = 1 + (ζ^2/Λ^2) s_k^2 / (4 M_k)     [diagonal, = 1 + ε/(4M_k)]
K_{M^i_i, M^j_j} = (ζ^2/Λ^2) s_i s_j / (4√(M_i M_j))  [off-diagonal, i ≠ j]
K_{X,X̄}          = 1 − |X|^2/(3μ^2)
K_{B,B̄}          = K_{B̃,B̃̄} = K_{Hu,Hū} = K_{Hd,Hd̄} = 1
```

Off-diagonal meson pairs M^i_j (i≠j) and all other cross-terms are zero.
Fractional corrections at bloom vacuum: O(ζ^2) ≈ 1.5×10^{-2} for diagonals, smaller off-diagonal.

---

## 3. Superpotential (complete)

In LaTeX form suitable for a paper:

```latex
W = \underbrace{\sum_{i=u,d,s} m_i \, M^i_{\;i}
    + X\!\left(\det M - B\tilde{B} - \Lambda^6\right)}_{W_{\text{Seiberg}}}
  + \underbrace{y_c \, H_u^0 \, M^d_{\;d}
    + y_b \, H_d^0 \, M^s_{\;s}
    + y_t \, H_u^0 \, t \bar{t}}_{W_{\text{Yukawa}}}
  + \underbrace{\lambda \, X \!\left(H_u^0 H_d^0 - H_u^+ H_d^-\right)}_{W_{\text{NMSSM}}}
  + \underbrace{c_3 \frac{(\det M)^3}{\Lambda^{18}}}_{W_{\text{instanton}}}
```

With explicit SU(2)_L × U(1)_Y structure:
```latex
W_{\text{NMSSM}} = \lambda \, X \, (H_u \cdot H_d)
  = \lambda \, X \, \epsilon^{\alpha\beta} H_{u,\alpha} H_{d,\beta}
  = \lambda \, X \, (H_u^0 H_d^0 - H_u^+ H_d^-)
```

### Term-by-term status

**W_Seiberg — DETERMINED** by Seiberg duality for N_f = N_c = 3 SQCD in the
confining phase. This is exact: the Seiberg constraint det M − BB̃ − Λ^6 = 0
is the defining algebraic variety of the moduli space. The quark mass matrix
Σ_i m_i M^i_i is required by holomorphy and the explicit breaking of the
SU(N_f) × SU(N_f) chiral symmetry. No free parameters beyond the input masses m_i.

**W_Yukawa (charm, bottom) — PARTIALLY DETERMINED**:
- y_c = 2 m_c / v = 1.032×10^{-2}: FITTED to m_c = 1270 MeV and v = 246220 MeV
- y_b = 2 m_b / v = 3.395×10^{-2}: FITTED to m_b = 4180 MeV and v = 246220 MeV
- The factor 2 (not √2) arises from the convention H_u^0 → v/√2 so that the
  mass is m_f = y_f × (v/√2) × √2 = y_f × v, consistent with y = 2m/v.
- The specific flavor assignments (charm couples to M^d_d, bottom to M^s_s)
  are ASSUMED based on the identification of ISS seesaw VEVs M^d_d ~ 1/m_d
  and M^s_s ~ 1/m_s. This identification is not derived from first principles.

**W_Yukawa (top) — ASSUMED**:
- y_t = √2 m_t / v = 0.9934: standard top Yukawa. Conventional factor reflects
  that the top couples as a weak singlet, not through a composite meson.
- The top quark is EXTERNAL to the SQCD sector (non-composite). The term y_t H_u^0 t t̄
  treats t and t̄ as elementary chiral superfields not in the meson matrix M.
- How the top couples to the meson sector: it does not appear in M^i_j directly.
  The SU(5) flavor symmetry (if present) would need to embed t within a larger
  representation, but this embedding is NOT specified in this Lagrangian.
  **This is a gap** (see Section 6).

**W_NMSSM — ASSUMED** (λ fitted):
- λ = 0.72: FITTED to reproduce m_h = 125.25 GeV at tree level with tan β = 1
- The identification X = S (Seiberg Lagrange multiplier = NMSSM singlet) is
  the central structural assumption of the model. It is self-consistent (both
  fields are SU(2)_L × U(1)_Y singlets with R-charge 2) but not derived from
  a UV completion.
- λ = √2 m_h / v = √2 × 125250/246220 = 0.7186: the formula is determined,
  the value is fitted to the Higgs mass.

**W_instanton — ASSUMED** (c_3 undetermined):
- This is a multi-instanton correction to the Seiberg superpotential.
  For N_f = N_c = 3, the standard one-instanton contribution is already the
  ADS superpotential. The "three-instanton" term (det M)^3/Λ^{18} would be
  a higher-order correction.
- c_3 = O(1): undetermined. Its effect is to shift the vacuum value of X by
  3 c_3 / Λ^6 without changing the seesaw structure for M_i. It does not
  affect any of the mass predictions.
- **WARNING**: The inclusion of this term is not standard in SQCD. It should
  be derived from an instanton calculation or removed if not needed.

### Parameter values (explicit)

| Parameter | Value | Determination |
|-----------|-------|---------------|
| m_u | 2.16 MeV | FREE (PDG input) |
| m_d | 4.67 MeV | FREE (PDG input, but constrained by Oakes: tan θ_C = √(m_d/m_s)) |
| m_s | 93.4 MeV | FREE (overall scale, PDG input) |
| Λ | 300 MeV | ASSUMED (= f_π × √(4π) ≈ SQCD confinement scale) |
| y_c | 1.032×10^{-2} | FITTED from m_c = 1270 MeV |
| y_b | 3.395×10^{-2} | FITTED from m_b = 4180 MeV |
| y_t | 0.9934 | ASSUMED (SM top Yukawa) |
| λ | 0.72 | FITTED to m_h = 125.25 GeV |
| c_3 | O(1) | UNDETERMINED |

---

## 4. F-term Potential

The F-term potential is:
```
V_F = K^{IJ̄} F_I F̄_J
```

The F-terms ∂W/∂Φ_I are:

```
F_X      = det M − BB̃ − Λ^6
F_{M^i_i} = m_i + X · (det M)(M^{-1})^i_i
            + 3c_3 (det M)^2 (det M)(M^{-1})^i_i / Λ^{18}
            + δ_{i,d} y_c H_u^0
            + δ_{i,s} y_b H_d^0
F_{M^i_j} = X · ∂(det M)/∂M^i_j     (for i ≠ j)
F_B      = −X B̃
F_{B̃}    = −X B
F_{H_u^0} = y_c M^d_d + λ X H_d^0
F_{H_d^0} = y_b M^s_s + λ X H_u^0
F_{H_u^+} = −λ X H_d^-
F_{H_d^-} = −λ X H_u^+
```

At the diagonal seesaw vacuum (M_i = C/m_i, B = B̃ = 0, H_u^0 = H_d^0 = v/√2):

```
F_X           = 0                    [Seiberg constraint satisfied]
F_{M^u_u}     = m_u − m_u = 0       [seesaw]
F_{M^d_d}     = m_d − m_d + y_c v/√2 = y_c v/√2 = m_c  [Yukawa breaks SUSY]
F_{M^s_s}     = m_s − m_s + y_b v/√2 = y_b v/√2 = m_b  [Yukawa breaks SUSY]
F_{H_u^0}     = y_c M^d_d + λ X H_d^0 = 2C/v + λ X v/√2
F_{H_d^0}     = y_b M^s_s + λ X H_u^0 = 2C/v + λ X v/√2
```

Key identity (flavor-universal F-term cancellation, exact algebraically):
```
y_j ⟨M^j_j⟩ = (2m_j/v)(C/m_j) = 2C/v    [independent of flavor j]
```

The SUSY-breaking scale is:
```
F = C/v = Λ^2(m_u m_d m_s)^{1/3} / v = 882297 / 246220 MeV = 3.58 MeV
```

This is an extremely low SUSY-breaking scale (3.58 MeV). Whether this is compatible
with electroweak symmetry breaking at v = 246 GeV requires the large NMSSM λ to
provide the tree-level Higgs quartic (see Section 5).

### Inverse Kähler metric (used in V_F)

```
K^{M^k_k, M^k_k} = 1 − (ζ^2/Λ^2) s_k^2/(4M_k) + O(ζ^4)
K^{M^i_i, M^j_j} = −(ζ^2/Λ^2) s_i s_j/(4√(M_i M_j)) + O(ζ^4)  [i ≠ j]
K^{X,X̄}          = 1/(1 − |X|^2/(3μ^2))
K^{B,B̄}          = K^{B̃,B̃̄} = K^{Hu,Hū} = K^{Hd,Hd̄} = 1
```

Bion corrections to V_F are O(ζ^2) ≈ 1.5×10^{-2} and perturbatively small.
The X metric correction diverges at the Kähler pole |X| = √3 μ.

---

## 5. D-term Potential

```
V_D = (g^2/2)(H_u^† T^a H_u + H_d^† T^a H_d)^2
    + (g'^2/2)(Y_{H_u}|H_u|^2 + Y_{H_d}|H_d|^2 + Σ_{ij} Q_{EM}(M^i_j)|M^i_j|^2)^2
```

where T^a are SU(2)_L generators and Y is the hypercharge.

**Hypercharge assignments**:
```
Y(H_u) = +1/2,  Y(H_d) = −1/2
Q_EM(M^i_j) = Q(i) − Q(j)    with Q(u) = 2/3, Q(d) = Q(s) = −1/3
```

Nonzero EM charges of meson fields:
```
Q_EM(M^u_d) = Q_EM(M^u_s) = +1  (6 charged mesons, pairs)
Q_EM(M^d_u) = Q_EM(M^s_u) = −1  (6 charged mesons, conjugates)
```
The 5 diagonal plus off-diagonal (M^d_s, M^s_d) entries are electrically neutral.

**D-flat direction at tan β = 1**:
At H_u^0 = H_d^0 = v/√2, the SU(2)_L × U(1)_Y D-terms for the neutral Higgs vanish:
```
D^a = 0,  D' = 0    at |H_u^0|^2 = |H_d^0^2|
```
The quartic (g^2 + g'^2)/8 × (|H_u^0|^2 − |H_d^0|^2)^2 vanishes on the D-flat direction.
This is why the NMSSM singlet λ X H_u · H_d is required for m_h ≠ 0 at tree level.

**STATUS: DETERMINED** — D-term structure follows from SU(2)_L × U(1)_Y gauge invariance
with standard SM hypercharge assignments. No free parameters.

---

## 6. Soft SUSY-Breaking Potential

```
V_soft = f_π^2 Tr(M^† M)
```

where f_π = 92 MeV.

**Explicit form**:
```
V_soft = f_π^2 Σ_{i,j=u,d,s} |M^i_j|^2
```

This acts as a universal mass-squared term for all 9 meson scalar components:
m̃^2 = f_π^2 = 8464 MeV^2.

**STATUS: ASSUMED** with physical motivation. The identification of the soft SUSY-breaking
scale with f_π^2 is motivated by the pion decay constant setting the chiral symmetry
breaking scale. It is NOT derived from a mediation mechanism.

**What it determines**: The supertrace STr[M^2] = 2 × 9 × f_π^2 = 152352 MeV^2
(18 scalars each getting mass f_π^2 from V_soft). Verified to be unchanged by the
NMSSM coupling.

**What is NOT included in V_soft** (and why):
- No soft mass for X: The pseudo-modulus VEV is fixed by the Kähler pole mechanism,
  not by a soft mass. Adding m̃_X^2 |X|^2 would compete with this mechanism.
- No soft mass for B, B̃: Baryons are assumed to decouple.
- No soft mass for H_u, H_d: These are generated effectively by the EWSB mechanism
  through the NMSSM structure (m_{H_u}^2 = b − y_c^2, m_{H_d}^2 = b − y_b^2
  at tan β = 1, where b is the bilinear B_μ term).
- No A-terms (soft trilinear couplings): Not included. This is an assumption.
- No gaugino masses: The model is in the confined phase of SQCD, so SU(3)_c
  gaugino masses do not appear. SU(2)_L and U(1)_Y gaugino masses are not specified.

---

## 7. Complete Lagrangian (assembled)

The full scalar potential is:

```
V = V_F + V_D + V_soft
  = K^{IJ̄} F_I F̄_J
    + (g^2/2)(H_u^† T^a H_u + H_d^† T^a H_d)^2
    + (g'^2/2)(Y_I |Φ_I|^2)^2
    + f_π^2 Tr(M^†M)
```

The fermion mass matrix is:

```
W_{IJ} = ∂^2 W / ∂Φ_I ∂Φ_J
```

with block structure:

**Central 6×6 block** {M^u_u, M^d_d, M^s_s, X, H_u^0, H_d^0}:
```
W[M^u_u, M^d_d] = X M_s          (= −1.143×10^{-5} MeV at vacuum)
W[M^u_u, M^s_s] = X M_d          (= −2.287×10^{-4} MeV)
W[M^d_d, M^s_s] = X M_u          (= −4.944×10^{-4} MeV)
W[M^u_u, X]     = M_d M_s        (= +1.785×10^{9} MeV)
W[M^d_d, X]     = M_u M_s        (= +3.859×10^{9} MeV)
W[M^s_s, X]     = M_u M_d        (= +7.717×10^{10} MeV)
W[M^d_d, H_u^0] = y_c            (= 1.032×10^{-2} MeV)
W[M^s_s, H_d^0] = y_b            (= 3.395×10^{-2} MeV)
W[X, H_u^0]     = λ ⟨H_d^0⟩ = λv/√2  (= 1.254×10^5 MeV)
W[X, H_d^0]     = λ ⟨H_u^0⟩ = λv/√2  (= 1.254×10^5 MeV)
W[H_u^0, H_d^0] = λ X            (= −8.714×10^{-10} MeV)
```

**Three off-diagonal meson 2×2 blocks**:
```
W[M^u_d, M^d_u] = X M_s    (= −1.143×10^{-5} MeV)
W[M^u_s, M^s_u] = X M_d    (= −2.287×10^{-4} MeV)
W[M^d_s, M^s_d] = X M_u    (= −4.944×10^{-4} MeV)
```

**Baryon 2×2 block**:
```
W[B, B̃] = −X    (= +1.210×10^{-9} MeV^{-4})
```

**Charged Higgs 2×2 block**:
```
W[H_u^+, H_d^-] = −λ X    (= +8.714×10^{-10} MeV^{-4})
```

---

## 8. Free Parameter Count

### Free parameters (not determined by the theory):

| Parameter | Value | Role |
|-----------|-------|------|
| m_u | 2.16 MeV | Light quark mass (lightest meson VEV scale) |
| m_d | 4.67 MeV | Light quark mass (Cabibbo angle via Oakes: tan θ_C = √(m_d/m_s)) |
| m_s | 93.4 MeV | Overall quark mass scale (fixes C = Λ^2(m_u m_d m_s)^{1/3}) |
| Λ | 300 MeV | SQCD confinement scale |
| c_3 | O(1) | Three-instanton coefficient (shifts X vacuum, no observable effect) |

**Total free parameters in W**: 5 (or 4, if c_3 is dropped as unphysical).

### Fitted parameters (determined by matching to data):

| Parameter | Value | What it determines |
|-----------|-------|-------------------|
| y_c = 2m_c/v | 1.032×10^{-2} | Charm Yukawa coupling (m_c = 1270 MeV) |
| y_b = 2m_b/v | 3.395×10^{-2} | Bottom Yukawa coupling (m_b = 4180 MeV) |
| y_t = √2 m_t/v | 0.9934 | Top Yukawa coupling (m_t = 172760 MeV) |
| λ | 0.72 | NMSSM coupling (m_h = 125.25 GeV at tan β = 1) |

**Note on m_c and m_b**: In the sBootstrap program, m_c is predicted by the
O'Raifeartaigh–Koide mechanism (gv/m = √3, giving m_c/m_s = (2+√3)^2/(2−√3)^2 = 13.93,
so m_c = 93.4 × 13.93 = 1301 MeV vs PDG 1270 MeV, 2.4% deviation), and m_b is predicted
by v_0-doubling (m_b = (3√m_s + √m_c)^2 = 4177 MeV vs PDG 4180 MeV, 0.07% deviation).
If these predictions are accepted, then y_c and y_b are predicted, not fitted.
In this Lagrangian as stated, they are FITTED.

### Parameters determined by self-consistency:

| Parameter | Determination |
|-----------|---------------|
| μ (Kähler mass) | μ = C/(√3 Λ^6) [Kähler pole ↔ seesaw vacuum coincidence] |
| ζ (monopole fugacity) | ζ = exp(−π/(3α_s)) [from α_s at confinement scale] |
| tan β | = 1 [required by D-flat direction for m_h = 0 at MSSM tree level; lifted by λ] |
| B_μ bilinear | b = m_{H_u}^2 + y_c^2 = m_{H_d}^2 + y_b^2 [from EWSB at tan β = 1] |

### Assumed parameters (not derived):

| Parameter | Assumption |
|-----------|------------|
| f_π = 92 MeV | Soft mass scale = pion decay constant |
| Flavor assignment | y_c couples to M^d_d, y_b couples to M^s_s (not derived) |
| SU(2)_L × U(1)_Y gaugino masses | Not specified |
| A-terms | Set to zero |

---

## 9. Missing Terms and Gaps

### Gap 1: Lepton masses — ABSENT

Lepton masses are not generated by any term in this Lagrangian. The central-block
fermion eigenstates of the 6×6 {M^u_u, M^d_d, M^s_s, X, H_u^0, H_d^0} matrix
have masses ∼ 0.77 MeV, 0.78 MeV, 10.4 MeV at the seesaw vacuum.
These do NOT match (m_e, m_μ, m_τ) = (0.511, 105.7, 1777) MeV; the ratios
are wrong by factors of 10^2–10^5.

Possible insertion points:
- Direct lepton Yukawa terms W_lepton = y_e H_d^0 E_i^c L_i are the standard
  MSSM terms but require elementary slepton fields not present in this spectrum.
- The folding-pattern argument (lepton ↔ meson parity) from the sBootstrap program
  suggests leptons are related to the SU(5) flavor sector, but the coupling is
  not written down.
- **This is an explicit gap in the Lagrangian.**

### Gap 2: Top quark coupling to the meson sector — ABSENT

The top appears as y_t H_u^0 t t̄ with t an elementary chiral superfield.
The relationship between the elementary top and the meson matrix M is not specified.
In particular:
- Is there a "top meson" T = t t̄ / Λ in an extension of the matrix M to 5×5 or 6×6?
- If SU(5) flavor symmetry is present, the top should appear in a 5 × 5 meson matrix
  with the 5th flavor = top. This would require a 5×5 extension of the Seiberg dual
  description, which is a different model (N_f = 5, N_c = 3 SQCD).
- The "top is external" assumption prevents closing the model. It means the Yukawa
  coupling y_t H_u^0 t t̄ must be put in by hand without a compositeness relation.
- **This is an explicit gap.**

### Gap 3: μ-term for Higgs — PARTIALLY SPECIFIED

The μ-term arises from ⟨X⟩ through the NMSSM coupling:
```
μ_eff = λ ⟨X⟩ = λ × (−C/Λ^6)
```

Numerically:
```
⟨X⟩ = −C/Λ^6 = −882297/(3×10^{14}) = −2.94×10^{-9} MeV^{-4}
μ_eff = 0.72 × (−2.94×10^{-9}) = −2.12×10^{-9} MeV^{-4}
```

This value has dimension [mass]^{-4}, NOT [mass]^1. The issue is that [X] = [mass]^{-4}
in the seesaw description (X has dimension mass^{-4}). The standard μ term in the MSSM has
[μ] = mass. Therefore, the direct identification μ_eff = λ ⟨X⟩ gives a dimensionally
incorrect μ-term. Either:
- The Lagrange multiplier X needs to be rescaled to a dimensionless or mass-dimension-1
  field before identification as the NMSSM singlet, or
- The coupling constant λ absorbs the dimensional mismatch (λ has dimension [mass]^5).

This dimensional issue must be resolved for a consistent model. The NMSSM coupling
W_NMSSM = λ X (H_u · H_d) is written with λ dimensionless in the task specification,
which is inconsistent with [X] = mass^{-4} and [H_u · H_d] = mass^2. A consistent
coupling would require [λ] = mass^3, making λ a dimensionful parameter.

**Working resolution**: Treat X as dimensionlessly rescaled: X̂ = X / X_0 where
X_0 = ⟨X⟩ is the seesaw VEV. Then ⟨X̂⟩ = 1 and the NMSSM coupling becomes
W_NMSSM = λ' X̂ (H_u · H_d) with λ' = λ × X_0 having dimension [mass]^{-2},
and the μ-term μ_eff = λ' with appropriate dimensions. The parameter λ = 0.72
is then the dimensionless coupling in the rescaled field basis.
**This rescaling is not specified in the task statement and is assumed here.**

### Gap 4: Gauge kinetic terms — ASSUMED CANONICAL

The gauge kinetic function is specified as f_ab = δ_ab / g^2 (canonical). However:
- SU(3)_c is confined: the correct gauge kinetic term is the strong coupling α_s(μ)
  running to the confinement scale. No non-canonical gauge kinetic corrections are
  included, but threshold corrections at Λ are not computed.
- SU(2)_L × U(1)_Y kinetic terms are canonical: g = 0.653, g' = 0.350 at M_Z.
  This is consistent with the D-term calculation.
- **Assumed**: No gauge kinetic mixing between U(1)_Y and the meson flavor U(1)'s.

### Gap 5: Off-diagonal meson Yukawa couplings — ABSENT

The Yukawa sector couples only to diagonal mesons M^d_d and M^s_s.
The off-diagonal mesons M^i_j (i ≠ j) are not coupled to the Higgs.
This is an assumption. Flavor-changing Higgs couplings through off-diagonal mesons
would be:
```
y_cu H_u^0 M^u_d + y_cs H_u^0 M^s_d + ...   [if included]
```
These are absent. Without a derivation from first principles, their absence is assumed.

### Gap 6: Neutrino masses — ABSENT

No right-handed neutrino superfields N^c_i are present. No Majorana mass term or
seesaw mechanism for neutrinos is specified. This is a known omission.

### Gap 7: Higher-dimensional Kähler operators

The task specification mentions K_bion with c = -1/12 and p = 1, suggesting
K_bion ∝ c |Tr(M^†M)| / Λ^{2p-2}.

The physically derived bion Kähler potential (Section 2c) has the form
(ζ^2/Λ^2) |Σ_k s_k √M_k|^2, which is a non-holomorphic function involving
square roots of diagonal meson fields. This is NOT the same as a power of
Tr(M^†M). The two forms are different and the task specification's parameterization
is not the one derived from bion physics.

If c = -1/12 and p = 1 is intended as K_bion = (-1/12) Tr(M^†M) / Λ^0 = (-1/12) Tr(M^†M),
this would be a negative kinetic term correction to all mesons. This is different
from the bion correction. The physical bion correction (complete_kahler.md) produces
c_bion = -1/12 for the X quartic, not for Tr(M^†M).

---

## 10. Prediction Table

| Prediction | Formula | Value | PDG | Deviation | Source |
|------------|---------|-------|-----|-----------|--------|
| m_h | λ v / √2 | 125.35 GeV | 125.25 GeV | 0.08% | W_NMSSM at tan β = 1 |
| m_b | (3√m_s + √m_c)^2 | 4177 MeV | 4180 ± 30 MeV | 0.07% (0.1σ) | K_bion v_0-doubling |
| m_c | m_s (2+√3)^2/(2−√3)^2 | 1301 MeV | 1270 ± 20 MeV | 2.4% | O'R–Koide (K_X pole) |
| tan θ_C | √(m_d/m_s) = √(4.67/93.4) | 0.2236 | 0.2243 | 0.3% | W_Seiberg + Oakes |
| Q(e,μ,τ) | Koide quotient | 2/3 (0.91σ) | measured | 0.91σ | Koide constraint (not W) |
| Q(−s,c,b) | Koide quotient | 0.675 (1.24% from 2/3) | — | — | K_bion minimum |
| Q(c,b,t) | Koide quotient | 0.669 (0.41% from 2/3) | — | — | Overdetermination |
| Q_dual | Q(1/m_d, 1/m_s, 1/m_b) | 0.665 | — | 0.22% from 2/3 | Seiberg seesaw |
| STr[M^2] | 18 f_π^2 | 152352 MeV^2 | — | — | V_soft (exact) |
| F_{H_u} = F_{H_d} | 2C/v | 3.58 MeV each | — | — | Flavor universality |
| m_{H^+} | √(m_A^2 + m_W^2) | (input m_A) | — | — | D-term + V_soft |
| Charged Higgsino mass | λ|X| | 8.71×10^{-10} MeV | — | — | W_NMSSM |

**Not predicted** (required as inputs):
- m_u, m_d, m_s (light quark masses)
- Λ (SQCD confinement scale)
- m_t (top mass, external)
- g, g', m_Z (EW sector, external)
- Lepton masses (no mechanism)
- Neutrino masses (no mechanism)

---

## 11. Complete Superpotential in LaTeX Form

```latex
\begin{align}
W &= \underbrace{%
  \sum_{i \in \{u,d,s\}} m_i \, M^i_{\;i}
  + X\!\left(\det M - B\tilde{B} - \Lambda^6\right)%
}_{W_{\text{Seiberg}}}
\notag \\
&\quad + \underbrace{%
  y_c \, H_u^0 \, M^d_{\;d}
  + y_b \, H_d^0 \, M^s_{\;s}
  + y_t \, H_u^0 \, t\bar{t}%
}_{W_{\text{Yukawa}}}
\notag \\
&\quad + \underbrace{%
  \lambda \, X \!\left(H_u^0 H_d^0 - H_u^+ H_d^-\right)%
}_{W_{\text{NMSSM}}}
\notag \\
&\quad + \underbrace{%
  c_3 \, \frac{(\det M)^3}{\Lambda^{18}}%
}_{W_{\text{instanton}}}
\end{align}
```

with parameter values:
```latex
\begin{align}
y_c &= \frac{2m_c}{v} = 1.032 \times 10^{-2},
  &\quad m_c &= 1270~\text{MeV},~v = 246220~\text{MeV}
\\
y_b &= \frac{2m_b}{v} = 3.395 \times 10^{-2},
  &\quad m_b &= 4180~\text{MeV}
\\
y_t &= \frac{\sqrt{2}\, m_t}{v} = 0.9934,
  &\quad m_t &= 172760~\text{MeV}
\\
\lambda &= \frac{\sqrt{2}\, m_h}{v} = 0.7186 \approx 0.72,
  &\quad m_h &= 125250~\text{MeV}
\end{align}
```

and complete Kähler potential:
```latex
\begin{align}
K &= \operatorname{Tr}(M^\dagger M) + |X|^2 + |B|^2 + |\tilde{B}|^2
    + |H_u|^2 + |H_d|^2
\notag \\
&\quad - \frac{|X|^4}{12\mu^2}
  + \frac{\zeta^2}{\Lambda^2}
    \Bigl|\sum_{k} s_k \sqrt{M^k_{\;k}}\Bigr|^2
\end{align}
```

and soft-breaking potential:
```latex
V_{\text{soft}} = f_\pi^2 \operatorname{Tr}(M^\dagger M),
  \qquad f_\pi = 92~\text{MeV}
```

---

## 12. Honest Summary of What Is and Is Not Determined

**Determined (no free parameters, follows from structure)**:
- The form of W_Seiberg (Seiberg duality, exact)
- The form of W_NMSSM (gauge invariance + holomorphy)
- The D-term structure V_D (gauge group + charge assignments)
- The Kähler pole mechanism (self-consistency of seesaw + K_X)
- The flavor universality identity F_{H_u} = F_{H_d} = 2C/v (algebraically exact)
- The supertrace STr[M^2] = 18 f_π^2 (exact, from V_soft structure)
- The Higgs mass formula m_h = λv/√2 at tan β = 1 (algebraically exact)

**Fitted (parameters adjusted to match known data)**:
- y_c (adjusted to m_c = 1270 MeV), y_b (adjusted to m_b = 4180 MeV)
- y_t (adjusted to m_t = 172760 MeV)
- λ (adjusted to m_h = 125.25 GeV)
- m_u, m_d, m_s (PDG inputs, 3 free parameters)
- Λ = 300 MeV (chosen to equal ≈ f_π × √(4π); the specific choice is an input)

**Assumed without derivation**:
- The identification X (Seiberg Lagrange multiplier) = S (NMSSM singlet)
- The flavor assignment: which meson diagonal couples to which Higgs
- V_soft = f_π^2 Tr(M^†M) with coefficient f_π^2 (not derived from mediation)
- Absence of A-terms and gaugino masses
- Absence of off-diagonal meson Yukawa couplings
- c_3 = O(1) (three-instanton coefficient)

**Missing (absent from Lagrangian, gap)**:
- Lepton masses: no mechanism specified
- Top-meson sector coupling: top is treated as external elementary field
- Neutrino masses: absent
- Dimensional consistency of X: the field has dimension [mass]^{-4}, creating
  an inconsistency in W_NMSSM = λ X (H_u · H_d) if λ is dimensionless
- Non-canonical gauge kinetic corrections at confinement scale
- Derivation of the bion mechanism in 4D (not just R^3 × S^1)

---

*Generated: 2026-03-04*
*Source files: unified_superpotential.md, complete_kahler.md, nmssm_spectrum.md,*
*higgs_potential.md, sqcd_yukawa_spectrum.md, bion_kahler.md, v0_doubling.md,*
*bloom_mechanism.md, oraifeartaigh_koide.md, seesaw_mechanism.md,*
*vacuum_structure.md, iss_cw_koide.md, lepton_yukawa.md, dterm_fi.md*
