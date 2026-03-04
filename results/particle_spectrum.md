# Complete N=1 SUSY Particle Spectrum

## Theory Setup

SU(3) SQCD with N_f = N_c = 3 (u, d, s), confined at Λ = 300 MeV.
Coupled to NMSSM Higgs sector (H_u, H_d, S) with tan β = 1.

### Parameters

| Parameter | Value |
|-----------|-------|
| Λ | 300 MeV |
| v | 246.22 GeV |
| tan β | 1 (v_u = v_d = v/√2 = 174.10 GeV) |
| f_π | 92 MeV |
| λ_S | 0.72 |
| m_B | 300 MeV |
| m_u | 2.16 MeV |
| m_d | 4.67 MeV |
| m_s | 93.4 MeV |
| m_c | 1275 MeV |
| m_b | 4180 MeV |
| m_t | 172760 MeV |

### Derived Vacuum Quantities

- C = (Λ⁶ m_u m_d m_s)^(1/3) = 882297.43 MeV²
- M_uu = C/m_u = 408471.030 MeV²
- M_dd = C/m_d = 188928.785 MeV²
- M_ss = C/m_s = 9446.439 MeV²
- X₀ = -C/Λ⁶ = -1.2103e-09 (units: MeV⁻⁴ per problem statement)
- det(M_vev)/Λ⁶ = 1.000000000000 ≈ 1.0 ✓

### Key F-term W Matrix Entries at Vacuum

| Entry | Formula | Value |
|-------|---------|-------|
| W_{M_uu, M_dd} | X₀ M_ss | -1.1433e-05 |
| W_{M_uu, M_ss} | X₀ M_dd | -2.2866e-04 |
| W_{M_dd, M_ss} | X₀ M_uu | -4.9437e-04 |
| W_{M_ud, M_du} | X₀·(-M_ss) = +|X₀| M_ss | 1.1433e-05 |
| W_{M_us, M_su} | X₀·(-M_dd) = +|X₀| M_dd | 2.2866e-04 |
| W_{M_ds, M_sd} | X₀·(-M_uu) = +|X₀| M_uu | 4.9437e-04 |
| W_{M_uu, X} | Λ⁶/M_uu | 1.7847e+09 |
| W_{M_dd, X} | Λ⁶/M_dd | 3.8586e+09 |
| W_{M_ss, X} | Λ⁶/M_ss | 7.7172e+10 |
| W_{B, B̃} | m_B - X₀ | 300.000000 ≈ 300.0 MeV |

---

## Part (a): Scalar Meson Masses

### Method

Scalar mass² from F-terms (with X integrated out, contributing no scalar d.o.f.):

    (m²_scalar)_IJ = Σ_{K: kinetic term} W_IK W*_JK + f_π² δ_IJ (mesons only)

Since X has no kinetic term, W_{M,X} entries do NOT contribute to scalar masses.
The soft term V_soft = f_π² Tr(M†M) adds f_π² to each meson scalar mass².

### Off-Diagonal Scalar Mesons (6 complex = 12 real scalars)

For M^a_b with a ≠ b, the only relevant second derivative is W_{M_ab, M_ba}.
The F-term contribution comes from the kinetic fields only:

    m²(M^a_b) = |W_{M_ab, M_ba}|² + 2·f_π²
              = |X₀|² M_c² + 2·f_π²

where c is the spectator flavor index, and the factor 2 comes from both
the real and imaginary components of the complex scalar (no B-terms).

| Scalar | Spectator | |X₀|² M_c² | m² (MeV²) | m (MeV) |
|--------|-----------|------------|-----------|---------|
| M^u_d and M^d_u | M_ss | 1.3071e-10 | 16928.0000 | 130.1076 |
| M^u_s and M^s_u | M_dd | 5.2284e-08 | 16928.0000 | 130.1076 |
| M^d_s and M^s_d | M_uu | 2.4440e-07 | 16928.0000 | 130.1076 |

All off-diagonal scalar masses ≈ √(2) × f_π = 130.1076 MeV ≈ 130.1 MeV
(The |X₀|² M² correction is negligible: 1.31e-10 << 2f_π² = 16928)

Each complex scalar decomposes into 2 real scalars (CP-even + CP-odd) with equal masses.
Total: 12 real scalars at m ≈ 130.1 MeV

### Diagonal Scalar Mesons (3 complex = 6 real scalars)

The 3×3 mass² matrix for {M_uu, M_dd, M_ss}:

    (m²_F)_ij = X₀² (M_i M_j - M_i² δ_ij)  = X₀² (outer product - diagonal²)
    (m²_total)_ij = (m²_F)_ij + f_π² δ_ij

Explicitly:
    (m²_F)_uu,uu = X₀²(M_ss² + M_dd²) = 5.2415e-08
    (m²_F)_dd,dd = X₀²(M_ss² + M_uu²) = 2.4453e-07
    (m²_F)_ss,ss = X₀²(M_dd² + M_uu²) = 2.9668e-07
    (m²_F)_uu,dd = X₀² M_dd M_uu      = 1.1304e-07
    f_π² = 8464.0000 MeV²

All F-term entries ≪ f_π². The mass² matrix is dominated by f_π² I_3.

Eigenvalues of full 3×3 diagonal mass² matrix:
    λ_1 = 8464.000000 MeV² → m = 92.000000 MeV
    λ_2 = 8464.000000 MeV² → m = 92.000000 MeV
    λ_3 = 8464.000000 MeV² → m = 92.000000 MeV

All diagonal scalars: m ≈ f_π = 92.0 MeV (to <0.01% correction from F-terms)
Total: 6 real scalars at m ≈ 92 MeV (CP-even and CP-odd degenerate)

### Scalar Meson Summary

| Type | Count (real) | Mass (MeV) | Origin |
|------|-------------|------------|--------|
| Diagonal (Re part) | 3 | 92.0 | V_soft = f_π² |
| Diagonal (Im part) | 3 | 92.0 | V_soft = f_π² |
| Off-diagonal (Re) | 6 | 130.11 | 2f_π² (factor 2 from 2 F-terms) |
| Off-diagonal (Im) | 6 | 130.11 | 2f_π² |
| **Total** | **18** | | |

---

## Part (b): Fermion Masses (Mesinos)

### Dimensional Analysis and Kahler Normalization

The meson field M^i_j is a composite operator (q_i q̄_j) with canonical dimension [M] = mass².
From C = (Λ⁶ m_u m_d m_s)^(1/3), we have [C] = MeV³ and [M_vev] = MeV².
The Lagrange multiplier satisfies [X₀] = MeV⁻³ (from [W] = [X][Λ⁶] = mass³).

**The Kahler metric** for the confined meson is non-canonical:

    K = (1/Λ²) Tr(M† M)   ⟹   normalized field: M_norm = M/Λ,   [M_norm] = MeV

The **physical (Kahler-normalized) mass** requires rescaling by Λ²:

    m_physical = Λ² × W_{M_ab, M_ba}

This converts the raw W entries (in MeV⁻¹) to physical masses in MeV.

### Block Structure of the 12×12 Fermion Mass Matrix

After integrating out ψ_X (no kinetic term → auxiliary fermion):

**Off-diagonal mesino blocks** (three independent 2×2 Majorana blocks):

For pair {ψ_{M_ab}, ψ_{M_ba}} with a ≠ b:

    W_{M_ab, M_ba} = X₀ × (cofactor sign) × M_c
    m_physical = Λ² × |W_{M_ab, M_ba}| = Λ² × |X₀| × M_c = C²/(Λ⁴ × m_spectator)

At diagonal vacuum, the cofactor for off-diagonal elements carries a relative minus sign.
With X₀ < 0: W_{M_ud, M_du} = X₀ × (-M_ss) = |X₀| × M_ss  (positive entry).

Each 2×2 block [[0, μ],[μ, 0]] has eigenvalues ±μ (one Dirac fermion of mass μ).

**Raw W entries** (before Kahler normalization, in MeV⁻¹):

| Block | W_IJ = |X₀| × M_c | Physical mass = Λ² × W_IJ (MeV) |
|-------|-------------------|----------------------------------|
| ψ_ud / ψ_du | 1.1433e-05 | 1.029 MeV |
| ψ_us / ψ_su | 2.2866e-04 | 20.579 MeV |
| ψ_ds / ψ_sd | 4.9437e-04 | 44.493 MeV |

Closed-form expression: m(ψ_{ab}) = C²/(Λ⁴ m_c) where c is the spectator flavor.

**Diagonal mesino block** (after integrating out ψ_X):

The 4×4 block {ψ_uu, ψ_dd, ψ_ss, ψ_X} with W_{XX} = 0:

    W_{4×4} = [[0,        X₀M_s,   X₀M_d,   Λ⁶/M_u],
               [X₀M_s,   0,        X₀M_u,   Λ⁶/M_d],
               [X₀M_d,   X₀M_u,   0,        Λ⁶/M_s],
               [Λ⁶/M_u,  Λ⁶/M_d,  Λ⁶/M_s,  0      ]]

Since W_{XX} = 0, ψ_X acts as a Lagrange multiplier enforcing the constraint:

    (Λ⁶/M_u) ψ_uu + (Λ⁶/M_d) ψ_dd + (Λ⁶/M_s) ψ_ss = 0

equivalently (using Λ⁶/M_i = m_i Λ³/C^{2/3}):

    m_u ψ_uu + m_d ψ_dd + m_s ψ_ss = 0   (one direction eliminated)

This leaves 2 physical diagonal mesinos. The 4×4 eigenvalues are:

| n | Raw eigenvalue | Physical mass = Λ² × |eigenvalue| | Type |
|---|---------------|-----------------------------------|------|
| 1 | 7.5476e-06 | 0.679 MeV | light (X₀·M scale) |
| 2 | 5.2459e-05 | 4.721 MeV | light (X₀·M scale) |
| 3 | ±7.7289e+10 | heavy, ~Λ⁹/C² | heavy (Λ⁶/M scale, integrated out) |

The two heavy eigenvalues (±7.73×10¹⁰ MeV⁻¹, normalized to ~Λ⁹/C² ~ heavy scale)
correspond to the X-dominated modes, integrated out at low energies.

Projected 3×3 eigenvalues (constraint-satisfying subspace):

    0 (constraint direction, absorbed by ψ_X)
    7.461e-06 → physical mass: 0.672 MeV
    5.240e-05 → physical mass: 4.716 MeV

Note: the slight difference between the 4×4 light modes and the projected 3×3 values
reflects the approximate nature of the block decomposition at this order.

### Physical Mesino Spectrum Summary

Physical masses after Kahler normalization m = Λ² × W_entry:

| Mesino | Physical Mass (MeV) | Formula | Description |
|--------|--------------------|---------| ------------|
| ψ_ud = ψ_du (Dirac) | 1.029 | C²/(Λ⁴ m_s) | Off-diag, spectator s |
| ψ_us = ψ_su (Dirac) | 20.579 | C²/(Λ⁴ m_d) | Off-diag, spectator d |
| ψ_ds = ψ_sd (Dirac) | 44.493 | C²/(Λ⁴ m_u) | Off-diag, spectator u |
| Diagonal mesino 1 | ≈ 0.68 | from 4×4 block | Q=0, neutral |
| Diagonal mesino 2 | ≈ 4.72 | from 4×4 block | Q=0, neutral |
| Heavy X-mode 1 | Λ⁹/C² scale | integrated out | X-dominated |
| Heavy X-mode 2 | Λ⁹/C² scale | integrated out | X-dominated |

Note: all off-diagonal mesino masses satisfy m(ψ_{ab}) ∝ 1/m_spectator (dual mass relation).

---

## Part (c): Baryons

With m_B = Λ = 300 MeV and B = B̃ = 0 at vacuum.

**Baryon mass term in W:**

    W ⊃ -X(BB̃) + m_B BB̃ = (m_B - X)BB̃

**Baryonino** (Weyl fermion pair ψ_B, ψ_B̃ → 1 Dirac fermion):

    m_baryonino = |W_{B,B̃}| = |m_B - X₀| ≈ m_B = 300.000000 MeV

(The correction X₀ ≈ -1.21×10⁻⁹ is negligible.)

**Baryon scalars** (2 complex = 4 real, from B and B̃):

F-term contribution: from F_B = (m_B - X₀)B̃ and F_B̃ = (m_B - X₀)B:

    m²_scalar(B) = |m_B - X₀|² = 90000.000001 MeV²
    m_scalar(B) = 300.000000 MeV ≈ Λ = 300 MeV

No soft SUSY-breaking term for baryons (V_soft only for mesons in this model).

| Particle | Spin | Mass (MeV) | Count |
|----------|------|-----------|-------|
| Baryon B (scalar) | 0 | 300.00 | 2 real (Re B) |
| Anti-baryon B̃ (scalar) | 0 | 300.00 | 2 real (Re B̃) |
| Baryonino (ψ_B, ψ_B̃) | 1/2 | 300.0000 | 1 Dirac |

---

## Part (d): Higgs Sector (NMSSM)

**Superpotential:** W ⊃ λ_S S H_u·H_d + κ/3 S³

**VEVs:** ⟨H_u⁰⟩ = v_u = v/√2, ⟨H_d⁰⟩ = v_d = v/√2 (tan β = 1)

**Tree-level lightest CP-even Higgs** (tan β = 1 → cos 2β = 0):

    m_h² = m_Z² cos²(2β) + λ_S² v² sin²(2β)/2
          = 0 + λ_S² v²/2  (since sin 2β = 1 at tan β = 1)

    m_h = λ_S v/√2 = 0.72 × 246.22 GeV / √2 = 125.3548 GeV ≈ 125 GeV

**Yukawa couplings** (from W ⊃ y_c H_u M^c_c + y_b H_d M^b_b + y_t H_u Q^t Q̄_t):

| Quark | Yukawa y | Fermion mass (MeV) | Scalar partner (squark) |
|-------|----------|-------------------|------------------------|
| c | y_c = m_c/v_u = 0.007323 | 1275.0 MeV | mass from soft terms |
| b | y_b = m_b/v_d = 0.024009 | 4180.0 MeV | mass from soft terms |
| t | y_t = m_t/v_u = 0.992281 | 172760.0 MeV | mass from soft terms |

**NMSSM spectrum** (soft-dependent, not fully computed):

| Particle | Spin | Mass | Notes |
|----------|------|------|-------|
| h (CP-even Higgs) | 0 | 125.355 GeV | Tree-level |
| H (heavy CP-even) | 0 | ≫ m_h | Soft dependent |
| A (CP-odd) | 0 | Soft dependent | |
| H± (charged Higgs) | 0 | √(m_A²+m_W²) | Soft dependent |
| S (singlet scalar) | 0 | κ⟨S⟩ | κ, ⟨S⟩ free |
| H̃_u (Higgsino) | 1/2 | λ_S⟨S⟩ = μ_eff | |
| H̃_d (Higgsino) | 1/2 | λ_S⟨S⟩ = μ_eff | |
| S̃ (singlino) | 1/2 | 2κ⟨S⟩ | |

---

## Part (e): Complete Spectrum Table

### Group 1: Light Confined Sector (meson scalars + off-diagonal mesinos)

Scalar masses from V_soft; fermion masses from Kahler-normalized W entries (m = Λ² × W_entry).

| Particle | Spin | Q_em | SU(2)_L | Mass (MeV) | Origin |
|----------|------|------|---------|-----------|--------|
| M_uu scalar Re | 0 | 0 | singlet | 92.0 | f_π (dominant soft) |
| M_uu scalar Im | 0 | 0 | singlet | 92.0 | f_π |
| M_dd scalar Re | 0 | 0 | singlet | 92.0 | f_π |
| M_dd scalar Im | 0 | 0 | singlet | 92.0 | f_π |
| M_ss scalar Re | 0 | 0 | singlet | 92.0 | f_π |
| M_ss scalar Im | 0 | 0 | singlet | 92.0 | f_π |
| M_ud scalar (Re,Im) | 0 | +1 | singlet | 130.11 | √2 f_π |
| M_du scalar (Re,Im) | 0 | −1 | singlet | 130.11 | √2 f_π |
| M_us scalar (Re,Im) | 0 | +2/3 | singlet | 130.11 | √2 f_π |
| M_su scalar (Re,Im) | 0 | −2/3 | singlet | 130.11 | √2 f_π |
| M_ds scalar (Re,Im) | 0 | −1/3 | singlet | 130.11 | √2 f_π |
| M_sd scalar (Re,Im) | 0 | +1/3 | singlet | 130.11 | √2 f_π |
| Mesino ψ_ud/ψ_du (Dirac) | 1/2 | ±1 | singlet | **1.029** | Λ²|X₀|M_ss = C²/Λ⁴m_s |
| Mesino ψ_us/ψ_su (Dirac) | 1/2 | ±2/3 | singlet | **20.579** | Λ²|X₀|M_dd = C²/Λ⁴m_d |
| Mesino ψ_ds/ψ_sd (Dirac) | 1/2 | ±1/3 | singlet | **44.493** | Λ²|X₀|M_uu = C²/Λ⁴m_u |

### Group 2: Diagonal Mesinos + Heavy X-dominated Modes

| Particle | Spin | Q_em | Physical Mass (MeV) | Origin |
|----------|------|------|--------------------|----|
| Diagonal mesino 1 (Q=0) | 1/2 | 0 | **≈ 0.68** | 4×4 block, light mode, Kahler normalized |
| Diagonal mesino 2 (Q=0) | 1/2 | 0 | **≈ 4.72** | 4×4 block, light mode, Kahler normalized |
| X-dominated mesino 1 | 1/2 | 0 | heavy (Λ⁹/C² scale) | Λ⁶/M coupling, integrated out |
| X-dominated mesino 2 | 1/2 | 0 | heavy (Λ⁹/C² scale) | Λ⁶/M coupling, integrated out |
| X scalar | none | — | not propagating | No kinetic term |

### Group 3: Baryons

| Particle | Spin | Q_em | SU(2)_L | Mass (MeV) | Origin |
|----------|------|------|---------|-----------|--------|
| B (complex scalar) | 0 | 0 | singlet | 300.00 | F_B = m_B B̃ |
| B̃ (complex scalar) | 0 | 0 | singlet | 300.00 | F_B̃ = m_B B |
| ψ_B, ψ_B̃ (Dirac) | 1/2 | 0 | singlet | 300.0000 | W_{B,B̃} = m_B |

### Group 4: Elementary Quarks (above confinement threshold)

| Particle | Spin | Q_em | SU(3)_c | SU(2)_L | Mass (MeV) | Origin |
|----------|------|------|---------|---------|-----------|--------|
| c quark (Dirac) | 1/2 | +2/3 | 3 | doublet | 1275.0 | y_c H_u M^c_c |
| b quark (Dirac) | 1/2 | −1/3 | 3 | singlet | 4180.0 | y_b H_d M^b_b |
| t quark (Dirac) | 1/2 | +2/3 | 3 | doublet | 172760.0 | y_t H_u Q_t Q̄_t |
| c~ squark (scalar) | 0 | +2/3 | 3 | doublet | soft-breaking | Superpartner of c |
| b~ squark (scalar) | 0 | −1/3 | 3 | singlet | soft-breaking | Superpartner of b |
| t~ squark (scalar) | 0 | +2/3 | 3 | doublet | soft-breaking | Superpartner of t |

### Group 5: Higgs Sector

| Particle | Spin | Q_em | SU(2)_L | Mass | Origin |
|----------|------|------|---------|------|--------|
| h (CP-even Higgs) | 0 | 0 | singlet | 125.355 GeV | λ_S v/√2, tan β=1 |
| H (heavy CP-even) | 0 | 0 | singlet | soft-dep. | Decoupling limit |
| A (CP-odd Higgs) | 0 | 0 | singlet | soft-dep. | |
| H± (charged Higgs) | 0 | ±1 | doublet | soft-dep. | |
| S (singlet scalar) | 0 | 0 | singlet | κ⟨S⟩ | W ⊃ κ/3 S³ |
| H̃_u± (Higgsino) | 1/2 | ±1 | doublet | μ_eff=λ_S⟨S⟩ | |
| H̃_u⁰ (Higgsino) | 1/2 | 0 | doublet | μ_eff | |
| H̃_d± (Higgsino) | 1/2 | ±1 | doublet | μ_eff | |
| H̃_d⁰ (Higgsino) | 1/2 | 0 | doublet | μ_eff | |
| S̃ (singlino) | 1/2 | 0 | singlet | 2κ⟨S⟩ | W_SS = κS term |

### Group 6: Gauge Bosons (after EW breaking)

| Particle | Spin | Q_em | Mass (MeV) | Notes |
|----------|------|------|-----------|-------|
| Gluons (8) | 1 | 0 | — | Confined, no propagating gluons |
| W± bosons | 1 | ±1 | 80379 (≈80.4 GeV) | EW breaking |
| Z⁰ boson | 1 | 0 | 91187.6 (≈91.2 GeV) | EW breaking |
| Photon γ | 1 | 0 | 0 | Unbroken U(1)_em |
| Gluino g̃ (8) | 1/2 | 0 | soft-dep. | SU(3)_c adjoint |
| Wino W̃± | 1/2 | ±1 | soft-dep. | SU(2)_L adjoint |
| Wino W̃⁰ | 1/2 | 0 | soft-dep. | SU(2)_L adjoint |
| Bino B̃ | 1/2 | 0 | soft-dep. | U(1)_Y adjoint |

---

## Degree of Freedom Count

### Confined sector:
| Field | Scalars (real) | Weyl fermions |
|-------|---------------|--------------|
| M^i_j diagonal (3 complex) | 6 | 3 (+ 1 constrained away by ψ_X) |
| M^i_j off-diagonal (6 complex) | 12 | 6 (3 Dirac pairs) |
| X (no kinetic) | 0 | 0 (auxiliary → Lagrange multiplier) |
| B, B̃ | 4 | 2 (1 Dirac) |
| **Total confined** | **22 real scalars** | **11 Weyl = 5.5 Dirac** |

### Dimensional Analysis

The confined SQCD meson M^i_j is a composite field with [M] = mass². Key units:
- [C] = (Λ⁶ m³)^(1/3) = MeV³,  [M_vev] = MeV²,  [X₀] = MeV⁻³
- W entries: [W_{M_ab,M_ba}] = |X₀| × M_c ~ MeV⁻¹ (not MeV as for elementary fields)

The Kahler potential for mesons in the confined phase is K = (1/Λ²) Tr(M† M),
giving a non-canonical kinetic term. The physical (Kahler-normalized) masses are:

    m_physical = Λ² × (raw W entry)   [converts MeV⁻¹ → MeV]

Off-diagonal mesino physical masses:
- m(ψ_ud) = Λ² × |X₀| × M_ss = C²/(Λ⁴ m_s) = **1.029 MeV**
- m(ψ_us) = Λ² × |X₀| × M_dd = C²/(Λ⁴ m_d) = **20.579 MeV**
- m(ψ_ds) = Λ² × |X₀| × M_uu = C²/(Λ⁴ m_u) = **44.493 MeV**

Scalar masses (diagonal ≈ 92 MeV, off-diagonal ≈ 130 MeV) are already physical
because V_soft = f_π² Tr(M†M) contains the correct Kahler normalization implicitly:
the soft mass is defined relative to the canonical kinetic term for the physical scalar.

The dominant mass scale in the confined sector is f_π = 92 MeV (scalar) and
C²/Λ⁴m_i (fermion). The scalar-fermion splitting is large, indicating explicit
SUSY-breaking by V_soft is the dominant effect in this sector.

---

*Generated by particle_spectrum.py*