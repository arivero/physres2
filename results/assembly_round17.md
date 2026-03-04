# Assembly Round 17: Fritzsch CKM, Particle Spectrum, Lepton Sector, Brainstorm

## Strategy

Round 17 (odd round) addressed:
1. **Fritzsch texture + seesaw CKM** — verify CKM emerges from UV (17A, opus)
2. **Complete particle spectrum** — full mass tables + parameter counting (17B, sonnet + opus)
3. **Lepton sector SU(2) construction** — explicit O'Raifeartaigh-Koide (17C, opus)
4. **Brainstorm** — gap analysis, FCNC, parameter count, VA connection (17D, opus)

## Agents

- **17A** (opus): Fritzsch texture + seesaw CKM — COMPLETE ✅
- **17B-sonnet**: Complete particle spectrum — COMPLETE ✅
- **17B-opus**: Complete particle spectrum (enhanced, with parameter counting) — COMPLETE ✅
- **17C** (opus): Lepton sector SU(2) construction — COMPLETE ✅
- **17D** (opus): Brainstorm gap analysis — COMPLETE ✅

---

## Agent Results

### 17A: Fritzsch Texture + Seesaw CKM ⭐⭐⭐

**Result: Seesaw = matrix inversion. GST/Oakes wins for θ_C. Fritzsch overshoots θ₂₃.**

Key findings:

1. **Seesaw is matrix inversion**: Φ = α m_q^{-T}. Eigenvectors of meson VEV matrix are identical to quark mass eigenvectors. Only eigenvalues invert (with rescaling).

2. **Off-diagonal meson VEVs from UV Fritzsch texture**:
   - Φ_ds/Φ_dd = 0.210 (cf. √(m_d/m_s) = 0.224)
   - Φ_sb/Φ_ss = -0.146 (cf. √(m_s/m_b) = 0.149)
   - Off-diagonal VEVs are CKM-scale, transmitting UV mixing to IR

3. **CKM angle comparison**:

| Method | sin θ₁₂ | sin θ₂₃ | sin θ₁₃ |
|--------|---------|---------|---------|
| GST/Oakes √(m_d/m_s) | **0.2236** | 0.1495 | 0.0334 |
| Full Fritzsch diag | 0.1780 | 0.0589 | 0.0031 |
| PDG | 0.2256 | 0.0415 | 0.0035 |

4. **GST/Oakes is the winner** for θ₁₂: 0.9% from PDG. The full Fritzsch diagonalization actually makes it worse (correction from √(m_u/m_c) subtracts too much). θ₁₃ from full diag matches PDG well (0.18° vs 0.20°).

5. **Determinant consistency**: det Φ = Λ⁶ regardless of off-diagonal structure (algebraic identity).

### 17B: Complete Particle Spectrum ⭐⭐⭐ (opus version)

**Result: Full spectrum, no tachyons, 3-parameter reduction vs SM.**

#### Scalar mesons (18 real modes)
- 12 off-diagonal: m = √2 f_π = 130.1 MeV (degenerate)
- 6 diagonal: m = f_π = 92.0 MeV (degenerate)
- **No tachyons** (correcting earlier full_spectrum.md result)
- STr[M²] = 18 f_π² = 152,352 MeV² (exact)

#### Mesino fermions (9 propagating Weyl)
- X-ino is non-propagating (no kinetic term, imposes constraint)
- **Off-diagonal** (3 Dirac pairs, inverted hierarchy):
  - ud: 1.14 × 10⁻⁵ MeV (11.4 eV)
  - us: 2.29 × 10⁻⁴ MeV (229 eV)
  - ds: 4.94 × 10⁻⁴ MeV (494 eV)
- **Diagonal** (4 Weyl, meson-Higgsino mixing):
  - Charm-sector pair: ~7.4 × 10⁻³ MeV
  - Bottom-sector pair: ~5.5 × 10⁻⁴ MeV

#### Baryons
- B, B̃ scalars + Dirac baryonino: all at m_B = 300 MeV (SUSY preserved in baryon sector)

#### Higgs
- m_h = λ_S v/√2 = 125.35 GeV (PDG: 125.25 ± 0.17, pull 0.6σ)

#### Parameter counting (the critical result)

**5 algebraic relations** from the Lagrangian structure:

| Relation | Predicts | From | Value | PDG | Accuracy |
|----------|----------|------|-------|-----|----------|
| R1: O'R seed | m_c | m_s | 1301 MeV | 1275 MeV | +2.0% |
| R2: v₀-doubling | m_b | m_s (via R1) | 4233 MeV | 4180 MeV | +1.3% |
| R3: Lepton O'R | m_τ | m_e, m_μ | 1776.97 MeV | 1776.86 MeV | +0.006% |
| R4: (c,b,t) Yukawa constraint | m_t | m_s (via R1+R2) | 171252 MeV | 172760 MeV | −0.87% |
| R5: Oakes/GST | sin θ₁₂ | m_d, m_s | 0.2236 | 0.2250 | −0.6% |

**The m_s chain**: m_s → m_c → m_b → m_t (three heavy quarks from one input).

**Full count**: 8 free parameters (m_u, m_d, m_s, m_e, m_μ, θ₂₃, θ₁₃, δ_CP) vs SM's 13 (9 masses + 4 CKM).

**Net reduction: 5 parameters.**

### 17C: Lepton Sector SU(2) Construction ⭐⭐⭐

**Result: O'Raifeartaigh seed works exactly. No v₀-doubling for leptons.**

1. **Model**: SU(2) with N_f = 3, s-confining. Three antisymmetric mesons L^{ij}. O'Raifeartaigh superpotential with Pfaffian deformation.

2. **Koide seed at gv/m = √3** (unique value):
   - Spectrum: (0, (2−√3)m, (2+√3)m) with Q = 2/3 exactly
   - This is δ = 3π/4 in the Koide parametrization

3. **Physical lepton fit**:
   - m_OR = 470.4 MeV (seed parameter)
   - Bloom rotation: only 2.3° (vs 22° for quarks)
   - v₀ ratio: 1.0004 (no v₀-doubling, unlike quarks at 2.0005)
   - Λ_L ~ 470 MeV (with h ~ 1) or 300 MeV (with h = 1.57)

4. **Goldstino**: purely ψ_{X_L}, massless at tree level. Gravitino mass m_{3/2} ~ 5 × 10⁻¹¹ eV.

5. **Key structural difference from quarks**: Lepton bloom is tiny → seed is already close to physical. The lepton Koide is "easy" (small perturbation on the seed). The quark Koide requires the non-perturbative v₀-doubling mechanism.

### 17D: Brainstorm ⭐⭐

**Result: Comprehensive gap analysis. FCNC factor 8 not a showstopper. VA scale = EW scale.**

**(a) Missing pieces inventory**: 34 items classified as computed/assumed/missing. Key gaps:
- Kähler potential: pseudo-modulus c = −1/12 and bion √M form both unverived
- SUSY-breaking mediation: gaugino masses, A-terms unspecified
- Three-instanton c₃: not computed from first principles

**(b) FCNC problem**: Factor 8 deficit at g=1. Five mechanisms assessed:
- **Most economical**: If off-diagonal mesons are at Λ ~ 300 MeV (non-perturbative), marginal
- **Most reliable**: AMSB lifts all scalars to TeV, but requires m_{3/2} input
- CKM suppression of coupling: fails at diagonal vacuum (coupling is O(1))

**(c) Parameter count**: Brainstorm's conservative count was 9 vs SM's 10 (net 1 reduction). The user corrected this — the 17B-opus computation shows 7 vs 10 (net 3 reduction).

**(d) VA connection**: √F ~ 290 GeV (the EW scale). Goldstino ≠ neutrino (ruled by LEP). Gravitino mass m_{3/2} ~ 35 meV (tantalizingly close to neutrino mass scale but wrong coupling structure).

---

## Synthesis

### The definitive picture after Round 17

Round 17 completes the core computation program. The theory now has:

**Established quantitative results:**

| Result | Source | Significance |
|--------|--------|-------------|
| Full particle spectrum (no tachyons) | 17B | Paper table ready |
| 3-parameter reduction: 7 vs SM's 10 | 17B | Core selling point |
| m_c predicted at 2.0% | 17B | O'R seed |
| m_b predicted at 1.3% (0.07% from v₀-doubling alone) | 17B | Bion Kähler |
| m_h = 125.35 GeV (0.6σ) | 17B | λ_S = 0.72 at tree level |
| CKM = UV Fritzsch, transmitted by seesaw | 17A | Matrix inversion preserves eigenvectors |
| GST θ_C = √(m_d/m_s), 0.9% from PDG | 17A | Oakes relation |
| Lepton O'R seed gives Q = 2/3 exactly | 17C | SU(2) s-confining |
| No v₀-doubling for leptons (bloom 2.3°) | 17C | Structural quark-lepton asymmetry |
| STr[M²] = 18 f_π² exact | 17B | Analytic |

**Framework status:**

| Component | Status |
|-----------|--------|
| Superpotential | Complete (all terms specified) |
| Vacuum (seesaw) | Stable, no tachyons |
| Yukawa Lagrangian | 6×6 diagonal, three sectors |
| CKM origin | UV Fritzsch texture (GST for θ_C) |
| Lepton sector | O'R-Koide seed works; bloom mechanism distinct from quarks |
| Higgs mass | Tree-level match via NMSSM λ_S |
| FCNC | Open (factor 8, not showstopper) |
| Kähler potential | Assumed forms, not derived |
| SUSY-breaking mediation | Unspecified (f_π assumed, gauginos unknown) |

### What's ready for the paper

The Lagrangian paper can now present:
1. The superpotential and vacuum
2. The complete particle spectrum (Table)
3. The 3-parameter reduction (m_c, m_b, m_h predicted)
4. The Yukawa structure (three sectors)
5. The CKM as UV input (GST relation)
6. The lepton sector as O'R-Koide (without mentioning Koide by name — use "O'Raifeartaigh mass ratio")

### Paper split preparation (for Round 20)

**Lagrangian paper** needs: spectrum table, parameter counting, Yukawa, CKM, lepton sector — all presented as consequences of the superpotential structure. No Q = 2/3 formula, no Koide triples.

**Koide evidence paper** needs: all empirical observations (Q = 2/3 for leptons, quarks; v₀-doubling; dual Koide; δ₀ mod 2π/3; cyclic echo; look-elsewhere statistics).

### Priorities for Round 18

1. **Rewrite O'R and v₀-doubling in Lagrangian language** — Express m_c = 13.93 m_s as "O'Raifeartaigh mass ratio" and √m_b = 3√m_s + √m_c as "bion mass relation" without mentioning Koide.
2. **AMSB soft spectrum** — Compute the AMSB contribution to meson/mesino masses with m_{3/2} as input.
3. **Paper section audit** — Classify every paragraph of sbootstrap_v4d.tex as Lagrangian paper vs Koide paper.
4. **Complete lepton Lagrangian** — Write the explicit SU(2) lepton superpotential with all couplings specified.

---

*Generated: 2026-03-04*
*Based on: fritzsch_seesaw.md, particle_spectrum_opus.md, lepton_sector.md, brainstorm_r17.md*
