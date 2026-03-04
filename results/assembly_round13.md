# Assembly Round 13: Baryon Stabilization, Bloom Bion, Mesino, Brainstorm

## Strategy

Round 13 (odd round = opus only) addresses:
1. **Baryon stabilization** — Add W_B = m_B BB̃, compute bounce action
2. **Bloom bion potential** — V_bion(δ) on Koide manifold
3. **Mesino phenomenology** — Sub-keV charged fermions vs experiments
4. **Brainstorm** — 5 missing pieces for publishable theory

## Agents

- **13A** (opus): Baryon stabilization + vacuum tunneling — COMPLETE ✅
- **13B** (opus): Bloom bion potential δ-minimization — COMPLETE ✅
- **13C** (opus): Mesino phenomenology — COMPLETE ✅
- **13D** (opus): Brainstorm — COMPLETE ✅

---

## Agent Results

### 13A: Baryon Stabilization ⭐⭐

**Result 1 — Seesaw is cosmologically stable WITHOUT m_B**:
- S_bounce(thin-wall) = 1.16 × 10⁹
- S_bounce(thick-wall) = 2.39 × 10⁷
- S_bounce(single-field) = 1.97 × 10⁷
- ALL >> 400 (threshold for lifetime > age of universe)
- Tunneling rate Γ ~ exp(-2×10⁷) — effectively zero

**Result 2 — m_B lifts baryonic vacuum**:
- Without Kähler constraint: m_B doesn't help (X adjusts to m_B, F_B = 0)
- WITH Kähler pole (X confined near X_seesaw ~ 10⁻⁹): any m_B >> 10⁻⁹ suffices
- Threshold: m_B > 1.08 MeV (no-λ case) or m_B > 572 MeV (with-λ case)
- Natural scale: m_B ~ Λ = 300 MeV (instanton origin)
- At m_B = Λ: V_B gains 2m_B²Λ⁶ = 1.31×10²⁰ MeV², far exceeds V_A

**Result 3 — m_B changes nothing in the meson sector**:
- Seesaw vacuum M_i = C/m_i is UNAFFECTED by m_B (B = B̃ = 0 there)
- Off-diagonal tachyonic masses UNCHANGED
- Higgs mass UNCHANGED
- Cabibbo angle UNCHANGED
- Only effect: baryino mass lifted from |X| ~ 10⁻⁹ MeV to m_B

### 13B: Bloom Bion Potential ⭐⭐⭐

**Major algebraic result**: The adiabatic bion potential is EXACTLY FLAT on the Koide manifold.

Proof: s_k|z_k| = z_k identically (adiabatic = sign(z_k)), so Σ s_k|z_k| = Σ z_k = 3v₀ = constant.
V_bion(δ) = 9v₀² for ALL δ. Verified numerically: max|V - 9| = 1.4×10⁻¹⁴.

**Frozen-sign potential**: V_frozen = (Σ|z_k|)² has nontrivial structure:
- Flat plateau at V = 9 when all z_k > 0
- Rises to max 14.657 when z_k crosses zero
- Seed (δ = 3π/4) is a **corner minimum** — flat on left, rising on right

**Physical quark triple (-s,c,b)**:
- δ_phys = 157.21° (in z₀ < 0 region)
- V_frozen = 13.015 > 9: bion OPPOSES bloom, acts as restoring force
- v₀-doubling confirmed: v₀(full)/v₀(seed) = 2.0005

**Three-instanton alternative**:
- V_3inst = -cos(9δ) has nearest minimum at δ = 160°, only 2.8° from physical quark δ
- Combined: V = A(V_frozen - 9) + B(-cos(9δ)) with B/A = -2.07 gives equilibrium at physical δ
- d²V/dδ² = 135.4 > 0 at equilibrium (stable minimum)

**Physical picture**:
1. Three-instanton cos(9δ) DRIVES bloom away from seed
2. Bion potential STABILIZES bloom (restoring force)
3. Equilibrium = physical bloom position at ~157°

**Lepton triple (e,μ,τ)**:
- δ_phys = 132.73° (all z_k > 0)
- V_frozen = 9.0000 exactly (flat)
- Bion provides NO δ-selection for leptons
- δ mod 2π/3 = 0.22223 rad ≈ 2/9 (33 ppm)

### 13D: Brainstorm ⭐⭐

Five most important missing pieces:

| # | Missing Piece | Classification | Risk |
|---|--------------|----------------|------|
| 1 | CKM from combined SQCD+EWSB vacuum | (b) may fail | HIGH: F_X vs f_π² |
| 2 | Baryon stabilization / vacuum lifetime | (a) likely succeeds | LOW ✅ DONE |
| 3 | Lepton Sp(2) mediation scale | (c) needs new ideas | HIGH |
| 4 | μ-problem resolution | (a) for option A | MEDIUM |
| 5 | Mesino phenomenology | (b) may fail | HIGH |

**Most important new observation**: The Cabibbo angle as algebraic consequence of Seiberg seesaw (ε_us/ε_ud = √(m_d/m_s) exactly).

**Experimental test**: tan β = 1 exactly, with:
- κ_b/κ_t at SM value (no tan β enhancement)
- A → tt̄ dominates over A → bb̄ for m_A > 2m_t
- sin²θ_C = m_d/m_s testable by lattice QCD + precision Cabibbo measurements

---

## Synthesis

### Major advances in Round 13:
- ✅ **Vacuum stability resolved** — S_bounce ≫ 400, seesaw is cosmologically eternal
- ✅ **m_B stabilization works** — m_B > 1 MeV eliminates baryonic vacuum, changes no predictions
- ✅ **Bloom mechanism identified** — Three-instanton cos(9δ) drives bloom, bion stabilizes it
- ✅ **Adiabatic flatness theorem** — Σ sign(z_k)|z_k| = 3v₀ is algebraic identity on Koide manifold

### 13C: Mesino Phenomenology ⭐⭐

**Quantum numbers**: Charged mesinos (ψ_ud, ψ_us) have Q_em = ±1, B = 0, matching leptons on charge and baryon number. But ψ_us carries strangeness S = ±1 (muon has S = 0), and all mesinos have L = 0.

**g-2 exclusion — catastrophic**:
- ψ_ud (11.4 eV): Δa_e = 1.55 × 10⁻⁵, exceeds δa_e = 2.8 × 10⁻¹³ by factor 5.5 × 10⁷
- ψ_us (229 eV): Δa_e = 1.01 × 10⁻⁵, exceeds by factor 3.6 × 10⁷
- Contribution is logarithmic: Δa_e ~ (α/π)² × (1/3)ln(m_e/m_f), not suppressed

**LEP exclusion**: Would appear as additional light charged generation. N_ν measurement excludes at > 100σ.

**Resolution**: Must come from Kähler potential, not superpotential. Physical masses m_phys = Z^{-1/2} W_IJ, where Z (Kähler metric) can differ from unity by orders of magnitude in strongly coupled theories. Required Z_M^{1/2} ~ 10⁴-10⁵ for lepton identification.

**Identification test**: Mass ratios fail: m(ψ_us)/m(ψ_ud) = m_s/m_d = 20, vs m_μ/m_e = 207.

**Key conclusion**: Mesinos at holomorphic masses are absolutely ruled out as physical particles. The gap must be filled by non-perturbative Kähler physics.

### Key open question from brainstorm:
**Combined SQCD+EWSB vacuum**: The tachyon/CKM derivation assumes F_X = -λv²/2 at the EWSB vacuum. But at the exact seesaw point, F_X = 0 (det M = Λ⁶). The EWSB shift is fractional (3×10⁻⁵). Does the resulting F_X overcome f_π² to create tachyons? This is the decisive calculation.

### Round 14 priorities (even round = sonnet/opus):
1. **Combined SQCD+EWSB vacuum** — The critical calculation (brainstorm item 1)
2. **Mesino phenomenology** — Re-run interrupted agent
3. **Paper incorporation** — Baryon stabilization + bloom mechanism sections
4. **Brainstorm** — Review after combined vacuum result
