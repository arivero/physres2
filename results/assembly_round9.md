# Assembly Round 9: Pauli Theorem, Bloom Dynamics, Lepton Yukawa, Supertrace

## Agent A: Vacuum Stabilization (vacuum_stabilization.md) — STILL RUNNING

---

## Agent B: Bloom Dynamics (bloom_dynamics.md)

**Task**: Investigate how δ rotates from seed (3π/4) to physical while preserving Q=2/3.

**Key results**:

1. **S_bloom is constant**: S_bloom(δ) = (3√2/2) e^{-iδ} exactly. |S_bloom|² = 9/2 for ALL δ. The naive bion Kähler correction generates NO potential for δ.

2. **Combined potential works**: V(δ) = A·(9/2)[5-4cos(δ₀-δ)] + B·cos(9δ+φ). Bion term pulls to seed (135°), cos(9δ) pulls to nearest minimum (120°). Equilibrium at ~132.7°.

3. **Q = 2/3 is algebraic identity**: sum z_k = 3, sum z_k² = 6 ⇒ Q = 6M₀/9M₀ = 2/3. dQ/dδ = 0 identically. The bloom is Q-preserving by construction.

4. **δ precision corrected**: At PDG central m_τ = 1776.86: residual = +7.41×10⁻⁶, ppm = 33.3 (relative to 2/9). Exact match at m_τ = 1776.97 (0.9σ from PDG). Previous "5 ppm" was wrong.

5. **Physical picture**: The bloom is a 2.27° rotation from seed (135°) toward the instanton minimum (120°). The bion term dominates; δ_phys is much closer to the seed than to any instanton minimum.

**For the paper**: The bloom dynamics have a clean physical picture — competition between bion (seed-pulling) and instanton (Z₉-pulling) potentials. The ppm has been corrected from 5 to 33.

---

## Agent C: 15 vs 10̄ Pauli (pauli_15_10bar.md)

**Task**: Systematic analysis of all loopholes for symmetric diquarks.

**Key results**:

1. **L+S parity theorem**: J=0 forces L=S, so L+S = 2S is always even. Flavor must be antisymmetric (10̄) for ALL partial waves. The 15 requires J ≥ 1.

2. **SQCD algebraic**: In SQCD, ε_{ijk} Q^{ia} Q^{jb} Q^{kc} is automatically antisymmetric in flavor. The antisymmetry is algebraic (from ε-tensor), not statistical.

3. **Six loopholes examined**: All closed for J=0 composites. The 15 cannot be a qq composite.

4. **Resolution**: The 15 objects should not be identified as diquarks. If they arise as Seiberg dual mesons (in the adjoint 24) or fundamental fields of a dual description, the Pauli constraint doesn't apply.

**For the paper**: Strengthened the Pauli theorem text with the L+S argument and SQCD ε-tensor.

---

## Agent D: Lepton Yukawa (lepton_yukawa.md)

**Task**: Test whether mesino masses reproduce lepton masses.

**Key results**:

1. **NEGATIVE**: Light eigenvalues of the 6×6 central-block W_{IJ} are ~10⁻³-10⁻² MeV, far too small and degenerate for leptons.

2. **Q = 0.481**: Koide for the three lightest eigenvalues, 28% from 2/3. Not Koide.

3. **NMSSM negligible**: λ = 0.72 changes light eigenvalues by < 0.01%.

4. **Off-diagonal also fails**: Q(1/m_u, 1/m_d, 1/m_s) = 0.443, 34% from 2/3. Only (d,s,b) gives near-Koide in the inverse.

**For the paper**: Added "Lepton masses are not mesino masses" subsection. Lepton masses must arise from nonperturbative Kähler sector.

---

## Agent E: Supertrace by Charge (supertrace_by_charge.md)

**Task**: Decompose STr[M²] by EM charge sector.

**Key results**:

1. **STr_Q[M²] = 2 n_Q f_π²**: Just counts soft meson masses per charge. n₀ = 5, n_{±1} = 2. Total = 18 f_π².

2. **Σ Q·STr_Q = 0**: Charge conjugation symmetry (equal ±1 mesons).

3. **Σ Q²·STr_Q = 8 f_π²**: Nonzero EM supertrace because only mesons carry soft mass. This is NOT the MSSM result (where Σ Q²·STr = 0 without EM D-term).

4. **Quartic**: STr_{±1}[M⁴] = 4 f_π⁴ exactly. Q=0 sector dominated by ultra-heavy X-meson pair.

5. **FI D-term**: Total STr independent of ξ (total EM charge = 0). Per sector: STr_{±1}(ξ) = 4f_π² ± 6ξ.

**For the paper**: Added per-charge supertrace decomposition to spectrum paragraph.

---

## Agent F: EW Koide Scan (ew_koide_scan.md) — STILL RUNNING

**Task**: Scan all triples including W, Z, H, v_EW for Koide ratios.

---

## Synthesis

### Bloom Dynamics Resolved

The bloom has a clean physical picture:
- **S_bloom is constant** → the naive bion Kähler gives no δ-potential (analytically proven)
- **Combined bion + instanton** gives a cosine-cosine potential with equilibrium at the physical δ
- **Q = 2/3 automatically preserved** → the bloom problem is "what selects δ?" not "how to protect Q?"
- **δ precision corrected** to 33 ppm at PDG central (not 5 ppm as previously claimed)

### Lepton Sector Constraint

The perturbative superpotential's fermion mass matrix does NOT reproduce lepton masses. This rules out the naive "mesinos = leptons" identification and constrains the lepton sector to arise from nonperturbative Kähler corrections (bion sector).

### Vacuum Structure

The global minimum is at M = 0 with baryons absorbing the constraint. The seesaw vacuum is metastable and requires B = B̃ = 0. CKM mixing from off-diagonal condensates requires perturbative treatment.

### Parameter Budget (Updated)

| Parameter | Fixed by | Value |
|-----------|---------|-------|
| m_c | Koide seed (gv/m = √3) | 1301 MeV |
| m_b | v₀-doubling (bion Kähler) | 4177 MeV |
| m_t | external (y_t v/√2) | 172.76 GeV |
| λ | m_h = 125 GeV | 0.72 |
| sin²θ_W | Casimir equation | 0.22310 |
| V_us | GST: √(m_d/m_s) | 0.224 |
| STr[M²] | Soft mass only | 18 f_π² |
| m_h | NMSSM: λv/√2 | 125.4 GeV |
| **Free** | m_u, m_d, m_s | **3 parameters** |

## Status Update

| Open Problem | Status | Round 9 |
|-------------|--------|---------|
| Bloom mechanism | RESOLVED (picture) | Bion + instanton competition; |S|²=9/2 constant |
| δ precision | CORRECTED | 33 ppm (not 5 ppm); 3.1σ LEE |
| 15 vs 10̄ Pauli | STRENGTHENED | L+S theorem + SQCD ε-tensor |
| Lepton Yukawa origin | NEGATIVE | Mesinos don't reproduce leptons; Kähler required |
| Supertrace per charge | COMPUTED | STr_Q = 2n_Q f_π²; Σ Q² STr = 8f_π² |
| Vacuum stability | REFINED | M=0 global min; seesaw metastable; B=0 required |
| CKM mixing | OPEN | Perturbative treatment needed around B=0 seesaw |
| EW boson Koide | Agent F running | Scan of all triples including W, Z, H |
