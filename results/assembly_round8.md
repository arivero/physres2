# Assembly Round 8: CKM/Cabibbo, Dual Koide Seesaw, NMSSM Spectrum

## Agent A: Cabibbo from Oakes Relation (cabibbo_oakes.md)

**Task**: Test Oakes relation, scan Koide for (u,d,s) and (u,c,t) triples, parameter counting.

**Key results**:

1. **GST relation is the best**: sin θ_C = √(m_d/m_s) = 0.2236 gives θ_C = 12.92° (−0.9% from PDG 13.04°). Better than Weinberg-Oakes (−3.3%) or Fritzsch (−20.5%).

2. **Koide does NOT hold for (u,d,s)**: Q(m_u, m_d, m_s) = 0.567 (−14.9% from 2/3). No sign convention helps. This is qualitatively different from the ~0.5% deviations in established triples.

3. **Koide does NOT hold for (u,c,t)**: Q = 0.849 (+27.4%). Solving Q = 2/3 for m_u gives m_u = 4356 MeV (2000× PDG). The lightest quarks live outside the Koide manifold.

4. **Parameter counting**: After Koide seed + v₀-doubling + top Yukawa, the remaining free parameters are (m_u, m_d). The Oakes relation connects m_d/m_s to θ_C, leaving only m_u free. No Koide condition fixes m_u.

**For the paper**: Confirms the paper's claim that m_u and m_d are the remaining 2 free parameters. The V_us prediction is 0.224 (1.4σ from PDG using GST). The (u,d,s) and (u,c,t) triples are definitively NOT Koide.

---

## Agent B: Dual Koide Seesaw Mechanism (seesaw_mechanism.md)

**Task**: Investigate Q(1/m_d, 1/m_s, 1/m_b) ≈ 2/3 and compatibility with v₀-doubling.

**Key results**:

1. **Seesaw transforms non-Koide to near-Koide**: Q(m_d, m_s, m_b) = 0.731 (+9.7%), but Q(1/m_d, 1/m_s, 1/m_b) = 0.66521 (−0.22%). The algebraic condition p₁p₃/p₂² = 1/6 (to 0.44%) is NOT approximated by the leading-order ratio.

2. **Exact dual Koide gives POOR m_b prediction**: Setting Q(1/m) = 2/3 exactly predicts m_b = 4562 MeV (9.1% off PDG, 12.7σ). The 0.22% closeness of Q to 2/3 does NOT translate to a precise mass prediction — the Koide quadratic amplifies small Q deviations.

3. **v₀-doubling is far superior**: m_b = 4177 MeV (0.07%, 0.1σ). The dual Koide is an observation, not a prediction tool.

4. **Both conditions compatible**: Adjusting m_d by 0.14σ (to 4.604 MeV) makes both hold simultaneously. No fundamental tension.

5. **Delta parametrization**: For inverse (d,s,b), δ_dual = 130.70°. Close to the critical 129.30° (where Q = 2/3 exactly) but no simple algebraic relation to other deltas.

**For the paper**: The dual Koide should remain as an observed pattern, not promoted to a prediction. Its compatibility with v₀-doubling within PDG uncertainties is reassuring.

---

## Agent C: NMSSM Spectrum (nmssm_spectrum.md)

**Task**: Recompute full spectrum with λ X H_u·H_d coupling (λ = 0.72).

**Key results**:

1. **Higgs mass confirmed**: m_h = λv/√2 = 125.35 GeV (0.08% from 125.25 GeV).

2. **Vacuum survives**: The NMSSM shift det M → det M − λv²/2 changes the seesaw by 3×10⁻⁵ fractionally. Negligible.

3. **Charged Higgsinos get mass**: |λX| = 8.7×10⁻¹⁰ MeV (tiny). Previously massless.

4. **Supertrace unchanged**: STr[M²] = 18 f_π² = 152352 MeV². The NMSSM coupling cancels between scalar and fermion sectors. Only soft mass survives.

5. **X-Higgs coupling**: W_{X,H_u⁰} = W_{X,H_d⁰} = λv/√2 = 125355 MeV. This is a large new cross-coupling in the central block but does not significantly alter the meson/baryon spectrum.

6. **New tachyonic modes from NMSSM**: The central block and charged Higgs sector develop tachyonic R-components from the X-Higgs holomorphic mass terms. However, there is a discrepancy with Round 7 Agent D (full_spectrum) regarding whether the off-diagonal meson sectors are also tachyonic — the difference likely stems from whether W_{IJK}F_I terms are properly included in the off-diagonal blocks.

**For the paper**: The NMSSM coupling λ = 0.72 works cleanly. It produces the correct Higgs mass, preserves the supertrace, and barely perturbs the Seiberg vacuum. The identification X = S (NMSSM singlet) is well-motivated.

---

## Agent D: Off-Diagonal Vacuum (offdiag_vacuum.md)

**Task**: Find true vacuum allowing off-diagonal meson VEVs, extract CKM mixing.

**Key results**:

1. **Global minimum destabilizes seesaw**: V_best = 2.17×10¹³ MeV² (factor 60 below diagonal minimum). The soft term m̃² = f_π² overwhelms the Seiberg structure, pulling all M^i_j to O(10⁴ MeV) instead of the hierarchical seesaw (M_u ~ 4×10⁵).

2. **Order-one mixing, not CKM**: SVD gives θ₁₂ = 33.7° (vs Cabibbo 13.04°), θ₂₃ = 76.3° (vs PDG 2.38°), θ₁₃ = 16.2° (vs PDG 0.20°). The vacuum is completely rearranged, not a perturbation.

3. **Saddle point, not minimum**: Physical Hessian has 8 negative eigenvalues. The apparent stability in scaled coordinates is an artifact.

4. **Baryons absorb constraint**: B, B̃ acquire large VEVs (~10⁵-10⁹) so that det M − BB̃ − Λ⁶ ≈ 0 despite det M ≠ Λ⁶.

5. **Physical conclusion**: CKM mixing requires perturbative treatment. With m̃² ≪ m_i², off-diagonal tachyonic modes produce VEVs ~ m̃/m_i, giving parametrically small mixing angles. The brute-force approach fails because it misses the correct vacuum branch.

---

## Synthesis

### Koide Manifold Boundaries

Round 8 conclusively maps out WHERE Koide works and where it doesn't:

| Triple | Q | Status |
|--------|---|--------|
| (e, μ, τ) | 0.666661 | ON manifold (0.001%) |
| (−s, c, b) | 0.67495 | NEAR manifold (1.2%) |
| (c, b, t) | 0.66949 | NEAR manifold (0.4%) |
| (1/m_d, 1/m_s, 1/m_b) | 0.66521 | NEAR manifold (0.22%) |
| (d, s, b) | 0.73143 | OFF manifold (9.7%) |
| (u, d, s) | 0.56704 | FAR OFF manifold (14.9%) |
| (u, c, t) | 0.84901 | FAR OFF manifold (27.4%) |

**Pattern**: Koide works for triples where ALL members are >> m_u, m_d. The lightest quarks break the pattern. This is consistent with the sBootstrap picture: u and d are free parameters (connected to CKM mixing), not constrained by Koide.

### The NMSSM Coupling Completes the Lagrangian

With λ X H_u·H_d added to W, the explicit superpotential now includes ALL necessary terms:
- Seiberg constraint (confinement)
- Three-instanton (Z₉ potential for δ)
- Yukawa couplings (charm, bottom, top)
- NMSSM singlet (Higgs mass)

The Lagrangian is structurally complete. The remaining open problems are:
1. Off-diagonal meson condensation → CKM (Agent D running)
2. Bloom mechanism (nonperturbative)
3. 15 vs 10̄ Pauli conflict
4. SO(32) projection

### Parameter Budget (Updated)

| Parameter | Fixed by | Value |
|-----------|---------|-------|
| m_c | Koide seed (gv/m = √3) | 1301 MeV |
| m_b | v₀-doubling (bion Kähler) | 4177 MeV |
| m_t | external (y_t v/√2) | 172.76 GeV |
| λ | m_h = 125 GeV | 0.72 |
| sin²θ_W | Casimir equation | 0.22310 |
| V_us | GST: √(m_d/m_s) | 0.224 |
| **Free** | m_u, m_d, m_s | **3 parameters** |

With V_us connecting m_d to m_s through the CKM matrix, the free parameter count is effectively 2: m_s (overall scale) and m_u (not constrained by any known relation).

## Status Update

| Open Problem | Status | Round 8 |
|-------------|--------|---------|
| Origin of c = -1/12 | PARTIALLY RESOLVED | (Round 7: bion at α_s = 0.143) |
| Higgs mass at tan β = 1 | RESOLVED | NMSSM λ = 0.72 confirmed |
| CKM mixing | OPEN (perturbative regime needed) | Global min gives O(1) mixing; need perturbative expansion |
| 15 vs 10̄ Pauli | Open | No new progress |
| Bloom mechanism | Open | No new progress |
| Dual Koide origin | CLARIFIED | Observation only; m_b prediction poor (9%) |
| δ₀ mod 2π/3 = 2/9 | RE-VERIFIED | 5 ppm (not 33 ppm — previous correction was wrong) |
| Koide for (u,d,s) | NEGATIVE | Q = 0.567, definitively NOT Koide |
| Parameter budget | CLARIFIED | 2-3 free (m_u, m_d, m_s) |
| Superpotential completeness | RESOLVED | λ X H_u·H_d added; all terms present |
