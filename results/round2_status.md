# Verification Status — All Tasks Complete

## Gap Resolution Summary

### RESOLVED by computation (16/16)

| Item | Task | Result |
|-----|-------|--------|
| Error propagation on Q | Numerical verification | |Q−3/2| = 0.91σ (dominated by δm_τ) |
| ppm denominator (7 vs 33 ppm) | Numerical/Statistics | **35 ppm** relative to 2/9. Initial computation was wrong. |
| Chain stall value | Extended chain | 0.03465 MeV confirmed |
| Cyclic echo | Numerical verification | **Exact** to machine epsilon. Chain is cyclic. |
| M_W error propagation | Numerical verification | 80.374 ± 0.002 GeV, **0.39σ** from PDG |
| Full ranked triplet table | Extended chain | All 84 triplets ranked. Rank 1 gap: 450× |
| Second scaling test | Extended chain | **FAILS**. M₂=9M₀ gives (478, 92, 16377) MeV — no match |
| 4-branch catalog | Extended chain | Complete: Step 1 has 4 branches, physical chain confirmed |
| δ₀ ≈ 2/9 significance | Statistics | **3.1σ** (look-elsewhere), **4.3σ** (specific) |
| Look-elsewhere for Q | Statistics | p = 0.24% = **2.8σ** after scanning 84 triplets |
| Anomaly cancellation | Statistics | A(24⊕15⊕15̄) = 0 ✓. **Anomaly-free.** |
| SU(5) → Diophantine mapping | SU(5) accounting | Mapped. rs=6 ✓, r²+s²−1=12 ✓, r(r+1)/2=6 partial |
| Neutral census | SU(5) accounting | **12 neutrals** (not 16). Corrects initial error. |
| 12 vs 16 discrepancy | SU(5) accounting | Resolved: census was wrong, 12 is correct |
| Charge normalization | SU(5) accounting | Q = Y + T₃ with (α,β) = (1,1) confirmed |
| Quartic uniqueness | Quartic analysis | **NOT unique** in 4-param space (codim 1). Unique under 3 structural constraints. |

### CORRECTIONS to initial computation

1. **Charge census was WRONG**: |Q|=0: 12 (not 16), |Q|=1/3: 2×6 (not 2×4), |Q|=2/3: 2×9 (not 2×10), |Q|=1: 2×2 (not 2×4), |Q|=4/3: 2×3 (not 2×1). Missing |Q|=2: 2×1
2. **±4/3 states misidentified**: They live in (3,2)₋₅/₆ and (3̄,2)₊₅/₆ of the **24** (adjoint), not in (1,3)₁ of the 15
3. **"7 ppm" is actually 35 ppm** (relative to 2/9)
4. **Second scaling fails**: M₂ = 9M₀ does NOT generate (c,b,t). The scaling stops at one step.

### New results not in initial computation

1. **|Q−3/2| is < 1σ** from PDG uncertainties (entirely from δm_τ)
2. **M_W prediction matches to 0.39σ**; disfavors CDF-II at 6.2σ
3. **Cyclic echo is algebraically exact** (not approximate)
4. **SU(5) with 24⊕15⊕15̄ is asymptotically free** (b₀ = 31/3)
5. **SO(32) → SU(16) × U(1) → SU(5) × SU(3) × U(1)²** chain verified with full U(1) charges
6. **Statistical significance**: Koide coincidence = 2.8σ (look-elsewhere), δ₀≈2/9 = 3.1σ

### R from quartic: codimension analysis
In the 4-parameter family x² + f(C)x + g(C) = 0:
- R = 0.2231 is a codimension-1 hypersurface (NOT isolated)
- But the specific polynomial is pinned by 3 constraints: f₀=0, g₀=0, g₁=−f₁
- These constraints have physical content: (1) no spin-independent mass, (2) massless spin-0, (3) coefficient antisymmetry
- R is then a CONSEQUENCE, not a free parameter

---

## Remaining COORDINATOR gaps (user must address)

### Critical for referee defense
- ⚠️ RENORMALIZATION SCHEME: On-shell sin²θ_W = 0.22306 matches R. MS-bar = 0.23122 is 36σ off. Why is on-shell the natural comparison?
- ⚠️ DYNAMICAL LINK: No mechanism connecting SO(32) to the Casimir quartic. Frame as independent observations?

### Physics interpretation
- Physical derivation of Diophantine equations (why rs, not 2rs?)
- Factor-of-two conventions (real vs complex d.o.f.)
- Why 24 ⊕ 15 ⊕ 15̄ and not 24 ⊕ 10 ⊕ 10̄?
- Asymptotic freedom bound derivation
- 291→270 extra states from SO(32): projection mechanism?
- Chain uniqueness (which two chains?)
- (10,6) colour-sextet diquarks
- "Simplest polynomial" → derive from first principles (can use quartic structural constraints)
- 122.4 GeV vs M_H = 125.25 GeV spin mismatch
- Role of SU(3) colour in Casimir formula
- M_q = 3M_ℓ connection to Casimir scale m = 106.6 GeV
- Renormalization scale mixing in quark masses
- Chain A vs B discrepancy (m_c = 1357 vs 1360)
- PDG vs predicted masses distinction

### Housekeeping
- M₀ significant figures vs published EPJC paper
- PDG mass table with edition year
- de Vries reference citation
- Second solution (7,4,28) — mention or ignore?
- Angle between lepton and quark Koide vectors

---

## Tasks Completed

| Task | Description | Status | Key output |
|-------|------|--------|------------|
| Numerical verification | Error propagation, ppm, stall, echo, M_W | ✅ | Q < 1σ, M_W at 0.39σ |
| Extended mass chain | Full triplet table, scaling, catalog | ✅ | Full 84 table, second scaling fails |
| Statistical significance | Look-elsewhere, Monte Carlo | ✅ | 2.8σ (Q), 3.1σ (δ₀≈2/9) |
| SU(5) rep accounting | Census, anomaly, charge operator | ✅ | Census corrected, anomaly-free |
| SO(32) U(1) charges | Full decomposition | ✅ | Full 496 decomposition with charges |
| Quartic uniqueness | Codimension analysis | ✅ | Codim 1 in 4-param; unique under 3 constraints |

## Results Files

Initial computation: agent1–8, assembly_bootstrap, assembly_mass_formula, assembly_group_theory
Verification: r2a_numerical, r2b_extended_chain, r2c_statistics, r2d_SU5_accounting, r2e_SO32_charges, r2f_quartic_uniqueness
Plan: round2_plan.md (original), round2_status.md (this file)
