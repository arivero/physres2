# Phenomenologist — Mass Relations and Iterative Chains

## Status: COMPLETED

## What was written

Inserted ~315 lines into `sbootstrap_v4d.tex` at line 429 as `\section{Mass Relations and Iterative Chains}` (label: `sec:mass_chains`).

### Subsections (7 total):

1. **The lepton Koide formula: how special?** (line 444)
   - Q(e,μ,τ) = 1.500014 ± 0.000015, deviation 0.91σ
   - Table: exhaustive scan of all 84 SM triplets. Rank 1 gap: 450×
   - Monte Carlo look-elsewhere: 2.8σ

2. **The Koide angle parametrization** (line 499)
   - √m_k = √M₀(1 + √2 cos(2πk/3 + δ₀))
   - δ₀ mod 2π/3 = 0.22223 ≈ 2/9, residual 35 ppm (corrected from "7 ppm")
   - Significance: 3.1σ look-elsewhere, 4.3σ specific

3. **The scaling rule: M₁ = 3M₀, δ₁ = 3δ₀** (line 548)
   - Predicts s=92.3 (1.2% off), c=1360 (7.1% off), b=4197 (0.4% off)
   - Predicted triple satisfies Q=3/2 to machine precision with negative √m_s
   - Renormalization caveat included

4. **The second scaling fails** (line 602)
   - M₂=9M₀, δ₂=9δ₀ gives (478, 92.3, 16377) MeV — no SM match
   - Reported as negative result

5. **Iterative chain: descending from t and b** (line 624)
   - Step 1: (b,t) → c = 1357 MeV (4 branches)
   - Step 2: (c,b) → s = 92.17 MeV (negative √m_s natural)
   - Step 3: stalls at 0.035 MeV
   - Step 4: cyclic echo — stall,s → c = 1357 to machine epsilon

6. **Cross-checks and master comparison** (line 687)
   - Master table: Chain A (from t,b) vs Chain B (from leptons) vs PDG
   - Two chains agree on m_s ≈ 92, m_c ≈ 1357-1360 to 0.2%
   - Three negative results given equal prominence

### Features:
- 3 tables with experimental comparisons
- All deviations reported with % and σ values
- Convention footnote for Q=2/3 vs Q=3/2
- Cross-references to existing sections
