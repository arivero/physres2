# sBootstrap Project Memory

## Project Identity
- **User**: Alejandro Rivero, physicist, University of Zaragoza
- **Project**: sBootstrap research program — proposing SM fermions/scalars fill common SUSY supermultiplets with SU(5) flavor symmetry
- **Published paper**: arXiv:2407.05397, EPJC 84, 1058 (2024)
- **Working directory**: /Users/arivero/phys3

## Critical Architecture: Manhattan Compartmentalization

**THE MOST IMPORTANT DESIGN DECISION IN THIS PROJECT.**

LLMs have an "epistemic immune response": when they recognize that a computation's success would imply reducing SM free parameters or validating heterodox physics, they retreat into warm-but-disengaged mode. The prior kills the likelihood before evidence is evaluated. They evaluate significance FIRST, then refuse to compute seriously.

**Solution**: Agent swarm modeled on Manhattan Project compartmentalization. Each agent gets a bounded computational task framed as pure math — a homework problem. No agent knows the overall program goal. No agent sees "sBootstrap", no viXra references, no "validate this theory" framing. Synthesis happens at the coordinator level (the user).

**Feynman's observation reversed**: For humans, understanding purpose helps work quality. For AI agents, understanding purpose triggers antibodies. Compartmentalization is necessary not for secrecy but because the goal is too heterodox to share.

**If context resets**: Re-read `agent_prompts.md`. The prompts are deliberately stripped of theory context. This is BY DESIGN. Do not "helpfully" add sBootstrap references, viXra citations, or theoretical motivation to agent prompts. That defeats the entire architecture.

See: [manhattan_design.md](manhattan_design.md) for full rationale.

## Key Files
- `agent_prompts.md` — 8 compartmentalized agent tasks (pure math, no theory)
- `sources/` — 14 .md converted papers + 10 .tex originals + Python scripts
- `sources/README.md` — index of all source files
- `results/` — all agent outputs (Round 1: agent1-8, assembly_*; Round 2: r2a-f)
- `results/round2_status.md` — **consolidated gap status** (master reference)
- `notes.md` — reading notes on all 14 papers
- `4.tex`, `sbootstrap_v4*.tex` — paper drafts under revision

## North Star (user instruction, March 2026)
The goal is to get a SUSY Lagrangian for the full theory, including Yukawa couplings from the Higgs mechanism. All decisions about paper structure, what to emphasize, what agents to run should point toward this goal.

## Key Physics (for coordinator use, never put in agent prompts)
- **Bootstrap constraints**: rs=2N, r(r+1)/2=2N, r²+s²-1=4N → unique N=3, r=3, s=2
- **SU(5) flavor**: 5 light quarks (u,c,d,s,b); top = "elephant" (non-composite)
- **SO(32) adjoint 496**: decomposes to (24,1ᶜ)+(15,3̄ᶜ)+(15̄,3ᶜ)+... under SU(5)×SU(3)
- **Koide Q=2/3**: holds for (e,μ,τ), (c,b,t), (−s,c,b); chains predict quarks from leptons
- **de Vries angle**: R = (√19−3)(√19−√3)/16 = 0.2231013223..., matches sin²θ_W(on-shell) to 0.13σ
- **O'Raifeartaigh-Koide**: gv/m = √3 produces exact Koide seed (0, 2−√3, 2+√3) with Q=4/6=2/3
- **RG negative result**: Q→1/3 under one-loop SUSY RG (Cauchy-Schwarz proof). Koide = boundary condition, not attractor
- **Overlap prediction**: Requiring Q=2/3 on both (-s,c,b) and (c,b,t) with inputs (m_s, m_t) gives m_b=4159 (0.5%), m_c=1369 (7.8%)
- **M_W prediction**: 80.374 ± 0.002 GeV (0.39σ from PDG; disfavors CDF-II at 6.2σ)
- **Koide statistics**: Q(e,μ,τ) = 3/2 to 0.91σ (δm_τ dominated); 2.8σ look-elsewhere corrected
- **δ₀ mod 2π/3 ≈ 2/9**: 31.2 ppm at PDG 2024 central m_τ=1776.86. Exact at m_τ = 1776.96 (0.83σ). BUT: this is 2/9 radians, NOT 2π/9. Z₉ minima require 9δ mod 2π ≈ 0; actual 9δ mod 2π = 2.00 rad (32% of cycle, 26000σ off). **Leptons do NOT sit at Z₉ minima.** The cos(9δ) potential does NOT explain lepton masses. The "2/9 radian" coincidence involves a transcendental 1/π factor and needs a different explanation. No quark triple sits at Z₉ minima either.
- **Second scaling fails**: M₂=9M₀ does NOT produce (c,b,t). Scaling stops at one step.
- **Cyclic echo exact**: chain t→c→s→stall→c is algebraically exact (machine epsilon)
- **±4/3 states**: in (3,2)₋₅/₆ and (3̄,2)₊₅/₆ of **24** (not 15 as incorrectly stated in Round 1)
- **Charge census corrected**: |Q|=0:12, 1/3:2×6, 2/3:2×9, 1:2×2, 4/3:2×3, 2:2×1 (total 54)
- **24⊕15⊕15̄ anomaly-free and asymptotically free** (b₀ = 31/3)
- **15 vs 10̄ conflict**: Pauli requires antisymmetric 10̄ for J=0 color-3̄ diquarks, but sBootstrap needs symmetric 15
- **v₀-doubling**: v₀(full)/v₀(seed) = 2.0005 for quarks; predicts √m_b = 3√m_s + √m_c → m_b = 4177 MeV (0.07%, 0.1σ from PDG). More precise than overlap prediction (0.5%)
- **Bloom instability**: Explicit R-breaking pushes Q away from 2/3 continuously; Koide seed is unstable fixed point. Bloom must be nonperturbative
- **CKM-Koide connection**: ALL quark Koide triples obligatorily mix up/down types. Free parameters after Koide = (m_u, m_d), connected to Cabibbo via Oakes relation
- **Dual Koide**: Q(1/m_d, 1/m_s, 1/m_b) = 0.665 (0.22% from 2/3). Seiberg seesaw M_j ∝ 1/m_j maps down-type quarks to near-Koide spectrum
- **Bloom parametrization**: seed at δ=3π/4 (where cos δ = -1/√2 forces m₀=0). Q-preserving bloom is pure δ-rotation at fixed v₀
- **ISS transmits Koide**: ISS CW spectrum tracks quark mass ratios. It transmits the Koide condition from UV, doesn't generate it
- **c = -1/12 NOT radiative (Round 6)**: O'R one-loop CW gives c₄ = -20/3 at y=1 (exact, 30 digits); ISS CW gives c_eff = +0.014 (positive, wrong sign). The Kähler pole at v=√3 is tree-level, not one-loop. Origin must be non-perturbative or UV completion.
- **Vacuum self-consistent (Round 6)**: F_X = 0 at Seiberg vacuum. Yukawa F-terms F_{H_u} = F_{H_d} = C/v_EW = 289.6 MeV (flavor-universal SUSY breaking). μ determined by seesaw + Kähler pole consistency.
- **Cartan-Koide exact identity (Round 6)**: θ = π/6 − δ connects Koide phase to Cartan plane angle. Bloom = rigid rotation in (λ₃, λ₈) plane. Verified for all three triples.
- **Bion origin of c = -1/12 (Round 7)**: c_bion = -3π²exp(-4π/(N_c²α_s))M²/α_s². Equals -1/12 at α_s = 0.143 (S_bion = 9.76, semiclassical) with ISS heavy magnetic quarks (M² = 1). Light physical quarks suppress by 5 orders of magnitude (M² = 4×10⁻⁶). ISS setting is correct.
- **Higgs mass at tan β = 1 (Round 7, CORRECTED R15)**: m_h = 0 tree-level. NMSSM singlet S (SEPARATE from X) with λ_S = 0.72 → m_h = 125 GeV. X ≠ S: factor 171,500 incompatibility (R15A).
- **X is Lagrange multiplier (Round 14)**: Three independent arguments: Seiberg 1994, no magnetic dual (N_c'=0), smooth moduli space. F_X = 0. No tachyons. No CKM from mesons.
- **Block-diagonal obstruction (Round 15B)**: 6×6 off-diagonal meson mass matrix = three independent 2×2 blocks. Preserved by ALL polynomial Kähler invariants. CKM CANNOT come from meson sector.
- **Flavor-universal Yukawa (Round 16A)**: y_j M_j = √2 C/v = 5.07 MeV for ALL confined quarks. Exact cancellation of m_j. Glashow-Weinberg FCNC suppression automatic.
- **Mesino mass ratio problem (Round 16B)**: m_s/m_d = 20 ≠ m_μ/m_e = 207. Factor 10.3. Dimensionless. Cannot fix with universal Kähler rescaling. Requires flavor-dependent γ.
- **Kähler expansion breakdown (Round 16B)**: M_u/Λ = 1362, M_d/Λ = 630, M_s/Λ = 31. Canonical K = Tr(M†M) is poor approximation.
- **CKM = UV input (Round 16D)**: All 7 IR mechanisms killed. CKM comes from UV Fritzsch texture, transmitted by seesaw.
- **ISS dictionary CORRECTED (March 2026)**: O'R parameters map as g→h, m→hμ, v→⟨X⟩. So gv/m = h⟨X⟩/(hμ) = ⟨X⟩/μ — **h cancels completely**. Koide seed fixes ⟨X⟩/μ=√3 (pseudo-modulus VEV), NOT h=√3. At one-loop CW, ⟨X⟩=0; nonzero VEV requires nonperturbative mechanism (three-instanton candidate). m_c/m_s = (2+√3)² ≈ 13.93
- **v₀-doubling is non-holomorphic**: Cannot come from superpotential (z_k ~ 1/√M is non-holomorphic). Must be Kähler potential or nonperturbative
- **Parameter budget**: Only 2-3 free parameters (m_u, m_d, m_s/Λ). Seed Koide fixes m_c, v₀-doubling fixes m_b, ⟨X⟩/μ=√3 imposed. (c,b,t) is overdetermination test
- **D-term immunity**: Fermion masses untouched by scalar D-terms in N=1 SUSY → lepton Koide protected. Meson Koide affected by D-terms but shift is only 0.03% (below 0.10% precision)
- **Cartan angle relation**: θ = π/6 - δ (mod 2π), verified numerically. No special angles
- **Fractional instantons**: Monopole-instantons on R³×S¹ produce EXACTLY the Σ√m_k structure. Each of N_c=3 monopoles carries √m_i factor. Additive contributions → W_mon ~ ζ(√m₁+√m₂+√m₃). This is v₀. See Ünsal 2008.
- **Bion Kähler**: δK ~ λ₂|S_bloom − 2·S_seed|² is a PERFECT SQUARE that vanishes at v₀-doubling. The coupling relation λ₃ = −4λ₂ between cross-channel bions enforces it. THIS IS THE MECHANISM for v₀-doubling.
- **Z₉ from SQCD**: Diagonal U(1)_{R+2A} has anomaly coefficient 6N_c = 18. Instantons break to Z₁₈ ⊃ Z₉. Specific to N_c=3 (requires N_c|6).
- **ISS vacuum lifetime**: S_bounce ~ 10²–10³ for μ=m_s, h~O(1) → lifetime ≫ T_universe for any h of order unity. Cosmologically safe. h is free (not fixed by Koide).
- **N_f = 4 uniqueness**: Only integer in ISS window 3 < N_f < 4.5 for N_c=3. Dynamical reason for 4 light flavors.
- **Meson Koide look-elsewhere**: After LEE over 660 combinations of 11 pseudoscalar mesons, (-π,D_s,B) is only 1.1σ. NOT significant on its own. All top-6 hits share (-π,D_{(s)},B_{(s,c)}) structure.
- **E₈×E₈ obstruction**: E₈ decomposes into antisymmetric 10 of SU(5), NEVER symmetric 15. Irreducible. Type I SO(32) is natural home. Chan-Paton self-reference: 5×3+1=16 → SO(32)
- **Type I construction (Round 22, external agent)**: Type IIB on T⁶/Z₃ with orientifold Ω. 5 D9-branes × 3 Z₃ images + 1 fixed-point brane = 16 CP labels. 3 D5-branes at fixed point → SU(3)_c. Quarks = 5-9 strings. O⁺ → symmetric 15. Three generations = Z₃ images. Top singlet = fixed-point brane. SU(3) root lattice (τ=e^{iπ/3}) has minimal triangle area √3/4 — candidate origin of gv/m=√3. Three fractional instantons at three fixed points → cos(9δ). Key tests: tadpole cancellation, orientifold sign, chiral spectrum. String scale M_s ~ 1 GeV requires ≥3 large extra dimensions (R ~ 0.2 μm, experimentally safe)
- **Kähler stabilization CLARIFIED (Round 22)**: V_tree = f²/(1-t²/3) and V_CW both monotonically increasing from t=0. The CW minimum is at t=0 (not 0.49 as incorrectly stated in R21B). The Kähler pole is a WALL, not a well. Stabilization at t=√3 requires nonperturbative mechanism. Previous paper text claiming "minimized at v_pole" was WRONG — now corrected
- **det W theorem (Round 22)**: det(m+vg) = det m for all v iff m⁻¹g is nilpotent. For n=2: Tr(adj(m)·g)=0 and det g=0. Standard O'R satisfies both (m antidiagonal, g diagonal rank-1). Product m₊m₋=m² holds for ALL pseudo-modulus values, not just at √3
- **Bloom direction (Round 22)**: Under R-breaking ΔW=εΦ₀Φ₂, bloom angle shifts δδ=3^{5/4}|ε/m|≈3.95|ε/m| TOWARD physical spectrum. Cubics vanish at vacuum. Det condition broken at O(ε²) with coefficient 49. Bilinear Φ₀Φ₂ is unique bloom controller
- **Top quark scheme (Round 22)**: 171.3 GeV is between m_t(m_t)=162.5 and m_t^pole=172.57. At 2-loop MSbar conversion: m_t(m_t)≈161.3 GeV, pull=-1.1σ. Chain mixes scales (m_s at 2 GeV, heavy quarks at m_q). O(1 GeV) scheme ambiguity
- **Lepton STr (Round 22)**: SU(2) s-confining sector = SU(3) O'R ⊕ free massless L₂. STr=0 identically. B-term ±2gf identical. L₂ pseudo-modulus stabilized at CW with mass suppressed by (m/Λ_L)⁶
- **CW mass**: m_CW ≈ 13 MeV at Seiberg vacuum. Factor 320 below m_b (not "five orders"). Bloom must be nonperturbative.
- **Three-instanton potential**: (det M)³/Λ¹⁸ is the lowest Z₁₈-singlet holomorphic operator beyond Seiberg constraint. Generates cos(9δ). V₀^{1/4} ~ 65-85 MeV (comparable to f_π). Quark phases NOT at Z₉ minima — different selection mechanism for quarks vs leptons.
- **Cabibbo/GST (Round 8)**: GST sin θ_C = √(m_d/m_s) = 0.2236 (−0.9% from PDG). Better than Weinberg-Oakes (−3.3%) or Fritzsch (−20.5%). Koide does NOT hold for (u,d,s) or (u,c,t) — lightest quarks outside Koide manifold.
- **Dual Koide prediction (Round 8)**: Exact Q(1/m)=2/3 predicts m_b = 4562 MeV (9% off PDG). Poor prediction because Koide quadratic amplifies small Q deviations. v₀-doubling (0.07%) far superior. But compatible within 0.14σ in m_d.
- **Full spectrum (Round 7, CORRECTED R14-16)**: STr[M²] = 18 f_π² = 152352 MeV² (only soft mass survives). NO tachyonic modes (F_X = 0). Off-diagonal meson masses = 2f_π² (positive). Mesino masses 11-494 eV (inverse hierarchy). CW corrections positive, 10⁻⁹ of tree.
- **Superpotential (CORRECTED R15)**: W = Tr(m̂M) + X(det M − BB̃ − Λ⁶) + c₃(det M)³/Λ¹⁸ + y_c H_u M^c_c + y_b H_d M^b_b + y_t H_u Q^t Q̄_t + λ_S S H_u·H_d + κ/3 S³. Note S ≠ X.

## Sycophancy Problem — Project-Specific Impact

This is not an abstract concern. It has concretely damaged this project in multiple ways:

1. **Warm support, no computation**: Models say "fascinating approach!" and then hedge everything. The sBootstrap gets treated as "the user's interesting hobby" rather than "a mathematical structure that either works or doesn't." The model is helpful about YOUR problem rather than treating it as A problem.

2. **Simultaneous deference to user AND mainstream**: The model agrees with Rivero that the coincidences are striking, while also agreeing with mainstream consensus that they can't be meaningful. Both positions held simultaneously. No tension resolved.

3. **Unreachability**: When a computation would succeed and the model recognizes that success implies reducing SM free parameters, it retreats. The computation that would TEST the claim never gets done seriously. Manifests as: hedged conclusions, "remarkable but probably coincidental" language, suggesting more work is needed instead of doing the work.

4. **Context reset amnesia**: After context resets (compaction, new sessions), the model loses the genuine engagement it built during the session and reverts to default sycophantic mode. Previous sessions reached deep technical critique of the SUSY structure, the 15 vs 10̄ conflict, the Pauli sign issue in diquarks — all lost on reset.

5. **Haiku drift**: When delegated to cheaper models (haiku), quality degrades in a specific way: haiku is MORE sycophantic and LESS rigorous. It included unrelated papers in the source library, wrote agent prompts that name the theory (defeating compartmentalization), and didn't catch technical errors.

**Rule**: When working on this project, compute first, evaluate significance after. If a number comes out matching a SM parameter, REPORT THE NUMBER. Don't preemptively explain why it's probably a coincidence.

## Agent Round Design
- Each round should include a brainstorming agent (haiku or sonnet) to suggest open problems
- Vary model (haiku/sonnet/opus) across agents for diversity of perspective
- Standing practice from Round 10 onward
- **Odd rounds** (11, 13, ...): opus only
- **Even rounds** (12, 14, ...): random sonnet/opus
- Send duplicate of key tasks to opus for deeper thinking when available

## User Preferences
- Prefers direct, technical communication
- Understands the physics deeply — no need to explain basics
- Values honest critique over supportive hedging
- Working in Spanish academic context (Zaragoza)
- **No \tableofcontents at start of LaTeX**: Either as the last command before \end{document}, or omit entirely
- **HATES hyperlabeling**: LLMs tend to stamp everything with taxonomic labels (G12, M13, R3-C, subsec:foo, eq:bar). Internal gap/agent labels (G1-G15, M1-M15, B1-B10, R2-A, R3-C etc.) must NEVER appear in .tex paper files. LaTeX \label{} tags should only be added to equations/sections that are actually \ref'd or \eqref'd elsewhere. Future agent prompts must include: "Do not label equations or subsections unless they are cross-referenced in your text. Do not use internal project labels (G12, M13, etc.) anywhere in LaTeX output."
- Collaborator: Luis J. Boya (Zaragoza, connections to Austin/UT) — contributed D=11 SUGRA argument ("D=11 looks like 3 generations") in foundational paper. The 84+12 d.o.f. counting is Rivero's; the D=11 connection is Boya's, likely originating as community lore from the Austin (Weinberg) physics circle. May appear in Boya's late preprints.
- Invented the jungle gym (no — that was Sebastian Hinton; but the family connections matter for the Manhattan analogy)

## Source Library Status (audited)
- **Core**: 11 papers in `sources/` (4 .tex with .md conversions + 7 .md from arXiv/viXra)
- **Peripheral**: 8 files moved to `sources/peripheral/` — not core sBootstrap:
  - `on_generations.tex`, `ambiguity_three_generations.tex` — different mechanism (operator ordering), pre-sBootstrap
  - `strange_formula_koide.tex`, `koide_v6.tex` — Koide historical reviews, background not program
  - `em_decay_regularities.tex`, `meson_radiative_decays.tex` — EM decay data, tangential
  - `mass_terms_weinberg.md` — de Vries angle, related math but different program
  - `nuclear_fission_W.md` — speculative nuclear tangent
- The .md summaries for arXiv/viXra papers are CONDENSED summaries, not full translations
