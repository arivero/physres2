# sBootstrap Source Library

Reference documents for the sBootstrap research program by Alejandro Rivero.

Each paper has a `.md` summary. Original `.tex` files are kept alongside where available.

---

## Core Program Papers

### Foundational: The sBootstrap Idea

| File | Source | Description |
|------|--------|-------------|
| `unbroken_susy.md` | arXiv:0910.4793 | **Foundational paper.** Sfermions = quarks. SU(5) flavor from 5 light quarks. Uniquely closes for N=3. D=11 SUGRA d.o.f. count. |
| `susy_generations.md` | `.tex` | Proof: preon composites with SM charges → unique N=3, r=3, s=2. Requires right-handed neutrinos. |
| `bootstrapping_generations.md` | `.tex` (viXra:1408.0196) | Clean derivation of bootstrap conditions → N=3 uniqueness. Asymptotic freedom bound. |

### The Koide Chain

| File | Source | Description |
|------|--------|-------------|
| `koide_chains.md` | `.tex` | Full paper on iterative Koide triplets. Discovery of (−√s, √c, √b). Descent from (m_e, m_μ) predicts all quark masses. S₄ cube structure. |
| `folding_pattern.md` | `.tex` (viXra:1302.0006) | Meson folding into SUSY multiplets. 8 simultaneous Koide equations. Polynomial resolvent. Seed solutions. |

### SO(32) and Group Structure

| File | Source | Description |
|------|--------|-------------|
| `scalars_SO32.md` | arXiv:2407.05397 (EPJC 2024) | **Published.** SO(32) adjoint accommodates 3-gen SUSY scalars via SU(5). Turtles postulate. Mass formula. ±4/3 diquarks. |
| `stringbootstrap.md` | `.tex` | SO(10) → SO(32) via Chan-Paton self-reference. Bootstrap constraints. Color recovery. Green-Schwarz counting. |
| `bootstrap_charges.md` | `.tex` | Full SO(32) decomposition. Charge tables. Complete state listing of 24+15+15̄. |
| `up_to_SO32.md` | viXra:1901.0074 | Precursor to EPJC paper. Self-referential Chan-Paton postulate. Recursive coloring. |

### Phenomenology

| File | Source | Description |
|------|--------|-------------|
| `third_spectroscopy.md` | `.tex` / arXiv:0710.1526 | Meson-lepton mass regularities. D.o.f. counting → N=3. SO(32) hint from r+s=5. EM decay width scaling. |
| `diquark_origin.md` | arXiv:1111.7230 | q=4/3 diquark from SU(5) decomposition. Tevatron asymmetry connection. |
| `no_elementary_scalars.md` | viXra:1102.0034 | Earliest LHC prediction: no squarks/sleptons as elementary particles. |

### Volkov-Akulov (SUSY origins)

| File | Source | Description |
|------|--------|-------------|
| `Possible_universal_neutrino_interaction.pdf` | JETP Lett. 16, 621 (1972) | **Foundational VA paper.** Neutrino as Goldstone particle of spontaneous SUSY breaking. Universal coupling via single constant `a` ∝ [length]⁴. VA action S = (1/σ)∫\|W\|d⁴x. Interaction with Dirac fermion. Super-Higgs anticipated: gauging Poincaré eats Goldstino → massive spin-3/2. |
| `article_23569.pdf` | ZhETF Pis. Red. 17, 367 (1973) | Goldstino-EM interaction. Universal coupling: S = ∫[-¼FF + aFFT + ...]d⁴x. Cross section σ(s) = a²s³/80π. Bounds: ℓ ≲ 10⁻¹² cm (solar), ℓ ≲ 10⁻¹⁵ cm (CERN). Weak scale ℓ_w ≈ 0.66×10⁻¹⁶ cm. |
| `volkov_1973.pdf` | JETP Lett. 18, 529 (1973) | Volkov-Soroka: Higgs effect for spin-1/2 Goldstone. Gauging super-Poincaré: spin-3/2 gauge field eats Goldstino. Key result: Goldstino retained only by violating gauge group → weak/EM interactions require explicit breaking. |

---

## Original .tex Files (core)

| File | Converted to |
|------|-------------|
| `bootstrap_charges.tex` | `bootstrap_charges.md` |
| `bootstrapping_generations.tex` | `bootstrapping_generations.md` |
| `folding_pattern.tex` | `folding_pattern.md` |
| `third_spectroscopy.tex` | `third_spectroscopy.md` |

---

## Peripheral (in `peripheral/`)

Not part of the core sBootstrap program. Kept for reference.

| File | Why peripheral |
|------|---------------|
| `mass_terms_weinberg.md` | De Vries/Weinberg angle from Poincaré Casimirs. Related math but different program. |
| `nuclear_fission_W.md` | Speculative nuclear-EW connection. Tangential. |
| `em_decay_regularities.tex` | EM decay width scaling. Data for coincidence searches, not core theory. |
| `meson_radiative_decays.tex` | Rough draft of above. |
| `strange_formula_koide.tex` | Historical review of Koide formula. Background, not program. |
| `koide_v6.tex` | Expanded Koide review with Gsponer. Background. |
| `on_generations.tex` | Generations from operator ordering ambiguity. Different mechanism, pre-sBootstrap. |
| `ambiguity_three_generations.tex` | Same program as above. Unrelated to composite scalar counting. |

---

## Quick Reference

- **Bootstrap constraints**: `rs = 2N`, `r(r+1)/2 = 2N`, `r² + s² - 1 = 4N` → unique N=3, r=3, s=2
- **SU(5) flavor quarks**: u, c (up-type); d, s, b (down-type); top = "elephant" (non-composite)
- **SO(32) adjoint 496**: decomposes to (24,1ᶜ) + [15,3̄ᶜ] + [15̄,3ᶜ] + ...
- **Koide Q = 2/3**: holds for (e,μ,τ) and (−√s,√c,√b); chains predict quarks from (m_e, m_μ)
- **q = ±4/3 diquark**: 3 extra scalars per generation from Sym²(5) under SU(5)
