# Innovations for Opus

Proposed by a Fable subagent review of this repo.

## 1. Verified calculation library with a paper-claim regression harness

The `calculations/` directory holds 65+ standalone scripts that each hardcode PDG 2024 masses, re-implement the Koide/mass-charge function Q, and emit print-only output that gets hand-copied into the four LaTeX papers — a workflow that has already produced one retraction (the m_t×m_ν = f_π² unit error noted in project memory) and several "initial computation was wrong" corrections in `results/round2_status.md`, and `coordination2.md` explicitly flags "PDG vs predicted distinction: yep we should track that". Opus would refactor the shared pieces into a small package (one PDG constants module with uncertainties, common functions for Q, the O'Raifeartaigh spectrum, the bloom parametrization, chain iteration) and write a pytest suite that recomputes every headline number quoted in `paper_lagrangian.tex` and `paper_experimental.tex` (f_π^SUSY = 91.19 MeV, g_πNN = 13.12, tau partner at 1810 MeV, δ_phys = 159°, m_d = 4.60 MeV, ...) and fails whenever a paper claim and its computation diverge. Scope: refactor plus roughly 30-50 regression tests keyed to labeled equations/tables in the papers.

**Size:** medium

## 2. Complete the CKM matrix from the fixed down-quark texture

Project memory marks "CKM beyond Cabibbo" as the hardest open problem: the orbifold route gives ε₁³ = 0.066 vs PDG 0.042 (56% high), and the 5-zero texture "can fit but needs definite D_d". The repo already fixes the ingredients piecewise — Cabibbo as the 2×2 rotation inside the (u,d)/(c,s) O'Raifeartaigh block, V_cb as bion eigenvalue mixing (`calculations/bloom_vcb_connection.py`, `ckm_orbifold.py`, `full_ckm_from_lagrangian.py`) — so the natural next step is to assemble one definite 3×3 D_d from the bloom eigenvector structure, diagonalize, and produce all nine |V_ij| plus the Jarlskog invariant with a χ² against the PDG/UTfit global fit and proper error propagation. The outcome is decisive either way: a working texture closes the papers' most-flagged gap, and a failed fit sharpens exactly which element the framework cannot reproduce.

**Size:** large

## 3. Joint overdetermination solve for the bloom angle

`paper_lagrangian.tex` (Section "Open problems") shows the physical bloom angle δ_phys = 159° sitting at a saddle of the instanton-plus-bion potential and conjectures it "may be fixed by the overdetermination of the mass chain rather than by a simple potential minimum" — a hypothesis stated in prose that has never been tested numerically. The constraints already exist as separate scripts: the bion radial condition, the three-instanton cos(3nδ) potential (`cos9delta_derivation.py`, `z9_analysis*.py`), the dual self-consistency condition cos 3δ = 5√2/8, the R-breaking instability slope δδ ≈ 3.95|ε| (`bloom_dynamics.py`, `dual_bloom_constraint.py`), and the mass-chain equations. Opus would combine these into one coupled system in (z₀, δ, ε), solve it globally, propagate PDG mass uncertainties, and report whether the constraint intersection lands on 159° with a look-elsewhere-corrected significance — turning the paper's central conjecture into a computed result suitable for a new subsection.

**Size:** medium
