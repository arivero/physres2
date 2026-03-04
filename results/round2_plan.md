# Verification Plan — Gap Triage and Task Assignments

## Summary
40 gaps across 3 assembly documents. Sorted into:
- **COMPUTE** (16): New compartmentalized agents can solve
- **COORDINATOR** (18): User/author must address (physics interpretation, conventions, framing)
- **REFEREE-CRITICAL** (3): Must be airtight before submission; straddle compute + interpretation
- **BIBLIOGRAPHY** (3): Straightforward lookups

---

## COMPUTE — Verification Tasks

### Task A: Numerical Verification Batch
Pure arithmetic, no theory context needed.
- Verify 7 ppm vs 33 ppm: δ₀ mod 2π/3 = 0.222230, residual/value vs residual/2π/3
- Error propagation: PDG mass uncertainties → Q values for all triplets
- Exact value and branch ID for the 0.035 MeV chain stall
- Cyclic echo: does Step 4 branch (+,−) = 1357 MeV exactly equal Step 1 output?
- Error propagation: δM_Z = 0.0021 GeV → δM_W from Casimir formula

### Task B: Extended Mass Chain
- Full ranked table of ALL 84 triplets (ranks 2–8, not just rank 1)
- Second scaling: set M₂ = 3M₁ = 9M₀, δ₂ = 3δ₁ = 9δ₀. Does this produce (c,b,t)? If not, document the failure explicitly.
- Complete 4-branch catalog at each chain step (all sign combinations, not just Q=3/2 ones)

### Task C: Look-Elsewhere Statistics
- Statistical significance of δ₀ mod 2π/3 ≈ 2/9: Monte Carlo with random lepton masses
- Look-elsewhere correction for Koide scan: how often does a random set of 6 masses produce a triplet with Q this close to 3/2?
- Anomaly cancellation: does 24 ⊕ 15 ⊕ 15̄ of SU(5) have gauge anomalies? Compute A(R) = Tr[T^a {T^b, T^c}] for each rep.

### Task D: SU(5) Rep Accounting
- Explicit SU(5) to Diophantine mapping: which SU(5) → SM components correspond to each Diophantine equation
- Neutral census: 16 neutrals — decompose into contributions from 24 vs 15 ⊕ 15̄
- Neutral state discrepancy (12 vs 16): eq 3 gives 12, census gives 16. Account for the 4 extras (Cartan generators of 24?)
- Charge normalization: verify Q = Y/2 + T₃ normalization against (α,β) = (±1, ±1/3)

### Task E: SO(32) Completion
- Add U(1)² charges to the branching table (needed for electric charge verification)
- Tensor product embedding SU(5) × SU(3) ⊂ SO(32): cite Dynkin classification or construct explicit 32×32 generators

### Task F: Quartic Uniqueness (REFEREE-CRITICAL)
- The quartic M⁴ − M²C₂ + C₁C₂ = 0 — is it the UNIQUE polynomial in (C₁, C₂) of degree ≤ 4 in M that gives real mass eigenvalues with this ratio? Or does a family of quartics exist with different R values? This is pure polynomial algebra.

---

## COORDINATOR — User Must Address

### Physics interpretation (cannot be compartmentalized)
- Physical derivation of Diophantine equations: why rs and not 2rs? Why symmetric pairing? Origin of "−1" in eq 3?
- Factor-of-two conventions: real vs complex d.o.f. bridge
- Why 24 ⊕ 15 ⊕ 15̄ and not 24 ⊕ 10 ⊕ 10̄?
- Asymptotic freedom bound: which gauge group gives n_f < 33?
- 291 "extra" states from SO(32) → SM: orientifold projection? String-scale masses?
- Chain uniqueness: which two branching chains, agreement at what level?
- (10,6) colour-sextet diquarks: phenomenology or decoupling?
- "Simplest polynomial" — replace with derivation from first principles
- 122.4 vs 125.25 GeV eigenvalue: spin-0 Higgs vs spin-1/2 Casimir state. Flag or explain.
- Role of SU(3) colour in Casimir formula
- Does M_q = 3M_ℓ connect to Casimir internal scale m = 106.6 GeV?
- Renormalization scale mixing: m_s(2 GeV), m_c(m_c), m_b(m_b), m_t(pole). Common scale needed?
- Chain A (m_c=1.357) vs Chain B (m_c=1.360) discrepancy: explain as different inputs
- PDG vs predicted masses distinction: Rank 9 uses PDG masses (Q≈1.482) vs chain uses predicted (Q=3/2 exact). Distinguish clearly.
- Angle between lepton and quark Koide vectors — include or not?

### Referee-critical interpretation
- ⚠️ RENORMALIZATION SCHEME: on-shell sin²θ_W = 0.22306 matches R. MS-bar = 0.23122 is 36σ off. Must argue why on-shell is the natural comparison for a tree-level kinematic result.
- ⚠️ DYNAMICAL LINK: Is there ANY mechanism connecting SO(32) string background to the Casimir quartic? If not, paper must frame them as independent observations.

### Housekeeping
- M₀ significant figures: check consistency with published EPJC paper
- Explicit PDG mass table with edition year
- de Vries reference: find full citation
- Second solution (7,4,28): mention or ignore?

---

## Priority Order for Verification Tasks

1. **Quartic uniqueness** — referee-critical, blocks the Casimir section
2. **Numerical verification** — catches errors before they propagate
3. **SU(5) accounting** — resolves the neutral state discrepancy
4. **Statistics** — look-elsewhere is essential for credibility
5. **Extended chain** — completeness, not urgency
6. **SO(32) completion** — table polishing

## Resource Note
6 tasks. Can run 3–4 in parallel without melting the Mac.
