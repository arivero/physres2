# Extended Mass Chain — Results

## Task 1: Full Ranked Table (All 84 Triplets)

| Rank | Triplet | Signs | Q | |Q−3/2| |
|------|---------|-------|---|--------|
| 1 | (e, μ, τ) | +++ | 1.5000138 | 1.38e-05 |
| 2 | (c, b, t) | +++ | 1.4937749 | 6.23e-03 |
| 3 | (e, μ, c) | +−− | 1.4846664 | 1.53e-02 |
| 4 | (τ, u, s) | +++ | 1.5164821 | 1.65e-02 |
| 5 | (s, c, b) | +−− | 1.4815820 | 1.84e-02 |
| ... | ... | ... | ... | ... |
| 84 | (τ, c, b) | +−+ | 0.7008550 | 7.99e-01 |

**Gap between rank 1 and rank 2: factor 450×.** The (e,μ,τ) triplet is dramatically isolated.

Rank 5 is (s,c,b) with signs +−−, i.e., (√s, −√c, −√b). Note: the predicted chain uses (−√s, √c, √b) which by squaring symmetry gives the same Q. This is the "quark Koide" at |Q−3/2| = 0.018.

## Task 2: Second Scaling (M₂ = 9M₀, δ₂ = 9δ₀)

### Lepton fit parameters
M₀ = 313.838 MeV, δ₀ = 2.31663 rad
δ₀ mod 2π/3 = 0.22223 (cf. 2/9 = 0.22222)

### First scaling: M₁ = 3M₀, δ₁ = 3δ₀ → (s, c, b) quarks
(Assignment: k=0→b, k=1→s, k=2→c)

| k | Predicted | Particle | PDG | Dev |
|---|-----------|----------|-----|-----|
| 0 | 4197 MeV | b | 4180 MeV | 0.4% |
| 1 | 92.28 MeV | s | 93.4 MeV | 1.2% |
| 2 | 1360 MeV | c | 1270 MeV | 7.1% |

### Second scaling: M₂ = 9M₀, δ₂ = 9δ₀

| k | Mass | Nearest particle | Match? |
|---|------|-----------------|--------|
| 0 | 477.9 MeV | — | No clear match |
| 1 | 92.26 MeV | s quark (93.4 MeV) | Echoes m_s |
| 2 | 16377 MeV | — | Between b and t |

Q for this triplet: **1.500000000000001** (exact by construction)

**The second scaling FAILS to produce (c, b, t) or any other known triplet.** The k=1 mass echoes m_s again rather than advancing. Documented as a negative result.

## Task 3: Complete 4-Branch Catalog

### Step 1: Q(m_b, m_?, m_t) = 3/2

| Branch | m_? (MeV) | m_? (GeV) | Signs | Note |
|--------|-----------|-----------|-------|------|
| 1 | 1357.0 | 1.357 | +++ | **≈ m_c** |
| 2 | 60280 | 60.28 | +−− | Near M_W? |
| 3 | 1,341,214 | 1341 | +−− | ~1.3 TeV |
| 4 | 3,549,509 | 3550 | +++ | ~3.5 TeV |

### Step 2: From physical branch (m_c = 1357 → m_? → m_b)

| Branch | m_? (MeV) | Signs | Note |
|--------|-----------|-------|------|
| 1 | 92.17 | +−+ | **≈ m_s** (negative √m_s) |
| 2 | 172690 | +++ | **Echoes m_t** (exact) |

### Step 3: From physical chain (m_s = 92.17 → m_? → m_c = 1357)

| Branch | m_? (MeV) | Signs | Note |
|--------|-----------|-------|------|
| 1 | **0.03465** | +++ | **Chain stall** |
| 2 | 1962 | +−− | Echoes τ at 10.4% |
| 3 | 4180 | +−− | **Echoes m_b** (exact) |
| 4 | 34434 | +++ | ~34 GeV |

### Key findings:
- **Stall value confirmed**: 0.03465 MeV
- **m_b echo at Step 3**: the chain remembers its starting point (cyclic structure)
- **Branch 2 of Step 1**: 60.28 GeV is intriguingly close to M_W
- **τ echo**: 1962 MeV from Step 3 branch 2 (10% from m_τ)
- Every chain step produces an exact echo of an input mass, confirming algebraic consistency
