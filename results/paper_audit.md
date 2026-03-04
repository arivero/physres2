# Paper Split Audit: sbootstrap_v4d.tex

**Paper A** ("Lagrangian paper"): SUSY Lagrangian, vacuum structure, particle spectrum, mass predictions as consequences of the O'Raifeartaigh vacuum and bion Kahler corrections. Must NOT mention "Koide" by name, nor display Q = 2/3 as an empirical observation. Mass relations are derived from the superpotential structure.

**Paper B** ("Koide evidence paper"): Empirical mass relations, statistical significance, the Koide formula Q = 2/3 as an observed regularity, look-elsewhere corrections, the de Vries angle, historical observations, framing mass relations as unexplained empirical facts.

**Both**: Content needed by both papers (notation, basic setup, sBootstrap postulate).

---

## Section-by-Section Audit

### Abstract (lines 29--98)

| Classification | **SPLIT** |
|---|---|
| Reason | The abstract mixes Lagrangian results with Koide-empirical language throughout. |
| Modifications for A | Rewrite entirely. Paper A abstract should lead with the Lagrangian, the O'Raifeartaigh superpotential, the bion Kahler correction, the $v_0$-doubling as a vacuum condition, the prediction chain $m_s \to m_c \to m_b \to m_t$, the tan(beta)=1 prediction, and the Diophantine uniqueness. Replace every instance of "Koide ratio/triple/formula" with "mass-ratio condition" or "energy-balance condition" or "O'Raifeartaigh mass relation". Never display $Q = (\sum m)/(\sum\sqrt{m})^2 = 2/3$ as a stand-alone empirical formula. |
| Modifications for B | Paper B abstract should lead with the empirical Koide formula, its statistical significance, the look-elsewhere corrections, the seed triples as observed patterns, and frame the meson Koide triple, de Vries angle, and iterative chains as empirical evidence requiring explanation. |

---

### Section 1: Introduction (lines 103--177)

| Classification | **SPLIT** |
|---|---|
| Reason | Introduces both the sBootstrap postulate (Both/A) and the Koide ratio as an observed formula (B). The three "themes" are: (1) signed charges = Both/A; (2) three Koide tuples = B; (3) Koide condition constrains superpotential = A. |
| Modifications for A | Rewrite to lead with the Lagrangian motivation. Theme 1 (signed mass charges) can stay but reframe: "The mass spectrum of the O'Raifeartaigh vacuum is parametrized by signed mass charges..." Theme 3 (superpotential constraint) is the core of Paper A. Theme 2 (three Koide tuples) should become: "The superpotential produces three mass triples satisfying the energy-balance condition..." The "New results" paragraph should emphasize Lagrangian-derived results. Remove all mentions of "Koide" by name; use "energy-balance condition" or "mass-ratio constraint". |
| Modifications for B | Theme 2 is the core. Keep Koide language. Add the historical context that Paper A omits. |

---

### Section 2: The Koide Formula as Energy Balance (lines 180--228)

#### 2.1 Charges, not masses (lines 184--188)

| Classification | **Both** |
|---|---|
| Reason | Defines $\mathfrak{z} = z_0 + z_i$, $m = \mathfrak{z}^2$. Pure notation needed by both papers. |
| Modifications for A | None; notation is theory-neutral. |

#### 2.2 Derivation of $Q = 2/3$ (lines 189--221)

| Classification | **SPLIT** |
|---|---|
| Reason | The algebra itself (conditions 1 and 2) is needed by Paper A as the "energy-balance condition" derived from the vacuum structure. But the framing as "The Koide Formula" and the boxed $Q = 2/3$ equation presented as an empirical regularity is Paper B. |
| Modifications for A | Rename section: "The Energy-Balance Condition". Present the two conditions as properties of the O'Raifeartaigh fermionic mass matrix. Do not call it "the Koide ratio" or "the Koide formula". Use $Q$ as an internal variable if needed, but present the result as: "The perturbation energies match the base scale: $\langle z_k^2 \rangle = z_0^2$, a necessary consequence of the superpotential constraint $gv/m = \sqrt{3}$." |
| Modifications for B | Keep as is. This is the mathematical heart of the Koide formula. |

#### 2.3 Why one mass can vanish (lines 223--228)

| Classification | **Both** |
|---|---|
| Reason | Structural property of the energy-balance condition relevant to both the Lagrangian (flat direction) and the empirical pattern. |
| Modifications for A | Reframe in terms of flat directions of the O'Raifeartaigh model. |

---

### Section 3: Seed Triples (lines 230--349)

#### 3.1 Analytical condition (lines 234--256)

| Classification | **Both** |
|---|---|
| Reason | Derives $r = 2 - \sqrt{3}$ from the energy-balance condition + one zero. Needed by Paper A (O'Raifeartaigh spectrum) and Paper B (observed pattern). |
| Modifications for A | Present as "the fermionic mass ratio of the O'Raifeartaigh model at $gv/m = \sqrt{3}$" rather than "any Koide triple with one massless member." |

#### 3.2 The quark seed: $(0, s, c)$ (lines 258--293)

| Classification | **B** |
|---|---|
| Reason | Empirical observation that $\sqrt{m_s}/\sqrt{m_c} \approx 2 - \sqrt{3}$. The scale ambiguity discussion and historical Harari note are purely empirical/historical. |
| Modifications for A | Paper A should not present $(0,s,c)$ as an independent empirical observation. Instead, the seed structure should emerge from the superpotential analysis in Sec. 10. If referencing the numerical match, present it as "the predicted mass ratio matches PDG values." |

#### 3.3 The meson seed: $(0, \pi, D_s)$ (lines 295--306)

| Classification | **B** |
|---|---|
| Reason | Purely empirical observation about meson masses. |
| Modifications for A | Omit entirely from Paper A. |

#### 3.4 The lepton triple as an "almost-seed" (lines 308--316)

| Classification | **B** |
|---|---|
| Reason | Empirical observation about lepton masses. |
| Modifications for A | Omit or mention only in passing when discussing the separate lepton confining sector. |

#### 3.5 Scale coincidence (lines 318--349)

| Classification | **B** (primarily), with one result for **A** |
|---|---|
| Reason | The $v_0^2$ table is empirical data. But $M_0 \approx f_\pi^2/27$ is used in Paper A as a structural relation (lepton scale tied to pion decay constant). |
| Modifications for A | Extract $M_0 \approx f_\pi^2/27$ and present it as a scale relation in the Lagrangian context (the SU(2) confining sector sets its scale by $f_\pi$). |

---

### Section 4: A Tale of Two Tuples (lines 351--485)

#### 4.1 From seed to full triple (lines 355--419)

| Classification | **SPLIT** |
|---|---|
| Reason | The seed-to-bloom narrative is central to both papers but with different framings. The $v_0$-doubling relation $\sqrt{m_b} = 3\sqrt{m_s} + \sqrt{m_c}$ is a **Paper A** result (derived from bion Kahler minimum). The "tale of two tuples" literary framing is **Paper B**. |
| Modifications for A | Extract the $v_0$-doubling as a consequence of the bion Kahler correction. Present the seed-to-bloom transition as a vacuum transition, not as an empirical pattern. Remove "Koide triple" language; use "mass triple" or "O'Raifeartaigh triple". |
| Modifications for B | Keep the narrative framing. Reference the Lagrangian derivation from Paper A. |

#### 4.2 Precision summary (lines 420--456)

| Classification | **B** |
|---|---|
| Reason | Table of $Q$ values and deviations from $2/3$ is pure empirical data presented as evidence for the Koide pattern. |
| Modifications for A | Paper A should have its own prediction table (Table 8 / tab:predictions already serves this role). Do not reproduce this $Q$-deviation table. |

#### 4.3 Orthogonality of seed and full triples (lines 458--485)

| Classification | **B** |
|---|---|
| Reason | Geometric analysis of empirical charge-space patterns. |
| Modifications for A | Omit. |

---

### Section 5: The Meson Koide Triple (lines 488--534)

| Classification | **B** |
|---|---|
| Reason | Entirely empirical: the meson triple $(-\pi, D_s, B)$, its precision, and the look-elsewhere correction. |
| Modifications for A | Omit entirely. The meson Koide triple has no Lagrangian derivation in the current framework. |

---

### Section 6: The Third Tuple: $(c,b,t)$ (lines 536--582)

#### 6.1 A tuple without mesons (lines 540--551)

| Classification | **SPLIT** |
|---|---|
| Reason | The observation that $(c,b,t)$ has $Q \approx 2/3$ is Paper B. The fact that the top does not hadronize is both. |
| Modifications for A | Present as: "The heavy-sector mass triple $(c,b,t)$ provides an overdetermination test of the Lagrangian predictions." |

#### 6.2 Diquarks and the top (lines 553--567)

| Classification | **A** |
|---|---|
| Reason | Structural discussion of the diquark sector and how top diquarks fit in. |
| Modifications for A | Remove "Koide triple" language; reframe as "energy-balance constraint on the diquark sector." |

#### 6.3 The $(c,b,t)$ triple has no seed (lines 569--582)

| Classification | **B** (primarily) |
|---|---|
| Reason | Empirical observation that no seed structure exists for this triple. |
| Modifications for A | Mention briefly in the discussion of the heavy sector: "The $(c,b,t)$ mass triple is not generated by the O'Raifeartaigh seed mechanism but by the Yukawa sector." |

---

### Section 7: Mass Relations and Iterative Chains (lines 585--908)

#### 7.1 The lepton Koide formula: how special? (lines 600--660)

| Classification | **B** |
|---|---|
| Reason | Statistical significance of the lepton Koide formula, look-elsewhere corrections, Monte Carlo studies, extended meson scan. Core empirical-evidence section. |
| Modifications for A | Omit entirely. |

#### 7.2 The Koide angle parametrization (lines 662--712)

| Classification | **SPLIT** |
|---|---|
| Reason | The parametrization itself is needed by Paper A (it appears in the Cartan-angle analysis and the bloom mechanism). The phase $\delta_0 \bmod 2\pi/3 = 2/9$ observation and its statistical significance are Paper B. |
| Modifications for A | Present the parametrization as a mathematical tool: "Any triple satisfying the energy-balance condition is parametrized by $(M_0, \delta_0)$..." Omit the $2/9$ observation and its significance. |
| Modifications for B | Keep in full. |

#### 7.3 The scaling rule: $M_1 = 3M_0$, $\delta_1 = 3\delta_0$ (lines 714--769)

| Classification | **B** |
|---|---|
| Reason | Empirical scaling observation between lepton and quark sectors. |
| Modifications for A | Omit. The Lagrangian paper derives quark masses from $m_s$ via the O'Raifeartaigh mechanism, not via scaling from leptons. |

#### 7.4 The second scaling fails (lines 771--789)

| Classification | **B** |
|---|---|
| Reason | Negative empirical result constraining the scaling pattern. |
| Modifications for A | Omit. |

#### 7.5 Iterative chain (lines 791--852)

| Classification | **B** |
|---|---|
| Reason | Algebraic chain descending from $(t,b)$ via the Koide equation. The cyclic echo is an algebraic property of the Koide condition. |
| Modifications for A | Omit. Paper A derives masses from the superpotential, not from iterative Koide chains. |

#### 7.6 Cross-checks and master comparison (lines 854--908)

| Classification | **B** |
|---|---|
| Reason | Comparison table of predicted quark masses from different empirical routes. |
| Modifications for A | Replace with the prediction table from the Lagrangian (tab:predictions already exists). |

---

### Section 8: Composite Scalar Counting and Generation Uniqueness (lines 910--1221)

#### 8.1 Setup and Diophantine system (lines 922--958)

| Classification | **A** |
|---|---|
| Reason | Foundational for the Lagrangian: the Diophantine system that uniquely selects $N = 3$. |
| Modifications for A | None needed. Already in Lagrangian language. |

#### 8.2 The uniqueness theorem (lines 966--1011)

| Classification | **A** |
|---|---|
| Reason | Core mathematical result. |
| Modifications for A | None. |

#### 8.3 Partial solutions and hexagonal-number family (lines 1012--1033)

| Classification | **A** |
|---|---|
| Reason | Mathematical completeness of the uniqueness argument. |
| Modifications for A | None. |

#### 8.4 Identification with SU(5) representations (lines 1035--1088)

| Classification | **A** |
|---|---|
| Reason | Maps counting to SU(5) representations; forces symmetric 15. |
| Modifications for A | None. Already clean. |

#### 8.5 Charge census (lines 1090--1150)

| Classification | **A** |
|---|---|
| Reason | Structural content of the representation. |
| Modifications for A | None. |

#### 8.6 Anomaly cancellation and asymptotic freedom (lines 1152--1202)

| Classification | **A** |
|---|---|
| Reason | Consistency checks on the gauge theory. |
| Modifications for A | None. |

#### 8.7 Summary of the counting argument (lines 1204--1221)

| Classification | **A** |
|---|---|
| Reason | Logical summary. |
| Modifications for A | None. |

---

### Section 9: The Diquark Sector and SU(5) Flavor (lines 1223--1345)

#### 9.1 Counting partners (lines 1227--1241)

| Classification | **A** |
|---|---|
| Reason | Structural. |
| Modifications for A | None. |

#### 9.2 The symmetric representation (lines 1243--1270)

| Classification | **A** |
|---|---|
| Reason | How the SUSY generator produces the 15. |
| Modifications for A | None. |

#### 9.3 The Pauli theorem for scalar diquarks (lines 1272--1345)

| Classification | **A** |
|---|---|
| Reason | Key theoretical result for the Lagrangian paper: the diquark cannot be a $qq$ composite. |
| Modifications for A | None. Already clean. |

---

### Section 10: Two Algebraic Structures (lines 1348--1687)

#### 10.1 Part A: The SO(32) embedding chain (lines 1362--1389)

| Classification | **A** |
|---|---|
| Reason | UV completion of the Lagrangian. |
| Modifications for A | None. |

#### 10.2 Decomposition of the adjoint (lines 1391--1451)

| Classification | **A** |
|---|---|
| Reason | Technical group theory for the embedding. |
| Modifications for A | None. |

#### 10.3 Physical content and the projection problem (lines 1453--1517)

| Classification | **A** |
|---|---|
| Reason | Identifies sBootstrap representations in SO(32) adjoint. |
| Modifications for A | None. |

#### 10.4 Part B: The Casimir eigenvalue equation (lines 1519--1560)

| Classification | **SPLIT** |
|---|---|
| Reason | The equation itself could appear in either paper. The "structural constraints (codimension analysis)" is more Lagrangian-flavored. But the comparison with $\sin^2\theta_W$ is empirical. |
| Modifications for A | Present the equation and its structural constraints. State the predicted ratio $R = 0.22310$ as a consequence. Omit the de Vries historical framing. |
| Modifications for B | Keep the de Vries attribution and the scheme-dependence discussion as empirical evidence. |

#### 10.5 The electroweak ratio (lines 1561--1628)

| Classification | **SPLIT** |
|---|---|
| Reason | The algebraic derivation of $R$ is Paper A. The comparison with PDG, the $M_W$ prediction, and the CDF-II tension are empirical tests (shared with B). |
| Modifications for A | Present the $W$-mass prediction as a Lagrangian consequence. |
| Modifications for B | Present as empirical evidence for the algebraic structure. |

#### 10.6 Open issues for the Casimir equation (lines 1630--1661)

| Classification | **B** (primarily) |
|---|---|
| Reason | Scheme dependence, the fourth eigenvalue curiosity, and the de Vries historical note are empirical/historical. |
| Modifications for A | Move scheme-dependence caveat to Paper A as a brief footnote. The fourth eigenvalue and de Vries history go to Paper B. |

#### 10.7 Independence of the two structures (lines 1663--1687)

| Classification | **Both** |
|---|---|
| Reason | States that content (SO(32)) and dynamics (Casimir) are independent. Relevant to both. |
| Modifications for A | Keep a brief version. |

---

### Section 11: SUSY Breaking: From Seeds to Spectrum (lines 1690--3137)

This is the **core** of Paper A. Most of this section belongs there. Individual subsections:

#### 11.0 Section intro (lines 1690--1708)

| Classification | **A** |
|---|---|
| Reason | Sets up the SUSY-breaking story from the VA relation. |
| Modifications for A | Reframe: avoid "Koide triples"; use "mass triples satisfying the energy-balance condition." |

#### 11.1 The Seiberg framework (lines 1709--1732)

| Classification | **A** |
|---|---|
| Reason | SQCD framework for the Lagrangian. |
| Modifications for A | None. Already clean. |

#### 11.2 The three chiral superfields of each seed (lines 1734--1758)

| Classification | **A** |
|---|---|
| Reason | Maps seed triples onto O'Raifeartaigh superfields. |
| Modifications for A | Replace "Koide energy-balance condition" with "energy-balance condition" or "O'Raifeartaigh mass-ratio condition". |

#### 11.3 F-term breaking: O'Raifeartaigh mechanism (lines 1760--1960)

| Classification | **A** |
|---|---|
| Reason | Core Lagrangian content: the superpotential, $gv/m = \sqrt{3}$, the bloom as phase rotation, the bion Kahler mechanism for $v_0$-doubling. |
| Modifications for A | Systematic replacement: "Koide seed" -> "O'Raifeartaigh seed" or "mass-ratio seed". "Koide ratio" -> "energy-balance ratio". "Koide condition" -> "energy-balance condition" or "mass-ratio condition". "Koide parametrization" -> "angular parametrization". The paragraph "The bloom as a phase rotation" (line 1841) references eq:koide_angle which uses "Koide" in its label -- relabel. The $v_0$-doubling derivation (bion mechanism, lines 1876--1937) is fully Paper A material. |

#### 11.4 The Goldstino direction and VA realization (lines 1961--2061)

| Classification | **A** |
|---|---|
| Reason | Identifies the muon as pseudo-Goldstino, $f_{VA} = f_\pi$. |
| Modifications for A | Keep. Remove any "Koide" references. |

#### 11.5 D-term breaking (lines 2062--2107)

| Classification | **A** |
|---|---|
| Reason | Electromagnetic splittings from FI term. |
| Modifications for A | Clean as is. |

#### 11.6 Mapping the breaking to the tuples (lines 2108--2167)

| Classification | **SPLIT** |
|---|---|
| Reason | The three-layer breaking structure is Paper A. References to "Koide triples" are Paper B language. |
| Modifications for A | Replace "Koide triples" with "mass triples" throughout. |

#### 11.7 The O'Raifeartaigh structure in SQCD (lines 2168--2319)

| Classification | **A** |
|---|---|
| Reason | ISS mechanism, Kahler stabilization, pseudo-modulus, dual Koide, ISS vacuum lifetime. |
| Modifications for A | The "dual Koide" discussion (lines 2273--2319) needs renaming. Replace "dual Koide" with "Seiberg-seesaw mass condition" or "inverse-mass energy-balance condition". The algebraic identity $Q(a,b) = Q(1/a,1/b)$ is mathematical, keep it. |

#### 11.8 How the three supermultiplets join (lines 2320--2373)

| Classification | **A** |
|---|---|
| Reason | Embedding in the meson matrix. |
| Modifications for A | Replace "Koide condition" with "energy-balance condition" and "Koide triples" with "mass triples". |

#### 11.9 Energy balance is not an RG attractor (lines 2374--2396)

| Classification | **A** |
|---|---|
| Reason | Important negative result for the Lagrangian paper: the condition is a boundary condition, not an attractor. |
| Modifications for A | Replace "Koide preservation condition" with "energy-balance preservation condition". |

#### 11.10 The exact seed and what breaks it (lines 2397--2517)

| Classification | **SPLIT** |
|---|---|
| Reason | The two conditions for a seed (flat direction + mass ratio) are Paper A. The deviation table (lines 2421--2434) is empirical (B). The "remarkable feature of the VA relation" paragraph (lines 2492--2506) is both. |
| Modifications for A | Keep conditions 1 and 2, the discussion of electroweak running, charge running. Replace "Koide target" with "energy-balance target". Remove the empirical deviation table or present it as "predicted vs. observed" rather than as evidence for an empirical pattern. |

#### 11.11 What the $(c,b,t)$ triple tells us (lines 2518--2597)

| Classification | **SPLIT** |
|---|---|
| Reason | Interpretations (1)--(3) are Paper A (Lagrangian consequences). The CKM mixing paragraph (lines 2559--2582) is Paper A (structural consequence of mixed triples). The "empirical scan" confirmation (all triples near $Q = 2/3$ are mixed) is Paper B evidence. |
| Modifications for A | Keep the CKM-mixing observation but present it as "the mass triples required by the superpotential obligatorily mix up-type and down-type quarks" rather than "the Koide triples obligatorily mix..." |

#### 11.12 Toward a Superpotential (lines 2598--3137)

| Classification | **A** |
|---|---|
| Reason | **This is the heart of Paper A.** The full Lagrangian: Kahler potential, superpotential, D-terms, soft breaking, the prediction table, the parameter count. |
| Modifications for A | Systematic Koide->energy-balance replacement. Specific instances: line 2890 "seed Koide condition" -> "seed mass-ratio condition"; line 2988 "Seed Koide" column -> "Seed condition"; line 2997 "3-instanton" row references $\delta_\ell$ which is "the Koide phase" -> "the angular phase". The prediction table (tab:predictions, lines 3103--3136) is central to Paper A; keep it but replace "Koide" in the caption/source column with "energy-balance" or "O'Raifeartaigh". |

---

### Section 12: The Electron Mass and the Lepton Seed (lines 3139--3172)

| Classification | **B** (primarily) |
|---|---|
| Reason | Empirical discussion of the electron mass as a near-cancellation in $\mathfrak{z}_e$. The "Koide embedding" framing is Paper B. |
| Modifications for A | Paper A can include a brief remark that the energy-balance condition naturally accommodates one near-zero mass. Omit the table of $|\mathfrak{z}|/v_0$ values. |

---

### Section 13: Symmetry Restoration (lines 3175--3206)

| Classification | **Both** |
|---|---|
| Reason | The symmetry-restoration sequence is structural (Paper A: how the vacuum transitions work) and empirical (Paper B: what the observed spectrum looks like as symmetry is restored). |
| Modifications for A | Present as the vacuum structure of the Lagrangian under varying $m_q$ and $\alpha_{em}$. Replace "Koide splitting" with "mass splitting". |

---

### Section 14: Open Problems (lines 3209--3764)

#### 14.1 The diquark sector (lines 3213--3241)

| Classification | **A** |
|---|---|
| Reason | Structural open problem for the Lagrangian. |
| Modifications for A | None. |

#### 14.2 What determines the Koide phase? (lines 3242--3368)

| Classification | **SPLIT** |
|---|---|
| Reason | The dynamical origin of $\delta$ (three-instanton potential, $\mathbb{Z}_9$ symmetry) is **Paper A**. The empirical $\delta \bmod 2\pi/3 = 2/9$ observation and its statistical significance are **Paper B**. |
| Modifications for A | Present the $\cos(9\delta)$ potential as a consequence of the $\mathbb{Z}_{18}$ anomaly structure. State that the potential has minima that include the observed lepton phase. Do not present $2/9$ as an empirical coincidence; present it as a predicted minimum of the instanton potential. Replace "Koide phase" with "angular phase $\delta$" or "mass-spectrum phase". |
| Modifications for B | Keep the full empirical analysis: the $33$ ppm observation, the $3.1\sigma$ significance, the look-elsewhere correction. |

#### 14.3 Connecting the three tuples (lines 3369--3388)

| Classification | **A** |
|---|---|
| Reason | How the Lagrangian connects the sectors. |
| Modifications for A | Replace "Koide condition" with "energy-balance condition". |

#### 14.4 Lepton masses are not mesino masses (lines 3390--3448)

| Classification | **A** |
|---|---|
| Reason | Technical analysis of lepton mass generation: why the Seiberg vacuum mesino masses do not reproduce lepton masses, and the SU(2) confining sector proposal. |
| Modifications for A | Replace "Koide" references: "O'Raifeartaigh-Koide seed" -> "O'Raifeartaigh seed". |

#### 14.5 The baryon sector (lines 3450--3493)

| Classification | **A** |
|---|---|
| Reason | Structural: Miyazawa supercharge, baryon-meson connection. |
| Modifications for A | None needed. |

#### 14.6 The strong CP problem (lines 3495--3497)

| Classification | **Both** |
|---|---|
| Reason | Brief mention, relevant to both. |
| Modifications for A | Keep. |

#### 14.7 The kaon-muon puzzle (lines 3499--3511)

| Classification | **A** |
|---|---|
| Reason | Limitation of the universal soft-splitting prediction. |
| Modifications for A | None. Already in Lagrangian language. |

#### 14.8 Flavour-changing neutral currents (lines 3513--3528)

| Classification | **A** |
|---|---|
| Reason | Phenomenological constraint on the Lagrangian. |
| Modifications for A | None. |

#### 14.9 The Casimir eigenvalue equation (lines 3530--3537)

| Classification | **B** |
|---|---|
| Reason | Open issues with the empirical Casimir equation. |
| Modifications for A | Brief mention that the equation lacks a Lagrangian derivation. |

#### 14.10 The Higgs mass at tan(beta) = 1 (lines 3539--3566)

| Classification | **A** |
|---|---|
| Reason | NMSSM Higgs sector, S != X argument. |
| Modifications for A | None. Already clean. |

#### 14.11 Vacuum stability, CKM mixing, and the Oakes relation (lines 3568--3659)

| Classification | **A** |
|---|---|
| Reason | Vacuum stability, $F_X = 0$, off-diagonal meson masses, GST relation, baryon mass term. |
| Modifications for A | Replace "Koide condition" at line 3632 with "energy-balance condition". The paragraph about $Q(0, m_u+m_d, m_s) = 0.665$ (lines 3632--3639) uses Koide language; rewrite as: "The isospin-summed mass $m_u + m_d$ satisfies the seed mass-ratio condition to $0.27\%$, extending the prediction chain to the first generation." |

#### 14.12 The bloom mechanism (lines 3661--3745)

| Classification | **A** |
|---|---|
| Reason | Dynamical analysis of the bloom: adiabatic bion, frozen-sign bion, three-instanton driving force. |
| Modifications for A | Replace "Koide manifold" with "energy-balance manifold" or "Q = 2/3 surface". Replace "Koide parametrization" with "angular parametrization". Replace "Koide-preserving bloom" with "energy-balance-preserving bloom". |

#### 14.13 The SO(32) projection (lines 3747--3757)

| Classification | **A** |
|---|---|
| Reason | Open problem for the UV completion. |
| Modifications for A | None. |

#### 14.14 Generation uniqueness (lines 3759--3763)

| Classification | **A** |
|---|---|
| Reason | Open problem for the Diophantine system. |
| Modifications for A | None. |

---

### Section 15: Epistemological classification (lines 3766--3831)

| Classification | **SPLIT** |
|---|---|
| Reason | The table classifies claims as derived/observed/assumed. Paper A should have its own version focusing on derived results from the Lagrangian. Paper B should have its version focusing on empirical observations. |
| Modifications for A | Retain the "Derived" and "Assumed" rows. Move "Observed" rows to Paper B, except where the observation is also a Lagrangian prediction (e.g., $v_0$-doubling, which is "Observed" but also "Derived via bion Kahler"). |
| Modifications for B | Retain the "Observed" rows and the "Derived" rows that support the empirical analysis. |

---

### Section 16: Conclusions (lines 3833--3993)

| Classification | **SPLIT** |
|---|---|
| Reason | Seven principal results -- some are Lagrangian (A), some empirical (B). |

| Result | Classification |
|---|---|
| 1. "The Koide formula is energy equipartition" | **SPLIT**: The algebra is A; the name "Koide formula" and the meson triple observation are B |
| 2. "Seed triples set the stage for SUSY breaking" | **A** (rewrite to remove "Koide" name) |
| 3. "The seed condition is one superpotential constraint" | **A** |
| 4. "Two breaking mechanisms plus heavy sector" | **A** |
| 5. "The Diophantine counting uniquely selects N=3" | **A** |
| 6. "An algebraic electroweak ratio" | **SPLIT**: derivation is A, empirical comparison is B |
| 7. "An explicit Lagrangian with 2-3 free parameters" | **A** |

| Modifications for A | Rewrite conclusions around the Lagrangian: lead with the explicit Lagrangian (result 7), the Diophantine uniqueness (5), the O'Raifeartaigh superpotential constraint (3), the breaking mechanisms (4), the energy-balance condition as vacuum property (1), and the prediction table. Remove "Koide formula/triple/ratio" throughout; use "energy-balance condition" and "mass-ratio condition". |

---

### Bibliography (lines 3994--4089)

| Classification | **Both** |
|---|---|
| Reason | Each paper will cite different subsets. |
| Modifications for A | Omit Koide1983 if Paper A never mentions "Koide" by name. Keep Seiberg, O'Raifeartaigh, ISS, Unsal, GST, GOR, Miyazawa, Catto-Gursey, Fayet-Iliopoulos, Witten, DGMLY, Rivero. |
| Modifications for B | Keep Koide1983, Harari, Rivero2005b (de Vries), and all empirical references. |

---

## Summary: Content Allocation

### Paper A receives (in order of presentation):
1. Introduction (rewritten): sBootstrap postulate, Lagrangian motivation, signed mass charges
2. Energy-balance condition (algebra from Sec. 2, rewritten without "Koide" name)
3. O'Raifeartaigh seed (from Sec. 3.1 analytical condition + Sec. 11.3)
4. Composite scalar counting and generation uniqueness (Sec. 8, entire)
5. Diquark sector and SU(5) flavor (Sec. 9, entire)
6. SO(32) embedding (Sec. 10, Part A)
7. Casimir eigenvalue equation (Sec. 10, Part B, derivation only)
8. SUSY breaking (Sec. 11, essentially entire, with Koide->energy-balance replacements)
9. The explicit Lagrangian (Sec. 11.12, heart of Paper A)
10. Symmetry restoration (Sec. 13)
11. Open problems: diquark, bloom mechanism, lepton masses, baryon sector, kaon-muon puzzle, FCNC, Higgs mass, vacuum stability, CKM, SO(32) projection, generation uniqueness
12. Epistemological table (derived + assumed rows)
13. Conclusions (results 2--5, 7, plus rewritten 1 and 6)

### Paper B receives:
1. Introduction: Koide formula history, empirical observations
2. The Koide Formula as Energy Balance (Sec. 2, with "Koide" name)
3. Seed triples: empirical evidence (Sec. 3.2--3.5)
4. Scale coincidence (Sec. 3.5)
5. A Tale of Two Tuples (Sec. 4, narrative framing)
6. The Meson Koide Triple (Sec. 5, entire)
7. The Third Tuple as empirical observation (Sec. 6.1, 6.3)
8. Mass relations and iterative chains (Sec. 7, entire)
9. Casimir equation: empirical comparison (Sec. 10.5, W-mass, CDF-II, de Vries history)
10. Statistical significance: look-elsewhere corrections for all patterns
11. Epistemological table (observed rows)
12. Conclusions (results 1, 6 in empirical framing)

---

## Content Gaps for Paper A

The following content is needed for Paper A but is not currently in the paper (or is present only in "Koide language"):

1. **Self-contained derivation of the energy-balance condition from the superpotential.** Currently, Sec. 2 derives $Q = 2/3$ from the signed-charge formalism, and then Sec. 11.3 shows the O'Raifeartaigh model produces it. Paper A needs a unified presentation: start from the superpotential, derive the fermionic mass matrix, show that $gv/m = \sqrt{3}$ produces the energy-balance condition, then present the condition as a theorem about the vacuum, not as an empirical observation.

2. **The Kahler stabilization at $c = -1/12$ needs to be presented as a self-contained argument**, not as an aside in the "Koide phase" discussion. Currently it appears in Sec. 11.3 and again in Sec. 14.2. Paper A should have a dedicated subsection.

3. **The prediction chain $m_s \to m_c \to m_b \to m_t$ needs a clean, forward presentation.** Currently the chain is assembled from pieces scattered across Sections 4, 7, and 11.12. Paper A should present it linearly: given $m_s$, the seed condition gives $m_c$; the bion Kahler gives $m_b$; the overdetermination test gives $m_t$.

4. **The $\tan\beta = 1$ prediction needs more prominence.** It currently appears in Sec. 11.12 (lines 2724--2749) but is somewhat buried. In Paper A it should be a highlighted result.

5. **The SU(2) lepton confining sector** (lines 3420--3448) is currently a brief sketch in "Open Problems". Paper A should either expand it into a proper section or clearly flag it as a conjecture.

6. **RG analysis**: The negative result (Sec. 11.9, energy balance is not an RG attractor) should be more prominent in Paper A, as it establishes that the condition must be a boundary condition.

7. **Parameter count**: The explicit statement "only $m_u$, $m_d$, $m_s$ remain free" (currently at line 3055) should be a highlighted conclusion of Paper A, with a clear diagram showing which Lagrangian term fixes which mass.

8. **The "mass spectrum from the Lagrangian" table (lines 2978--3001) should be the centerpiece of Paper A**, not buried deep in Sec. 11.12.

---

## Paragraphs Currently in "Koide Language" That Can Be Rewritten for Paper A

| Location | Current Koide Language | Rewrite for Paper A |
|---|---|---|
| Line 122--127 | "The Koide ratio $Q = 2/3$ is then not numerology but the energy-balance condition" | "The fermionic mass matrix of the O'Raifeartaigh vacuum satisfies the energy-balance condition $\langle z_k^2 \rangle = z_0^2$ when $gv/m = \sqrt{3}$" |
| Line 129 | "the full sBootstrap has three Koide tuples" | "the superpotential produces three mass triples satisfying the energy-balance condition" |
| Line 143 | "the Koide condition constrains the superpotential" | "the mass-ratio condition constrains the superpotential" |
| Line 200 | $Q = 2/3 \Longleftrightarrow \langle z_k^2 \rangle = z_0^2$ (boxed) | Keep the boxed equation but present as "Theorem: The O'Raifeartaigh fermionic spectrum satisfies..." |
| Line 1751--1752 | "enforced by the Koide energy-balance condition applied to the seed" | "enforced by the energy-balance condition at the O'Raifeartaigh vacuum" |
| Line 1826--1827 | "The Koide condition on the seed is therefore one algebraic constraint" | "The mass-ratio condition on the seed is therefore one algebraic constraint" |
| Lines 2140, 2322, 2361, 2366 | "Koide triples" | "mass triples" or "energy-balance triples" |
| Line 2383--2387 | "Koide preservation condition" | "energy-balance preservation condition" |
| Line 2890 | "seed Koide condition" | "O'Raifeartaigh seed condition" |
| Line 2988 | "Seed Koide" in table | "Seed condition" |
| Lines 3663--3668 | "Koide parametrization" and "Koide manifold" | "angular parametrization" and "energy-balance manifold" |

---

## Sections That Should Be Split

| Section | Lines | What Goes to A | What Goes to B |
|---|---|---|---|
| Abstract | 29--98 | Lagrangian results, predictions | Empirical Koide observations, statistical significance |
| Sec. 1 (Introduction) | 103--177 | Themes 1 and 3 (charges, superpotential) | Theme 2 (three observed tuples) |
| Sec. 2.2 (Derivation of Q=2/3) | 189--221 | Energy-balance algebra as vacuum property | "The Koide Formula" as empirical observation |
| Sec. 4.1 (Seed to full triple) | 355--419 | $v_0$-doubling as bion Kahler consequence | Narrative "tale of two tuples" |
| Sec. 7.2 (Koide angle parametrization) | 662--712 | The $(M_0, \delta_0)$ parametrization as mathematical tool | The $\delta_0 \bmod 2\pi/3 = 2/9$ empirical observation |
| Sec. 10.4--10.5 (Casimir equation) | 1519--1628 | Derivation, structural constraints, predicted $R$ | Comparison with $\sin^2\theta_W$, $M_W$ prediction, CDF-II |
| Sec. 11.10 (Exact seed) | 2397--2517 | Conditions 1 and 2, running discussion | Empirical deviation table |
| Sec. 11.11 ($(c,b,t)$ triple) | 2518--2597 | Interpretations, CKM-mixing structure | Exhaustive scan confirming mixed triples |
| Sec. 14.2 (Koide phase) | 3242--3368 | $\cos(9\delta)$ from instanton potential | $2/9$ observation and its significance |
| Sec. 15 (Epistemology table) | 3766--3831 | Derived + Assumed rows | Observed rows |
| Sec. 16 (Conclusions) | 3833--3993 | Results 2--5, 7; rewritten 1, 6 | Results 1, 6 in empirical framing |
