# Agent Swarm: Computational Tasks

## Architecture Note (coordinator eyes only)

This document contains self-contained computational tasks for a compartmentalized agent swarm.

**Design principles:**
- Each agent gets a bounded mathematical problem. No theory, no motivation, no program name.
- No agent sees another agent's task. No cross-references.
- Agents work with numbers, formulas, and group theory — not with interpretations.
- References to viXra or heterodox frameworks are omitted entirely.
- Success is defined by computational output, not by agreement with any theoretical framework.
- Synthesis is the coordinator's job. The agents produce ingredients.

**Embedded self-verification (the Joan Hinton principle):**
Each task contains a built-in cross-check via an independent computational route. The agent verifies its own work without knowing it's being checked — the spy plays viola in the quartet, she doesn't audit the rehearsal from outside. This eliminates the need for external overseer agents that would see multiple compartments and potentially reconstruct the program.

**Parallelization:** All 8 tasks are independent. Launch simultaneously.

---

## Agent 1: Mass Triplet Formula — Exhaustive Scan

Write a Python script. Define the quantity:

$$Q(m_1, m_2, m_3) = \frac{(\sqrt{m_1} + \sqrt{m_2} + \sqrt{m_3})^2}{m_1 + m_2 + m_3}$$

Use PDG 2024 particle masses (pole masses where available, $\overline{MS}$ at self-scale for quarks):

| particle | mass |
|----------|------|
| e | 0.51099895 MeV |
| μ | 105.6583755 MeV |
| τ | 1776.86 MeV |
| u | 2.16 MeV |
| d | 4.67 MeV |
| s | 93.4 MeV |
| c | 1.27 GeV |
| b | 4.18 GeV |
| t | 172.69 GeV |

For all $\binom{9}{3} = 84$ possible triplets drawn from these nine masses:

1. Compute $Q$ with all positive square roots.
2. Compute $Q$ with each of the 7 possible sign-flip combinations (one negative, two negative, all three negative).
3. For each triplet and sign combination, record $|Q - 3/2|$.

**Self-check:** Also compute $Q$ via the equivalent angle formula. Define:

$$\cos\theta = \frac{\sqrt{m_1} + \sqrt{m_2} + \sqrt{m_3}}{\sqrt{3(m_1 + m_2 + m_3)}}$$

Then $Q = 3\cos^2\theta$, so $Q = 3/2$ iff $\theta = 45°$. For each triplet, compute $Q$ both ways. Report any discrepancy $> 10^{-10}$ between the two formulas.

**Output:** A ranked table of all triplets with $|Q - 3/2| < 0.05$, showing the triplet, the sign combination used, the value of $Q$, the deviation from $3/2$, and the angle $\theta$ in degrees. Also compute $Q$ for the charged lepton triplet $(e, \mu, \tau)$ to 10 significant figures and report $\theta$ to 8 significant figures.

---

## Agent 2: Iterative Mass Chain from Quadratic Constraint

Write a Python script. Given the constraint $Q(m_1, m_2, m_3) = 3/2$ (formula defined above), solve for $m_3$ as a function of $m_1, m_2$:

$$\sqrt{m_3} = \left(\sqrt{m_1} + \sqrt{m_2}\right) \left(2 - \sqrt{3 + \frac{6\sqrt{m_1 m_2}}{(\sqrt{m_1} + \sqrt{m_2})^2}}\right)$$

This gives $\sqrt{m_3}$ (which may be positive, negative, or zero). The full solution has four branches from sign choices in both the outer and inner square roots. For each step, evaluate all four branches and report them.

**Chain A (descending, empirical masses):**
Start from $m_t = 172.69$ GeV, $m_b = 4.18$ GeV.
- Step 1: Solve $Q(m_b, m_?, m_t) = 3/2$ for $m_?$. Report all four branches. Identify the branch nearest to $m_c = 1.27$ GeV.
- Step 2: Using the result from Step 1 as $m_c$, solve $Q(m_c, m_?, m_b) = 3/2$. The physically relevant branch requires $\sqrt{m_?}$ to be **negative**. Report all four branches. Identify the one nearest to $m_s = 93.4$ MeV.
- Step 3: Continue with $(m_s, m_?, m_c)$. Identify the branch nearest to $m_u$ or $m_d$.
- Step 4: Continue one more step. Report.

**Chain B (ascending, lepton-seeded):**
Fit the charged lepton masses $(m_e, m_\mu, m_\tau)$ to the parametrization:

$$\sqrt{m_k} = \sqrt{M_0} \left(1 + \sqrt{2} \cos\left(\frac{2\pi k}{3} + \delta_0\right)\right), \quad k = 0, 1, 2$$

Extract $M_0$ and $\delta_0$. Then construct a second triplet with parameters $M_1 = 3 M_0$ and $\delta_1 = 3 \delta_0$. Compute the three masses of this second triplet. Compare to PDG values of $s$, $c$, $b$ quarks.

**Self-check:** For every predicted mass at each step, substitute it back into $Q(m_1, m_2, m_3)$ and verify $Q = 3/2$. Report the back-substitution residual $|Q - 3/2|$ at each step. Any residual $> 10^{-8}$ indicates a computational error.

**Output:** Full table of all chains, all branches at each step, back-substitution residuals, and comparison to PDG masses. Tabulate: particle, predicted mass, PDG mass, % deviation.

---

## Agent 3: Simultaneous Constraint on a 3×2 Mass Grid

Write a Python/Mathematica script. Place 6 non-negative real numbers in a $3 \times 2$ grid:

$$\begin{pmatrix} a_1 & b_1 \\ a_2 & b_2 \\ a_3 & b_3 \end{pmatrix}$$

A "triplet" is formed by choosing one element from each row. There are $2^3 = 8$ such triplets. Define $Q$ as in Agent 1.

**Problem:** Find all configurations of 6 non-negative reals such that all 8 triplets simultaneously satisfy $Q = 3/2$, allowing sign flips on $\sqrt{m}$.

Approach:
1. Let $k^+ = 2 + \sqrt{3}$, $k^- = 2 - \sqrt{3}$. Note $k^+ k^- = 1$, and the triplet $(k^-, 1, k^+)$ satisfies $Q = 3/2$ with all positive roots.
2. Start from the fully degenerate case: all $a_i = b_i$ (column 1 = column 2). Characterize all such degenerate solutions.
3. Look for solutions where exactly one pair $a_j \neq b_j$. How far can they differ while all 8 constraints hold?
4. Classify solutions by how many pairs are non-degenerate.
5. For each solution found, scale $M_0$ so that the column-1 masses approximate $(0.122, 1.70, 3.64)$ GeV. Report what column-2 masses become under this scaling.

**Self-check:** For each solution, verify all 8 triplets by direct substitution of $Q$. Also check independently: does the 3×2 grid possess any symmetry (row permutations, column swaps, joint scaling $m \to \lambda m$) that was not imposed by construction? Report all discovered symmetries.

**Output:** Complete catalog of solution types. For each, the 3×2 grid of values (in units of $M_0$), which sign combinations are needed, all 8 $Q$ values by direct substitution, and any emergent symmetries.

---

## Agent 4: Diophantine System — Integer Solutions

Solve the following system of Diophantine equations over positive integers $(N, r, s)$:

**System A (three equations):**

$$rs = 2N \tag{1}$$
$$\frac{r(r+1)}{2} = 2N \tag{2}$$
$$r^2 + s^2 - 1 = 4N \tag{3}$$

1. From (1) and (2), eliminate $N$ to get a relation between $r$ and $s$.
2. Substitute into (3). Find all positive integer solutions.
3. List all solutions with $N \leq 1000$.

**System B (two equations only, drop equation 3):**

Using only equations (1) and (2):
1. Show that $2N$ must be a triangular number $T_r = r(r+1)/2$.
2. List all solutions with $N \leq 200$.
3. For each solution, compute $r^2 + s^2 - 1$ and compare to both $4N$ and $2N$. Which equation (if either) is satisfied?

**System C (modified equation 3):**

Replace equation (3) with $r^2 + s^2 - 1 = 2N$. Find all positive integer solutions of the system (1), (2), (C). Are there any?

**Bonus:** For each solution of System B, compute $n_f = r + s$. Note which solutions have $2 n_f < 33$.

**Self-check:** Solve System A by two independent methods: (i) algebraic elimination as above, (ii) brute-force enumeration of all $(r, s)$ pairs with $1 \leq r, s \leq 100$, testing all three equations. Compare results. Any disagreement indicates an error.

**Output:** Complete solution tables for Systems A, B, C. For System A, prove whether the solution is unique. Report whether algebraic and brute-force methods agree.

---

## Agent 5: Lie Algebra Branching Rules — SO(32) Adjoint

Compute the branching rules for the adjoint representation of $SO(32)$ along **two independent chains** and verify they produce the same final result.

The adjoint of $SO(32)$ has dimension $\frac{32 \times 31}{2} = 496$.

**Chain 1:**

$$SO(32) \supset SU(16) \times U(1) \supset SU(15) \times U(1)^2 \supset SU(5) \times SU(3) \times U(1)^2$$

**Chain 2:**

$$SO(32) \supset SO(30) \times U(1) \supset SO(10) \times SO(20) \times U(1)$$

then branch the $SO(10)$ piece to $SU(5) \times U(1)$ and the $SO(20)$ piece to $SU(3) \times \ldots$

For each chain:
1. List all irreps that appear at each step, with dimensions and $U(1)$ charges.
2. Verify that dimensions sum to 496 at every level.
3. At the final level $SU(5) \times SU(3) \times U(1)^2$, present the full decomposition as a table: $(SU(5)_{\text{irrep}}, SU(3)_{\text{irrep}})_{Y_1, Y_2}$ with multiplicity and dimension.

**Self-check:** The two chains must produce the same final decomposition table. Any discrepancy indicates an error in one or both chains. Report whether the results agree.

If a computer algebra system (LiE, SageMath, or similar) is available, use it as a third verification. Otherwise, use Slansky's "Group Theory for Unified Model Building" (Physics Reports 79, 1981).

**Output:** The complete decomposition table from each chain, dimension checks, and explicit comparison of the two final results.

---

## Agent 6: SU(5) Representation Decomposition and Charge Assignment

Decompose the following representations of $SU(5)$ under the subgroup $SU(3) \times SU(2) \times U(1)$:

$$\mathbf{24}, \quad \mathbf{15}, \quad \mathbf{\overline{15}}, \quad \mathbf{10}, \quad \mathbf{\overline{10}}, \quad \mathbf{5}, \quad \mathbf{\bar{5}}$$

For each decomposition, the $SU(3) \times SU(2)$ content is standard. The task concerns the $U(1)$ charge assignment.

**Problem:** Find all $U(1)$ charge assignments (normalizations) such that, when combined into a single electric charge formula:

$$Q = \alpha Y_1 + \beta Y_2$$

(where $Y_1$ is the $U(1)$ charge from $SU(5) \supset SU(3) \times SU(2) \times U(1)$ and $Y_2$ is an independent $U(1)$), the states in $\mathbf{24} \oplus \mathbf{15} \oplus \mathbf{\overline{15}}$ produce electric charges exclusively from the set $\{0, \pm 1/3, \pm 2/3, \pm 1, \pm 4/3\}$.

For each valid charge assignment:
1. Count how many states have each charge value.
2. List the charge content as a table: $(SU(3), SU(2))_{Y_1}$ irrep, charge $Q$, multiplicity.
3. Separately count the states with $Q = \pm 4/3$. What $SU(3) \times SU(2)$ irrep do they live in?

**Self-check:** For each decomposition, verify dimensions by summing: the components of $\mathbf{24}$ must sum to 24, of $\mathbf{15}$ to 15, etc. Also verify independently via the tensor product route: $\mathbf{24} = \mathbf{5} \otimes \mathbf{\bar{5}} - \mathbf{1}$ and $\mathbf{15} = \text{Sym}^2(\mathbf{5})$. Decompose the tensor products under $SU(3) \times SU(2)$ and confirm agreement with the direct branching rule result.

**Output:** All valid charge assignments with the resulting charge spectrum table. Dimension checks and tensor product cross-verification for each irrep.

---

## Agent 7: Mass Splitting from Poincaré Casimir Invariants

The Poincaré group has two Casimir invariants:

$$C_1 = P_\mu P^\mu = m^2$$
$$C_2 = W_\mu W^\mu = -m^2 s(s+1)$$

where $s$ is the spin. Consider a mass operator $M_s^2$ that depends on both Casimirs and satisfies:

- $\lim_{s \to \infty} M_s^2 / m^2 = 1$ (asymptotic Regge behavior)
- $M_s^2$ is a solution of a polynomial equation in $C_1$ and $C_2$

The simplest such equation is:

$$M^4 - M^2 C_2 + C_1 C_2 = 0$$

which gives:

$$M_s^2 = \frac{1}{2}\left(C_2 \pm \sqrt{C_2^2 - 4 C_1 C_2}\right)$$

Substituting $C_1 = m^2$ and $C_2 = m^2 s(s+1)$:

$$M_s^2 = \frac{m^2 s(s+1)}{2}\left(1 \pm \sqrt{1 - \frac{4}{s(s+1)}}\right)$$

**Task:**

1. For $s = 1/2$: compute both roots $M_{1/2,+}^2$ and $M_{1/2,-}^2$ in terms of $m$.
2. For $s = 1$: compute both roots $M_{1,+}^2$ and $M_{1,-}^2$ in terms of $m$.
3. Compute the ratio:

$$R = 1 - \frac{M_{1/2,+}^2}{M_{1,+}^2}$$

Give the exact algebraic form and a numerical value to 12 significant digits.

4. Compare $R$ to the experimental quantity $\sin^2\theta_W = 1 - M_W^2/M_Z^2 = 0.22306 \pm 0.00033$ (from the on-shell definition using $M_W = 80.3692 \pm 0.0133$ GeV and $M_Z = 91.1876 \pm 0.0021$ GeV).
5. If $m$ is set so that $M_{1,+} = M_Z = 91.1876$ GeV, compute $M_{1/2,+}$, $M_{1,-}$, and $M_{1/2,-}$ in GeV.
6. The electroweak vacuum expectation value is $\langle v \rangle = (\sqrt{2} G_F)^{-1/2} = 246.22$ GeV, giving $\langle v \rangle / \sqrt{2} = 174.10$ GeV. Compare $|M_{1,-}|$ to this value.

**Self-check via independent classical derivation:** Consider the relativistic equation for a de Broglie standing wave on a circular orbit of period $T_r = h/(m_0 c^2)$. With the Landé-Pauli angular momentum substitution $L \to \sqrt{j(j+1)}\hbar$, the velocity condition becomes:

$$\frac{\beta^2}{\sqrt{1 - \beta^2}} = \sqrt{j(j+1)}$$

Solve for $\beta(j)$ at $j = 1/2$ and $j = 1$. Compute:

$$R_{\text{classical}} = 1 - \frac{\beta_{1/2}^2}{\beta_1^2}$$

This should equal $R$ from the Casimir calculation. Report both values and their difference.

**Output:** Exact algebraic expressions, numerical values, comparisons to experimental quantities, and the cross-check between quantum (Casimir) and classical (de Broglie) derivations.

---

## Agent 8: Numerical Coincidences in Pseudoscalar Meson and Lepton Masses

Using PDG 2024 values for charged lepton masses and pseudoscalar meson masses:

**Leptons:** $m_e = 0.51099895$ MeV, $m_\mu = 105.6583755$ MeV, $m_\tau = 1776.86$ MeV.

**Pseudoscalar mesons:** $m_{\pi^0} = 134.9768$ MeV, $m_{\pi^\pm} = 139.57039$ MeV, $m_{K^0} = 497.611$ MeV, $m_{K^\pm} = 493.677$ MeV, $m_\eta = 547.862$ MeV, $m_{\eta'} = 957.78$ MeV.

**Also define:** $m_{\eta_8} = \sqrt{(4 m_{K^0}^2 - m_{\pi^0}^2)/3}$ (Gell-Mann–Okubo prediction for the octet $\eta$).

Compute the following quantities and report each to 6 significant figures:

**Set A (mixed lepton products):**
For all $\{a, b, c\} \subset \{e, \mu, \tau\}$ with $m_b > m_c$, compute:
- $\sqrt{m_a} \cdot \sqrt{m_b - m_c}$

**Set B (meson-lepton differences):**
- $m_{\pi^0} - m_\mu$
- $m_{\pi^\pm} - m_\mu$
- $m_{\eta_8} - m_{\pi^0}$
- $m_\eta - m_{\pi^0}$

**Set C (de Vries' isospin relation):**
- Compute $\left|\frac{m_{\pi^\pm}}{m_{\pi^0}} - 1\right|^2$ and compare to $\frac{m_\mu}{m_Z}$, where $m_Z = 91187.6$ MeV.

**Set D (electromagnetic decay width scaling):**
For each of $\pi^0, \eta, \omega(782), J/\psi, Z^0$, using PDG total and electromagnetic partial widths, compute $(\Gamma/M^3)^{-1/2}$ in GeV. Report all values in a table.

**Task:** Scan all pairs of quantities from Sets A and B. For each pair $(x, y)$, compute $|x - y| / \max(x, y)$. Report all pairs agreeing to better than 1%.

**Self-check (look-elsewhere correction):** Sets A and B together contain $N$ computed quantities. The number of pairwise comparisons is $\binom{N}{2}$. If masses were random draws from a log-uniform distribution over $[0.1, 2000]$ MeV, estimate by Monte Carlo simulation (10,000 trials) the probability of finding $k$ or more coincidences at the 1% level from $\binom{N}{2}$ comparisons. Report the look-elsewhere-corrected significance of the actual number of coincidences found.

**Output:** Table of all coincidences found, sorted by precision of agreement (best first). For each, give the two expressions, their numerical values, and the percentage discrepancy. Then the look-elsewhere analysis.

---

## Coordinator Notes (not for agents)

**Synthesis map** — what the coordinator reconstructs from the 8 outputs:

| Agent | Provides | Connects to |
|-------|----------|-------------|
| 1 | Which mass triplets satisfy Q=3/2 | Identifies the empirical Koide triplets |
| 2 | Iterative chain predictions | Shows whether 7 masses are determined from 2 inputs |
| 3 | Simultaneous 8-Koide solutions | Establishes the folding/resolvent structure |
| 4 | Unique (N,r,s) = (3,3,2) | Proves 3 generations is forced |
| 5 | SO(32) → SU(5) × SU(3) decomposition | Shows 496 accommodates 54 scalar states |
| 6 | Charge spectrum of 24⊕15⊕15̄ | Identifies the ±4/3 diquark content |
| 7 | R = 0.22310... ≈ sin²θ_W | Links Poincaré structure to EW scale |
| 8 | Meson-lepton mass coincidences | Documents the empirical regularities |

No agent sees this table. Each computes; the coordinator integrates.

**Self-verification design:** Each agent cross-checks its own results via an independent route embedded in the task itself. Agent 5 uses two different branching chains. Agent 7 uses both a quantum (Casimir) and classical (de Broglie) derivation. Agent 4 uses algebraic + brute-force. The verification is structural — it's part of the homework problem, not an external audit. The spy plays in the string quartet.

**Deployment:** All 8 tasks are independent. Launch in parallel via Claude API, separate chat sessions, or any multi-agent framework. Each agent writes code and reports numerical results. No agent needs context beyond its own prompt.
