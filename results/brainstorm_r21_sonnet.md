# Brainstorm: Proposed Calculations for the sBootstrap Program
**Date:** 2026-03-04

---

## Summary of the Current Impasse

The program has three well-grounded computational results (Koide seed from O'Raifeartaigh at $t = \sqrt{3}$, $v_0$-doubling giving $m_b$ to 0.1$\sigma$, SO(32) bootstrap self-reference) and four unresolved structural problems (pseudo-modulus stabilization, $\mathbf{15}$ vs $\overline{\mathbf{10}}$ Pauli conflict, SO(32)$\to$54-state projection, and the bloom mechanism). The proposals below are ordered by what is both tractable and highest-leverage.

---

## Proposal 1: The $\mathfrak{z}^2 = m$ Condition as a Jordan Algebra Constraint

**Priority: High**

### What to compute

The mass parametrization $m_k = (v_0 + z_k)^2$ where $\sum z_k = 0$ and the Koide condition $\langle z_k^2 \rangle = v_0^2$ defines a real quadratic surface in the space of the $z_k$. The coordinates $\mathfrak{z}_k = \sqrt{m_k}$ are the Jordan square roots of the masses, and the constraint $m_k = \mathfrak{z}_k^2$ means the masses live in a Jordan algebra $J$ where $\mathfrak{z}_k \circ \mathfrak{z}_k = m_k$ under the symmetric product.

The specific task: identify whether the seed triple $(0, \mathfrak{z}_-, \mathfrak{z}_+) = (0, (\sqrt{6}-\sqrt{2})/2, (\sqrt{6}+\sqrt{2})/2)$ can be embedded as a minimal idempotent decomposition $e_0 + e_- + e_+ = 1$ in one of the four exceptional Jordan algebras $J_3(\mathbb{O}), J_3(\mathbb{H}), J_3(\mathbb{C}), J_3(\mathbb{R})$.

The critical check: the exceptional Jordan algebra $J_3(\mathbb{O})$ (Albert algebra, dimension 27) has an idempotent spectral theorem — every element has a unique decomposition into primitive idempotents with real eigenvalues. The Freudenthal-Tits construction of $J_3(\mathbb{O})$ produces eigenvalues from the characteristic polynomial $\lambda^3 - \text{tr}(A)\lambda^2 + \frac{1}{2}((\text{tr}A)^2 - \text{tr}(A^2))\lambda - \det A = 0$. Ask: for what trace and determinant values does this characteristic polynomial have roots $(0, 2-\sqrt{3}, 2+\sqrt{3})$?

Answer by direct substitution: $\text{tr} = 4$, discriminant condition from $\lambda^3 - 4\lambda^2 + (2+\sqrt{3})(2-\sqrt{3})\lambda = \lambda^3 - 4\lambda^2 + \lambda = 0$, so $\text{tr}(A^2) = 4^2 - 2\cdot 1 = 14$ and $\det A = 0$. The constraint $\det A = 0$ in $J_3(\mathbb{O})$ defines the "rank-2 boundary" of the Jordan cone — which is precisely the locus of elements where one eigenvalue is zero, i.e., the Koide seed boundary.

**Tools:** Representation theory of Jordan algebras. No computer needed for this part; it is algebraic verification.

**What constitutes a positive result:** The seed triple sits on the rank-2 boundary of the Albert algebra's positive cone, and the bloom (moving $\delta$ away from $3\pi/4$) moves off this boundary into the interior. This gives a geometric reason why the seed has $m_0 = 0$ and the bloom breaks it: the seed is a degenerate element of $J_3(\mathbb{O})$, and bloom = deformation into the generic stratum.

**What constitutes a negative result:** The eigenvalues $(0, 2-\sqrt{3}, 2+\sqrt{3})$ do not correspond to any natural element of $J_3(\mathbb{O})$, and the identification is numerological. This would close the Jordan algebra avenue.

**Why this matters:** $J_3(\mathbb{O})$ appears in the Freudenthal magic square with $E_8$, and the $E_8 \times E_8$ heterotic string is already known to be blocked. If $J_3(\mathbb{O})$ is also blocked (as expected from the $\overline{\mathbf{10}}$ issue), but $J_3(\mathbb{H})$ (dimension 15, connected to $Sp(6)$) or $J_3(\mathbb{C})$ (dimension 9, connected to $SU(3)$) is relevant, this may clarify the $\mathbf{15}$ vs $\overline{\mathbf{10}}$ conflict. The symmetric $\mathbf{15}$ is precisely the Jordan product $\mathbf{5} \circ \mathbf{5}$ in $J_3(\mathbb{C})$, not the antisymmetric exterior product. This reframes the "Pauli conflict" as a question about which product structure is appropriate — Jordan (symmetric, gives $\mathbf{15}$) or Grassmann (antisymmetric, gives $\overline{\mathbf{10}}$). In SUSY, chiral superfields commute under the Jordan product in the chiral ring $\mathcal{R} = \mathbb{C}[\phi_i]/(\partial W)$.

---

## Proposal 2: Stabilize the Pseudo-Modulus with a Competing Sector Superpotential

**Priority: High**

### What to compute

The Coleman-Weinberg calculation (pseudomodulus_cw.md) establishes unambiguously that the one-loop minimum sits at $t \approx 0.49$, far below the Koide point $t = \sqrt{3} \approx 1.73$. The Kahler correction with $c = -1/12$ places a pole at $t = \sqrt{3}$ but does not create a minimum there — it creates a wall. No known perturbative mechanism fixes $t = \sqrt{3}$.

The proposal is to determine whether a competing sector generates an effective potential that has a minimum at $t = \sqrt{3}$ via the following mechanism: introduce a second O'Raifeartaigh model (the lepton sector, proposed as $Sp(2)$) coupled to the first through a messenger field $S$. The superpotential coupling is $W \supset \kappa S (t - \sqrt{3}) = \kappa S \cdot X / (m/g) - \kappa S \sqrt{3}$ (schematically). The $F_S$ equation then enforces $\langle X \rangle = \sqrt{3} m/g$ at the minimum. This is a NMSSM-type stabilization but with the stabilization point determined by the Koide condition algebraically.

**Explicit computation:**

Write the coupled superpotential:
$$W = fX + m\phi\tilde\phi + gX\phi^2 + \kappa S(gX - \sqrt{3}m) + \mu^2 S$$

The $F_S = 0$ equation gives $\kappa(gX - \sqrt{3}m) + \mu^2 = 0$, so $gX = \sqrt{3}m - \mu^2/\kappa$. In the limit $\mu^2 \to 0$, this enforces $X = \sqrt{3}m/g$ exactly. The $F_X$ equation: $F_X = f + g\phi^2 + \kappa g S = 0$. At $\phi = 0$ and $S = -f/(\kappa g)$, this is satisfied. But then $F_S = \kappa g X - \kappa \sqrt{3} m + \mu^2 = \mu^2 \neq 0$ — this breaks SUSY, with $\langle F_S \rangle = \mu^2$.

Check the scalar potential at the proposed vacuum $(X = \sqrt{3}m/g, \phi = \tilde\phi = 0, S = -f/(\kappa g))$:
$$V = |F_X|^2 + |F_\phi|^2 + |F_S|^2 = 0 + 0 + \mu^4$$

This is SUSY-breaking with $V_0 = \mu^4$, and the O'Raifeartaigh F-term $f$ does not contribute directly to the vacuum energy (it is cancelled at $X = \sqrt{3}m/g$ when $\phi = 0$). Compute the mass spectrum at this vacuum and check whether the Koide condition $Q = 2/3$ is preserved for the $(\phi, \tilde\phi)$ sector fermions.

**Tools:** Standard superpotential algebra. F-term equations are linear in the fields. One-loop corrections should be computed but the tree-level stabilization is tractable analytically.

**Positive result:** A two-coupling system $(f, \mu^2, \kappa)$ with $\mu^2 \ll \kappa^2 m$ produces $X = \sqrt{3}m/g + O(\mu^2/\kappa)$ with the Koide seed in the fermion spectrum. The bloom is then the deviation from this vacuum when $\mu^2 \neq 0$, controlled by $\mu^2/\kappa \ll m/g$.

**Negative result:** The $S$-field coupling generically lifts the goldstino (Goldstino is eaten, i.e., there is no SUSY breaking in the $\phi, \tilde\phi$ sector anymore) or the mass spectrum at the new vacuum no longer has $Q = 2/3$.

---

## Proposal 3: The SO(32) $\to$ 54 Projection via Wilson Lines

**Priority: High**

### What to compute

The SO(32) adjoint $\mathbf{496}$ decomposes under $SU(5)_f \times SU(3)_c$ as:
$$\mathbf{496} = (\mathbf{24},\mathbf{1}) + (\mathbf{15},\bar{\mathbf{3}}) + (\overline{\mathbf{15}},\mathbf{3}) + (\mathbf{1},\mathbf{8}) + (\mathbf{24},\mathbf{8}) + (\mathbf{10},\mathbf{6}) + (\overline{\mathbf{10}},\bar{\mathbf{6}}) + (\mathbf{5},\mathbf{3}) + (\bar{\mathbf{5}},\bar{\mathbf{3}}) + \ldots$$

The physical content is the first four terms (122 states). The remaining 374 states are "heavy" in some sense. The question is: what geometric or algebraic mechanism selects exactly 24 + 15 + 15-bar + 8 = 62 representations from 496?

The computation: in the Type I string context (where SO(32) arises as open string Chan-Paton gauge group), a Wilson line $A_\mu = \text{diag}(a_1 I_5, a_2 I_3, a_3 I_3, a_4 I_5)$ in the compact dimensions with $5+3+3+5 = 16$ Chan-Paton labels breaks $SO(32) \to SU(5) \times SU(3) \times SU(3) \times SU(5) \times U(1)^k$. Compute the explicit branching under this Wilson line for the adjoint representation and identify which components survive as massless in the low-energy theory (i.e., which components have vanishing Wilson line contribution to their mass).

For the meson sector to emerge, the relevant massless fields are those where the Wilson line acts trivially — this selects the block-diagonal generators $(\mathbf{24}, \mathbf{1})$ and $(\mathbf{1}, \mathbf{8})$ — plus the off-diagonal generators where both SU(5) and SU(3) generators contribute with opposite Wilson line phases. The diquarks $(\mathbf{15}, \bar{\mathbf{3}})$ live in the $5 \times 3 = 15$ off-diagonal block connecting the first SU(5) with SU(3).

**Explicit task:** Write the $16 \times 16$ matrix representation of a generic SO(32) generator in the block form $\begin{pmatrix} A_5 & C \\ -C^T & B_3 \end{pmatrix}$ where $A_5$ is $5\times 5$ antisymmetric, $B_3$ is $3\times 3$ antisymmetric, and $C$ is $5\times 3$ arbitrary. Count: $\dim A_5 = 10$, $\dim B_3 = 3$, $\dim C = 15$, total $= 28$. This accounts for $SO(16) \subset SO(32)$ — not the full 496. The remaining generators connect the two SO(16) blocks.

Redo with the full $32 \times 32$ antisymmetric matrix interpretation. The $16 + 16 = 32$ Chan-Paton labels split as $(5_f + 3_c) + (5_f + 3_c) = \{q, \bar q\}$ labels. The adjoint antisymmetric $32\times 32$ matrix $M_{IJ}$ decomposes by the $(i,\alpha) = (\text{flavor},\text{color})$ index structure:

$$M_{(i\alpha)(j\beta)} - M_{(j\beta)(i\alpha)}: \quad (ij,\alpha\beta) \in \mathbf{5}\otimes\mathbf{5} \times \mathbf{3}\otimes\mathbf{3}$$

The flavor-symmetric, color-antisymmetric part $\mathbf{15} \times \bar{\mathbf{3}}$ is exactly the diquark. The flavor-antisymmetric part $\overline{\mathbf{10}} \times \bar{\mathbf{3}}$ also appears. Identify the mass of each component from the Wilson line at a given moduli point and determine what value of the Wilson line modulus projects out the $\overline{\mathbf{10}} \times \bar{\mathbf{3}}$ while keeping $\mathbf{15} \times \bar{\mathbf{3}}$ massless.

**Tools:** Standard Lie algebra branching rules, Wilson line calculation in Type I string theory, or equivalently, the orbifold/orientifold projection mechanism.

**Positive result:** There exists a discrete Wilson line (e.g., $\mathbb{Z}_2$ or $\mathbb{Z}_3$ orbifold) that projects out exactly $(\mathbf{24},\mathbf{1}) + (\mathbf{15},\bar{\mathbf{3}}) + (\overline{\mathbf{15}},\mathbf{3}) + (\mathbf{1},\mathbf{8})$ and no others. This simultaneously resolves the $\mathbf{15}$ vs $\overline{\mathbf{10}}$ conflict: the Wilson line projection selects the symmetric block $\mathbf{15}$ and kills $\overline{\mathbf{10}}$ by a phase argument, not by Pauli statistics.

**Negative result:** Any Wilson line that kills $\overline{\mathbf{10}}$ also kills $\mathbf{15}$, since they appear in the same $\mathbf{5} \otimes \mathbf{5}$ decomposition. The 15-vs-10-bar conflict is then fundamental and unresolvable within the Type I framework.

---

## Proposal 4: Fritzsch Texture CKM Calculation with Koide-Determined Masses

**Priority: High**

### What to compute

The quark masses from the sBootstrap program are:
- $m_s = 93.4$ MeV (input)
- $m_c = 13.93 \, m_s = 1301$ MeV (Koide seed prediction)
- $m_b = 4177$ MeV ($v_0$-doubling prediction)
- $m_t = 172760$ MeV (external)
- $m_u = 2.16$ MeV, $m_d = 4.67$ MeV (free)

The task: write the Fritzsch-type up and down mass matrices
$$M_u = \begin{pmatrix} 0 & A_u & 0 \\ A_u^* & 0 & B_u \\ 0 & B_u^* & m_t \end{pmatrix}, \quad M_d = \begin{pmatrix} 0 & A_d & 0 \\ A_d^* & 0 & B_d \\ 0 & B_d^* & m_b \end{pmatrix}$$

with $|A_u|^2 \approx m_u m_c$, $|B_u|^2 \approx m_c m_t$, $|A_d|^2 \approx m_d m_s$, $|B_d|^2 \approx m_s m_b$.

Diagonalize $M_u M_u^\dagger$ and $M_d M_d^\dagger$ to get the left-diagonalizing unitaries $U_L^u$ and $U_L^d$. The CKM matrix is $V = (U_L^u)^\dagger U_L^d$. Extract $|V_{us}|, |V_{cb}|, |V_{ub}|$ and compare to PDG.

This is a calculation that can be done analytically to leading order in the small ratios $\epsilon_u = m_u/m_c$ and $\epsilon_d = m_d/m_s$:
$$|V_{us}| \approx \sqrt{m_d/m_s} - \sqrt{m_u/m_c} \cos\delta + O(\epsilon^{3/2})$$
$$|V_{cb}| \approx \sqrt{m_s/m_b} - \sqrt{m_c/m_t} \cos\delta' + O(\epsilon^{3/2})$$

Evaluate numerically:
- $\sqrt{m_d/m_s} = \sqrt{4.67/93.4} = 0.2236$ (Oakes prediction for $|V_{us}|$, matches PDG $0.2243$)
- $\sqrt{m_s/m_b} = \sqrt{93.4/4177} = 0.1495$
- PDG $|V_{cb}| = 0.0408$

There is a discrepancy of factor 3.7 in $|V_{cb}|$ at leading order. The question is whether the subleading term (with phase $\delta'$) can bring it into agreement.

**Tools:** $3 \times 3$ eigenvalue and eigenvector calculation. Numerically trivial; analytically tractable to two orders.

**Positive result:** The CKM matrix from the Fritzsch texture with sBootstrap masses reproduces $|V_{us}|$ exactly (Oakes), $|V_{cb}|$ within 50%, and gives a rough prediction for $|V_{ub}|$ and $\delta_{CP}$. This establishes the framework as a predictive model for flavor.

**Negative result:** The Fritzsch texture is ruled out — for example, $|V_{cb}|$ cannot be brought to $0.041$ for any choice of phases $\delta, \delta'$. This directs attention to alternative textures (e.g., Georgi-Jarlskog, or anarchic).

**Note on the leading-order discrepancy:** The ratio $|V_{cb}|/\sqrt{m_s/m_b} = 0.041/0.150 \approx 0.27$ suggests there is destructive interference at leading order, which is possible if $\cos\delta' \approx 0.73$. This should be checked explicitly before concluding the texture is ruled out.

---

## Proposal 5: Compute the Bloom Potential from the Instanton-Induced Superpotential

**Priority: Medium**

### What to compute

The bloom mechanism requires moving the vacuum from the seed $\delta = 3\pi/4$ to the physical $\delta \approx 157°$ while preserving $Q \approx 2/3$ and doubling $v_0$. The bloom is argued to be nonperturbative because the Coleman-Weinberg mass is 5 orders of magnitude too small.

The specific proposal: compute the effective potential for $\delta$ from the Affleck-Dine-Seiberg (ADS) superpotential in the SU(3) SQCD electric theory. The ADS superpotential for $N_f = N_c = 3$ is:
$$W_\text{ADS} = c_3 \left(\frac{\det M}{\Lambda^6}\right)^{1/(N_f - N_c)} \to W_{N_f = N_c} = c_3 \frac{(\det M)^n}{\Lambda^{6n}}$$

For $N_f = N_c$ (the relevant case), the ADS superpotential does not exist as a power law; instead, there is a quantum-deformed constraint $\det M - B\bar B = \Lambda^6$. The F-term of the Lagrange multiplier $X$ enforces this constraint, and the effective superpotential on the constrained surface is the meson superpotential $W = \text{Tr}(\hat m M)$ restricted to $\det M = \Lambda^6$.

The bloom in the Koide parametrization corresponds to moving $\delta$ in the parametrization $M_{kk} = (v_0 + r\cos(2\pi k/3 + \delta))^2$. The determinant constraint $\det M = \Lambda^6$ with $M = \text{diag}(M_1, M_2, M_3)$ gives:
$$M_1 M_2 M_3 = \Lambda^6$$

In the Koide parametrization:
$$\prod_k (v_0 + r\cos(2\pi k/3 + \delta))^2 = \Lambda^6$$

This is a constraint on $(v_0, r, \delta)$ that the bloom must satisfy. Compute the determinant constraint explicitly as a function of $\delta$ at fixed $Q = 2/3$ (so $r = \sqrt{2} v_0$) and determine how $v_0(\delta)$ must change to maintain $\det M = \Lambda^6$. Check whether the doubling $v_0(157°)/v_0(135°) \approx 2$ is consistent with the determinant constraint.

**Tools:** Trigonometric algebra, numerical verification.

**Positive result:** The determinant constraint $\prod_k m_k = \Lambda^6$ (or equivalently $m_s m_c m_b = \Lambda_{\rm QCD}^6$) is approximately satisfied with PDG masses ($93.4 \times 1270 \times 4180 = 4.95 \times 10^8$ MeV$^3$ vs $\Lambda^6 = (300 \text{ MeV})^6 = 7.3 \times 10^{14}$ MeV$^6$... these are off by a factor). Check with the predicted masses $(93.4, 1301, 4177)$ and a renormalization-group corrected $\Lambda$.

**Note on the numbers:** The ratio $m_s m_c m_b / \Lambda_{\rm QCD}^6$ involves a mismatch in dimensions (or equivalently, a mismatch in the normalization of $\Lambda$). The ADS constraint $\det M = \Lambda^{N_f}$ in the magnetic dual theory has $\Lambda$ in the magnetic frame. The mapping between electric and magnetic $\Lambda$ via the Seiberg duality relation $\Lambda_e^{N_c} \Lambda_m^{N_f - N_c} = (-1)^{N_f - N_c} \mu^{N_f}$ is the key numerical check.

**Negative result:** The determinant constraint is incompatible with both Koide relations being satisfied simultaneously, ruling out the confining SQCD interpretation of the bloom.

---

## Proposal 6: The Dual Koide as a Seiberg Duality Prediction

**Priority: Medium**

### What to compute

The observed dual Koide $Q(1/m_d, 1/m_s, 1/m_b) = 0.665$ (0.22% from 2/3) is currently an observation without a mechanism. The Seiberg seesaw gives $\langle M_{jj} \rangle \propto 1/m_j$ (meson VEV is inverse proportional to quark mass). If the meson VEVs form a Koide triple, then the dual Koide is a prediction of the Seiberg vacuum.

**Explicit computation:**

Let $\tilde m_j = \langle M_{jj} \rangle = C/m_j$ with $C = \Lambda^2/m_s^{1/3} m_d^{1/3} m_b^{1/3}$ (from the Seiberg constraint $\det M = \Lambda^6$ at the seesaw vacuum). Compute:
$$Q(\tilde m_d, \tilde m_s, \tilde m_b) = Q(1/m_d, 1/m_s, 1/m_b)$$

This is exactly the dual Koide. Now compute: is $Q(1/m_j)$ uniquely close to 2/3 for the down-type quarks, or does it also hold for $(u, c, t)$?

Check $(1/m_u, 1/m_c, 1/m_t)$ with PDG masses. Then compute whether the dual Koide condition $Q(1/m_d, 1/m_s, 1/m_b) = 2/3$ is equivalent (via the Seiberg map $m \to 1/m$) to the original Koide $Q(m_d, m_s, m_b) = 2/3$, or whether it is a genuinely independent condition. The answer:

The Koide condition $Q = 2/3$ is NOT invariant under $m_j \to 1/m_j$ in general. However, the specific seed $(0, m_s, m_c)$ is self-dual because $Q(0, m_s, m_c) = Q(0, 1/m_s, 1/m_c)$ — both equal $\infty/(1+\infty)$ ... actually this needs careful handling of the $m_0 = 0$ limit. Compute the self-duality condition explicitly and check whether the self-dual seed $(0, m_s, m_c)$ is what generates the dual Koide for the full triple.

**Tools:** Algebra, numerical calculation.

**Positive result:** The self-duality of the seed $(0, m_s, m_c)$ under Seiberg duality implies the dual Koide is a theorem (given the seed Koide), not a separate coincidence. This eliminates the dual Koide as a free observation and absorbs it into the existing structure.

**Negative result:** The dual Koide is not derivable from the seed self-duality — it is an independent coincidence, weakening the structural case.

---

## Proposal 7: Compute What the Bloom Does to the Supertrace

**Priority: Medium**

### What to compute

The supertrace $\text{STr}(M^2) = 0$ is an exact identity for O'Raifeartaigh models with canonical Kahler potential. The bloom deforms the mass spectrum from $(0, m_s, m_c)$ to $(-m_s, m_c, m_b)$, changing the supertrace.

Before bloom (seed): $\text{STr} = |0|^2 + m_s^2 + m_c^2 - 2(\text{fermionic masses})^2$. After bloom: $\text{STr} = m_s^2 + m_c^2 + m_b^2 - 2(\text{fermionic masses})^2$.

The question: what is the supertrace of the full physical spectrum in the sBootstrap (after bloom)? Is it zero, or is the nonzero supertrace a measure of the non-holomorphic (Kahler-sector) origin of the bloom?

**Compute explicitly:**

Using the fermion spectrum from or_full_spectrum.md at $t = \sqrt{3}$, the fermionic masses are $(0, m(2-\sqrt{3}), m(2+\sqrt{3}))$. The bosonic spectrum at $t = \sqrt{3}$, $y$ small is approximately $(0, 0, m\sqrt{7-6}, m\sqrt{7+6})$ (from the $\alpha_\pm$ formulas). Compute:
$$\text{STr}(M^4) = \sum_\text{bosons} M_b^4 - 2\sum_\text{fermions} m_f^4$$

This equals $8y^2 m^4$ (constant in $t$, from the pseudomodulus_cw.md calculation). After the bloom, the mass spectrum shifts. Compute $\text{STr}(M^4)$ for the bloomed spectrum $(-m_s, m_c, m_b)$ and determine whether the bloom can be described as a shift in $y$ (the SUSY-breaking parameter) or whether it requires a genuine modification of the Kahler potential.

**Positive result:** The bloom corresponds to a specific shift in the O'Raifeartaigh parameter $y$, with $y_\text{bloom} - y_\text{seed}$ determined by the mass change $\delta m = m_b - m_c$. This means the bloom is a deformation within the O'Raifeartaigh moduli space and does not require a new mechanism.

**Negative result:** The bloomed supertrace is incompatible with any $y$ in the O'Raifeartaigh family — the bloom requires genuinely modifying the Kahler potential, confirming the nonperturbative interpretation.

---

## Proposal 8: The $\mathbf{15}$ vs $\overline{\mathbf{10}}$ Conflict Resolution via the Chiral Ring

**Priority: Medium**

### What to compute

The Pauli theorem is exact for ordinary quantum mechanics but may have a loophole in SUSY: chiral superfields satisfy the chiral ring relation $\phi_i \phi_j = \phi_j \phi_i$ (they commute, not anticommute). The epsilon tensor $\epsilon_{ijk}$ that antisymmetrizes color indices in the SUSY composite $B^{ab} = \epsilon_{ij} Q^{ia} Q^{jb}$ makes the composite antisymmetric in color. But for $N_c = 3$, $\epsilon_{ijk} Q^{ia} Q^{jb} Q^{kc}$ is a baryon, not a diquark — it is symmetric in flavor if $Q$ commutes.

The key question: in the $N_c = 2$ SQCD diquark $B^{ab} = \epsilon_{ij} Q^{ia} Q^{jb}$, is the composition symmetric or antisymmetric in the flavor indices $(a, b)$?

Since $\epsilon_{ij} = -\epsilon_{ji}$ (antisymmetric) and $Q^{ia} Q^{jb} = Q^{jb} Q^{ia}$ (chiral superfields commute), we have:
$$B^{ab} = \epsilon_{ij} Q^{ia} Q^{jb} = -\epsilon_{ji} Q^{ia} Q^{jb} = -\epsilon_{ji} Q^{jb} Q^{ia} = -B^{ba}$$

Wait — this gives $B^{ab} = -B^{ba}$, i.e., antisymmetric in flavor, i.e., $\overline{\mathbf{10}}$. But this is the same result as the fermionic case. The commutativity of chiral superfields does NOT save the $\mathbf{15}$.

**But:** for $N_c = 3$ with $N_f = 5$, the relevant composite is not a diquark but a meson $M^a_b = \tilde Q_a Q^b$ and a baryon $B^{abc} = \epsilon_{ijk} Q^{ia} Q^{jb} Q^{kc}$. The meson $M^a_b$ lives in $\mathbf{5} \otimes \bar{\mathbf{5}} = \mathbf{24} \oplus \mathbf{1}$, not in $\mathbf{5} \otimes \mathbf{5}$. The $\mathbf{15}$ does not appear in the meson sector at all.

**So where does the $\mathbf{15}$ come from in the sBootstrap?** The SO(32) embedding (so32_embedding.md) puts $(\mathbf{15}, \bar{\mathbf{3}})$ in the adjoint as an off-diagonal block between two different $\mathbf{5}_f$ Chan-Paton factors. This is not a diquark composite of the $SU(3)_c$ SQCD fields — it is a state in the string spectrum associated to open strings stretching between two different sets of D-branes.

**The computation:** Identify exactly what field in the string/SQCD dual description corresponds to the $(\mathbf{15}, \bar{\mathbf{3}})$ states. Options:
1. They are fundamental fields of the magnetic dual (not composites, no Pauli issue)
2. They are baryon operators in a different SQCD with $N_c = 2$, $N_f = 5$
3. They arise from the $Sp$ sector (relevant since $Sp(2N)$ theories naturally produce symmetric composites)

**For option 3:** The $Sp(2)$ SQCD with $N_f$ fundamentals has baryons $B^{ab} = \Omega_{ij} Q^{ia} Q^{jb}$ where $\Omega_{ij}$ is the symplectic form (symmetric — $\Omega_{ij} = +\Omega_{ji}$ in appropriate conventions). This gives a baryon symmetric in flavor, i.e., the $\mathbf{15}$ of SU(5). Check: $\Omega_{ij}$ for $Sp(2)$ has the property $\Omega_{12} = +1, \Omega_{21} = +1$, and $B^{ab} = \Omega_{ij} Q^{ia} Q^{jb} = +\Omega_{ji} Q^{ia} Q^{jb} = +\Omega_{ji} Q^{jb} Q^{ia} = +B^{ba}$.

**This is the resolution:** The $\mathbf{15}$ composites are $Sp(2)$ baryons (with symplectic color), not $SU(3)$ diquarks. The $Sp(2)$ confinement produces symmetric flavor composites, naturally giving the $\mathbf{15}$ without violating any Pauli principle for the $SU(3)$ sector.

**Task:** Verify this explicitly. Write the $Sp(2)$ SQCD superpotential for $N_f = 5$, compute the baryon operators, and show they transform as $\mathbf{15}$ of $SU(5)_f$. Identify what gauge group they are color-neutral under and how they couple to the $SU(3)_c$ sector.

**Positive result:** $Sp(2)$ baryons with $N_f = 5$ give the $\mathbf{15}$ automatically. The sBootstrap needs both an $SU(3)_c$ and an $Sp(2)_c$ confining sector. The $Sp(2)$ sector is the same one proposed for the lepton sector — unifying two open problems.

**Negative result:** The $Sp(2)$ baryon is antisymmetric (I made an error in the symplectic form convention above — need to check), or the $Sp(2)$ sector decouples from the $SU(5)_f$ flavor symmetry. This would rule out the $Sp(2)$ resolution and force option 1 (fundamental field in the magnetic dual).

---

## Summary Table

| # | Proposal | Opens | Priority | Method |
|---|----------|-------|----------|--------|
| 1 | Jordan algebra structure of $m_k = \mathfrak{z}_k^2$ | Resolves $\mathbf{15}$ vs $\overline{\mathbf{10}}$ via product structure | High | Algebraic |
| 2 | Two-sector stabilization of pseudo-modulus | Stabilizes $t = \sqrt{3}$ without ad hoc Kahler | High | SUSY algebra |
| 3 | SO(32) $\to$ 54 via Wilson line | Projects physical spectrum, resolves $\overline{\mathbf{10}}$ | High | String/branching |
| 4 | Fritzsch texture + CKM | Direct prediction for CKM from Koide masses | High | Numerical algebra |
| 5 | ADS determinant constraint vs bloom | Tests SQCD interpretation of bloom | Medium | Algebra + numerics |
| 6 | Dual Koide as seed self-duality theorem | Absorbs dual Koide into existing structure | Medium | Algebra |
| 7 | Supertrace after bloom | Tests holomorphic vs non-holomorphic origin | Medium | Algebra |
| 8 | $Sp(2)$ baryons as $\mathbf{15}$ of $SU(5)$ | Unifies diquark problem with lepton sector | Medium | SUSY algebra |

**Recommended execution order:** Proposals 4 and 6 are numerically trivial and should be done first to establish what the program already predicts. Proposals 1 and 8 are closely related (Jordan product vs Grassmann exterior product) and could be combined. Proposal 3 is the highest-stakes structural question. Proposal 2 is the most important if the goal is a complete Lagrangian.

---

*Generated: 2026-03-04*
