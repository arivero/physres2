# Round 21: Proving the Spectrum — Agent Compartmentalization

## Goal

Strengthen `paper_lagrangian.tex` by:
1. Building a **layered spectrum table** — the centerpiece — showing what each Lagrangian term predicts
2. Fixing the **Kähler pole formulation** (user criticism: "not correctly formulated")
3. Providing an **explicit F-term SUSY breaking demonstration** (user criticism: SQCD section doesn't demonstrate it)
4. Expressing the spectrum in **abelian charge language** ($m_k = \mathfrak{z}_k^2$) throughout
5. Addressing the **soft mass scale** — why $\tilde{m}^2 \sim f_\pi^2$ and what this means

All agents get pure math/physics tasks. No agent sees "sBootstrap", no viXra, no "validate this theory."

---

## Agent 21A (opus): Complete O'Raifeartaigh Scalar + Fermion Spectrum

**Task**: Given the O'Raifeartaigh superpotential
$$W = f\Phi_0 + m\Phi_1\Phi_2 + g\Phi_0\Phi_1^2$$
compute the **complete** spectrum at the SUSY-breaking vacuum ($\phi_1 = 0$, $\phi_0 = v$) as a function of the dimensionless parameter $t = gv/m$.

Deliver:
1. The 6 real scalar masses (3 complex fields → 6 real components). Express the scalar mass-squared matrix and diagonalize it.
2. The 3 fermion masses (eigenvalues of $W_{ij}$).
3. Verify $\text{STr}[M^2] = -2f^2g^2/(m^2 + g^2v^2)$ (the field-independent supertrace).
4. **Specialize to $t = \sqrt{3}$**: give all 9 masses (6 scalar + 3 fermion) numerically in units of $m$.
5. Express the fermion spectrum in **mass-charge parametrization**: define $\mathfrak{z}_k$ such that $m_k = \mathfrak{z}_k^2$, with $\sum \mathfrak{z}_k = $ const. Compute $Q = \sum m_k / (\sum |\mathfrak{z}_k|)^2$.
6. Present results as a LaTeX table.

**Do not**: add labels/tags, reference any specific theory, or speculate about physical applications.

---

## Agent 21B (opus): Kähler Pole Effective Potential — Precise Derivation

**Task**: Consider a single pseudo-modulus $X$ with:
- Non-canonical Kähler: $K(X, \bar{X}) = |X|^2 - \frac{|X|^4}{12\mu^2}$
- Tree-level F-term: $F_X = f$ (constant, from an O'Raifeartaigh model)
- One-loop Coleman-Weinberg contribution from a fermion spectrum with masses $m_\pm = (\sqrt{t^2+1} \pm t)m$, where $t = g|X|/m$.

Derive:
1. $K_{X\bar{X}} = 1 - |X|^2/(3\mu^2)$. The pole is at $|X|_{\text{pole}} = \sqrt{3}\,\mu$.
2. The tree-level potential $V_{\text{tree}}(v) = |f|^2 / K_{X\bar{X}}$ as function of $v = |X|$.
3. The one-loop CW potential $V_{\text{CW}}(v)$ from the O'R spectrum.
4. The total $V_{\text{eff}}(v) = V_{\text{tree}} + V_{\text{CW}}$ and its behavior:
   - At $v = 0$: $V = |f|^2$ (canonical value)
   - As $v \to \sqrt{3}\mu$: $V \to ?$ (does it diverge or not?)
   - Location and nature of the minimum
5. **Key question**: Is the potential minimized just below the pole, or does it diverge at the pole? Under what conditions is there a local minimum?
6. Set $\mu = m/g$. Then the pole is at $v = \sqrt{3}m/g$, i.e., $t = \sqrt{3}$. Verify consistency.
7. Expand $V_{\text{eff}}$ near the pole: $V \approx V_0 + a(v_{\text{pole}} - v)^{-1} + \ldots$

**Do not**: Use "Koide" anywhere. Do not speculate about physical applications.

---

## Agent 21C (opus): Multi-Stage Spectrum Table

**Task**: A theory has 13 free parameters in the fermion sector (9 masses + 4 mixing angles). The following mass relations can be "turned on" one at a time, each reducing the parameter count:

| Level | Relation | Formula | Determines | Remaining free |
|-------|----------|---------|------------|----------------|
| 0 | None | — | — | 13 |
| 1 | O'Raifeartaigh mass ratio | $m_c = (2+\sqrt{3})^2 \times m_s$ | $m_c$ from $m_s$ | 12 |
| 2 | Bion mass relation | $\sqrt{m_b} = 3\sqrt{m_s} + \sqrt{m_c}$ | $m_b$ from $m_s, m_c$ | 11 |
| 3 | Yukawa eigenvalue constraint | $Q(m_c, m_b, m_t) = 2/3$ | $m_t$ from $m_c, m_b$ | 10 |
| 4 | Lepton energy balance | $Q(m_e, m_\mu, m_\tau) = 2/3$ | $m_\tau$ from $m_e, m_\mu$ | 9 |
| 5 | Seesaw-Fritzsch (GST) | $\sin\theta_C = \sqrt{m_d/m_s}$ | $\theta_{12}$ from $m_d, m_s$ | 8 |

where $Q(a,b,c) = (a+b+c)/(\sqrt{a}+\sqrt{b}+\sqrt{c})^2$.

For each level:
1. Compute the predicted values using SM inputs: $m_s = 93.4$ MeV, $m_e = 0.511$ MeV, $m_\mu = 105.66$ MeV, $m_d = 4.67$ MeV.
2. Compare to PDG values. Compute pulls.
3. Express each predicted mass in the **mass-charge parametrization**: $m = \mathfrak{z}^2$, where $\mathfrak{z} = z_0 + z_k$ are abelian charges with $\sum z_k = 0$.
4. Show explicitly how the charges $z_k$ change at each level.
5. Present as a clean LaTeX table suitable for a paper. Two formats: (a) compact summary table, (b) expanded table showing charges.

Use Python to verify all numerical values.

**Do not**: Use "Koide" anywhere. Call $Q = 2/3$ the "energy-balance condition." Do not use internal labels.

---

## Agent 21D (opus): F-term Breaking in ISS Meta-Stable Vacuum

**Task**: In the Intriligator-Seiberg-Shih (ISS) model of meta-stable SUSY breaking:
- $SU(N_c)$ SQCD with $N_f$ flavors, $N_c < N_f < 3N_c/2$
- Magnetic dual: $SU(N_f - N_c)$ with $N_f$ flavors $q, \tilde{q}$ and gauge-singlet meson $\Phi$
- Superpotential: $W = h\,\text{Tr}(q\Phi\tilde{q}) - h\mu^2\,\text{Tr}(\Phi)$

Demonstrate explicitly:
1. The F-term equations and why they cannot all vanish simultaneously (the rank condition).
2. The meta-stable vacuum: $\langle q\rangle = \langle\tilde{q}\rangle = \mu\,\mathbb{1}_{N_c}$, $\langle\Phi\rangle = 0$.
3. The nonzero F-terms: $F_\Phi \neq 0$ for the $(N_f - N_c) \times (N_f - N_c)$ block.
4. The vacuum energy: $V_0 = (N_f - N_c)|h^2\mu^4|$.
5. The pseudo-modulus and its Coleman-Weinberg stabilization.
6. **Specialize to $N_f = N_c = 3$**: the electric theory is at the boundary of the conformal window. Use Seiberg's quantum-modified constraint instead of the ISS magnetic dual. Show that adding quark masses $W_{\text{tree}} = \sum m_i M^i{}_i$ leads to the seesaw vacuum $\langle M^i{}_i \rangle \propto 1/m_i$, and that the Lagrange-multiplier $X$ has $F_X \neq 0$ when the O'Raifeartaigh deformation is added.
7. Present the result as a self-contained derivation suitable for a paper section.

**Do not**: Use "Koide" or "sBootstrap." Use "energy-balance condition" for $Q = 2/3$.

---

## Synthesis Plan

After agents complete:
1. Insert 21C's layered spectrum table into Section 5 of `paper_lagrangian.tex` (or create a new subsection)
2. Replace Section 4.4 (Kähler stabilization) with 21B's precise derivation
3. Strengthen Section 4.5 (O'R structure in SQCD) with 21D's explicit F-term demonstration
4. Add 21A's complete spectrum to Section 8 (Complete particle spectrum)
5. Strengthen mass-charge language throughout using 21A and 21C outputs

## What this round does NOT do
- No updates to `paper_koide.tex`
- No new empirical observations
- No RG analysis (already done, negative result stands)
- No bloom mechanism derivation (remains open problem)
