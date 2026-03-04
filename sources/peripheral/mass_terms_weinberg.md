# Mass terms to break susy-like degeneration

**Author:** Alejandro Rivero (EUPT, University of Zaragoza)
**arXiv:** hep-ph/0606171
**Submitted:** June 15, 2006

## Abstract

A very simple but general operator to break mass degeneration between representations of the Poincare group having spin 1 and 1/2. During this process, a quantity resembling Weinberg's angle emerges, matching experimental values within 0.13 standard deviations.

---

## Mass Terms from Casimir Invariants

Mass operator $M_s^2$ is constructed from Poincaré Casimir invariants $\mathcal{C}_1$ and $\mathcal{C}_2$ with eigenvalues $c_1 = m^2$ and $c_2 = -m^2 s(s+1)$.

Asymptotic constraint (preserving Regge trajectories): $\lim_{s \to \infty} m_s^2 = m$.

The simplest solution:

$$M_s^4 - M_s^2 \mathcal{C}_2 + \mathcal{C}_1 \mathcal{C}_2 = 0$$

Using Pauli matrices to avoid square roots:

$$M_s^2 = \sigma^+ \otimes \mathcal{C}_1 \mathcal{C}_2 + \sigma^- \otimes \mathbf{I} + \frac{\mathbf{I}-\sigma_z}{2} \otimes \mathcal{C}_2$$

---

## The de Vries' Angle (Central Finding)

Hans de Vries discovered that the positive eigenvalues for $s=1/2$ and $s=1$ yield:

$$s_{dV}^2 \equiv 1 - \frac{m_{s=1/2,+}^2}{m_{s=1,+}^2} = 0.22310132\ldots$$

The electroweak Weinberg angle:

$$s_W^2 = 1 - \frac{M_W^2}{M_Z^2} = 0.22306 \pm 0.00033$$

Agreement: $s_{dV}^2 / s_W^2 = 0.9998 \pm 0.0015$ (**0.13 sigma**).

---

## Breaking SUSY Multiplet Degeneracy

For degenerate mass $(m, s_i)$ SUSY multiplet representations, the operator:

$$M_{(s)}^2 = \frac{1}{2}\left(\mathcal{C}_2 + \sqrt{(\mathcal{C}_2)^2 - 4\mathcal{C}_1 \mathcal{C}_2}\right)$$

breaks degeneracy while maintaining asymptotic Regge behavior.

Using $M_Z^2 = (91.1874 \text{ GeV})^2$ as input:
- $m_{s=1/2,+}^2 = (80.3717 \text{ GeV})^2$ (near $M_W$)
- $m_{s=1,-}^2 = -(176.154 \text{ GeV})^2$
- $m_{s=1/2,-}^2 = -(122.384 \text{ GeV})^2$

The vacuum expectation value: $\langle v \rangle / \sqrt{2} = 174.1042 \pm 0.00075$ GeV, with coupling $\lambda_h \approx 1$.

---

## De Broglie / De Vries Classical Connection

De Vries combined de Broglie's relativistic orbit rule with the Landé-Pauli angular momentum substitution. Setting orbital period $T_r = h/(m_0 c^2)$:

$$\frac{\beta^2}{\sqrt{1-\beta^2}} = \sqrt{j(j+1)}$$

Solutions:

$$\beta_{1/2} = \sqrt{\frac{3}{8}(\sqrt{19/3}-1)} = 0.7541414\ldots$$

$$\beta_1 = \sqrt{\sqrt{3}-1} = 0.8555996\ldots$$

$$s_{dV}^2 \equiv 1 - \left(\frac{\beta_{1/2}}{\beta_1}\right)^2 = 0.22310132\ldots$$

---

## Additional Numerical Coincidences

$$\cos\theta = \sinh^{-1}(1) - s_W^2 = 0.2231806$$

$$\sin\theta/\cos\theta = \beta_1^4 - s_W^2 = 0.223112151$$

---

## Main Results

- A Casimir-based mass operator naturally splits spin-1 from spin-1/2 masses.
- The ratio of the resulting masses reproduces the Weinberg angle $\sin^2\theta_W$ to **0.13 sigma** from pure Poincaré group theory.
- The splitting has a classical relativistic mechanics interpretation via de Broglie/de Vries.
- The top Yukawa coupling $\lambda_t = 0.991 \pm 0.013$ aligns with these calculations.
- The electroweak scale may be defined as the point where the running Weinberg angle reaches de Vries' value.
