# An interpretation of scalars in SO(32)

**Author:** Alejandro Rivero
**arXiv:** 2407.05397 [hep-ph]
**Published:** European Physical Journal C 84, 1058 (2024)

## Abstract

The paper proposes an interpretation for the adjoint representation of the SO(32) group to classify the scalars of a generic Supersymmetric Standard Model having just three generations of particles, via a flavour group SU(5). The interpretation emerges from a self-consistency principle for scalar composites. The model uses only color and electric charge, and introduces an additional chiral +4/3 quark per generation.

---

## The Flavour Group in SO(32)

For SO(2n) groups with SU(3) color, the adjoint of SO(32) decomposes into 496 components. Three terms provide scalars matching a minimal three-generation SUSY model:

$$496 = \mathbf{(1,24,1^c)} + \mathbf{[1,15,\bar{3}^c]} + \mathbf{[1,\bar{15},3^c]} + (1,24,8^c) + \ldots$$

SU(5) flavor decomposes under $SU(3) \times SU(2) \times U_2(1)$:

$$15 = (1,3)_{-6} + (3,2)_{-1} + (6,1)_4$$

$$24 = (1,1)_0 + (1,3)_0 + (3,2)_5 + (\bar{3},2)_{-5} + (8,1)_0$$

Electric charge formula:

$$Q = \frac{1}{5}\left(\frac{2}{3}Y_1 - Y_2\right)$$

Alternative decomposition chain:

$$SO(32) \supset SU(16) \times U(1) \supset SU(15) \times U(1) \supset SU(5) \times SU(3)$$

---

## Section 3.1: Turtles and Elephants (Self-Consistency Postulate)

Scalars are treated as composites of quark pairs. With $N$ generations, $k_u$ up-type ("turtle") quarks and $k_d$ down-type ("turtle") quarks, the consistency conditions for squarks are:

$$2N = k_u k_d \quad (8)$$

$$2N = \frac{k_d(k_d+1)}{2} \quad (9)$$

Adding slepton conditions:

$$4N = k_u^2 + k_d^2 - 1 \quad (11)$$

**Unique solution:** $N = 3$, $k_u = 2$, $k_d = 3$

Five "turtle" quarks $(u, c, d, s, b)$ + one "elephant" quark (top) that does not participate in composites.

Three extra scalar quarks of charge $\pm 4/3$ emerge from all turtle combinations.

---

## Discussion on Related Groups

- **SU(15):** Does not produce the same clean 3-generation structure.
- **SU(8):** Does not achieve the same elegance.
- **E8 × E8:** Does not provide the correct decomposition.

Only SO(32) gives the clean $\mathbf{(1,24,1^c)} + \mathbf{[1,15,\bar{3}^c]} + \mathbf{[1,\bar{15},3^c]}$ structure.

---

## The Mass Formula

Classical preon model energy for pair $(q_a, q_b)$:

$$E(q_a, q_b) = (q_a + q_b)^2 K_\Omega \quad (14)$$

For three-preon composites, two conditions:

$$z_1 + z_2 + z_3 = 0 \quad (15)$$

$$z_0^2 = \frac{z_1^2 + z_2^2 + z_3^2}{3} \quad (16)$$

Angular parametrization of abelian charges from SU(3):

$$z_1 = \tfrac{1}{2}\cos\alpha + \tfrac{1}{2\sqrt{3}}\sin\alpha$$
$$z_2 = -\tfrac{1}{\sqrt{3}}\sin\alpha$$
$$z_3 = -\tfrac{1}{2}\cos\alpha + \tfrac{1}{2\sqrt{3}}\sin\alpha$$
$$z_0 = \pm\tfrac{1}{\sqrt{6}}$$

For $\alpha = 0.745821$, the formula recovers exactly $m_e, m_\mu, m_\tau$ using $k = m_e + m_\mu + m_\tau$.

At $\alpha = 0$: exact mass pairs throughout all representations (unbroken SUSY condition).

Squark mass analysis gives:

| quark | mass (MeV) |
|-------|-----------|
| u = c | 313.85 |
| s | 0 |
| d = b | 470.8 |

Producing: $(u,c,t) = (0, 470.8, 1883)$ MeV and $(d,s,b) = (15.8, 313.85, 1553.4)$ MeV.

---

## Conclusions

1. SO(32)'s adjoint representation naturally accommodates three-generation SUSY scalar content.
2. A self-consistency ("turtles all the way down") postulate uniquely determines $N=3$.
3. The decomposition identifies five "light" turtle quarks and one heavy elephant (top).
4. Classical mass formulas produce equal-mass pairs required by unbroken SUSY.
5. Each generation produces two extra $\pm 4/3$ scalar quarks — phenomenologically exotic but theoretically necessary.
6. **"Supersymmetry could be hiding in plain sight, not broken but distorted."**
7. Related groups SO(30), SU(15), SU(8), E8 × E8 are discussed but none achieve the same elegance.
