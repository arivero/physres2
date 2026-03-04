# Folding a pattern

**Author:** Alejandro Rivero

## Abstract

We propose a reorganisation of the standard model and their mesons in order to build supersymmetric multiplets. The presentation is open to improvements to choose the adequate candidates in each recombination.

---

## Koide Chain Mass Predictions

A chain of four Koide equations does a good work to predict SM masses. Gray cells are inputs. All triplets in the sequence meet Koide equation. For $(bcs)$, sign of $\sqrt{m_s}$ is minus.

| particle | PDG 2012 | exact ($y_t=1, y_u=0$) | rotated ($m_e, m_\mu$ inputs) |
|----------|---------|--------|---------|
| t | $173.5 \pm 1.0$ | 174.10 (input) | 173.26 |
| b | $4.18 \pm 0.03$ | 3.64 | 4.197 |
| c | $1.275 \pm 0.025$ | 1.698 | 1.359 |
| $\tau$ | $1.77682(16)$ | 1.698 | 1.776968 |
| s | $95 \pm 5$ | 121.95 | 92.275 |
| $\mu$ | 105.65837 | 121.95 | 105.6584 (input) |
| d | $\sim 4.8$ | 8.75 | 5.32 |
| u | $\sim 2.3$ | 0 (input) | 0.03564 |
| e | 0.5109989 | 0 | 0.5109989 (input) |

---

## Folding the Standard Model into SUSY Multiplets

### Step 1: SM as it is

For each QCD string, consider its fundamental state:

$$\begin{array}{lllllll}
&\nu_1, t_{rgb}& & & & \\
&\nu_2, b_{rgb}& B^+,B_c^+ & bu, bc & bb, bs, bd & \eta_b, B^0, B^0_s, \bar B^0, \bar B^0_s \\
cc,cu & \tau, c_{rgb} & D^+, D_s^+ & sc,dc & & \eta_c, D^0, \bar{D^0} \\
uu & \mu, s_{rgb} & \pi^+, K^+ & su, du & ss, sd, dd & \eta_8, \pi^0, K^0, \bar K^0 \\
&\nu_3, d_{rgb} \\
&e, u_{rgb}
\end{array}$$

### Step 2: Fold t and u (identify by Koide mass level)

$$\begin{array}{lllllll}
&\nu_2, b_{rgb}, e, u_{rgb}& B^+,B_c^+ & bu, bc & bb, bs, bd & \eta_b, B^0, B^0_s, \bar B^0, \bar B^0_s \\
cc,cu & \tau, c_{rgb} & D^+, D_s^+ & sc,dc & & \eta_c, D^0, \bar{D^0} \\
uu & \mu, s_{rgb}, \nu_1, t_{rgb} & \pi^+, K^+ & su, du & ss, sd, dd & \eta_8, \pi^0, K^0, \bar K^0 \\
&\nu_3, d_{rgb}
\end{array}$$

### Steps 3 & 4: Fold quarks and assign neutral mesons

The $\pm 4/3$ diquarks are related to electroweak symmetry breaking via some condensation.

**Final folded pattern (with antiparticles):**

$$\begin{array}{||l|l|llll||}
\hline
cc & \nu_2, b_{rgb}, e, u_{rgb} & B^\pm,B_c^\pm & bu, bc & bs, bd & B^0, B^0_s, \bar B^0, \bar B^0_s \\
cu & \tau, c_{rgb}, \nu_3, d_{rgb} & D^\pm, D_s^\pm & sc,dc & bb,dd & \eta_b, \eta_c, D^0, \bar{D^0} \\
uu & \mu, s_{rgb}, \nu_1, t_{rgb} & \pi^\pm, K^\pm & su, du & ss, sd & \eta_8, \pi^0, K^0, \bar K^0 \\
\hline
\end{array}$$

The point of the grouping on the fermion side is that we were looking at the Koide triplets of the waterfall; any triplet is built by taking a fermion from each of the three lines.

---

## The Eight Simultaneous Koide Equations

A lateral question: can mass values be built so that any of the eight possible triplets is a solution?

### Degenerate seed solutions ($m_u = 0$)

Let $k^+ = 2 + \sqrt{3}$, $k^- = 2 - \sqrt{3}$, $k^+k^- = 1$.

The triplets $(k^+, 8, k^-)$ and $(k^+, 0, k^-)$ are Koide triplets. Using them to build "mass spectra" that saturate all equations.

From the analysis, there are only **three possibilities**:
1. Keep fully degenerate with the trivial triplet $(0, 2-\sqrt{3}, 2+\sqrt{3})$
2. Change one zero to a nonzero from the complementary triplets (8 or 96)
3. Multiply one of the $2-\sqrt{3}$ masses with $k^4$, scaling it to $k^3 = 26+15\sqrt{3}$

In the former case, the $u,b$ pair gains mass (b quark above the other two). In the latter, it should be the $t,s$ pair, separating the top quark — the most allowed by the set of eight simultaneous Koide equations with a massless quark still present.

With mass scale in GeV, the likeness with the SM is apparent:

| column 1 | column 2 | column 3 | column 4 |
|----------|----------|----------|---------|
| $0 \mid 43.6918$ | $0 \mid 3.64098$ | $0 \mid 0$ | $0 \mid 0$ |
| $0.12195 \mid 0.12195$ | $0.12195 \mid 0.12195$ | $0.12195 \mid 0.12195$ | $1.69854 \mid 0.008755$ |
| $1.69854 \mid 1.69854$ | $1.69854 \mid 1.69854$ | $1.69854 \mid 1.69854$ | $0.12195 \mid 0.12195$ |

---

## Resolvent Structure

From the polynomial resolvent (via Mathematica/Maxima), the polynomials produce the four numerical solutions, each with extra real triplets that form zeroless but degenerate solutions. These are inside the continuous spectra of answers, but specially signaled in the resolvent, and connected to the zeroed solutions.

In particular, the strange-charm-bottom triplet:

$$[2-\sqrt 3, 1, 2\sqrt 3-2]$$

and the charm-bottom-top triplet:

$$[1, 2\sqrt 3-2, 7\sqrt 3-2]$$

appear explicitly. So on one side a zeroed solution, and on the other, from the resolvent, a particular zeroless solution that happens to be related to our folding:

$$\begin{array}{|ll|}
\hline
3.64098 & 0 \\
1.69854 & 1.69854 \\
0.12195 & 0.12195 \\
\hline
\end{array}
\quad \Longrightarrow \quad
\begin{array}{|ll|}
\hline
b & u \\
d & c \\
s & t \\
\hline
\end{array}
\quad \Longrightarrow \quad
\begin{array}{|ll|}
\hline
3.640 & 3.640 \\
1.698 & 1.698 \\
0.1219 & 174.1 \\
\hline
\end{array}$$

Note that some of the zeroed solutions were also able to produce good values for the down quark. Perhaps the right symmetry group is beyond $S_4$.
