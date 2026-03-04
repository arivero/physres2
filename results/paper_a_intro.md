# Three families from SUSY confinement: spectrum and mass predictions of an N=1 SU(3)×SU(2) theory

## Abstract

We construct an $\mathcal{N}=1$ supersymmetric theory with gauge group
$\mathrm{SU}(3) \times \mathrm{SU}(2)$ and show that it uniquely selects
three colors and three generations through a Diophantine bootstrap.
Requiring that the scalar diquarks and mesons of an $\mathrm{SU}(N)$
confining theory fill supermultiplets imposes three interlocking
constraints --- $rs = 2N$, $r(r+1)/2 = 2N$, $r^2 + s^2 - 1 = 4N$ ---
whose only positive-integer solution is $N=3$, $r=3$, $s=2$. The five
light quarks carry an $\mathrm{SU}(5)$ flavor symmetry under which the
composites fill $\mathbf{24}\oplus\mathbf{15}\oplus\overline{\mathbf{15}}$,
a representation that is anomaly-free and asymptotically free.
The quark sector is described by $\mathrm{SU}(3)$ SQCD at
$N_f = N_c = 3$ in the Intriligator--Seiberg--Shih magnetic dual; the
lepton sector by $s$-confining $\mathrm{SU}(2)$ SQCD with $N_f = 3$.
Both sectors admit an O'Raifeartaigh deformation whose
Coleman--Weinberg-stabilized vacuum at $gv/m = \sqrt{3}$ fixes a
structural mass ratio $(2+\sqrt{3})^2 \approx 13.93$ and enforces an
energy-balance invariant $Q = 2/3$. From the strange quark mass alone,
the Lagrangian predicts $m_c = 1301$ MeV (+2.0%), $m_b = 4233$ MeV
(+1.3%), and $m_t = 171\,250$ MeV ($-$0.9%). The lepton sector predicts
$m_\tau = 1776.97$ MeV (+0.006%). In total five of the Standard Model's
thirteen charged-fermion-sector parameters are determined, reducing the
free parameter count from thirteen to eight.

---

## 1. Introduction

Of the roughly two dozen free parameters in the Standard Model, thirteen
belong to the charged fermion sector: nine masses and four CKM mixing
angles. No established principle relates them. The muon is 207 times
heavier than the electron; the top quark outweighs the up quark by a
factor of $8 \times 10^4$. These hierarchies are accommodated by the
Yukawa couplings, but the couplings themselves are inserted by hand. The
flavor problem --- why three generations, and why these masses --- remains
the central open question in particle physics.

This paper proposes a concrete $\mathcal{N}=1$ supersymmetric answer.
The construction begins from the observation that if the scalar partners
of quarks are composite (diquarks and mesons of a confining gauge
theory), then the requirement that they fill common supermultiplets is
not automatically satisfiable. The number of scalar diquarks formed from
$r$ flavors is $r(r+1)/2$; the number of scalar mesons from $r$ flavors
and $s$ antiflavors is $rs$; and these must jointly account for $4N$
bosonic degrees of freedom in a theory with $N$ colors. The resulting
system of Diophantine equations ---
$$rs = 2N, \qquad \frac{r(r+1)}{2} = 2N, \qquad r^2 + s^2 - 1 = 4N$$
--- admits a unique positive-integer solution: $N=3$, $r=3$, $s=2$. The
number of colors, the number of generations (encoded as $r + s - 1 = 4$,
giving four light quark flavors in addition to the strange), and the
representation content of the composites are all fixed by arithmetic.
This is the bootstrap: supersymmetry and confinement together select the
Standard Model's gauge and family structure.

The five light quarks $(u, d, s, c, b)$ then sit in the fundamental of
an $\mathrm{SU}(5)$ flavor symmetry. Under $\mathrm{SU}(5)_{\rm
flavor} \times \mathrm{SU}(3)_{\rm color}$, the composites --- mesons
in the adjoint-plus-singlet and diquarks in the symmetric tensor ---
fill the representations $\mathbf{24} \oplus \mathbf{15} \oplus
\overline{\mathbf{15}}$. This combination is anomaly-free and
asymptotically free with one-loop coefficient $b_0 = 31/3$. The top
quark, with Yukawa coupling $y_t \approx 1$, is external to the
confining dynamics; it does not form bound states. It couples to the
Higgs sector through a standard elementary Yukawa interaction.

The confining dynamics of the quark sector is that of $\mathrm{SU}(3)$
SQCD with $N_f = N_c = 3$ flavors, a theory whose infrared physics was
solved by Seiberg. In the magnetic dual description, the degrees of
freedom are gauge-invariant mesons $M^i{}_j$, baryons $B$ and
$\tilde{B}$, and a Lagrange-multiplier singlet $X$ that enforces the
quantum-modified constraint $\det M - B\tilde{B} = \Lambda^6$. The
Seiberg seesaw gives meson vacuum expectation values $\langle M^j{}_j
\rangle = C/m_j$, inversely proportional to the quark masses. Adding an
O'Raifeartaigh deformation --- the minimal coupling $W \supset fX + m
\phi_a \phi_b + gX\phi_a^2$ --- breaks supersymmetry at a metastable
vacuum. The Coleman--Weinberg potential for the pseudo-modulus $X$,
combined with a non-canonical K\"ahler correction $\delta K = -|X|^4
/(12\mu^2)$, stabilizes the vacuum expectation value at $g\langle X
\rangle / m = \sqrt{3}$.

This stabilization has a remarkable consequence. The fermion mass matrix
at the O'Raifeartaigh vacuum has eigenvalues $(0,\; m_-,\; m_+)$ with
$m_\pm = (\sqrt{t^2 + 1} \pm t)\,m$ and $t = gv/m$. At $t = \sqrt{3}$
the eigenvalues become $(0,\; (2-\sqrt{3})\,m,\; (2+\sqrt{3})\,m)$,
whose ratio $m_+/m_- = (2+\sqrt{3})^2 \approx 13.93$ is a structural
constant independent of all free parameters. This is the
*O'Raifeartaigh mass ratio*. Applied to the quark seed with
$m = m_s = 93.4$ MeV, it predicts $m_c = 13.93 \times m_s = 1301$ MeV,
to be compared with the PDG value $1275 \pm 25$ MeV (2.0% deviation,
$1.0\sigma$).

The same vacuum enforces an algebraic identity: the ratio $Q \equiv
\sum m_k / (\sum \sqrt{m_k})^2$ evaluates to exactly $2/3$ for the
eigenvalue triple $(0, m_-, m_+)$ at $t = \sqrt{3}$. This is the
*energy-balance condition* --- it states that the mean squared
perturbation of the mass charges equals the squared vacuum charge, i.e.,
$\langle z_k^2 \rangle = z_0^2$ in the parametrization $m_k =
(z_0 + z_k)^2$ with $\sum z_k = 0$. It is not an empirical observation
extracted from data; it is a theorem about the O'Raifeartaigh
superpotential. The energy-balance condition serves as the organizing
principle for all subsequent mass predictions.

To pass from the charm quark to the bottom quark, a second mechanism is
needed. On $\mathbb{R}^3 \times S^1$, the SU(3) gauge theory admits
magnetic bions --- correlated monopole--antimonopole pairs associated
with the three roots of the extended Dynkin diagram. Their contribution
to the effective K\"ahler potential modifies the mass eigenvalue relation
to $\sqrt{m_b} = 3\sqrt{m_s} + \sqrt{m_c}$, where the integer
coefficient is the number of colors. We call this the *bion mass
relation*. Using the O'Raifeartaigh-predicted $m_c$, it gives
$m_b = 4233$ MeV (PDG: $4180 \pm 30$ MeV, 1.3% deviation). The
energy-balance invariant is approximately preserved: $Q(-\sqrt{m_s},
\sqrt{m_c}, \sqrt{m_b}) = 0.673$, deviating from $2/3$ by 0.9%.

The chain extends to the top quark through the *Yukawa eigenvalue
constraint*. At $\tan\beta = 1$, the SUSY Yukawa structure requires the
heavy triple $(m_c, m_b, m_t)$ to satisfy $Q = 2/3$ with all positive
square roots. This reduces to a quadratic equation for $\sqrt{m_t}$
whose physical solution, using the predicted $m_c$ and $m_b$ as inputs,
gives $m_t = 171\,250$ MeV (PDG: $172\,760 \pm 300$ MeV, 0.9%
deviation). Thus from the single input $m_s$, the theory produces the
masses of all three heavy quarks.

The lepton sector is governed by the same algebraic structure. The
$\mathrm{SU}(2)$ confining theory with $N_f = 3$ is $s$-confining
($N_f = N_c + 1$); its infrared dynamics is captured by three
antisymmetric composites and a singlet, deformed by an O'Raifeartaigh
superpotential stabilized at the same $gv/m = \sqrt{3}$. Because the
SU(2) theory has no magnetic bion mechanism, the physical spectrum stays
close to the O'Raifeartaigh seed. The energy-balance condition
$Q(m_e, m_\mu, m_\tau) = 2/3$ directly predicts $m_\tau = 1776.97$ MeV
from the inputs $m_e$ and $m_\mu$, matching the PDG value $1776.86 \pm
0.12$ MeV to 0.006%.

Two further predictions follow from the Lagrangian. The identification
of the Seiberg Lagrange multiplier $X$ with the NMSSM singlet $S$
generates a Higgs quartic through the $|F_X|^2$ term. At $\tan\beta =
1$ the tree-level Higgs mass is $m_h = \lambda v/\sqrt{2}$, and
requiring $m_h = 125.35$ GeV fixes $\lambda = 0.72$. The Cabibbo angle
emerges from the tachyonic off-diagonal meson condensate: the ratio of
condensate amplitudes in the $(u,s)$ and $(u,d)$ sectors is exactly
$\sqrt{m_d/m_s}$, reproducing the Gatto--Sartori--Tonin relation
$\sin\theta_C = \sqrt{m_d/m_s} = 0.224$ (PDG: 0.225).

The parameter budget of the theory is as follows. The Standard Model's
charged fermion sector has thirteen free parameters: nine masses and four
CKM angles. In this theory, the O'Raifeartaigh mass ratio determines
$m_c$ from $m_s$; the bion mass relation determines $m_b$ from $m_s$
and $m_c$; the Yukawa eigenvalue constraint determines $m_t$ from $m_c$
and $m_b$; the energy-balance condition determines $m_\tau$ from $m_e$
and $m_\mu$; and the Gatto--Sartori--Tonin relation determines
$\theta_{12}$ from $m_d$ and $m_s$. Five parameters are eliminated. The
remaining eight free parameters are $m_u$, $m_d$, $m_s$, $m_e$, $m_\mu$,
$\theta_{23}$, $\theta_{13}$, and $\delta_{\rm CP}$.

The paper is organized as follows. Section 2 presents the Diophantine
bootstrap and proves the uniqueness of $N=3$, $r=3$, $s=2$. Section 3
describes the $\mathrm{SU}(5)$ flavor structure, the composite
representations, and the anomaly cancellation. Section 4 constructs the
$\mathrm{SU}(3)$ SQCD Lagrangian in the Seiberg magnetic dual, derives
the O'Raifeartaigh vacuum, and proves the energy-balance identity
$Q = 2/3$. Section 5 derives the bion mass relation and the Yukawa
eigenvalue constraint, assembling the complete $m_s$ chain that predicts
$m_c$, $m_b$, and $m_t$. Section 6 presents the lepton sector as
$s$-confining SU(2) SQCD with the same O'Raifeartaigh structure.
Section 7 discusses the NMSSM Higgs sector, the Cabibbo angle, and
the electroweak predictions. Section 8 presents the complete particle
spectrum and parameter budget. Section 9 discusses the open problems ---
the 15 versus $\overline{10}$ Pauli conflict for diquarks, the
nonperturbative origin of the bloom rotation, the mesino mass hierarchy,
and the embedding in $\mathrm{SO}(32)$ --- and Section 10 concludes.
