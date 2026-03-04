# The Koide formula and quark-lepton mass relations: evidence and statistical assessment

## Abstract

We compile and statistically evaluate a family of mass relations among Standard Model
fermions built on the Koide ratio $Q(m_1,m_2,m_3)\equiv(m_1+m_2+m_3)/(\sqrt{m_1}+\sqrt{m_2}+\sqrt{m_3})^2$.
The original charged-lepton relation $Q(m_e,m_\mu,m_\tau)=2/3$ holds to $0.001\%$
($0.91\sigma$ in $m_\tau$); a Monte Carlo look-elsewhere analysis over all fermion triples
gives $2.8\sigma$.  Extending the formula to signed mass charges
$\mathfrak{z}_k=\pm\sqrt{m_k}$, we find the quark triples
$(-m_s,m_c,m_b)$ and $(m_c,m_b,m_t)$ satisfy $Q\approx 2/3$ to $1.2\%$ and $0.4\%$
respectively, both obligatorily mixing up- and down-type quarks.
A seed triple $(0,m_s,m_c)$ with the O'Raifeartaigh mass ratio
$\sqrt{m_s}/\sqrt{m_c}=2-\sqrt{3}$ is confirmed to $1.2\%$.
The vacuum-charge doubling relation $\sqrt{m_b}=3\sqrt{m_s}+\sqrt{m_c}$ predicts
$m_b=4177$ MeV ($0.1\sigma$ from PDG).
The Koide phase parameter satisfies
$\delta_0\bmod 2\pi/3 \approx 2/9$ to 35 ppm ($3.1\sigma$ with look-elsewhere, $4.3\sigma$ specific).
A meson triple $(-\pi,D_s,B)$ extends $Q\approx 2/3$ to composite hadrons at $0.1\%$.
The dual Koide $Q(1/m_d,1/m_s,1/m_b)=0.665$ ($0.22\%$ from $2/3$) suggests a
Seiberg-dual structure.
The de Vries angle $R=(\sqrt{19}-3)(\sqrt{19}-\sqrt{3})/16=0.22310\ldots$ matches
the on-shell $\sin^2\theta_W$ to $0.13\sigma$,
and the Gatto--Sartori--Tonin relation $\sin\theta_C=\sqrt{m_d/m_s}$ holds to $0.9\%$.
Taken individually, several of these relations have modest significance; taken together,
their interconnections---shared algebraic structure, overlapping mass triplets, and
common Koide phase---point to a single organizing principle.
The companion paper [I] constructs the $\mathcal{N}=1$ Lagrangian that realizes this
principle via an O'Raifeartaigh superpotential with structural invariant $Q=2/3$.

---

## 1. Introduction

### 1.1 The Koide formula

In 1981, Koide derived a relation among the charged lepton masses from a composite model
with $S_3$ permutation symmetry acting on preon constituents [1].  Stripped of its
model-dependent origin, the relation states that the ratio
$$
Q(m_1,m_2,m_3) \;\equiv\; \frac{m_1+m_2+m_3}{\bigl(\sqrt{m_1}+\sqrt{m_2}+\sqrt{m_3}\bigr)^2}
$$
equals $2/3$ when evaluated on the electron, muon, and tau masses.  With current PDG
values [2],
$$
Q(m_e,\,m_\mu,\,m_\tau) \;=\; 0.666661 \pm 0.000007\,,
$$
where the uncertainty is dominated by $\delta m_\tau = 0.12$ MeV.  If the formula is
exact, it predicts $m_\tau = 1776.97$ MeV, consistent with the experimental value
$1776.86 \pm 0.12$ MeV to $0.91\sigma$.

Foot reformulated this as a geometric condition: the vector of square roots of masses
makes a $45^\circ$ angle with the democratic vector $(1,1,1)$ [3].  Equivalently, in the
angular parametrization
$$
\sqrt{m_k} \;=\; \sqrt{M_0}\,\Bigl(1 + \sqrt{2}\,\cos\!\Bigl(\tfrac{2\pi k}{3}+\delta_0\Bigr)\Bigr)\,,
\qquad k=0,1,2\,,
$$
the condition $Q=2/3$ is satisfied identically for any $M_0$ and $\delta_0$, and the
three masses are determined by two parameters [4].  The phase for charged leptons is
$\delta_0 = 0.2222$ (in natural units of $2\pi/3$).

The precision of the lepton relation---better than one part in $10^5$---has prompted
four decades of theoretical attempts to explain it, ranging from discrete symmetries
($S_3$, $A_4$) acting on flavon fields [5,6] to radiative mechanisms that protect
the formula against renormalization group running [7].  Yet no consensus model exists.
The formula has been called an "unexplained coincidence" [8] and a "tantalizing hint"
[9] in roughly equal measure.

### 1.2 The skeptic's objection: look-elsewhere and scheme dependence

Two objections are standard.

**Look-elsewhere.**  With six quarks and three charged leptons forming
$\binom{9}{3} = 84$ triples (or more, allowing signed square roots), a near-hit at
$Q=2/3$ somewhere in the data is not implausible by chance.  We address this head-on.
A Monte Carlo simulation drawing nine masses log-uniformly over the range
$[0.1,\,200{,}000]$ MeV and scanning all 84 triples times 8 sign assignments finds that
a match as precise as the charged leptons occurs in only $0.24\%$ of random mass
spectra.  This is a look-elsewhere corrected significance of $2.8\sigma$ [10].  Not
overwhelming---but the lepton Koide is only the first entry in a longer list.

**Renormalization scheme.**  The Koide formula is stated in terms of pole masses (or
$\overline{\text{MS}}$ masses evaluated at the mass itself, $m_q(m_q)$).  Running to the
GUT scale shifts $Q$ by a few percent for quarks and by less than $0.1\%$ for leptons
[7,11].  This is sometimes used to argue that the relation is an accident of the
infrared.  But it is equally true that many well-established mass relations---the
Gell-Mann--Okubo formula, the Gatto--Sartori--Tonin (GST) relation---are statements about
low-energy masses, not running couplings at a unification scale.  If the Koide relation
originates from boundary conditions on an effective superpotential (as argued in [I]),
then the infrared is exactly where it should hold.

### 1.3 Beyond the leptons: a web of mass relations

The real case for the Koide formula rests not on a single relation but on a web of
interlocking observations.  In 2011, Rodejohann and Zhang [12] noted that the triple
$(m_c, m_b, m_t)$ satisfies $Q = 0.6693$, within $0.4\%$ of $2/3$.  This was surprising
because it mixes quark generations and charge sectors.  In [13], one of us showed that
allowing signed square roots reveals a further triple $(-m_s, m_c, m_b)$ with $Q = 0.675$
($1.2\%$ from $2/3$), and that these triples chain together: the output of one is the
input of the next, covering five of the six quark masses.

The present paper collects all known mass relations of this type---ten observations
spanning charged leptons, quarks, and pseudoscalar mesons---and evaluates each one
statistically.  The observations, in order of discovery, are:

1. **The charged-lepton Koide** [1]:
   $Q(m_e,m_\mu,m_\tau) = 2/3$ to $0.001\%$.  Significance: $0.91\sigma$ in $m_\tau$;
   $2.8\sigma$ after look-elsewhere.

2. **The quark seed** [13]:
   $Q(0,\,m_s,\,m_c) \approx 2/3$, meaning $\sqrt{m_s}/\sqrt{m_c} \approx 2-\sqrt{3}$.
   Observed ratio $0.271$, theoretical $0.268$; deviation $1.2\%$.

3. **The signed quark triple** [13]:
   $Q(-m_s,\,m_c,\,m_b) = 0.675$, deviation $1.24\%$ from $2/3$.

4. **The heavy triple** [12]:
   $Q(m_c,\,m_b,\,m_t) = 0.6693$, deviation $0.40\%$.

5. **The meson triple** [14]:
   $Q(-m_\pi,\,m_{D_s},\,m_B) = 0.6674$, deviation $0.10\%$.

6. **$v_0$-doubling** [this work]:
   $\sqrt{m_b} = 3\sqrt{m_s} + \sqrt{m_c}$ to $0.07\%$, predicting $m_b = 4177$ MeV
   (PDG: $4180 \pm 30$ MeV, $0.1\sigma$).  In isolation, the look-elsewhere probability
   is $24.9\%$ (not significant); the significance comes from its structural role within
   the triple framework.

7. **The phase relation** [15,16]:
   $\delta_0 \bmod 2\pi/3 \approx 2/9$ to 35 ppm.  Significance: $3.1\sigma$ with
   look-elsewhere over 128 simple fractions; $4.3\sigma$ for the specific value $2/9$.

8. **The dual Koide** [this work]:
   $Q(1/m_d,\,1/m_s,\,1/m_b) = 0.665$, $0.22\%$ from $2/3$.
   Suggests a Seiberg-dual seesaw structure $M_j \propto 1/m_j$.

9. **The cyclic echo** [this work]:
   Iterating the Koide equation as $t \to c \to s \to \text{stall} \to c$ produces an
   algebraically exact cycle (agreement to machine epsilon).

10. **The second scaling fails** [this work]:
    $M_2 = 9M_0$, $\delta_2 = 9\delta_0$ does *not* produce $(c,b,t)$ or any physical
    triple.  The pattern stops at one step.  This negative result constrains the
    mechanism.

Two additional relations, while not Koide triples, belong to the same family:

- **The de Vries angle** [17]:
  $R = (\sqrt{19}-3)(\sqrt{19}-\sqrt{3})/16 = 0.22310\ldots$ matches
  $\sin^2\theta_W(\text{on-shell})$ to $0.13\sigma$ in $M_W$.

- **The GST/Oakes relation** [18,19]:
  $\sin\theta_C = \sqrt{m_d/m_s} = 0.224$, matching $|V_{us}| = 0.2243$ to $0.9\%$.

### 1.4 What this paper adds

Previous discussions of these relations appear scattered across individual papers,
conference proceedings, and blog posts.  No systematic statistical assessment exists.
This paper provides three things:

**First**, a self-contained compilation of all the observations, with current PDG values
and error bars, in a single reference.  We state every number with its uncertainty and
pull, not just a percentage deviation.

**Second**, a rigorous look-elsewhere analysis.  For each observation we define the trial
space---how many similar relations could have been examined---and compute corrected
$p$-values.  We are explicit about which observations are significant in isolation
($2.8\sigma$ for the lepton Koide, $3.1\sigma$ for the phase relation) and which derive
their significance only from the structural context ($v_0$-doubling, dual Koide).

**Third**, an argument that the observations are not independent.  The quark seed, the
signed triple, the heavy triple, and the $v_0$-doubling relation share three of the same
masses ($m_s$, $m_c$, $m_b$) and are linked by a common algebraic structure.
An exhaustive scan of all $\binom{6}{3} = 20$ quark triples reveals that *every* triple
close to $Q = 2/3$ obligatorily mixes up-type and down-type quarks [13]; pure
same-charge triples $(d,s,b)$ and $(u,c,t)$ give $Q = 0.73$ and $Q = 0.85$
respectively, both far from $2/3$.  This is not a selection effect: while $90\%$ of
random triples from six quarks are mixed by combinatorics, the $100\%$ rate among
Koide-satisfying triples is a structural feature pointing toward CKM-type mixing.

The companion paper [I] constructs the $\mathcal{N}=1$ supersymmetric Lagrangian that
produces these relations.  The key result there is that the O'Raifeartaigh superpotential
$$
W = f\,\Phi_0 + m\,\Phi_a\Phi_b + g\,\Phi_0\Phi_a^2
$$
produces a Koide seed triple $(0,\,m_-,\,m_+)$ with $Q = 2/3$ exactly when
$gv/m = \sqrt{3}$, where $v = \langle\Phi_0\rangle$ is the pseudo-modulus vacuum
expectation value.  The mass hierarchy $m_s/m_c = (2-\sqrt{3})^2 \approx 0.072$ emerges
from the *minimal* model without additional dynamics: it is the unique mass ratio
consistent with a massless goldstino (the seed zero), the Koide condition, and the single
algebraic constraint on the superpotential couplings.

The present paper is deliberately agnostic about the Lagrangian.  We state the empirical
relations, assess their statistical weight, and leave the dynamical interpretation to [I].
The question this paper asks is precise: *given the data, is there evidence for a common
algebraic structure in the fermion mass spectrum, or can the observations be dismissed as
look-elsewhere artifacts?*

### 1.5 Conventions and mass inputs

Throughout this paper, quark masses refer to $\overline{\text{MS}}$ values at the
self-energy scale $m_q(m_q)$ for $c$ and $b$, and the pole mass for the top quark:
$m_s = 93.4 \pm 5.8$ MeV, $m_c = 1270 \pm 20$ MeV, $m_b = 4180 \pm 30$ MeV,
$m_t = 172.76 \pm 0.30$ GeV [2].  Light quark masses are $m_u = 2.16 \pm 0.49$ MeV,
$m_d = 4.67 \pm 0.48$ MeV.  Lepton masses are pole masses: $m_e = 0.51100$ MeV,
$m_\mu = 105.658$ MeV, $m_\tau = 1776.86 \pm 0.12$ MeV.  Meson masses are PDG
averages: $m_\pi = 139.57$ MeV, $m_{D_s} = 1968.35$ MeV, $m_B = 5279.34$ MeV.

We use natural units $c = \hbar = 1$.  The Koide ratio is denoted $Q$.  We write
$\mathfrak{z}_k = \pm\sqrt{m_k}$ for signed mass charges, with the sign conventions
specified for each triple.  Pulls are quoted as
$(x_{\text{pred}} - x_{\text{obs}})/\sigma_{\text{obs}}$ in units of $\sigma$.

### 1.6 Outline

Section 2 reviews the Koide parametrization and the seed structure.
Section 3 presents the quark chain: the seed, the signed triple, and the heavy triple.
Section 4 discusses the $v_0$-doubling and the overlap prediction.
Section 5 extends the analysis to pseudoscalar mesons.
Section 6 covers the phase relation $\delta_0 \bmod 2\pi/3 \approx 2/9$ and the dual
Koide.
Section 7 presents the de Vries angle and the GST relation.
Section 8 gives the full statistical assessment, including combined significance.
Section 9 discusses scheme dependence, renormalization stability, and negative results.
Section 10 summarizes the conclusions and the connection to [I].

---

### References (partial, to be completed)

[1] Y. Koide, Phys. Rev. Lett. 47, 1241 (1981).

[2] Particle Data Group, R.L. Workman et al., Prog. Theor. Exp. Phys. 2022, 083C01.

[3] R. Foot, arXiv:hep-ph/9402242.

[4] D. Look and M.D. Sheppeard (see [13] and references therein).

[5] Y. Koide, Mod. Phys. Lett. A 5, 2319 (1990).

[6] E. Ma, Phys. Rev. D 73, 077301 (2006).

[7] Z.-z. Xing, H. Zhang, and S. Zhou, Phys. Lett. B 666, 462 (2008).

[8] M. Ratz, Lecture notes (2013).

[9] F. Wilczek, Nature 439, 547 (2006).

[10] Computed in this work; Monte Carlo with $10^4$ trials, see Section 8.

[11] S. Esposito and P. Santorelli, Mod. Phys. Lett. A 10, 3077 (1995).

[12] W. Rodejohann and H. Zhang, Phys. Lett. B 698, 152 (2011).

[13] A. Rivero, arXiv:1111.7232 [hep-ph].

[14] A. Rivero, arXiv:0710.1526 [hep-ph].

[15] P. Zenczykowski, Phys. Rev. D 86, 117303 (2012).

[16] M.D. Sheppeard, PhD thesis (2007).

[17] H. de Vries, unpublished; discussed in A. Rivero, arXiv:hep-ph/0606171.

[18] R. Gatto, G. Sartori, and M. Tonin, Phys. Lett. B 28, 128 (1968).

[19] R.J. Oakes, Phys. Lett. B 29, 683 (1969).
