# Third Spectroscopy with a hint of superstrings

**Author:** A. Rivero
**Affiliation:** Institute for Biocomputation and Physics of Complex Systems, University of Zaragoza
**arXiv:** 0710.1526

## Abstract

Some new regularities found in the pseudoscalar meson octet are reported. They invite to reconsider models where elementary fermions and composite QCD open strings can be grouped in common supermultiplets.

---

## Introduction

Weisskopf coined the term "third spectroscopy" to refer to the subnuclear energy range, stressing the empirical parallel with atomic and nuclear spectroscopies. It is particularly intriguing that the general organization of the QCD sector is patterned in the same four sectors as the GWS sector. The proximity of the masses of the muon and pion — one elementary, one composite — is still unexplained.

---

## Mass Regularities in the Pseudoscalar Octet

With the breaking of flavour SU(3) to SU(2) isospin plus U(1) hypercharge, the pseudoscalar octet has three mass parameters linked via a Gell-Mann-Okubo equation:

$$m_{\eta_8}^2 = \frac{1}{3}(4m_K^2 - m_\pi^2)$$

Using measured values $m_\pi = 134.9766 \pm 0.0006$ MeV and $m_K = 497.648 \pm 0.022$ MeV, we get $m_{\eta_8} = 569.326 \pm 0.026$ MeV.

### New Regularities (main results)

$$\begin{align}
\sqrt{m_e}\sqrt{m_\tau - m_\mu} &= 29.2233 \pm 0.0017 \\
m_\pi - m_\mu &= 29.3182 \pm 0.0006
\end{align}$$

Agreement: **99.6%**

$$\begin{align}
\sqrt{m_\tau}\sqrt{m_\mu - m_e} &= 432.246 \pm 0.024 \\
\sqrt{m_\mu}\sqrt{m_\tau - m_e} &= 433.232 \pm 0.024 \\
m_{\eta_8} - m_\pi &= 434.349 \pm 0.026
\end{align}$$

Agreement: **99.7%**

The first pair was noticed by T.A. Mir (Taarik) during research on mass multiplicities. The last three were suggested by the author.

### De Vries' Isospin Relation

Hans de Vries noticed:

$$\left\vert {m_{\pi^+} \over m_\pi} -1 \right\vert^2 = {m_\mu \over m_Z}$$

With $m_{\pi^+}=139.57018 \pm 0.00035$ MeV and $m_Z=91.1876 \pm 0.0021$ GeV:
- LHS: $0.00115821 \pm .00000049$
- RHS: $0.00115869 \pm .00000003$

Agreement: **99.95%** (discrepancy of order of experimental error).

Foot's rewriting of Koide's relation:

$$\mathbf{(1,1,1)} \measuredangle (\sqrt{m_e},\sqrt{m_\mu},\sqrt{m_\tau}) = 45°$$

---

## Boson/Fermion Degrees of Freedom and the String

### Oriented string sector (mesons / sleptons)

Consider $n$ generations with $r$ quarks of type D and $s$ quarks of type U binding into mesons, classified by $SU(N)$, $N=r+s$. The coincidence of degrees of freedom for fermions and bosons requires:

$$\begin{align}
2n &= r \cdot s \\
4n &= r^2+s^2-1
\end{align}$$

This constrains $r = s \pm 1$ and then $n$ must be a **triangular number**. The extra constraint $n \geq r,s$ forbids the trivial solution. Thus **the minimal solution asks explicitly for three generations**, with exactly one unmatched family ($r = s \pm 1$).

The best scenario presents a **24** representation decomposing in three [super]multiplets, each with an octet of scalars.

### Unoriented string sector (diquarks / squarks)

Two more equations, one for U type and one for D type quarks (DD and UD terminated strings):

$$\begin{align}
2n &= r(r+1)/2 \\
2n &= r \cdot s
\end{align}$$

These ask $2n$ to be in the series of **hexagonal numbers**; the lowest even number is $6$, thus the simplest possibility is $n=3$.

Furthermore, this second pair brings $2s = r+1$, which jointly with $r = s \pm 1$ **fixes uniquely**:

$$s = 2, \quad r = 3, \quad n = 3$$

which is the **exact content of the standard model**: three generations with three D-type and two U-type quarks in the low mass sector.

### The ±4/3 diquarks

This sector produces three extra diquarks (type UU) with charge $\pm 4/3$, living in the **15** representation from $\mathbf{5} \otimes \mathbf{5}$. Some extra argument is needed to dispose of them — e.g., impossibility of grouping them in a doublet of the SU(2) electroweak group.

### SO(32) connection

$r+s=5$ is encouraging because the quantization of such string asks for a $SO(2^5)$ group [Marcus-Sagnotti], and this $SO(32)$ is the unique group we can use in a theory of Type I superstrings.

---

## Caveats

Finding a mechanism to give mass to supermultiplets and then break them so finely is a serious challenge, because almost every work on SUSY starts from the principle of a hard separation between fermions and their spartners. The muon is still near a whole octet of spin 0 particles.

The hint of a SO(32) group raises the doubt about if we must pursue a 4D or a 10D string formalism. A related theme is the finding by Connes et al. of a hidden topological dimension in the structure of the SM fields, so that in some aspect they already inhabit a 10 (mod 8) dimensional space.

---

## Electromagnetic Decay Widths

Another intriguing regularity: all known complete electromagnetic decay widths of neutral particles, with the exception of the Upsilon meson, coincide when scaled with the cube of the mass of the decaying particle. And they coincide with the total decay width of $Z^0$, if scaled with the same rule:

| particle | $(\Gamma_\mathrm{TOT}/M^3)^{-1/2}$ (GeV) | $(\Gamma_\mathrm{EM}/M^3)^{-1/2}$ (GeV) |
|---------|---|---|
| $\pi^0$ | 561 | — |
| $\eta$ | 357 | 548 |
| $\omega$ | — | 591 |
| $\Sigma^0$ | 138 | — |
| $J/\Psi$ | 571 | — |
| $Z^0$ | 551 | — |

The current uncertainty in the half-life of $\pi^0$ allows predicting it via this naive scaling from $Z^0$, the result lying well inside the 1-sigma error bar.
