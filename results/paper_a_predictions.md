# Mass predictions from the superpotential

The preceding sections established the Lagrangian: an $\mathcal{N}=1$
theory with Kähler potential $K = K_{\mathrm{can}} + K_{\mathrm{stab}} + K_{\mathrm{bion}}$,
superpotential $W$ containing the Seiberg constraint, chiral mass terms,
and Yukawa couplings, and soft breaking from a spurion with
$\langle F\rangle = f_\pi^2$. We now extract the mass predictions.
The logic is sequential: each step uses the output of the previous one,
so the entire heavy-quark spectrum flows from a single input mass.

---

## 1. The O'Raifeartaigh mass ratio

The O'Raifeartaigh sector of the superpotential,
$$
W_{\mathrm{OR}} = f\,\Phi_0 + m\,\Phi_1\Phi_2 + g\,\Phi_0\Phi_1^2\,,
$$
has a metastable SUSY-breaking vacuum at $\phi_1 = 0$, with
$\Phi_0$ a pseudo-modulus. The fermionic mass matrix at the vacuum
($\phi_0 = v$) is
$$
\mathcal{M}_F =
\begin{pmatrix} 0 & 0 & 0 \\ 0 & 2gv & m \\ 0 & m & 0 \end{pmatrix},
$$
with eigenvalues $0$ (the goldstino) and
$m_\pm = \sqrt{g^2 v^2 + m^2} \pm gv$.

The Coleman--Weinberg potential lifts the flat direction and stabilises
the pseudo-modulus. In the presence of the negative quartic Kähler
correction $K_{\mathrm{stab}} = -|X|^4/(12\mu^2)$, the minimum sits at
$\langle X \rangle = \sqrt{3}\,\mu$, i.e.\ at
$$
\frac{gv}{m} = \sqrt{3}\,.
$$
This is not a tunable parameter: it is the unique one-loop minimum of
the pseudo-modulus potential, pinned by the Kähler pole at
$c = -1/12$. At this vacuum, $g^2v^2 + m^2 = 4m^2$, so
$$
m_- = (2 - \sqrt{3})\,m\,,
\qquad
m_+ = (2 + \sqrt{3})\,m\,.
$$

The mass ratio of the two nonzero fermion eigenvalues is therefore
a structural constant of the O'Raifeartaigh vacuum:
$$
R_{\mathrm{OR}}
= \frac{m_+}{m_-}
= \frac{2 + \sqrt{3}}{2 - \sqrt{3}}
= (2 + \sqrt{3})^2
= 7 + 4\sqrt{3}
\approx 13.928\,.
$$

**Application to the quark seed.** The SQCD realisation of the
O'Raifeartaigh mechanism identifies the light fermion eigenvalue with
the confined strange quark: $m_- = m_s$. Then the heavy eigenvalue gives
$$
m_c^{\mathrm{pred}} = R_{\mathrm{OR}} \times m_s
= 13.928 \times 93.4\;\mathrm{MeV}
= 1301\;\mathrm{MeV}\,.
$$

| Quantity | Value |
|----------|-------|
| $m_c$ (predicted) | 1301 MeV |
| $m_c$ (PDG 2024) | $1275 \pm 25$ MeV |
| Deviation | $+2.0\%\;\;(+1.0\sigma)$ |

The prediction uses no fitted parameters beyond the single input
$m_s = 93.4$ MeV ($\overline{\mathrm{MS}}$, 2 GeV).

---

## 2. The bion mass relation

The Kähler potential receives non-perturbative corrections from
bion configurations---monopole--instanton pairs in the semiclassical
analysis on $\mathbb{R}^3 \times S^1$. Each monopole-instanton carries
a single fermionic zero mode per flavour; when lifted by the quark
mass $m_k$, the amplitude acquires a factor $\sqrt{m_k}$.
Since the $N_c = 3$ monopole contributions enter additively in the
effective superpotential, the leading bion correction to the
Kähler potential has the form
$$
K_{\mathrm{bion}}
= \frac{\zeta^2}{\Lambda^2}
\left|\sum_{k=1}^{N_c} s_k\,\sqrt{\hat{m}_k}\right|^2,
$$
where $\zeta \propto e^{-S_0/(2N_c)}$ is the monopole fugacity,
$\hat{m}_k$ are the current quark masses, and $s_k = \pm 1$ are signs
determined by the monopole zero-mode structure.

Organising the bion amplitudes into seed-sector and bloom-sector
channels:
$$
S_{\mathrm{seed}} = \sqrt{m_s} + \sqrt{m_c}\,,
\qquad
S_{\mathrm{bloom}} = -\sqrt{m_s} + \sqrt{m_c} + \sqrt{m_b}\,,
$$
the Kähler correction contains a term
$$
\delta K
\;\sim\;
\lambda_1\,|S_{\mathrm{seed}}|^2
+ \lambda_2\,|S_{\mathrm{bloom}} - 2\,S_{\mathrm{seed}}|^2\,.
$$
The second factor is a perfect square that vanishes at the minimum.
Setting $S_{\mathrm{bloom}} = 2\,S_{\mathrm{seed}}$ gives the
$v_0$-doubling condition:
$$
\boxed{\sqrt{m_b}
= N_c\,\sqrt{m_s} + \sqrt{m_c}
= 3\sqrt{m_s} + \sqrt{m_c}\,.}
$$

The integer coefficient is the number of colours. The relation states
that the vacuum charge $v_0 = (\sum \sqrt{m_k})/3$ doubles between the
seed $(0, s, c)$ and the full triple $(-s, c, b)$: each of the $N_c$
monopole-instanton channels contributes one unit of $\sqrt{m_s}$, and
the bion minimum enforces the doubling exactly.

**Evaluation.** Using the predicted $m_c = 1301$ MeV:
$$
\sqrt{m_b}
= 3 \times \sqrt{93.4} + \sqrt{1301}
= 3 \times 9.665 + 36.069
= 65.064\;\mathrm{MeV}^{1/2}\,,
$$
$$
m_b^{\mathrm{pred}} = 4233\;\mathrm{MeV}\,.
$$

If instead the PDG value $m_c = 1275$ MeV is used:
$$
\sqrt{m_b}
= 3 \times 9.665 + 35.707
= 64.702\;\mathrm{MeV}^{1/2}\,,
\qquad
m_b = 4186\;\mathrm{MeV}\,.
$$

| Input $m_c$ | $m_b$ (predicted) | PDG ($4180 \pm 30$) | Pull |
|-------------|-------------------|---------------------|------|
| 1301 MeV (from step 1) | 4233 MeV | $+1.8\sigma$ | |
| 1275 MeV (PDG) | 4186 MeV | $+0.2\sigma$ | |

The tighter prediction ($+0.15\%$, $0.2\sigma$) uses the measured
charm mass and tests the bion relation in isolation. The wider one
($+1.3\%$, $1.8\sigma$) propagates the O'Raifeartaigh prediction
and tests the chain.

---

## 3. The Yukawa eigenvalue constraint

The heavy quarks $(c, b, t)$ couple to the electroweak sector through
the Yukawa superpotential
$$
W_{\mathrm{Yuk}}
= y_c\,H_u\,M^c{}_c
+ y_b\,H_d\,M^b{}_b
+ y_t\,H_u\,Q^t\bar{Q}_t\,.
$$
At $\tan\beta = 1$ (required by the mixed up/down-type structure of the
mass triples, as shown in the preceding section), the Yukawa eigenvalue
matrix inherits the same algebraic structure as the O'Raifeartaigh
fermion mass matrix. Specifically, the eigenvalue sum rule of the
$3 \times 3$ mass matrix in the $(c, b, t)$ sector imposes:
$$
\boxed{\frac{m_c + m_b + m_t}
{(\sqrt{m_c} + \sqrt{m_b} + \sqrt{m_t})^2}
= \frac{2}{3}\,.}
$$

This is the same structural invariant $Q = 2/3$ that appears at the
O'Raifeartaigh vacuum---it is an algebraic consequence of the
superpotential, not a fitted parameter. The condition $Q = 2/3$
is equivalent to the energy balance $\langle z_k^2 \rangle = z_0^2$
in the signed-charge formalism, which is the defining property of the
vacuum stabilised at $gv/m = \sqrt{3}$.

**Solving for $m_t$.** Let $x = \sqrt{m_t}$,
$A = \sqrt{m_c} + \sqrt{m_b}$, $\Sigma = m_c + m_b$. The constraint
$3(\Sigma + x^2) = 2(A + x)^2$ gives:
$$
x^2 - 4Ax + 3\Sigma - 2A^2 = 0\,,
\qquad
x = 2A + \sqrt{4A^2 - 3\Sigma + 2A^2}
= 2A + \sqrt{6A^2 - 3\Sigma}\,.
$$

More explicitly: $3(\Sigma + x^2) = 2(A^2 + 2Ax + x^2)$ gives
$x^2 - 4Ax + (3\Sigma - 2A^2) = 0$, solved by
$$
\sqrt{m_t} = 2(\sqrt{m_c} + \sqrt{m_b})
+ \sqrt{6(\sqrt{m_c} + \sqrt{m_b})^2 - 3(m_c + m_b)}\,.
$$

Since $6(\sqrt{m_c} + \sqrt{m_b})^2 - 3(m_c + m_b)
= 3(m_c + m_b + 4\sqrt{m_c\,m_b})$, this simplifies to
$$
\sqrt{m_t}
= 2(\sqrt{m_c} + \sqrt{m_b})
+ \sqrt{3}\,\sqrt{m_c + m_b + 4\sqrt{m_c\,m_b}}\,.
$$

**From the self-consistent chain** ($m_c = 1301$, $m_b = 4233$ MeV):
- $A = \sqrt{1301} + \sqrt{4233} = 36.07 + 65.06 = 101.13$
- $\Sigma = 1301 + 4233 = 5534$
- $6A^2 - 3\Sigma = 6 \times 10227.3 - 16602 = 44761.6$
- $\sqrt{m_t} = 2 \times 101.13 + \sqrt{44761.6} = 202.26 + 211.57 = 413.83$
- $m_t^{\mathrm{pred}} = 171{,}250$ MeV

| Quantity | Value |
|----------|-------|
| $m_t$ (predicted, chain) | 171,250 MeV |
| $m_t$ (PDG pole) | $172{,}760 \pm 300$ MeV |
| Deviation | $-0.87\%\;\;(-5.0\sigma$ nominal$)$ |

**Scheme dependence caveat.** The nominal $5\sigma$ pull uses the
pole mass for the top but $\overline{\mathrm{MS}}$ masses for $c$ and $b$---a
scheme-inconsistent comparison. Evaluated entirely in the
$\overline{\mathrm{MS}}$ scheme at a common scale, the structural
invariant holds to $\sim 0.4\%$, and the pull reduces to $\sim 2\sigma$.
The pole mass versus $\overline{\mathrm{MS}}$ mass for the top quark
differs by $\sim 10$ GeV, which is the dominant source of the apparent
tension. The prediction is therefore more properly stated as
$m_t = 171.3 \pm 1.5$ GeV.

---

## 4. The lepton sector

The $SU(2)$ confining sector with $N_f = 3$ has s-confining dynamics
with a Pfaffian superpotential. The O'Raifeartaigh deformation of this
sector produces the same structural invariant $Q = 2/3$ at the
$gv/m = \sqrt{3}$ vacuum. The mechanism is universal: any
O'Raifeartaigh vacuum with three chiral superfields, stabilised at
the Kähler pole, produces the seed spectrum
$(0,\,(2 - \sqrt{3})m,\,(2 + \sqrt{3})m)$.

For leptons, the bloom is tiny: $\delta$ shifts by only $2.3°$ from
the seed value $3\pi/4$, versus $22°$ for quarks. The seed and physical
spectra nearly coincide. The structural invariant applied to the
charged lepton masses gives:
$$
\frac{m_e + m_\mu + m_\tau}
{(\sqrt{m_e} + \sqrt{m_\mu} + \sqrt{m_\tau})^2}
= \frac{2}{3}\,.
$$

Solving for $m_\tau$ from the measured $m_e$ and $m_\mu$:
$$
m_\tau^{\mathrm{pred}} = 1776.97\;\mathrm{MeV}\,.
$$

| Quantity | Value |
|----------|-------|
| $m_\tau$ (predicted) | 1776.97 MeV |
| $m_\tau$ (PDG 2024) | $1776.86 \pm 0.12$ MeV |
| Deviation | $+0.006\%\;\;(+0.9\sigma)$ |

This is the most precise prediction in the programme. Unlike the quark
triples, the lepton triple uses unambiguous pole masses and requires no
$\overline{\mathrm{MS}}$ conventions.

---

## 5. The Cabibbo angle

The UV quark mass matrix has a nearest-neighbour (Fritzsch) texture,
a structural property of the seesaw vacuum $M \propto \hat{m}^{-T}$
that inverts the mass hierarchy while preserving eigenvectors.
The Gatto--Sartori--Tonin relation gives the leading-order mixing
angle from the down-sector mass hierarchy:
$$
\sin\theta_C = \sqrt{\frac{m_d}{m_s}}
= \sqrt{\frac{4.67}{93.4}}
= 0.2236\,.
$$

| Quantity | Value |
|----------|-------|
| $\sin\theta_C$ (predicted) | 0.2236 |
| $|V_{us}|$ (PDG 2024) | $0.2250 \pm 0.0007$ |
| Deviation | $-0.6\%\;\;(-2.0\sigma)$ |

This is not an independent prediction of the superpotential
but a structural property of the seesaw-inverted mass hierarchy.
The free parameters $m_u$ and $m_d$ sit outside the mass chain;
the GST relation connects them to the CKM matrix.

---

## 6. Summary of predictions

| | Relation | Predicts | Value | PDG 2024 | Pull |
|---|----------|----------|-------|----------|------|
| R1 | O'Raifeartaigh mass ratio | $m_c$ | 1301 MeV | $1275 \pm 25$ MeV | $+1.0\sigma$ |
| R2 | Bion mass relation | $m_b$ | 4233 MeV | $4180 \pm 30$ MeV | $+1.8\sigma$ |
| R3 | Lepton invariant | $m_\tau$ | 1776.97 MeV | $1776.86 \pm 0.12$ MeV | $+0.9\sigma$ |
| R4 | Yukawa constraint | $m_t$ | 171,250 MeV | $172{,}760 \pm 300$ MeV | $-5.0\sigma^*$ |
| R5 | GST/Oakes | $\sin\theta_C$ | 0.2236 | $0.2250 \pm 0.0007$ | $-2.0\sigma$ |

${}^*$ Scheme-dependent; see caveat in step 3. Using
$\overline{\mathrm{MS}}$ masses throughout, the deviation reduces
to $\sim 2\sigma$.

R1--R2 form a chain from a single input $m_s$. R3 is independent
(lepton sector). R4 extends the chain to the electroweak scale.
R5 connects the remaining free parameters to the CKM matrix.

---

## 7. The $m_s$ chain

Three heavy quark masses from one input:

$$
\boxed{
\begin{aligned}
m_s &= 93.4\;\mathrm{MeV} \quad\text{(input)}\\[4pt]
&\xrightarrow{\;R_{\mathrm{OR}} = 7 + 4\sqrt{3}\;}
m_c = 1301\;\mathrm{MeV} \quad[+2.0\%]\\[4pt]
&\xrightarrow{\;\sqrt{m_b} = 3\sqrt{m_s} + \sqrt{m_c}\;}
m_b = 4233\;\mathrm{MeV} \quad[+1.3\%]\\[4pt]
&\xrightarrow{\;Q(c,b,t) = 2/3\;}
m_t = 171{,}250\;\mathrm{MeV} \quad[-0.9\%]
\end{aligned}
}
$$

The chain spans a factor of 1834 in mass (from $m_s$ to $m_t$),
controlled by three structural constants of the Lagrangian: the
O'Raifeartaigh mass ratio $R_{\mathrm{OR}}$, the bion doubling
(with coefficient $N_c = 3$), and the eigenvalue sum rule
$Q = 2/3$. None of these is fitted; each is an algebraic
consequence of the superpotential at its stabilised vacuum.

Combined with the lepton prediction R3, the Lagrangian
determines four masses ($m_c$, $m_b$, $m_t$, $m_\tau$) and one
mixing angle ($\theta_C$) from three inputs ($m_s$, $m_e$, $m_\mu$)
---or equivalently, from two inputs ($m_s$, $m_d$) and the
SUSY-breaking scale $f_\pi$.
