# Abelian Charges and the Mass-Charge Identity

## LaTeX section for paper

```latex
\subsection{Abelian charges and the mass-charge identity}

Write each fermion mass as the square of an abelian charge,
\begin{equation}
  m_k = \mathfrak{z}_k^2, \qquad \mathfrak{z}_k = z_0 + z_k, \qquad k=1,2,3,
\end{equation}
with $z_0 \equiv \bar{\mathfrak{z}} = \tfrac{1}{3}\sum_k \mathfrak{z}_k$ and the tracelessness
constraint $\sum_k z_k = 0$.
Then the Koide relation $Q = 2/3$ takes the form of a
\emph{charge-variance identity}: expanding
\begin{equation}
  \frac{1}{3}\sum_k z_k^2
  = \frac{1}{3}\sum_k \mathfrak{z}_k^2 - z_0^2
  = \frac{\sum m_k}{3} - \frac{(\sum\sqrt{m_k})^2}{9} ,
\end{equation}
so that
$\langle z_k^2 \rangle = z_0^2$ is equivalent to
$\tfrac{1}{3}\sum m_k = 2z_0^2$, which rearranges to
\begin{equation}
  \frac{\sum m_k}{\bigl(\sum\sqrt{m_k}\bigr)^2} = \frac{2}{3}.
\end{equation}
In statistical language the identity states that the population variance of
the charges equals the squared mean:
$\operatorname{Var}(\mathfrak{z}) = \bar{\mathfrak{z}}^{\,2}$.

\paragraph{Angular parametrisation.}
The constraint $\sum z_k = 0$ restricts $(z_1,z_2,z_3)$ to a
two-dimensional hyperplane in $\mathbb{R}^3$, parametrised as
\begin{equation}\label{eq:zk-angular}
  z_k = r\cos\!\bigl(\delta + \tfrac{2\pi k}{3}\bigr), \qquad k=1,2,3.
\end{equation}
For any $\delta$ the tracelessness $\sum\cos(\delta+2\pi k/3)=0$ is
automatic, and the identity
$\sum\cos^2(\delta+2\pi k/3) = 3/2$ gives
$\langle z_k^2\rangle = r^2/2$.
The charge-variance identity $\langle z_k^2\rangle = z_0^2$ then fixes
\begin{equation}
  r = \sqrt{2}\,z_0,
\end{equation}
and the angle $\delta$---the \emph{bloom angle}---parametrises the family
of mass triples at fixed vacuum charge $z_0$.

Substituting back, the full charges read
$\mathfrak{z}_k = z_0\bigl[1 + \sqrt{2}\cos(\delta + 2\pi k/3)\bigr]$,
which reproduces the standard Koide parametrisation.

\paragraph{The seed angle.}
At the O'Raifeartaigh seed, $m_1=0$ forces $\mathfrak{z}_1=0$, hence
$z_1 = -z_0$.  From (\ref{eq:zk-angular}),
\begin{equation}
  r\cos\delta = -z_0 = -\frac{r}{\sqrt{2}}
  \quad\Longrightarrow\quad
  \cos\delta = -\frac{1}{\sqrt{2}},
  \qquad \delta = \frac{3\pi}{4}.
\end{equation}
Every O'Raifeartaigh seed sits at $\delta=3\pi/4$.
The remaining two eigenvalues satisfy
$t^2 - 4mt + m^2 = 0$
with $z_0^2 = 2m/3$, giving $m_\pm = (2\pm\sqrt{3})\,m$, which is the
spectrum at the vacuum $gv/m=\sqrt{3}$.  Their product
$m_+m_- = m^2$ is the \emph{determinant condition}: it equals
$\det W_{ij}$ restricted to the massive $2\times 2$ block.
The proof that this follows from $\langle z_k^2\rangle = z_0^2$
is direct: at the seed, $z_2+z_3=z_0$ and
$z_2 z_3 = -z_0^2/2$, so
\begin{equation}
  m_+ m_- = \bigl[(z_0{+}z_2)(z_0{+}z_3)\bigr]^2
  = \Bigl[z_0^2 + z_0^2 - \tfrac{z_0^2}{2}\Bigr]^2
  = \tfrac{9z_0^4}{4}
  = m^2.
\end{equation}
Conversely, requiring $m_+m_-=m^2$ with $m_++m_-=4m$ forces
$Q_{\text{charge}}=1$, i.e.\ the Koide identity.
Thus the abelian-charge identity at the seed is the determinant of the
superpotential mass matrix.

\paragraph{Signed charges and the quark bloom.}
The charge $\mathfrak{z}_k$ need not be positive: a negative sign
$\mathfrak{z}_1 = -\sqrt{m_1}$ still gives $m_1 = \mathfrak{z}_1^2>0$
while changing the sign of the field's mass contribution in the
superpotential.  The quark bloom triple $(-\sqrt{m_s},\sqrt{m_c},\sqrt{m_b})$
has vacuum charge
\begin{equation}
  v_0^{\text{bloom}}
  = \frac{-\sqrt{m_s}+\sqrt{m_c}+\sqrt{m_b}}{3}\,.
\end{equation}

\paragraph{Vacuum-charge doubling.}
The seed triple $(0,\sqrt{m_s},\sqrt{m_c})$ has
$v_0^{\text{seed}} = (\sqrt{m_s}+\sqrt{m_c})/3$.
If $\sqrt{m_b} = 3\sqrt{m_s}+\sqrt{m_c}$, substitution gives
\begin{equation}
  v_0^{\text{bloom}}
  = \frac{-\sqrt{m_s}+\sqrt{m_c}+3\sqrt{m_s}+\sqrt{m_c}}{3}
  = \frac{2(\sqrt{m_s}+\sqrt{m_c})}{3}
  = 2\,v_0^{\text{seed}}.
\end{equation}
Conversely, $v_0^{\text{bloom}}=2v_0^{\text{seed}}$ requires
$\sqrt{m_b}=3\sqrt{m_s}+\sqrt{m_c}$,
predicting $m_b = 4177$\,MeV (PDG: $4180\pm30$, 0.1$\sigma$).
The bion relation and the vacuum-charge doubling are
algebraically equivalent statements: one is the other rewritten in the
charge language.
```

## Numerical verification (Python)

All numerical values verified with PDG inputs (MeV):
- m_s = 93.4, m_c = 1270, m_b = 4180, m_t = 172760
- m_e = 0.511, m_mu = 105.658, m_tau = 1776.86

| Quantity | Value | Expected |
|----------|-------|----------|
| Q_Koide(e,mu,tau) | 0.66666082 | 2/3 = 0.66666667 |
| Q_charge leptons | 0.999982 | 1 |
| Q_signed(-s,c,b) | 0.67495 | 2/3 |
| Q_charge quarks | 1.02486 | 1 |
| v0_bloom/v0_seed | 2.0005 | 2 |
| Bion pred m_b | 4177 MeV | 4180 +/- 30 |
| O'R seed angle | 135.00 deg | 3pi/4 = 135 deg |
| O'R det m_+*m_- | 1.000000 | m^2 = 1 |
| <cos^2(theta+2kpi/3)> | 0.500000 | 1/2 (any theta) |
| Bloom angle delta | 158.96 deg | -- |

### Key relations verified:
1. **Koide <=> charge identity**: Q = (1 + Q_charge)/3, so Q = 2/3 iff Q_charge = 1
2. **Variance = squared mean**: Var(fz) = mean(fz)^2 for leptons (to 0.002%)
3. **Bion <=> v0-doubling**: algebraically equivalent (exact, confirmed symbolically)
4. **Seed at 3pi/4**: Any (0, a, b) with Koide has cos(delta) = -1/sqrt(2), exactly
5. **Determinant condition**: m_+*m_- = m^2 requires Q_charge = 1 (unique physical solution)
6. **Q(0,a,b) = 2/3 requires sqrt(a/b) = 2 +/- sqrt(3)**: the O'Raifeartaigh eigenvalue ratio
