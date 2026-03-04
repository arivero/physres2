# Complete Spectrum of the O'Raifeartaigh Model

We consider the O'Raifeartaigh superpotential with three chiral superfields $\Phi_0$, $\Phi_1$, $\Phi_2$:
$$W = f\,\Phi_0 + m\,\Phi_1\Phi_2 + g\,\Phi_0\Phi_1^2$$
with real positive parameters $f$, $m$, $g$. The dimensionless coupling ratio is $t \equiv gv/m$, where $v = \langle\phi_0\rangle$ is the pseudo-modulus VEV. A second dimensionless parameter $y \equiv gf/m^2$ controls the SUSY-breaking scale relative to the supersymmetric mass.

The SUSY-breaking vacuum has $\phi_1 = 0$, $\phi_0 = v$ (pseudo-modulus), $\phi_2 = 0$ (flat direction set to zero without loss of generality). The nonvanishing $F$-term is $F_0^* = f$, giving vacuum energy $V_0 = |f|^2$.

## Fermion masses

The fermion mass matrix $\mathcal{W}_{ij} = \partial^2 W/\partial\phi_i\partial\phi_j$ evaluated at the vacuum is
$$\mathcal{W} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 2gv & m \\ 0 & m & 0 \end{pmatrix}.$$
The field $\Phi_0$ decouples completely: its row and column vanish, producing a massless goldstino (the fermionic partner of the SUSY-breaking field). The remaining $2\times 2$ block has characteristic equation
$$\lambda^2 - 2gv\,\lambda - m^2 = 0,$$
with roots
$$m_\pm = \bigl(\sqrt{t^2+1} \pm t\bigr)\,m.$$
The physical (positive) fermion masses are therefore

| State | Mass |
|:---:|:---:|
| Goldstino | $0$ |
| $\psi_-$ | $m_- = (\sqrt{t^2+1} - t)\,m$ |
| $\psi_+$ | $m_+ = (\sqrt{t^2+1} + t)\,m$ |

Two useful identities:
$$m_+\,m_- = m^2, \qquad m_+ + m_- = 2\sqrt{t^2+1}\;m.$$

## Scalar masses

The scalar potential $V = \sum_i |\partial_i W|^2$ gives rise to a $6\times 6$ real mass-squared matrix. Writing each complex scalar as $\phi_j = (\alpha_j + i\beta_j)/\sqrt{2}$, the matrix decomposes into two $3\times 3$ blocks:

$$\mathcal{M}^2_{\text{Re}} = \mathcal{W}^2 + \mathcal{B}, \qquad \mathcal{M}^2_{\text{Im}} = \mathcal{W}^2 - \mathcal{B},$$

where $\mathcal{W}^2$ is the fermion mass-squared matrix
$$\mathcal{W}^2 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 4g^2v^2+m^2 & 2gmv \\ 0 & 2gmv & m^2 \end{pmatrix}$$
and $\mathcal{B}$ is the SUSY-breaking contribution. The only nonvanishing third derivative of the superpotential is $W_{011} = W_{101} = W_{110} = 2g$, so $\mathcal{B}_{ij} = \sum_k F_k^* W_{kij}$ gives
$$\mathcal{B} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 2gf & 0 \\ 0 & 0 & 0 \end{pmatrix}.$$

In each sector, $\Phi_0$ again decouples with eigenvalue zero. The $\Phi_1$--$\Phi_2$ block is a $2\times 2$ matrix whose eigenvalues are found from the trace and determinant. In units of $m^2$, defining $y = gf/m^2$:

$$\frac{M^2_{R,\pm}}{m^2} = (1+2t^2) + y \pm \sqrt{(y+2t^2)^2 + 4t^2}\,,$$
$$\frac{M^2_{I,\pm}}{m^2} = (1+2t^2) - y \pm \sqrt{(y-2t^2)^2 + 4t^2}\,.$$

The full spectrum of six real scalars:

| Real scalar | $M^2/m^2$ |
|:---:|:---:|
| $\alpha_0$ (pseudo-modulus) | $0$ |
| $\beta_0$ (R-axion) | $0$ |
| $\alpha_{-}$ | $(1{+}2t^2)+y - \sqrt{(y{+}2t^2)^2+4t^2}$ |
| $\alpha_{+}$ | $(1{+}2t^2)+y + \sqrt{(y{+}2t^2)^2+4t^2}$ |
| $\beta_{-}$ | $(1{+}2t^2)-y - \sqrt{(y{-}2t^2)^2+4t^2}$ |
| $\beta_{+}$ | $(1{+}2t^2)-y + \sqrt{(y{-}2t^2)^2+4t^2}$ |

The pseudo-modulus $\alpha_0 = \mathrm{Re}(\phi_0)$ and R-axion $\beta_0 = \mathrm{Im}(\phi_0)$ are massless at tree level. The former acquires a one-loop Coleman--Weinberg mass; the latter remains massless unless $R$-symmetry is explicitly broken.

The tree-level vacuum is stable (no tachyonic eigenvalue) provided all $M^2 \geq 0$. The most dangerous eigenvalue is $M^2_{I,-}$; the stability bound constrains the allowed range of $y$ for given $t$.

## Supertrace

The supertrace is
$$\mathrm{STr}\,\mathcal{M}^2 \;=\; \sum_{\text{scalars}} M_s^2 \;-\; 2\sum_{\text{fermions}} m_f^2.$$
The scalar trace receives contributions
$$\mathrm{Tr}(\mathcal{M}^2_{\text{Re}}) + \mathrm{Tr}(\mathcal{M}^2_{\text{Im}}) = \mathrm{Tr}(\mathcal{W}^2+\mathcal{B}) + \mathrm{Tr}(\mathcal{W}^2-\mathcal{B}) = 2\,\mathrm{Tr}(\mathcal{W}^2).$$
The fermion trace is
$$2\sum_f m_f^2 = 2\,\mathrm{Tr}(\mathcal{W}^2).$$
Therefore
$$\boxed{\mathrm{STr}\,\mathcal{M}^2 = 0}$$
identically, independent of $t$, $y$, and all parameters. This is a standard result for O'Raifeartaigh models without gauge interactions: the $\mathcal{B}$-matrix is traceless in the sum over real and imaginary parts, and the $\mathcal{W}^2$ contributions cancel exactly between bosons and fermions.

## Spectrum at $t = \sqrt{3}$

At $t = \sqrt{3}$ we have $\sqrt{t^2+1} = 2$, and the fermion masses become

| State | Mass (units of $m$) | Numerical |
|:---:|:---:|:---:|
| Goldstino | $0$ | $0$ |
| $\psi_-$ | $(2-\sqrt{3})\,m$ | $0.26795\,m$ |
| $\psi_+$ | $(2+\sqrt{3})\,m$ | $3.73205\,m$ |

The four nonzero scalar masses-squared at $t = \sqrt{3}$ (in units of $m^2$, parametric in $y$):

| Scalar | $M^2/m^2$ |
|:---:|:---:|
| $\alpha_0$ | $0$ |
| $\beta_0$ | $0$ |
| $\alpha_\pm$ | $7 + y \pm \sqrt{(y+6)^2 + 12}$ |
| $\beta_\pm$ | $7 - y \pm \sqrt{(y-6)^2 + 12}$ |

The supertrace remains $\mathrm{STr}\,\mathcal{M}^2 = 0$.

## Mass-charge parametrization at $t = \sqrt{3}$

The fermion masses $(0,\; (2{-}\sqrt{3})m,\; (2{+}\sqrt{3})m)$ admit a clean square-root decomposition. Define mass charges $\mathfrak{z}_k = \sqrt{m_k}$ (in units of $\sqrt{m}$):

$$\mathfrak{z}_0 = 0, \qquad \mathfrak{z}_- = \frac{\sqrt{6}-\sqrt{2}}{2}, \qquad \mathfrak{z}_+ = \frac{\sqrt{6}+\sqrt{2}}{2}.$$

One verifies $\mathfrak{z}_-^2 = 2-\sqrt{3}$ and $\mathfrak{z}_+^2 = 2+\sqrt{3}$, and $\mathfrak{z}_-\,\mathfrak{z}_+ = 1$.

The mean mass charge is
$$\bar{z} \equiv \frac{\mathfrak{z}_0 + \mathfrak{z}_- + \mathfrak{z}_+}{3} = \frac{\sqrt{6}}{3}\,.$$

Define deviations $z_k \equiv \mathfrak{z}_k - \bar{z}$, so that $\sum_k z_k = 0$:
$$z_0 = -\frac{\sqrt{6}}{3}, \qquad z_- = \frac{\sqrt{6}-3\sqrt{2}}{6}, \qquad z_+ = \frac{\sqrt{6}+3\sqrt{2}}{6}.$$

The energy-balance parameter is
$$Q \;\equiv\; \frac{\sum_k m_k}{\bigl(\sum_k \sqrt{m_k}\bigr)^2} \;=\; \frac{\sum_k m_k}{(3\bar{z})^2}\,.$$
Substituting the masses:
$$Q = \frac{0 + (2-\sqrt{3}) + (2+\sqrt{3})}{(\sqrt{6})^2} = \frac{4}{6} = \frac{2}{3}\,.$$

**Equivalence**: The condition $Q = 2/3$ is equivalent to $\langle z_k^2\rangle = \bar{z}^2$. The proof is immediate from expanding the variance:
$$\langle z_k^2\rangle = \frac{\sum m_k - 3\bar{z}^2}{3}$$
so $\langle z_k^2\rangle = \bar{z}^2$ iff $\sum m_k = 6\bar{z}^2 = \tfrac{2}{3}(\sum\sqrt{m_k})^2$ iff $Q = 2/3$.

At $t = \sqrt{3}$, we verify: $\langle z_k^2\rangle = \tfrac{2}{3}$ and $\bar{z}^2 = \tfrac{6}{9} = \tfrac{2}{3}$.

**Uniqueness**: The energy-balance parameter for general $t$ is
$$Q(t) = \frac{\sqrt{t^2+1}}{\sqrt{t^2+1}+1}\,,$$
derived from $\sum m_k = 2\sqrt{t^2+1}\,m$ and $(\sum\sqrt{m_k})^2 = 2(\sqrt{t^2+1}+1)\,m$. Setting $Q = 2/3$ requires $\sqrt{t^2+1} = 2$, i.e., $t = \sqrt{3}$. This value is unique.

## Comparison table

### Fermion sector

$$\begin{array}{l|ccc}
 & t = 0 & t = 1 & t = \sqrt{3} \\ \hline
m_0/m & 0 & 0 & 0 \\
m_-/m & 1 & \sqrt{2}-1 \approx 0.4142 & 2-\sqrt{3} \approx 0.2679 \\
m_+/m & 1 & \sqrt{2}+1 \approx 2.4142 & 2+\sqrt{3} \approx 3.7321 \\
m_+\cdot m_-/m^2 & 1 & 1 & 1 \\
\sum m_k^2/m^2 & 2 & 6 & 14 \\
Q & 1/2 & 2-\sqrt{2} \approx 0.5858 & 2/3
\end{array}$$

### Scalar sector (in units of $m^2$, parametric in $y = gf/m^2$)

$$\begin{array}{l|ccc}
 & t=0 & t=1 & t=\sqrt{3} \\ \hline
M^2_{\alpha_0}/m^2 & 0 & 0 & 0 \\
M^2_{\beta_0}/m^2 & 0 & 0 & 0 \\
M^2_{\alpha_-}/m^2 & 1 & 3{+}y{-}\sqrt{(y{+}2)^2{+}4} & 7{+}y{-}\sqrt{(y{+}6)^2{+}12} \\
M^2_{\alpha_+}/m^2 & 1{+}2y & 3{+}y{+}\sqrt{(y{+}2)^2{+}4} & 7{+}y{+}\sqrt{(y{+}6)^2{+}12} \\
M^2_{\beta_-}/m^2 & 1{-}2y & 3{-}y{-}\sqrt{(y{-}2)^2{+}4} & 7{-}y{-}\sqrt{(y{-}6)^2{+}12} \\
M^2_{\beta_+}/m^2 & 1 & 3{-}y{+}\sqrt{(y{-}2)^2{+}4} & 7{-}y{+}\sqrt{(y{-}6)^2{+}12} \\
\mathrm{STr}\,\mathcal{M}^2 & 0 & 0 & 0
\end{array}$$

At $t = 0$ the discriminants collapse: $\sqrt{y^2} = y$ (for $y > 0$), giving the familiar degenerate spectrum with eigenvalues $1$, $1+2y$, $1$, $1-2y$ (plus two zeros). Tree-level stability at $t = 0$ requires $y < 1/2$.

### Summary of $Q$ values

$$\begin{array}{c|ccc}
t & 0 & 1 & \sqrt{3} \\ \hline
Q(t) & 1/2 & 2{-}\sqrt{2} & 2/3
\end{array}$$

The energy-balance condition $Q = 2/3$ selects uniquely $t = \sqrt{3}$, equivalently $gv/m = \sqrt{3}$.

## LaTeX-ready tables

### Fermion spectrum

```latex
\begin{table}[h]
\centering
\begin{tabular}{lccc}
\toprule
& $t=0$ & $t=1$ & $t=\sqrt{3}$ \\
\midrule
$m_0/m$ & $0$ & $0$ & $0$ \\[2pt]
$m_{-}/m$ & $1$ & $\sqrt{2}-1$ & $2-\sqrt{3}$ \\[2pt]
$m_{+}/m$ & $1$ & $\sqrt{2}+1$ & $2+\sqrt{3}$ \\[2pt]
$m_+m_-/m^2$ & $1$ & $1$ & $1$ \\[2pt]
$Q$ & $\tfrac{1}{2}$ & $2-\sqrt{2}$ & $\tfrac{2}{3}$ \\
\bottomrule
\end{tabular}
\caption{Fermion masses in the O'Raifeartaigh model at three values of $t = gv/m$.
The energy-balance parameter $Q = \sqrt{t^2+1}/(\sqrt{t^2+1}+1)$ equals $2/3$
uniquely at $t = \sqrt{3}$.}
\end{table}
```

### Scalar spectrum

```latex
\begin{table}[h]
\centering
\begin{tabular}{lc}
\toprule
Real scalar & $M^2/m^2$ \\
\midrule
$\alpha_0$ (pseudo-modulus) & $0$ \\[2pt]
$\beta_0$ (R-axion) & $0$ \\[2pt]
$\alpha_{\pm}$ & $(1{+}2t^2) + y \pm \sqrt{(y{+}2t^2)^2 + 4t^2}$ \\[2pt]
$\beta_{\pm}$ & $(1{+}2t^2) - y \pm \sqrt{(y{-}2t^2)^2 + 4t^2}$ \\
\bottomrule
\end{tabular}
\caption{Tree-level scalar mass-squared spectrum of the O'Raifeartaigh model,
where $t = gv/m$ and $y = gf/m^2$.}
\end{table}
```

### Mass charges at $t = \sqrt{3}$

```latex
\begin{table}[h]
\centering
\begin{tabular}{lccc}
\toprule
& $k = 0$ & $k = -$ & $k = +$ \\
\midrule
$m_k/m$ & $0$ & $2{-}\sqrt{3}$ & $2{+}\sqrt{3}$ \\[3pt]
$\mathfrak{z}_k/\sqrt{m}$ & $0$ & $\frac{\sqrt{6}-\sqrt{2}}{2}$
  & $\frac{\sqrt{6}+\sqrt{2}}{2}$ \\[5pt]
$z_k/\sqrt{m}$ & $-\frac{\sqrt{6}}{3}$ & $\frac{\sqrt{6}-3\sqrt{2}}{6}$
  & $\frac{\sqrt{6}+3\sqrt{2}}{6}$ \\[5pt]
$z_k^2/m$ & $\frac{2}{3}$ & $\frac{2-\sqrt{3}}{3}$
  & $\frac{2+\sqrt{3}}{3}$ \\
\bottomrule
\end{tabular}
\caption{Mass charges and deviations at $t = \sqrt{3}$. The mean $\bar{z} = \sqrt{6}/3$
satisfies $\bar{z}^2 = \langle z_k^2\rangle = 2/3$, equivalent to the energy-balance
condition $Q = 2/3$.}
\end{table}
```

## Computational verification

All results were verified numerically using Python (NumPy + SymPy):

- Fermion eigenvalues of $\mathcal{W}$ match $0, (\sqrt{t^2+1}\pm t)m$ at $t = 0, 0.5, 1, \sqrt{3}, 2, 3, 5, 10$.
- Scalar eigenvalues of $\mathcal{W}^2 \pm \mathcal{B}$ match closed-form expressions at multiple $(t, y)$ values.
- $\mathrm{STr}\,\mathcal{M}^2 = 0$ verified to machine precision ($< 10^{-14}\,m^2$) at all test points.
- $Q(t) = \sqrt{t^2+1}/(\sqrt{t^2+1}+1)$ matches direct ratio $\sum m_k / (\sum\sqrt{m_k})^2$ at all test points.
- $\mathfrak{z}_\pm^2 = 2 \pm \sqrt{3}$ verified symbolically (exact).
- $\langle z_k^2\rangle = \bar{z}^2 = 2/3$ verified both symbolically and numerically.
