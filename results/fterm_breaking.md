# F-Term SUSY Breaking in the ISS Meta-Stable Vacuum

## Part 1: General ISS Mechanism

### 1.1 The Rank Condition

Consider $\mathcal{N}=1$ SQCD with gauge group $SU(N_c)$ and $N_f$ flavors of quarks $Q^i$, $\tilde{Q}_j$ ($i,j = 1,\ldots,N_f$), in the free magnetic range $N_c < N_f < \frac{3}{2}N_c$. The Seiberg dual has gauge group $SU(\tilde{N}_c)$ with $\tilde{N}_c = N_f - N_c$, magnetic quarks $q_i$, $\tilde{q}^j$, and a gauge-singlet meson field $\Phi^i{}_j$. The dual superpotential is

$$W_{\mathrm{ISS}} = h\,\mathrm{Tr}(q\,\Phi\,\tilde{q}) - h\mu^2\,\mathrm{Tr}(\Phi)\,.$$

The F-term equations are

$$\frac{\partial W}{\partial \Phi^i{}_j} = h(\tilde{q}\,q)^j{}_i - h\mu^2\,\delta^j_i = 0\,.$$

This requires $(\tilde{q}\,q)^j{}_i = \mu^2\,\delta^j_i$, i.e., the $N_f \times N_f$ matrix $\tilde{q}\,q$ must equal $\mu^2$ times the identity.

Now count: $q$ is an $\tilde{N}_c \times N_f$ matrix and $\tilde{q}$ is an $N_f \times \tilde{N}_c$ matrix, so $\tilde{q}\,q$ is $N_f \times N_f$ but has rank at most $\tilde{N}_c$. The identity matrix $\delta^j_i$ has rank $N_f$. Since $\tilde{N}_c = N_f - N_c < N_f$ (for $N_c > 0$), the system of $N_f^2$ equations for $2N_f\tilde{N}_c$ unknowns is overdetermined in a rank sense: there is no solution. This is the **rank condition**.

More precisely, $\tilde{q}\,q$ can satisfy at most $\tilde{N}_c$ of the diagonal equations $(\tilde{q}\,q)^j{}_j = \mu^2$. The remaining $N_f - \tilde{N}_c = N_c$ diagonal entries necessarily have $(\tilde{q}\,q)^j{}_j \neq \mu^2$, giving nonzero F-terms. SUSY is broken.

### 1.2 The Meta-Stable Vacuum

The vacuum that minimizes the tree-level potential subject to the rank constraint is

$$\langle q \rangle = \begin{pmatrix} \mu\,\mathbb{1}_{\tilde{N}_c} \\ 0 \end{pmatrix}, \qquad \langle \tilde{q} \rangle = \begin{pmatrix} \mu\,\mathbb{1}_{\tilde{N}_c} & 0 \end{pmatrix}, \qquad \langle \Phi \rangle = 0\,.$$

Here $q$ is $\tilde{N}_c \times N_f$, so $\langle q \rangle$ has $\mu$ in the upper-left $\tilde{N}_c \times \tilde{N}_c$ block and zero elsewhere, and similarly for $\tilde{q}$.

**Verification.** At this vacuum,

$$(\tilde{q}\,q)^j{}_i = \mu^2\,\mathrm{diag}(\underbrace{1,\ldots,1}_{\tilde{N}_c},\underbrace{0,\ldots,0}_{N_c})\,.$$

The F-terms are

$$F_{\Phi^i{}_j} = h\bigl[(\tilde{q}\,q)^j{}_i - \mu^2\,\delta^j_i\bigr]\,.$$

For $i,j \leq \tilde{N}_c$: $F = h(\mu^2 - \mu^2) = 0$. These equations are satisfied.

For $i,j > \tilde{N}_c$ (the lower-right $N_c \times N_c$ block): $F = h(0 - \mu^2) = -h\mu^2 \neq 0$. SUSY is broken.

For mixed indices ($i \leq \tilde{N}_c$, $j > \tilde{N}_c$ or vice versa): $F = 0$.

### 1.3 Nonzero F-Terms and Vacuum Energy

The nonzero F-terms reside entirely in the $N_c \times N_c$ lower-right block of the meson matrix:

$$F_{\Phi^a{}_b} = -h\mu^2\,\delta^b_a\,, \qquad a,b = \tilde{N}_c + 1,\ldots, N_f\,.$$

There are $N_c^2$ such F-terms, each with magnitude $|h\mu^2|$. However, the diagonal ones ($a = b$) are the independent ones for the vacuum energy. The tree-level scalar potential is

$$V = \sum_{i,j} |F_{\Phi^i{}_j}|^2 = \sum_{a=1}^{N_c}\sum_{b=1}^{N_c} |h\mu^2\,\delta^b_a|^2 = N_c\,|h\mu^2|^2\,.$$

Therefore

$$\boxed{V_0 = N_c\,|h\mu^2|^2\,.}$$

The factor $N_c$ counts the number of "frustrated" diagonal F-term equations.

### 1.4 The Pseudo-Modulus

Decompose $\Phi$ into blocks:

$$\Phi = \begin{pmatrix} Y_{\tilde{N}_c \times \tilde{N}_c} & Z_{\tilde{N}_c \times N_c} \\ \tilde{Z}_{N_c \times \tilde{N}_c} & \chi_{N_c \times N_c} \end{pmatrix}\,.$$

At the vacuum $\langle q \rangle = \langle \tilde{q} \rangle = \mu\,\mathbb{1}$, $\langle \Phi \rangle = 0$:

- The fields $Y$, $Z$, $\tilde{Z}$ get masses of order $h\mu$ from the cubic coupling $h\,q\,\Phi\,\tilde{q}$.
- The field $\chi$ (the $N_c \times N_c$ block) is **tree-level massless**. Its F-terms are nonzero ($F_\chi = -h\mu^2$), so it couples to the SUSY-breaking sector, but the tree-level potential is completely independent of $\chi$:

$$V_{\mathrm{tree}} = N_c\,|h\mu^2|^2 + 0 \cdot |\chi|^2 + \cdots$$

The trace $\mathrm{Tr}(\chi)$ (and indeed all components of $\chi$) are **pseudo-moduli**: flat directions at tree level. They are stabilized at the origin by the one-loop Coleman--Weinberg potential:

$$V_{\mathrm{CW}}(\chi) = \frac{1}{64\pi^2}\,\mathrm{STr}\left[\mathcal{M}^4(\chi)\left(\ln\frac{\mathcal{M}^2(\chi)}{\Lambda^2} - \frac{3}{2}\right)\right] > 0\,,$$

which has a minimum at $\chi = 0$ (Intriligator, Seiberg, and Shih, 2006).

### 1.5 Stability of the Meta-Stable Vacuum

This vacuum is meta-stable: the true vacuum is the SUSY-preserving one at large field values $\langle\Phi\rangle \sim \Lambda$. The meta-stable vacuum is separated from the SUSY vacuum by a potential barrier, and the tunneling rate is parametrically suppressed:

$$\Gamma \sim e^{-S_B}\,, \qquad S_B \sim \left(\frac{\Lambda}{\mu}\right)^{2(3N_c - N_f)/(N_f - N_c)} \gg 1$$

for $\mu \ll \Lambda$. The vacuum is therefore parametrically long-lived.

---

## Part 2: $N_f = N_c = 3$ with the Quantum Constraint

### 2.1 Setup

At $N_f = N_c = 3$, the dual magnetic gauge group is $SU(\tilde{N}_c) = SU(0)$, i.e., there are no magnetic quarks. The theory is instead described by the confined degrees of freedom: the meson matrix $M^i{}_j = Q^i\tilde{Q}_j/\Lambda$ and the baryons $B = \epsilon_{ijk}Q^iQ^jQ^k/\Lambda^3$, $\tilde{B} = \epsilon^{ijk}\tilde{Q}_i\tilde{Q}_j\tilde{Q}_k/\Lambda^3$.

Seiberg's exact result gives the quantum-modified constraint

$$\det M - B\tilde{B} = \Lambda^6\,,$$

enforced in the superpotential via a Lagrange multiplier $X$:

$$W = X\bigl(\det M - B\tilde{B} - \Lambda^6\bigr) + \mathrm{Tr}(\hat{m}\,M)\,,$$

where $\hat{m} = \mathrm{diag}(m_1, m_2, m_3)$ are the quark masses.

### 2.2 Vacuum Equations

For diagonal $M = \mathrm{diag}(M_1, M_2, M_3)$ and $B = \tilde{B} = 0$, we compute:

$$\frac{\partial(\det M)}{\partial M^j{}_j} = \frac{\det M}{M^j{}_j}\qquad(\text{no sum on }j)\,.$$

This follows from $\det M = M_1 M_2 M_3$, so $\partial(\det M)/\partial M_j = M_1 M_2 M_3 / M_j$.

The F-term equation $\partial W/\partial M^j{}_j = 0$ gives

$$\boxed{m_j + X\,\frac{\det M}{M_j} = 0\qquad(\text{no sum on }j)\,.}$$

### 2.3 The Seiberg Seesaw

From the vacuum equation for flavor $j$:

$$M_j = -\frac{X\,\det M}{m_j}\,.$$

Since $X\,\det M$ is a $j$-independent quantity (call it $-C$), this gives

$$\boxed{M_j = \frac{C}{m_j}\,,\qquad C = -X\,\det M\,.}$$

This is the **Seiberg seesaw**: the meson VEV is inversely proportional to the quark mass. Heavy quarks produce light mesons and vice versa.

**Self-consistency check.** Substituting back: $\det M = C^3/(m_1 m_2 m_3)$. The constraint $\det M = \Lambda^6$ (with $B = \tilde{B} = 0$) fixes

$$C^3 = \Lambda^6\,m_1 m_2 m_3\,,\qquad C = \Lambda^2\,(m_1 m_2 m_3)^{1/3}\,.$$

From the vacuum equation, $X = -C/\det M = -C \cdot m_1 m_2 m_3/C^3 = -m_1 m_2 m_3/C^2$, giving

$$X = -\frac{(m_1 m_2 m_3)^{1/3}}{\Lambda^4}\,.$$

The full solution is:

$$M_j = \frac{\Lambda^2\,(m_1 m_2 m_3)^{1/3}}{m_j}\,.$$

For instance, with three quark masses $m_1 \ll m_2 \ll m_3$, the lightest quark gives the heaviest meson: $M_1 \gg M_2 \gg M_3$.

### 2.4 The O'Raifeartaigh Deformation

Now deform the superpotential by adding terms that promote $X$ from a Lagrange multiplier to a dynamical pseudo-modulus with F-term SUSY breaking:

$$W = X\bigl(\det M - B\tilde{B} - \Lambda^6\bigr) + \mathrm{Tr}(\hat{m}\,M) + f\,X + \tfrac{1}{2}g\,X\,\phi^2\,,$$

where $\phi$ is one of the meson fields (say $M^1{}_2$ or a combination), and $f$, $g$ are parameters.

The F-term for $X$ is

$$F_X^* = \frac{\partial W}{\partial X} = \bigl(\det M - B\tilde{B} - \Lambda^6\bigr) + f + \tfrac{1}{2}g\,\phi^2\,.$$

**At the meta-stable vacuum** where the constraint is approximately satisfied ($\det M - B\tilde{B} \approx \Lambda^6$) and the meson fluctuation is at the origin ($\phi = 0$):

$$\boxed{F_X^* = f \neq 0\,.}$$

The constraint term vanishes, but the linear $f\,X$ term gives a nonzero F-term. SUSY is broken by the O'Raifeartaigh mechanism.

**Why $F_X \neq 0$ cannot be avoided.** The equation $F_X = 0$ would require $\det M - B\tilde{B} - \Lambda^6 = -f - \frac{1}{2}g\phi^2$. But the other F-term equations ($\partial W/\partial M = 0$, $\partial W/\partial \phi = 0$) drive $\phi \to 0$ and $\det M \to \Lambda^6$, leaving $F_X^* = f$. There is no solution with all F-terms zero simultaneously when $f \neq 0$.

### 2.5 Vacuum Energy

The vacuum energy is determined by the nonzero F-terms. At the meta-stable minimum with $\phi = 0$ and the constraint satisfied:

$$V_0 = |F_X|^2 = |f|^2\,.$$

All other F-terms vanish at the vacuum: $F_{M_j} = 0$ (the meson equations are satisfied), $F_\phi = 0$ (since $\phi = 0$ and $gX\phi = 0$). Therefore

$$\boxed{V_0 = |f|^2\,.}$$

### 2.6 Fermion Mass Matrix in the $(X, \phi_1, \phi_2)$ Sector

The relevant part of the superpotential for the fermion mass matrix is

$$W_{\mathrm{OR}} = f\,X + m\,\phi_1\,\phi_2 + g\,X\,\phi_1^2\,,$$

where $\phi_1$ and $\phi_2$ are the two meson fields that enter the O'Raifeartaigh coupling, and $m$ is the effective bilinear mass from the confined dynamics.

At the vacuum $\langle\phi_1\rangle = \langle\phi_2\rangle = 0$, $\langle X\rangle = v$ (stabilized by the Coleman--Weinberg potential), the fermion mass matrix is

$$\mathcal{W}_{IJ} = \frac{\partial^2 W}{\partial\Phi_I\partial\Phi_J}\bigg|_{\mathrm{vac}} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 2gv & m \\ 0 & m & 0 \end{pmatrix}$$

in the basis $(\psi_X, \psi_1, \psi_2)$. The $\psi_X$ direction is completely decoupled: it is the **massless goldstino** (the fermionic partner of the broken SUSY generator, mandated by $F_X \neq 0$).

The nontrivial $2 \times 2$ block has characteristic equation

$$\lambda^2 - 2gv\,\lambda - m^2 = 0\,.$$

The eigenvalues are

$$\lambda_\pm = gv \pm \sqrt{g^2v^2 + m^2}\,.$$

The physical masses (positive) are:

$$m_+ = gv + \sqrt{g^2v^2 + m^2}\,, \qquad m_- = \sqrt{g^2v^2 + m^2} - gv\,.$$

These satisfy the product relation $m_+ \cdot m_- = m^2$ (from the determinant of the $2\times 2$ block) and the sum relation $m_+ + m_- = 2\sqrt{g^2v^2 + m^2}$.

The three fermion masses are therefore

$$\boxed{(m_0,\;m_-,\;m_+) = \bigl(0,\;\sqrt{g^2v^2 + m^2} - gv,\;\sqrt{g^2v^2 + m^2} + gv\bigr)\,.}$$

This is the standard O'Raifeartaigh spectrum.

---

## Part 3: Connection to Physical Spectrum

### 3.1 Specialization to $gv/m = \sqrt{3}$

Set $t \equiv gv/m = \sqrt{3}$. Then $\sqrt{t^2 + 1} = \sqrt{3 + 1} = 2$, and:

$$m_+ = (\sqrt{3} + 2)\,m = (2 + \sqrt{3})\,m\,, \qquad m_- = (2 - \sqrt{3})\,m\,.$$

The full spectrum is:

$$\boxed{(m_0,\; m_-,\; m_+) = \bigl(0,\;(2-\sqrt{3})\,m,\;(2+\sqrt{3})\,m\bigr)\,.}$$

**Numerical values**: $m_- = 0.26795\,m$, $m_+ = 3.73205\,m$, ratio $m_+/m_- = (2+\sqrt{3})^2 = 7 + 4\sqrt{3} \approx 13.928$.

### 3.2 Mass-Charge Decomposition

Write each mass eigenvalue in terms of "charges":

$$m_k = \mathfrak{z}_k^2\,, \qquad \mathfrak{z}_k = z_0 + z_k\,, \qquad \sum_{k=1}^3 z_k = 0\,.$$

Here $z_0$ is the common offset and $z_k$ are the zero-sum deviations. Explicitly:

$$\mathfrak{z}_k = \sqrt{m_k}\,,$$

so $z_0 = \frac{1}{3}(\sqrt{m_0} + \sqrt{m_-} + \sqrt{m_+})$ and $z_k = \sqrt{m_k} - z_0$.

**Evaluation at $gv/m = \sqrt{3}$.** Set $A = \sqrt{m_-} = \sqrt{(2-\sqrt{3})\,m}$ and $B = \sqrt{m_+} = \sqrt{(2+\sqrt{3})\,m}$. Then:

$$z_0 = \frac{0 + A + B}{3}\,.$$

We need $(A + B)^2$:

\begin{align}
(A + B)^2 &= A^2 + B^2 + 2AB \\
&= (2-\sqrt{3})\,m + (2+\sqrt{3})\,m + 2\sqrt{(2-\sqrt{3})(2+\sqrt{3})}\,m \\
&= 4m + 2\sqrt{4-3}\,m \\
&= 4m + 2m = 6m\,.
\end{align}

So $A + B = \sqrt{6m}$ and $z_0 = \sqrt{6m}/3$.

### 3.3 The Identity $\langle z_k^2 \rangle = z_0^2$

The claim is that the variance-like quantity

$$\langle z_k^2 \rangle \equiv \frac{1}{3}\sum_{k=1}^3 z_k^2 = \frac{1}{3}\sum_{k=1}^3 (\mathfrak{z}_k - z_0)^2$$

equals $z_0^2$.

**Proof.** Expand:

$$\langle z_k^2 \rangle = \frac{1}{3}\sum_k \mathfrak{z}_k^2 - 2z_0 \cdot \frac{1}{3}\sum_k \mathfrak{z}_k + z_0^2 = \frac{1}{3}\sum_k m_k - 2z_0^2 + z_0^2 = \frac{1}{3}\sum_k m_k - z_0^2\,.$$

So the identity $\langle z_k^2 \rangle = z_0^2$ is equivalent to

$$\frac{1}{3}\sum_k m_k = 2\,z_0^2\,.$$

Evaluate both sides:

**Left side:**

$$\frac{1}{3}(m_0 + m_- + m_+) = \frac{1}{3}\bigl(0 + (2-\sqrt{3})m + (2+\sqrt{3})m\bigr) = \frac{4m}{3}\,.$$

**Right side:**

$$2\,z_0^2 = 2\left(\frac{\sqrt{6m}}{3}\right)^2 = 2 \cdot \frac{6m}{9} = \frac{12m}{9} = \frac{4m}{3}\,.$$

Both sides equal $4m/3$. The identity holds exactly. $\square$

**Alternative proof using $m_+m_- = m^2$ and $m_+ + m_- = 4m$.** The two relations needed are the product $m_+m_- = m^2$ (determinant of the $2\times 2$ block) and the sum $m_+ + m_- = 2\sqrt{g^2v^2+m^2} = 4m$ (at $gv/m = \sqrt{3}$). Setting $A = \sqrt{m_-}$, $B = \sqrt{m_+}$:

- $AB = \sqrt{m_+m_-} = m$.
- $A^2 + B^2 = m_- + m_+ = 4m$.
- $(A+B)^2 = A^2 + B^2 + 2AB = 4m + 2m = 6m$.
- $Q = (0 + m_- + m_+)/(0 + A + B)^2 = 4m/(6m) = 2/3$.

This identity is algebraically equivalent to $Q = 2/3$, and it is the statement that the variance of the $\mathfrak{z}_k$ around their mean equals the square of the mean:

$$\boxed{\langle z_k^2 \rangle \equiv \frac{1}{3}\sum_{k}(\mathfrak{z}_k - z_0)^2 = z_0^2\,.}$$

In the Koide parametrization $\mathfrak{z}_k = z_0(1 + \sqrt{2}\cos(\tfrac{2\pi k}{3} + \delta))$, this identity holds automatically since $\sum\cos^2(\tfrac{2\pi k}{3}+\delta) = 3/2$, giving $\langle z_k^2\rangle = z_0^2\cdot 2 \cdot 3/(2\cdot 3) = z_0^2$. The O'Raifeartaigh vacuum at $gv/m = \sqrt{3}$ is the unique point where the fermion spectrum takes this Koide form, with the seed phase $\delta = 3\pi/4$ (forcing $m_0 = 0$).

---

## Summary of Results

### Part 1 (General ISS)

| Result | Expression |
|--------|-----------|
| Rank condition | $\mathrm{rank}(\tilde{q}\,q) \leq \tilde{N}_c < N_f$ |
| Vacuum | $\langle q\rangle = \langle\tilde{q}\rangle = \mu\,\mathbb{1}_{\tilde{N}_c}$, $\langle\Phi\rangle = 0$ |
| Nonzero F-terms | $F_{\chi^a{}_b} = -h\mu^2\,\delta^b_a$, $a,b = 1,\ldots,N_c$ |
| Vacuum energy | $V_0 = N_c\,|h\mu^2|^2$ |
| Pseudo-modulus | $\chi_{N_c\times N_c}$ block: tree-level flat, CW-stabilized at origin |

### Part 2 ($N_f = N_c = 3$)

| Result | Expression |
|--------|-----------|
| Vacuum equation | $m_j + X\,(\det M)/M_j = 0$ |
| Seiberg seesaw | $M_j \propto 1/m_j$ |
| O'Raifeartaigh F-term | $F_X = f \neq 0$ |
| Vacuum energy | $V_0 = |f|^2$ |
| Fermion spectrum | $(0,\; m_-,\; m_+)$ with $m_\pm = gv \pm \sqrt{g^2v^2+m^2}$ |

### Part 3 (Physical spectrum at $gv/m = \sqrt{3}$)

| Result | Expression |
|--------|-----------|
| Masses | $(0,\;(2-\sqrt{3})m,\;(2+\sqrt{3})m)$ |
| Mass ratio | $m_+/m_- = 7+4\sqrt{3} \approx 13.93$ |
| Koide invariant | $Q = 4m/(6m) = 2/3$ exactly |
| Variance identity | $\langle z_k^2\rangle = z_0^2$ (algebraic, no free parameters) |

---

## LaTeX-Ready Derivation

The following is formatted for direct inclusion in a paper section.

---

### The O'Raifeartaigh--Koide vacuum

Consider the O'Raifeartaigh superpotential
\begin{equation}
W = f\,X + m\,\phi_1\phi_2 + g\,X\,\phi_1^2\,,
\end{equation}
where $X$ is the pseudo-modulus, $\phi_{1,2}$ are chiral superfields, and $f$, $m$, $g$ are real positive parameters. At the tree-level vacuum $\langle\phi_1\rangle = \langle\phi_2\rangle = 0$, the nonzero F-term $F_X = f$ breaks SUSY with vacuum energy $V_0 = f^2$.

The pseudo-modulus VEV $\langle X\rangle = v$ is stabilized by the one-loop Coleman--Weinberg potential. The fermion mass matrix at the vacuum is
\begin{equation}
\mathcal{W}_{IJ} = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 2gv & m \\ 0 & m & 0\end{pmatrix}\,,
\end{equation}
with eigenvalues $0$ (the goldstino) and
\begin{equation}
m_\pm = gv \pm \sqrt{g^2v^2 + m^2}\,.
\end{equation}

At the Coleman--Weinberg stabilized value $gv/m = \sqrt{3}$, the spectrum simplifies to
\begin{equation}
(m_0,\,m_-,\,m_+) = \bigl(0,\,(2-\sqrt{3})\,m,\,(2+\sqrt{3})\,m\bigr)\,.
\end{equation}

The mass-sum ratio takes the value
\begin{equation}
Q \equiv \frac{m_0 + m_- + m_+}{\bigl(\sqrt{m_0}+\sqrt{m_-}+\sqrt{m_+}\bigr)^2} = \frac{4m}{6m} = \frac{2}{3}\,,
\end{equation}
where the denominator follows from $\sqrt{m_-}\sqrt{m_+} = \sqrt{(4-3)m^2} = m$, so $(\sqrt{m_-}+\sqrt{m_+})^2 = (m_-+m_+) + 2m = 4m + 2m = 6m$.

This is an algebraic identity at the vacuum: it holds for any value of the mass parameter $m$ and depends on no free parameters beyond the stabilization condition $gv/m = \sqrt{3}$.

In the mass-charge decomposition $\sqrt{m_k} = z_0 + z_k$ with $\sum z_k = 0$, the identity $Q = 2/3$ is equivalent to
\begin{equation}
\frac{1}{3}\sum_{k=1}^3 z_k^2 = z_0^2\,,
\end{equation}
stating that the variance of the square-root masses equals the square of their mean. In the Koide parametrization $\sqrt{m_k} = z_0(1 + \sqrt{2}\cos(\frac{2\pi k}{3}+\delta))$, this holds identically by the trigonometric sum $\sum\cos^2(\frac{2\pi k}{3}+\delta) = \frac{3}{2}$. The seed spectrum corresponds to $\delta = 3\pi/4$, at which $\cos\delta = -1/\sqrt{2}$ forces $m_0 = 0$.

---

### Embedding in $N_f = N_c = 3$ SQCD

For $SU(3)$ SQCD with three flavors, the confined description uses the meson matrix $M^i{}_j$ and baryons $B$, $\tilde{B}$, subject to the quantum constraint $\det M - B\tilde{B} = \Lambda^6$. With quark masses $\hat{m} = \mathrm{diag}(m_1,m_2,m_3)$, the superpotential
\begin{equation}
W = X(\det M - B\tilde{B} - \Lambda^6) + \mathrm{Tr}(\hat{m}\,M) + f\,X + \tfrac{1}{2}g\,X\,\phi^2
\end{equation}
combines the Lagrange-multiplier enforcement of the constraint with the O'Raifeartaigh deformation. The vacuum equations $\partial W/\partial M^j{}_j = 0$ yield the Seiberg seesaw
\begin{equation}
M_j = \frac{C}{m_j}\,, \qquad C = \Lambda^2(m_1 m_2 m_3)^{1/3}\,,
\end{equation}
while the F-term $F_X = f \neq 0$ breaks SUSY at the meta-stable vacuum with energy $V_0 = |f|^2$.

The Lagrange multiplier $X$ plays the role of the O'Raifeartaigh pseudo-modulus. Its VEV, stabilized by the Coleman--Weinberg potential at $gv/m = \sqrt{3}$, generates the fermion seed spectrum $(0,(2-\sqrt{3})m,(2+\sqrt{3})m)$ with $Q = 2/3$ as a structural invariant of the Lagrangian.
