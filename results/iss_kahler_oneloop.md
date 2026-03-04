# One-Loop Kahler Corrections for the ISS Pseudo-Modulus

## Setup

ISS model: $SU(N_c)$ magnetic theory with $N_c = 3$, $N_f = 4$.

Fields:
- Magnetic quarks $q^a_i$, $\tilde{q}^a_i$ ($a=1\ldots 3$ color, $i=1\ldots 4$ flavor)
- Meson singlets $\Phi^i_j$ ($4\times 4$ matrix)

Superpotential:
$$W = h\, q^a_i \Phi^i_j \tilde{q}^a_j - h\mu^2 \operatorname{Tr}\Phi$$

Metastable vacuum: $\langle q^a_i\rangle = \langle\tilde{q}^a_i\rangle = \mu\,\delta^a_i$ for $a=1\ldots 3$, $i=1\ldots 3$.

SUSY breaking by the rank condition: $F_{\Phi^4_4} = -h\mu^2 \neq 0$.

The pseudo-modulus is $X = \Phi^4_4$, the component not stabilized at tree level.

| Parameter | Value |
|-----------|-------|
| $h$ | 1.0 |
| $\mu$ | 300 MeV |
| $N_c$ | 3 |
| $N_f$ | 4 |

## 1. Tree-Level Mass Spectrum

The pseudo-modulus $X$ couples to $\chi_a = q^a_4$ and $\tilde\chi_a = \tilde{q}^a_4$ ($a=1\ldots N_c$) through the effective superpotential $W_{\text{eff}} = h X\, \chi_a \tilde\chi_a - h\mu^2 X$. This is the standard O'Raifeartaigh structure.

Setting $\langle X\rangle = v$ and defining the dimensionless variable $x = hv/\mu$:

| Mode | Mass$^2$ (units of $(h\mu)^2$) | Multiplicity | Type |
|------|--------------------------------|--------------|------|
| Scalar $+$ | $x^2 + 1$ | $N_c$ complex | boson |
| Scalar $-$ | $x^2 - 1$ | $N_c$ complex | boson |
| Fermion | $x^2$ | $N_c$ Dirac | fermion |
| Goldstino | 0 | 1 Weyl | fermion |
| $X$ (pseudo-modulus) | 0 (tree-level) | 1 complex | boson |

At the vacuum ($x = 0$): $m_+ = h\mu = 300$ MeV, $m^2_- = -(h\mu)^2$ (tachyonic), $m_f = 0$.

The supertrace identities hold exactly at all $x$:
$$\operatorname{STr}[M^2] = N_c[(x^2+1) + (x^2-1) - 2x^2](h\mu)^2 = 0$$
$$\operatorname{STr}[M^4] = 2N_c(h\mu)^4 = 4.86\times 10^{10}\ \text{MeV}^4 \quad\text{(constant)}$$

The constant $\operatorname{STr}[M^4]$ drives the nontrivial (logarithmic) shape of the CW potential.

### Full field count

| Sector | Count | Status |
|--------|-------|--------|
| $q^a_i$, $\tilde{q}^a_i$ ($i \le N_c$) | $2 N_c^2 = 18$ complex | eaten by Higgsing $SU(N_c)$ |
| $\chi_a$, $\tilde\chi_a$ ($i = N_f$) | $2 N_c = 6$ complex | X-dependent masses |
| $\Phi^i_j$ ($i,j \le N_c$) | $N_c^2 = 9$ complex | massive from $q$-VEV |
| $\Phi^i_4$, $\Phi^4_i$ ($i \le N_c$) | $2N_c = 6$ complex | link fields |
| $\Phi^4_4 = X$ | 1 complex | pseudo-modulus |

## 2. Coleman-Weinberg Potential

$$V_{\text{CW}} = \frac{N_c}{64\pi^2}\operatorname{STr}\!\left[M^4\!\left(\ln\frac{M^2}{\Lambda^2} - \frac{3}{2}\right)\right]$$

In dimensionless form, $V_{\text{CW}} = N_c\frac{(h\mu)^4}{64\pi^2}\,f(x)$ with $f(x)$ per color:

$$f(x) = (x^2+1)^2\!\left[\ln(x^2+1) - \tfrac{3}{2}\right] + (x^2-1)^2\!\left[\ln|x^2-1| - \tfrac{3}{2}\right] - 2x^4\!\left[\ln x^2 - \tfrac{3}{2}\right]$$

| $x$ | $f(x)$ | $f'(x)$ |
|-----|--------|---------|
| 0 | $-3.000$ | 0 |
| 0.1 | $-2.999$ | 0.045 |
| 0.5 | $-2.640$ | 2.376 |
| 1.0 | $-0.227$ | 5.545 |
| 1.5 | $+1.588$ | 2.762 |
| 2.0 | $+2.762$ | 2.021 |
| 5.0 | $+6.437$ | 0.800 |

**Monotonicity**: $V_{\text{CW}}$ is strictly increasing for all $x > 0$. The global minimum is at $x = 0$ (confirmed numerically over $x \in [0, 20]$).

## 3. Expansion Near $X = 0$

The potential is non-analytic at $x = 0$ because the fermion mass $m_f = h\mu x$ vanishes there, producing a $x^4\ln x^2$ term. With $u = x^2$:

$$f(u) = -3 + 3u^2 - 2u^2\ln u - \frac{u^4}{6} + O(u^5)$$

Equivalently:

$$f = -3 + x^4\!\left[3 + 2\ln\frac{1}{x^2}\right] - \frac{x^8}{6} + O(x^{10})$$

**Derivation**: Expanding $(1+u)^2\ln(1+u) + (1-u)^2\ln(1-u) = 3u^2 - u^4/6 + \ldots$ (the odd-power terms cancel by the symmetry of the sum), and combining with the constant $-3(1+u^2)$ and the fermion term $-2u^2\ln u + 3u^2$.

Numerical verification shows agreement to machine precision for $x \le 0.1$ and to $10^{-8}$ at $x = 0.3$.

The effective mass-squared of $X$ at finite $v$:

| $v$ (MeV) | $x = hv/\mu$ | $m_X^2$ (MeV$^2$) | $m_X$ (MeV) |
|-----------|--------------|-------------------|-------------|
| 1 | 0.003 | 1.3 | 1.2 |
| 10 | 0.033 | 81 | 9.0 |
| 50 | 0.167 | 1116 | 33 |
| 100 | 0.333 | 2879 | 54 |
| 200 | 0.667 | 4829 | 69 |

The mass grows logarithmically: $m_X^2(v) \sim N_c h^4\mu^2/(8\pi^2)\,[\text{const} + \ln(\mu^2/h^2 v^2)]$.

## 4. Effective Kahler Correction

From the expansion $f = -3 + x^4[3 + 2\ln(1/x^2)] + \ldots$, the quartic coefficient at scale $x_0$ is:

$$c_{\text{quartic}}(x_0) = 3 + 2\ln\frac{1}{x_0^2}$$

The physical effective coupling:

$$c_{\text{eff}}(x_0) = \frac{N_c h^8}{64\pi^2}\left[3 + 2\ln\frac{1}{x_0^2}\right]$$

| $x_0$ | $c_{\text{quartic}}$ | $c_{\text{eff}}$ | Sign |
|--------|---------------------|-------------------|------|
| 0.01 | 21.4 | $1.0\times 10^{-1}$ | $+$ |
| 0.1 | 12.2 | $5.8\times 10^{-2}$ | $+$ |
| 0.5 | 5.8 | $2.7\times 10^{-2}$ | $+$ |
| 1.0 | 3.0 | $1.4\times 10^{-2}$ | $+$ |
| $\sqrt{3}$ | 0.80 | $3.8\times 10^{-3}$ | $+$ |
| $e^{3/4} = 2.12$ | 0 | 0 | 0 |

$c_{\text{quartic}}$ changes sign at $x_0 = e^{3/4} \approx 2.117$, but all physically relevant values ($x < 2$) give **positive** $c_{\text{eff}}$.

**Answer to the key question**: $c_{\text{eff}}$ is NOT negative.

At $x = 1$: $c_{\text{eff}} = +0.014$ vs. the target $-1/12 = -0.083$. The ISS one-loop potential produces a quartic correction with the **opposite sign** and roughly 6 times smaller magnitude.

The one-loop Kahler metric correction:

$$\delta Z(x) = -\frac{N_c h^2}{16\pi^2}\ln\frac{(x^2+1)|x^2-1|}{x^4}$$

This diverges logarithmically as $x \to 0$ (the standard log-enhanced pseudo-modulus mass) and approaches zero for $x \gg 1$ (decoupling).

## 5. External Correction to Shift the Minimum

Since the CW potential is monotonically increasing, the minimum can only be shifted away from the origin by adding an **external** contribution with the opposite sign.

### Pure quartic subtraction

Adding $V_{\text{ext}} = c_{\text{pot}}\,x^4$ to $f(x)$ and requiring $f'(\sqrt{3}) + 4c_{\text{pot}}(\sqrt{3})^3 = 0$ gives:

$$c_{\text{pot}} = -0.113$$

However, the second derivative at $x = \sqrt{3}$ is:

$$f''(\sqrt{3}) + 12 c_{\text{pot}}\cdot 3 = -1.47 + 12(-0.113)(3) = -5.5 < 0$$

This is a **maximum**, not a minimum. The CW potential is concave at $x = \sqrt{3}$ ($f'' < 0$ there), so a simple quartic subtraction creates a saddle point.

### Combined correction

A combined $c_2 x^2 + c_4 x^4$ is needed. The conditions at $x = \sqrt{3}$ are:

$$c_2 + 6c_4 = -f'(\sqrt{3})/(2\sqrt{3}) = -0.680 \qquad \text{(stationarity)}$$
$$2c_2 + 36 c_4 > -f''(\sqrt{3}) = 1.467 \qquad \text{(local minimum)}$$

With $c_4 = -1/12$: $c_2 = -0.180$, giving $d^2V/dx^2 = -4.8 < 0$ (still a maximum).

For a true minimum one needs $c_4 > 0.118$, requiring $c_4$ to be **positive** (the opposite sign from $-1/12$). With $c_2 = -0.680 - 6c_4$, the quadratic term must be sufficiently negative to create the stationary point while the positive quartic provides the upward curvature.

### Kahler metric route

With $K = |X|^2(1 + c|X|^2/\Lambda_K^2)$ and $\Lambda_K = \mu/h$:
$$V_{\text{tree}} = \frac{h^2\mu^4}{1 + 4cx^2}$$

Three solutions exist for $dV_{\text{total}}/dx = 0$ at $x = \sqrt{3}$:

| $c$ | $d^2V/dx^2$ | Nature | Global? |
|-----|-------------|--------|---------|
| $-0.084$ | $0$ | inflection | no |
| $+0.003$ | $-0.038$ | maximum | no |
| $+2.70$ | $+0.035$ | minimum | yes (global) |

The third solution ($c = 2.70$) produces a global minimum at $x = \sqrt{3}$, but requires an unnaturally large Kahler coefficient.

## 6. Summary

| Result | Value |
|--------|-------|
| $f(0)$ (CW vacuum energy per color) | $-3.000$ |
| $f'(x) > 0$ for all $x > 0$ | confirmed |
| $c_{\text{quartic}}$ at $x = 1$ | $+3.0$ |
| $c_{\text{eff}}$ at $x = 1$ | $+0.014$ |
| Sign of $c_{\text{eff}}$ | **positive** |
| $-1/12$ target | $-0.083$ |
| Pure $c_4$ for stationary point at $\sqrt{3}$ | $-0.113$ (gives maximum) |
| Minimum at $\sqrt{3}$ via Kahler | requires $c = +2.70$ |

The ISS one-loop CW potential for the pseudo-modulus:

1. Is monotonically increasing for $x > 0$, selecting $X = 0$ as the global minimum.
2. Produces a **positive** effective quartic Kahler correction $c_{\text{eff}} > 0$, not a negative one.
3. Cannot be overcome by a simple quartic subtraction; the concavity of $f$ at $x = \sqrt{3}$ means any quartic-only correction creates a maximum.
4. Moving the minimum to $x = \sqrt{3}$ through the Kahler metric requires $c \approx +2.7$, a large non-perturbative correction.

The perturbative ISS CW mechanism robustly stabilizes $X = 0$. An external non-perturbative contribution (e.g., the bion Kahler potential) is needed to compete with the one-loop stabilization and shift the vacuum away from the origin.

---

*Generated by iss_kahler_oneloop.py*
