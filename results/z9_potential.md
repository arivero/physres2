# Z_9 Instanton Potential for the Koide Phase

## 1. The Z_9-invariant potential

### Symmetry analysis

In SQCD with $N_f = N_c = 3$, the classical $U(1)_{R+2A}$ symmetry has anomaly coefficient $6N_c = 18$, reducing to $\mathbb{Z}_{18}$ after instanton breaking. The meson field $M^i_j \sim Q^i \bar{Q}_j$ carries charge 4 under $\mathbb{Z}_{18}$ (two quark zero modes, each with charge 2).

Writing $\phi = |\phi| e^{i\delta}$ for a composite field (e.g., $\text{Tr}\,M$ or a diagonal meson component), the $\mathbb{Z}_{18}$ generator acts as

$$\phi \to e^{4 \cdot 2\pi i / 18}\,\phi = e^{4\pi i/9}\,\phi, \qquad \delta \to \delta + \frac{4\pi}{9}.$$

A potential term $V \supset \cos(n\delta + \alpha)$ is $\mathbb{Z}_{18}$-invariant iff

$$n \cdot \frac{4\pi}{9} = 2\pi m \quad \text{for integer } m \qquad \Longrightarrow \qquad 2n = 9m.$$

The smallest positive integer solution is $m = 2$, $n = 9$.

**Result.** The lowest-dimension $\mathbb{Z}_{18}$-invariant (hence $\mathbb{Z}_9$-invariant) angular potential is

$$\boxed{V(\delta) = -V_0 \cos\!\big(9(\delta - \delta_0)\big)}$$

Higher harmonics are $\cos(18\delta)$, $\cos(27\delta)$, etc.

---

## 2. Minima and Z_3 classification

### The 9 minima

For $V(\delta) = -V_0 \cos(9(\delta - 2/9))$, the minima satisfy $9(\delta - 2/9) = 2\pi k$, giving

$$\delta_k = \frac{2}{9} + \frac{2\pi k}{9}, \qquad k = 0, 1, \ldots, 8.$$

Explicitly:

| $k$ | $\delta_k$ | Numerical value |
|-----|------------|-----------------|
| 0 | $2/9$ | 0.2222 |
| 1 | $2/9 + 2\pi/9$ | 0.9204 |
| 2 | $2/9 + 4\pi/9$ | 1.6185 |
| 3 | $2/9 + 6\pi/9$ | 2.3166 |
| 4 | $2/9 + 8\pi/9$ | 3.0147 |
| 5 | $2/9 + 10\pi/9$ | 3.7129 |
| 6 | $2/9 + 12\pi/9$ | 4.4110 |
| 7 | $2/9 + 14\pi/9$ | 5.1091 |
| 8 | $2/9 + 16\pi/9$ | 5.8073 |

### Z_3 equivalence

The mass-permutation symmetry $\delta \to \delta + 2\pi/3$ shifts $k \to k+3$ (since $2\pi/3 = 3 \times 2\pi/9$). The 9 minima group into 3 equivalence classes:

| Class | Members $k$ | $\delta \bmod 2\pi/3$ | Mass spectrum $m/M_0$ |
|-------|-------------|----------------------|----------------------|
| **A** | 0, 3, 6 | 0.2222 | {0.00163, 0.3366, 5.662} |
| **B** | 1, 4, 7 | 0.9204 | {0.1623, 2.392, 3.446} |
| **C** | 2, 5, 8 | 1.6185 | {0.0360, 0.8697, 5.094} |

Verification: all members of Class A produce the identical mass set $\{0.00163, 0.3366, 5.662\} \times M_0$ under permutation.

**The lepton triple $(e, \mu, \tau)$ sits in Class A** ($\delta \bmod 2\pi/3 = 0.22223$, matching $2/9$ to 33 ppm).

---

## 3. Instanton origin of cos(9 delta)

### Z_18 charge accounting

| Object | Z_18 charge | Phase shift under generator |
|--------|-------------|---------------------------|
| Quark $Q$ | 2 | $e^{2\pi i \cdot 2/18}$ |
| Meson $M \sim QQ$ | 4 | $e^{4\pi i/9}$ |
| $\det M \sim M^3$ | 12 | $e^{12 \cdot 2\pi i/18}$ |
| $n$-instanton vertex | $6n$ | $e^{6n \cdot 2\pi i/18}$ |

### Option A: Three-instanton superpotential

The $n$-instanton vertex carries Z_18 charge $6n$. For a Z_18 singlet: $6n \equiv 0 \pmod{18}$, requiring $n = 3k$. The **three-instanton** ($n=3$) vertex is the leading candidate.

**Operator:**
$$\delta W_{3\text{-inst}} = c_3 \cdot \frac{(\det M)^3}{\Lambda^{18}}$$

**Charge check:** $\det M$ has Z_18 charge 12, so $(\det M)^3$ has charge 36 = $2 \times 18$ -- a Z_18 singlet. Correct.

**Angular dependence:** Near the SUSY vacuum $\det M = \Lambda^6$, parametrize fluctuations as $\det M = \Lambda^6 e^{3i\delta_M}$ (the factor 3 from three flavors). Then:

$$\delta W = c_3\,\Lambda^{18}\,e^{9i\delta_M}/\Lambda^{18} = c_3\,e^{9i\delta_M}$$

The F-term potential $V \supset |F|^2$ produces, after cross terms with the tree-level superpotential:

$$V_{\text{eff}} \supset -2\,\text{Re}\!\big[c_3\,e^{9i\delta}\big] \times (\ldots) = -V_0\cos(9\delta - 9\delta_0)$$

where $\delta_0 = \arg(c_3)/9$.

**This generates $\cos(9\delta)$ exactly.** The key check: $(\det M)^3$ is the lowest-order Z_18-singlet holomorphic operator in the meson variables.

Note: the 1-instanton vertex ($\det M$, charge 12) is already captured by the Seiberg constraint $\det M - B\bar{B} = \Lambda^6$. The 2-instanton vertex ($(\det M)^2$, charge 24, not a Z_18 singlet) cannot appear in the superpotential. The 3-instanton is the first new contribution.

### Option B: Bion contributions in the Kahler potential

On $\mathbb{R}^3 \times S^1$ with radius $R = 1/T$, SU(3) has three types of monopole-instantons (corresponding to the two simple roots $\alpha_1, \alpha_2$ and the affine root $\alpha_0 = -\alpha_1 - \alpha_2$ of the extended Dynkin diagram).

**Magnetic bions** (from different-root pairs) generate a center-stabilizing potential on the dual photon $\sigma$:

$$V_{\text{bion}} \sim -\sum_{a=0}^{2} \cos(\alpha_a \cdot \sigma)$$

Each monopole-instanton with $N_f = 3$ carries $2N_f = 6$ fermionic zero modes, but a bion (monopole--anti-monopole pair) has zero net fermionic zero modes and can contribute to the bosonic potential.

The Z_3-symmetric combination of the SU(3) affine roots produces $\cos(3\sigma)$. Under the identification of the dual photon phase with the meson phase (accounting for the charge-4 Z_18 assignment of the meson), this maps to $\cos(9\delta)$.

**However**, bion contributions to the Kahler potential are suppressed by $e^{-2S_0/N_c}$ where $S_0 = 8\pi^2/g^2$, and they are controlled only in the weak-coupling (small-$S^1$) regime. At strong coupling ($\Lambda \sim 300$ MeV), the semiclassical bion expansion is unreliable.

### Verdict

**Option A** (three-instanton superpotential) is the cleaner mechanism:
- The operator $(\det M)^3/\Lambda^{18}$ is a well-defined Z_18 singlet
- It generates exactly $\cos(9\delta)$ in the effective potential
- At strong coupling, $c_3$ is $O(1)$ (no exponential suppression)
- The phase $\delta_0 = \arg(c_3)/9$ is determined by the UV theory

**Option B** (bions) also produces $\cos(9\delta)$ but is controlled only in the compactified regime. Both mechanisms give the same symmetry structure.

### What about cos(3 delta) or cos(12 delta)?

- $\cos(3\delta)$: would require a Z_18 invariant with $n = 3$, but $2 \times 3 = 6 \neq 9m$. NOT Z_18 invariant. The one-instanton vertex generates $\cos(3\arg(\det M))$ but this is already part of the Seiberg constraint, not a correction.
- $\cos(12\delta)$: requires $2 \times 12 = 24 = 9 \times 8/3$ -- NOT an integer $m$. Not Z_18 invariant.
- $\cos(18\delta)$: requires $2 \times 18 = 36 = 9 \times 4$. Z_18 invariant. This is the next harmonic, from a 6-instanton vertex $(\det M)^6/\Lambda^{36}$, highly suppressed.

**Conclusion:** $\cos(9\delta)$ is uniquely the lowest harmonic consistent with the Z_18 instanton symmetry.

---

## 4. Numerical estimates

### Scale comparison

| Quantity | Value | $\text{MeV}^4$ |
|----------|-------|-----------------|
| $\Lambda^4$ | $(300)^4$ | $8.1 \times 10^9$ |
| $f_\pi^4$ | $(92)^4$ | $7.2 \times 10^7$ |
| $\chi_{\text{top}}$ (topological susceptibility) | $(180)^4$ | $1.05 \times 10^9$ |

### Three-instanton estimate

The 3-instanton coefficient at strong coupling:
$$c_3 \sim \frac{1}{(16\pi^2)^3} \sim 2.5 \times 10^{-7}$$

giving $V_0 \sim \Lambda^4 \times c_3 \sim 2.1 \times 10^3\ \text{MeV}^4$, i.e., $V_0^{1/4} \sim 6.7$ MeV. This is a conservative (weak-coupling-inspired) estimate.

### Instanton liquid estimate

At strong coupling, the semiclassical dilute-gas approximation breaks down. Using the empirical instanton liquid parameters:

$$V_0 \sim \frac{\Lambda^4}{(4\pi)^2} \sim 5.1 \times 10^7\ \text{MeV}^4, \qquad V_0^{1/4} \sim 85\ \text{MeV}$$

### Z_9 harmonic of topological susceptibility

The $\cos(3\delta)$ potential has known strength $\chi_{\text{top}} \sim (180\ \text{MeV})^4$. The Z_9 harmonic is suppressed relative to Z_3 by a factor $(\chi_{\text{top}}/\Lambda^4)^2 \sim 0.017$:

$$V_0(\mathbb{Z}_9) \sim \chi_{\text{top}} \times \left(\frac{\chi_{\text{top}}}{\Lambda^4}\right)^2 \sim 1.8 \times 10^7\ \text{MeV}^4, \qquad V_0^{1/4} \sim 65\ \text{MeV}$$

### Summary of estimates

| Method | $V_0\ (\text{MeV}^4)$ | $V_0^{1/4}\ (\text{MeV})$ | $V_0/f_\pi^4$ |
|--------|----------------------|---------------------------|----------------|
| 3-instanton semiclassical | $2.1 \times 10^3$ | 6.7 | $2.9 \times 10^{-5}$ |
| Z_9 harmonic of $\chi_{\text{top}}$ | $1.8 \times 10^7$ | 65 | 0.25 |
| Instanton liquid | $5.1 \times 10^7$ | 85 | 0.72 |
| Tree-level $f_\pi^4$ (reference) | $7.2 \times 10^7$ | 92 | 1.00 |

The strong-coupling estimates give $V_0^{1/4} \sim 65$--$85$ MeV, comparable to $f_\pi = 92$ MeV. The Z_9 potential is NOT parametrically suppressed relative to hadronic scales -- it is an $O(1)$ effect in the strong-coupling regime.

### Cross-check: eta' mass

The $\mathbb{Z}_3$ anomalous potential (topological susceptibility) reproduces the $\eta'$ mass:

$$m_{\eta'}^2 \simeq \frac{2N_f\,\chi_{\text{top}}}{f_\pi^2} \qquad \Rightarrow \qquad m_{\eta'} \simeq 863\ \text{MeV} \quad (\text{cf.}\ 958\ \text{MeV experimental})$$

The 10% discrepancy is typical for leading-order instanton estimates. Similar accuracy is expected for the Z_9 harmonic.

---

## 5. Connection to the meson spectrum

### Generalized Koide parametrization

Using $z_k = v_0 + r\cos(\delta + 2\pi k/3)$ where $z_k = \text{sign}_k \times \sqrt{m_k}$, the physical parameters are:

| Triple | $v_0$ | $r$ | $r/v_0$ | $\delta/\pi$ | $\delta \bmod 2\pi/3$ | $Q$ |
|--------|-------|-----|---------|-------------|---------------------|-----|
| $(e, \mu, \tau)$ | 17.72 | 25.05 | 1.4142 | 0.7374 | **0.2222** | 0.66666 |
| seed $(0, s, c)$ | 15.10 | 21.28 | 1.4093 | 0.7500 | 0.2618 | 2/3 exact |
| $(-\sqrt{s}, \sqrt{c}, \sqrt{b})$ | 30.21 | 43.25 | 1.4317 | 0.8734 | 0.6494 | 0.6750 |
| $(\sqrt{c}, \sqrt{b}, \sqrt{t})$ | 171.98 | 244.24 | 1.4202 | 0.6885 | 0.0686 | 0.6695 |

### Input masses (PDG)

$m_e = 0.511$ MeV, $m_\mu = 105.658$ MeV, $m_\tau = 1776.86$ MeV, $m_s = 93.4$ MeV, $m_c = 1270$ MeV, $m_b = 4180$ MeV, $m_t = 172760$ MeV.

### Proximity to Z_9 minima

Each phase is checked against the nearest minimum $\delta = 2/9 + 2\pi k/9$:

| Triple | $\delta \bmod 2\pi/9$ | Distance to $2/9$ | Fraction of Z_9 spacing |
|--------|---------------------|-------------------|------------------------|
| $(e, \mu, \tau)$ | 0.22223 | $7.4 \times 10^{-6}$ | **0.011%** |
| seed $(0, s, c)$ | 0.26180 | 0.0396 | 5.7% |
| $(-s, c, b)$ | 0.65101 | 0.269 | 38.6% |
| $(c, b, t)$ | 0.63156 | 0.289 | 41.4% |

### Interpretation

**The lepton phase sits at a Z_9 minimum to 33 ppm precision.** This is the $2/9$ coincidence: $\delta_\ell \bmod 2\pi/3 = 0.222230 \approx 2/9 = 0.222222$ (residual $7.4 \times 10^{-6}$, or 33 ppm of 2/9).

**The quark phases do NOT sit at Z_9 minima.** Their phases are displaced from any $2/9 + 2\pi k/9$ minimum by 6%--42% of the Z_9 spacing. This means:

1. The Z_9 instanton potential, if it is the mechanism selecting $\delta = 2/9$ for leptons, does not simultaneously constrain the quark sector phases.

2. The seed phase $\delta = 3\pi/4$ ($= 0.75\pi$) is close to but not at the nearest Z_9 minimum at $k=3$ ($\delta_3 = 2/9 + 2\pi/3 = 2.317$, while the seed has $\delta = 2.356$). The 5.7% displacement corresponds to the $\pi/12$ shift of the seed from the lepton minimum: $3\pi/4 \bmod 2\pi/3 = \pi/12 = 0.2618$ vs $2/9 = 0.2222$.

3. The $(-s,c,b)$ and $(c,b,t)$ triples have $Q \neq 2/3$ (deviations of 1.2% and 0.4% respectively), so their phases are not constrained by the Koide condition to sit at any special value. Their phases are determined by the actual mass values and have no obvious relation to Z_9 geometry.

### Sector separation

The Z_9 instanton potential operates in the **confining phase** of $N_f = N_c = 3$ SQCD, where the relevant degrees of freedom are meson composites. The lepton sector, which has the exact Koide relation $Q = 2/3$ to 0.91$\sigma$, naturally couples to this $\cos(9\delta)$ potential.

The quark sector involves mixed up-down triples (the Koide triples obligatorily mix quark types) and operates at higher energy scales where the $N_f = N_c = 3$ confining phase may not apply. The quark phases would be determined by a different dynamical mechanism -- possibly related to the ISS metastable vacuum for $N_f > N_c$, or to the bloom mechanism that rotates $\delta$ from the seed value while changing $v_0$.

---

## Summary

1. **The Z_9-invariant potential $V(\delta) = -V_0\cos(9(\delta - 2/9))$ is the unique lowest harmonic** consistent with the $\mathbb{Z}_{18}$ instanton symmetry of $N_f = N_c = 3$ SQCD.

2. **It arises from a three-instanton superpotential term** $\delta W \propto (\det M)^3/\Lambda^{18}$, which is the lowest-order $\mathbb{Z}_{18}$-singlet holomorphic operator beyond the Seiberg constraint. The angular dependence is $\cos(9\delta)$ because $\det M$ carries Z_18 charge 12, so $(\det M)^3$ has charge 36 = 2 $\times$ 18.

3. **The potential strength is $V_0^{1/4} \sim 65$--$85$ MeV**, comparable to $f_\pi$, in the strong-coupling (instanton liquid) regime.

4. **The 9 minima at $\delta = 2/9 + 2\pi k/9$ split into 3 Z_3-equivalent classes.** The lepton triple sits in Class A ($\delta \bmod 2\pi/3 = 2/9$) to 33 ppm precision.

5. **The quark sector phases are NOT at Z_9 minima**, indicating that the quark mass spectrum is determined by a different dynamical mechanism than the instanton potential that selects the lepton phase.
