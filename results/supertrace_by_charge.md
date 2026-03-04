# Supertrace Decomposed by EM Charge Sector

## Setup

**16 chiral superfields** with NMSSM coupling $\lambda X H_u \cdot H_d$ ($\lambda = 0.72$).

**Superpotential:**
$$W = \sum_i m_i M^i_i + X(\det M - B\tilde{B} - \Lambda^6) + y_c H_u^0 M^d_d + y_b H_d^0 M^s_s + \lambda X (H_u^0 H_d^0 - H_u^+ H_d^-)$$

**Soft SUSY breaking:** $V_{\text{soft}} = f_\pi^2 \text{Tr}(M^\dagger M)$, applied to the 9 meson fields only.

**EM charges:** $Q(u)=2/3$, $Q(d)=Q(s)=-1/3$. Mesons: $Q(M^i_j) = Q_i - Q_j$. Baryons: $Q(B)=Q(\tilde{B})=0$.

---

## Charge Sector Decomposition

| Sector | Fields | $n_{\text{fields}}$ | $n_{\text{meson}}$ (soft) |
|--------|--------|----------|-------------------|
| $Q = -1$ | Mdu, Msu, Hd- | 3 | 2 |
| $Q = +0$ | Muu, Mdd, Mds, Msd, Mss, X, B, Bt, Hu0, Hd0 | 10 | 5 |
| $Q = +1$ | Mud, Mus, Hu+ | 3 | 2 |

---

## Fermion Spectra by Charge Sector

### $Q = -1$ sector

Fields: Mdu, Msu, Hd-

**Fermion eigenvalues:**

| # | Eigenvalue (MeV) |
|---|------------------|
| 1 | +0.0000000000e+00 |
| 2 | +0.0000000000e+00 |
| 3 | +0.0000000000e+00 |

**Scalar mass-squared eigenvalues:**

| # | $m^2$ (MeV$^2$) | Status |
|---|------------------|--------|
| 1 | +0.0000000000e+00 | OK |
| 2 | +0.0000000000e+00 | OK |
| 3 | +8.4640000000e+03 | OK |
| 4 | +8.4640000000e+03 | OK |
| 5 | +8.4640000000e+03 | OK |
| 6 | +8.4640000000e+03 | OK |

### $Q = +0$ sector

Fields: Muu, Mdd, Mds, Msd, Mss, X, B, Bt, Hu0, Hd0

**Fermion eigenvalues:**

| # | Eigenvalue (MeV) |
|---|------------------|
| 1 | -7.7287405924e+10 |
| 2 | -1.0417580065e-02 |
| 3 | -7.6947261131e-04 |
| 4 | -4.9437110449e-04 |
| 5 | -1.2103086908e-09 |
| 6 | +1.2103086908e-09 |
| 7 | +4.9437110449e-04 |
| 8 | +7.7964377699e-04 |
| 9 | +1.0467155278e-02 |
| 10 | +7.7287405924e+10 |

**Scalar mass-squared eigenvalues:**

| # | $m^2$ (MeV$^2$) | Status |
|---|------------------|--------|
| 1 | -3.1536890747e-07 | OK |
| 2 | +1.4648471270e-18 | OK |
| 3 | +1.4648471270e-18 | OK |
| 4 | +1.4648471270e-18 | OK |
| 5 | +1.4648471270e-18 | OK |
| 6 | +1.4490912609e-06 | OK |
| 7 | +1.0823747261e-04 | OK |
| 8 | +1.1000200264e-04 | OK |
| 9 | +8.4624418225e+03 | OK |
| 10 | +8.4624418225e+03 | OK |
| 11 | +8.4638338065e+03 | OK |
| 12 | +8.4639764515e+03 | OK |
| 13 | +8.4640234947e+03 | OK |
| 14 | +8.4641641627e+03 | OK |
| 15 | +8.4655581780e+03 | OK |
| 16 | +8.4655581780e+03 | OK |
| 17 | +5.9733431145e+21 | OK |
| 18 | +5.9733431145e+21 | OK |
| 19 | +5.9733431145e+21 | OK |
| 20 | +5.9733431145e+21 | OK |

### $Q = +1$ sector

Fields: Mud, Mus, Hu+

**Fermion eigenvalues:**

| # | Eigenvalue (MeV) |
|---|------------------|
| 1 | +0.0000000000e+00 |
| 2 | +0.0000000000e+00 |
| 3 | +0.0000000000e+00 |

**Scalar mass-squared eigenvalues:**

| # | $m^2$ (MeV$^2$) | Status |
|---|------------------|--------|
| 1 | +0.0000000000e+00 | OK |
| 2 | +0.0000000000e+00 | OK |
| 3 | +8.4640000000e+03 | OK |
| 4 | +8.4640000000e+03 | OK |
| 5 | +8.4640000000e+03 | OK |
| 6 | +8.4640000000e+03 | OK |

---

## Quadratic Supertrace per Charge Sector

**Analytic result (exact):**

$$\text{STr}_Q[M^2] = 2 \, \text{Tr}(M^2_{\text{soft}}|_Q) = 2 \, n_{\text{meson}}(Q) \, f_\pi^2$$

This is exact because in $\text{STr} = \text{Tr}(m^2_R + m^2_I) - 2\text{Tr}(m^2_f)$,
the $W^2$ contributions cancel between scalars and fermions, and the holomorphic
contribution $M^2_{\text{hol}}$ cancels between R and I scalar components.
The soft mass is the only surviving term.

| $Q_{\text{em}}$ | $n_{\text{meson}}$ | $\text{STr}_Q / f_\pi^2$ | Value (MeV$^2$) |
|---------|-----------|--------------------------|-----------------|
| $-1$ | 2 | 4 | 33856.0 |
| $+0$ | 5 | 10 | 84640.0 |
| $+1$ | 2 | 4 | 33856.0 |
| **Total** | **9** | **18** | **152352.0** |

---

## EM-Weighted Supertraces

$$\sum_Q Q_{\text{em}} \cdot \text{STr}_Q[M^2] = (+1)(4) + (-1)(4) + (0)(10) = 0$$

$$\sum_Q Q_{\text{em}}^2 \cdot \text{STr}_Q[M^2] = (1)(4) + (1)(4) + (0)(10) = 8 f_\pi^2$$

Numerically: $\sum Q^2 \cdot \text{STr}_Q = 67712.0$ MeV$^2 = 8 f_\pi^2$.

The first sum vanishes by charge conjugation symmetry of the soft sector
(equal number of $Q=+1$ and $Q=-1$ mesons with equal soft masses).

The second sum gives $8 f_\pi^2$ --- a nonzero electromagnetic supertrace.
In the MSSM, $\sum Q^2 \cdot \text{STr}[M^2] = 0$ when there is no EM D-term.
Here the nonzero value reflects the fact that soft breaking is applied
only to mesons (which carry EM charge), not to the Higgs fields or X, B, $\tilde{B}$.

---

## Quartic Supertrace $\text{STr}[M^4]$

$$\text{STr}[M^4] = 4 \text{Tr}(W^2 M_{\text{soft}}) + 2 \text{Tr}(M_{\text{hol}}^2) + 2 \text{Tr}(M_{\text{soft}}^2)$$

| Term | Value (MeV$^4$) |
|------|-----------------|
| $4\text{Tr}(W^2 M_{\text{soft}})$ | +2.022335e+26 |
| $2\text{Tr}(M_{\text{hol}}^2)$ | +3.061626e+19 |
| $2\text{Tr}(M_{\text{soft}}^2)$ | +1.289507e+09 |
| **Total** | **+2.022335e+26** |

Per-sector:

| $Q_{\text{em}}$ | $\text{STr}_Q[M^4]$ (MeV$^4$) | $\text{STr}_Q[M^4] / f_\pi^4$ |
|---------|-------------------------------|--------------------------------|
| $-1$ | +2.865572e+08 | +4.000000 |
| $+0$ | +2.022335e+26 | +2822941407737450496.000000 |
| $+1$ | +2.865572e+08 | +4.000000 |

**Total:** $\text{STr}[M^4] = +2.022335e+26$ MeV$^4$

Unlike $\text{STr}[M^2]$, the quartic supertrace does **not** reduce to a pure soft-mass expression
in general. It receives contributions from $W^2 M_{\text{soft}}$ (fermion masses crossed with soft breaking)
and $M_{\text{hol}}^2$ (Yukawa-induced F-term splittings).

**Charged sectors are special:** For $Q = \pm 1$, all fermion eigenvalues vanish
and all holomorphic splittings vanish (because the off-diagonal meson pairs $M^a_c, M^c_a$
carry opposite charges and are split across different sectors). The only nonzero scalar
masses are exactly $m^2 = f_\pi^2$. Consequently, for the charged sectors:
$$\text{STr}_{Q=\pm 1}[M^{2n}] = 4 f_\pi^{2n} \quad \text{for all } n \geq 1$$

The quartic supertrace is dominated by the $Q=0$ sector, where the $4\text{Tr}(W^2 M_{\text{soft}})$
term is enormous ($\sim 2 \times 10^{26}$ MeV$^4$) due to the ultra-heavy X--meson modes
at $\sim 7.7 \times 10^{10}$ MeV.

---

## Fayet-Iliopoulos D-term

Adding an EM FI term $\xi$, each scalar gets $m^2 \to m^2 + Q_{\text{em}} \xi$.

Per-sector shift:
$$\text{STr}_Q(\xi) = \text{STr}_Q(0) + 2 n_Q \cdot Q \cdot \xi$$

where $n_Q$ is the number of complex fields with charge $Q$.

| $Q$ | $n_Q$ | $2 n_Q Q$ | Effect |
|-----|-------|-----------|--------|
| $-1$ | 3 | -6 | $+(-6)\xi$ |
| $+0$ | 10 | +0 | unchanged |
| $+1$ | 3 | +6 | $+(+6)\xi$ |

**Total:** $\sum_Q n_Q Q = +0$ --- the total supertrace is **independent of $\xi$**,
because the total EM charge of all chiral fields vanishes.

Per-sector, the FI term breaks the degeneracy between $Q=+1$ and $Q=-1$ sectors:

$$\text{STr}_{+1}(\xi) = 4 f_\pi^2 + 6\xi, \qquad \text{STr}_{-1}(\xi) = 4 f_\pi^2 - 6\xi$$

The EM-weighted sums:
$$\sum_Q Q \cdot \text{STr}_Q(\xi) = 12\xi$$
$$\sum_Q Q^2 \cdot \text{STr}_Q(\xi) = 8 f_\pi^2$$

The $Q^2$-weighted sum is independent of $\xi$ because $\sum Q^3 n_Q = 0$.

---

*Generated by supertrace_by_charge.py*