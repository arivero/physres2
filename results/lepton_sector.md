# Lepton Sector: SU(2) SQCD with O'Raifeartaigh Deformation

## Model Specification

**Gauge group**: SU(2) with N_f = 3 fundamental flavors (s-confining, N_f = N_c + 1).

**UV fields**: Psi^i_a (i = 1,2,3 flavor; a = 1,2 color).

**Confined degrees of freedom**: Three antisymmetric mesons

    L^{ij} = epsilon^{ab} Psi^i_a Psi^j_b,   i < j

with L^{12} (= L1), L^{13} (= L2), L^{23} (= L3).

**Superpotential** (canonical normalization, all fields dimension 1):

    W = f X_L + m L1 L3 + g X_L L1^2 + eps L1 L2 L3

where X_L is a Lagrange-multiplier superfield. Dimensions: [f] = mass^2, [m] = mass, [g] = dimensionless, [eps] = dimensionless (the Pfaffian coupling in canonical normalization).

The last term is the dynamical s-confining superpotential (Pfaffian). In the original composite variables, W_dyn = L^{12} L^{13} L^{23} / Lambda_L^3, where Lambda_L^{b_0} with b_0 = 3N_c - N_f = 3.

---

## Part (a): F-Term Equations and SUSY Breaking

The F-term conditions dW/dPhi_I = 0:

    (1)  F_{X_L} = f + g L1^2
    (2)  F_{L1}  = m L3 + 2g X_L L1 + eps L2 L3
    (3)  F_{L2}  = eps L1 L3
    (4)  F_{L3}  = m L1 + eps L1 L2

### Metastable vacuum (eps small)

From (3): eps L1 L3 = 0. For eps nonzero, either L1 = 0 or L3 = 0.

**Taking L1 = 0** (the metastable branch):

- (1) gives F_{X_L} = f. **SUSY is broken**: F_{X_L} != 0.
- (2) gives m L3 = 0, so L3 = 0.
- (3) and (4) are automatically satisfied.
- X_L and L2 remain undetermined: **pseudo-moduli**, lifted at one loop by the Coleman-Weinberg potential.

**Vacuum expectation values**:

    <L1> = <L3> = 0,    <X_L> = v (pseudo-modulus),    <L2> = arbitrary

### Distant SUSY-preserving vacuum

When L1 != 0 and L3 = 0, equation (4) gives L2 = -m/eps (parametrically large for eps << 1), and (1) gives L1^2 = -f/g. This is the Pfaffian-restored vacuum, exponentially far in field space. The metastable vacuum is long-lived for eps << 1.

---

## Part (b): Fermion Mass Matrix W_{IJ}

The full 4x4 matrix of second derivatives (basis: X_L, L1, L2, L3):

    W_{IJ} = | 0       2g L1          0                0               |
             | 2g L1   2g X_L         eps L3           m + eps L2      |
             | 0       eps L3         0                eps L1          |
             | 0       m + eps L2     eps L1           0               |

At the metastable vacuum (<L1> = <L3> = 0, <X_L> = v, eps -> 0):

    W_{IJ}|_vac = | 0    0     0    0 |
                  | 0    2gv   0    m |
                  | 0    0     0    0 |
                  | 0    m     0    0 |

---

## Part (c): Fermion Mass Spectrum

X_L decouples (entire row and column vanish) -- this is the **goldstino** direction.

L2 decouples (entire row and column vanish) -- the **pseudo-modulus partner**.

The massive 2x2 block in the (L1, L3) sector:

    M_2x2 = | 2gv   m |
            | m     0 |

Characteristic polynomial: lambda^2 - 2gv lambda - m^2 = 0.

    lambda_pm = gv +/- sqrt(g^2 v^2 + m^2)

Physical masses (absolute values):

    m_+ = gv + sqrt(g^2 v^2 + m^2)
    m_- = sqrt(g^2 v^2 + m^2) - gv

### The Koide seed at gv/m = sqrt(3)

Setting t = gv/m = sqrt(3):

    lambda_+ = m(sqrt(3) + 2) = (2 + sqrt(3))m
    lambda_- = m(sqrt(3) - 2)  =>  |lambda_-| = (2 - sqrt(3))m

The full L-sector spectrum is:

    (m_0, m_-, m_+) = (0, (2-sqrt(3))m, (2+sqrt(3))m)

**Numerical verification** (m = 1):

    M_3x3 eigenvalues:  [-0.26794919,  0.0,  3.73205081]
    Physical masses:     [0, 0.26794919, 3.73205081]
    Analytic:            [0, 0.26794919, 3.73205081]
    All match to machine precision.

### Koide Q = 2/3 exactly

    Q = (m_0 + m_- + m_+) / (sqrt(m_0) + sqrt(m_-) + sqrt(m_+))^2

Numerator: 0 + (2-sqrt(3))m + (2+sqrt(3))m = 4m.

Denominator: Let A = sqrt(2-sqrt(3)), B = sqrt(2+sqrt(3)). Then:
- A*B = sqrt((2-sqrt(3))(2+sqrt(3))) = sqrt(4-3) = 1
- A^2 + B^2 = (2-sqrt(3)) + (2+sqrt(3)) = 4
- (A + B)^2 = A^2 + B^2 + 2AB = 4 + 2 = 6

Therefore (sqrt(m_-) + sqrt(m_+))^2 = m(A+B)^2 = 6m.

    Q = 4m / 6m = 2/3    [exact]

Numerical: Q = 0.666666666666667, |Q - 2/3| = 1.1 x 10^{-16}.

This is the **exact Koide seed** at delta = 3pi/4.

---

## Part (d): Mediation and Physical Lepton Mass Scale

### Koide parametrization

The seed spectrum (0, (2-sqrt(3))m, (2+sqrt(3))m) corresponds to the Koide parametrization

    sigma_k = v_0 (1 + sqrt(2) cos(delta + 2pi k/3))

at delta = 3pi/4, where cos(3pi/4) = -1/sqrt(2) forces sigma_1 = 0 (the massless state).

The total mass at the seed: M_0^{seed} = sum m_k = 4m_OR.

The relationship M_0 = 6 v_0^2 gives v_0^{seed} = sqrt(2m_OR/3).

### Physical lepton masses

    m_e = 0.511 MeV,   m_mu = 105.658 MeV,   m_tau = 1776.86 MeV
    M_0 = m_e + m_mu + m_tau = 1883.03 MeV
    v_0 = (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))/3 = 17.716 MeV^{1/2}
    Q(e,mu,tau) = 0.66666 (deviation -0.0009% from 2/3, within 0.91 sigma)

### Bloom for leptons

    delta_seed    = 135.000 deg
    delta_phys    = 132.733 deg
    Bloom rotation = -2.267 deg   (very small)
    v_0^{phys} / v_0^{seed} = 1.0004   (no v_0-doubling for leptons)

The lepton bloom is a tiny delta-rotation of 2.3 degrees, unlike the quark sector (22 degrees). The v_0 ratio is essentially 1, confirming that **v_0-doubling does not apply to leptons** (consistent with bloom_mechanism.md).

### O'Raifeartaigh mass parameter

Since the bloom is negligible, the seed and physical scales are nearly identical:

    m_OR = M_0^{seed} / 4 = 470.4 MeV

The seed spectrum:

    m_0 = 0
    m_- = (2-sqrt(3)) x 470.4 = 126.0 MeV   (blooms to m_e = 0.511 MeV)
    m_+ = (2+sqrt(3)) x 470.4 = 1755.5 MeV   (blooms to m_tau = 1776.9 MeV)

After the 2.3-degree bloom rotation, the reconstructed masses are:

    m_1 = 0.511 MeV     (physical m_e = 0.511 MeV)
    m_2 = 105.66 MeV    (physical m_mu = 105.66 MeV)
    m_3 = 1776.88 MeV   (physical m_tau = 1776.86 MeV)

The bloom converts the seed massless state into m_e and redistributes mass from m_- to m_mu, while leaving m_+ nearly unchanged as m_tau.

---

## Part (e): Confinement Scale

For SU(2) with N_f = 3, the one-loop beta function coefficient is b_0 = 3N_c - N_f = 3. The dynamical scale Lambda_L^3 = mu^3 exp(-8pi^2/g^2(mu)).

In canonical normalization, the O'Raifeartaigh parameters relate to the confinement scale through the ISS-type mapping:

    m = h Lambda_L     (magnetic Yukawa x confinement scale)
    f = h Lambda_L^2   (SUSY-breaking scale)
    g = h               (magnetic Yukawa coupling)

where h is the (dimensionless) magnetic Yukawa coupling.

The Koide condition gv/m = sqrt(3) translates to <X_L>/Lambda_L = sqrt(3) (independent of h, as shown in or_iss_homework.py).

**Scale determination**: With h ~ O(1):

    Lambda_L ~ m_OR = 470 MeV

This is larger than Lambda_QCD ~ 300 MeV but of the same order. For the lepton confining sector to be separate from QCD, this SU(2) theory must be a distinct gauge group.

If Lambda_L = 300 MeV (as suggested in some earlier estimates), then h = m_OR/Lambda_L = 470/300 = 1.57, which is perturbative and acceptable.

### Dimensional analysis check

The Pfaffian coupling in canonical variables is eps = O(1) (dimensionless). The metastable vacuum lifetime scales as exp(m^2/(eps^2 Lambda_L^2)). For eps << 1, the vacuum is parametrically long-lived.

---

## Part (f): The Goldstino

At the SUSY-breaking vacuum, the F-terms are:

    F_{X_L} = f != 0
    F_{L1} = F_{L2} = F_{L3} = 0

The goldstino direction in field space is:

    psi_goldstino = (sum_I F_I psi_I) / sqrt(sum_I |F_I|^2) = psi_{X_L}

The goldstino is **purely the fermionic component of X_L**.

This is algebraically consistent: the X_L row and column of W_{IJ}|_vac are identically zero (since <L1> = 0), so psi_{X_L} is a zero eigenvector of the fermion mass matrix.

### Gravitino mass (if SUGRA is included)

    m_{3/2} = F / (sqrt(3) M_Pl) = f / (sqrt(3) M_Pl)

With f ~ Lambda_L^2 ~ (470 MeV)^2:

    m_{3/2} ~ 2.2 x 10^5 MeV^2 / (1.73 x 2.4 x 10^{21} MeV)
            ~ 5.3 x 10^{-17} MeV
            ~ 5.3 x 10^{-11} eV

This is an extraordinarily small gravitino mass, reflecting the low SUSY-breaking scale (sqrt(f) ~ 470 MeV << sqrt(f)_gauge-mediation ~ 10^5 GeV).

---

## Summary Table

| Quantity | Value | Notes |
|----------|-------|-------|
| Gauge group | SU(2), N_f = 3 | s-confining |
| Confinement scale Lambda_L | ~ 470 MeV | from m_OR with h ~ 1 |
| O'R mass parameter m_OR | 470.4 MeV | = M_0^{lepton}/4 |
| gv/m (Koide condition) | sqrt(3) | fixes pseudo-modulus VEV |
| Seed spectrum | (0, 126.0, 1755.5) MeV | delta = 3pi/4 |
| Koide Q (seed) | 2/3 exactly | algebraic identity |
| Bloom rotation | -2.27 deg | seed -> physical |
| v_0 ratio (phys/seed) | 1.0004 | no v_0-doubling |
| Physical spectrum | (0.511, 105.66, 1776.88) MeV | matches (m_e, m_mu, m_tau) |
| Q(e,mu,tau) | 0.66666 | 0.0009% from 2/3 |
| Goldstino | psi_{X_L} | purely X_L fermion |
| SUSY-breaking F-term | f ~ (470 MeV)^2 | low-scale breaking |
| Gravitino mass | ~ 5 x 10^{-11} eV | if SUGRA included |

---

## Key Structural Observations

1. **The Koide seed requires exactly gv/m = sqrt(3).** For the triple (0, m_-, m_+) with m_+m_- = m^2, the Koide ratio is Q = 2 sqrt(t^2+1) / (2 sqrt(t^2+1) + 2) where t = gv/m. Setting Q = 2/3 gives sqrt(t^2+1) = 2, hence t = sqrt(3). This is the UNIQUE value of the pseudo-modulus ratio that produces the Koide seed. A scan over t confirms: Q(t=0.5) = 0.528, Q(t=1) = 0.586, Q(sqrt(3)) = 2/3 exactly, Q(t=2) = 0.691. (Note: an earlier computation in or_iss_homework.py incorrectly stated Q = 2/3 for any r; this is wrong.)

2. **The seed mass ratio is universal.** At the Koide seed (delta = 3pi/4), the ratio c_3^2/c_2^2 = (2+sqrt(3))/(2-sqrt(3)) = (2+sqrt(3))^2 = 13.93, which is exactly the O'R ratio m_+/m_-. The Koide parametrization and the O'R mass matrix give identical structure.

3. **The lepton bloom is tiny.** Unlike the quark sector (bloom = 22 deg, v_0-doubling = 2.0005), the lepton bloom is only 2.3 degrees and the v_0 ratio is 1.0004. The lepton spectrum is very close to the seed — the electron mass is a small perturbation on the (0, m_mu-like, m_tau-like) seed.

4. **L^{13} is the pseudo-modulus partner.** The field L^{13} has zero mass at tree level and is not fixed by the F-terms. Its scalar VEV is lifted at one loop by the Coleman-Weinberg potential. The fermionic component of L^{13} remains massless at tree level (second massless state besides the goldstino). At one loop, with the lepton-sector parameters (f ~ (470 MeV)^2, m ~ 470 MeV, g ~ 1), the CW-induced mass is parametrically m_CW ~ g^2 f / (16 pi^2 m) ~ 3 MeV. This is an order of magnitude above m_e = 0.511 MeV. (The earlier estimate of 0.07 MeV in bloom_mechanism.md used quark-sector parameters, not the lepton-sector ones.) Whether the CW mass for L^{13} can be identified with the electron mass after the bloom requires a more careful one-loop computation.

5. **The Pfaffian restores SUSY at a distant vacuum.** The eps-deformation ensures that a SUSY-preserving vacuum exists at L2 ~ m/eps >> Lambda_L. The metastable vacuum lifetime is exponentially long for eps << 1, making the SUSY-breaking vacuum cosmologically stable.
