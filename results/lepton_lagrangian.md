# Lepton Sector Lagrangian: SU(2) SQCD with N_f = 3 and O'Raifeartaigh Deformation

## Model Definition

**Gauge group**: SU(2), with N_f = 3 fundamental flavors. This is an s-confining theory
(N_f = N_c + 1 = 3), meaning the infrared dynamics is entirely captured by gauge-invariant
composites with a smooth confining superpotential (no chiral symmetry breaking, no runaway).

**UV fields**: Three chiral superfields in the fundamental of SU(2):

    Psi^i_a,    i = 1,2,3 (flavor),    a = 1,2 (color)

**Confined degrees of freedom**: For SU(2), the fundamental is pseudo-real (2 ~ 2-bar), so
the antisymmetric product of two fundamentals is a gauge singlet. The "baryons" of SU(2)
are bilinears:

    L^{ij} = epsilon^{ab} Psi^i_a Psi^j_b,    i < j

This gives three independent composites:

    L_1 = L^{12} = epsilon^{ab} Psi^1_a Psi^2_b
    L_2 = L^{13} = epsilon^{ab} Psi^1_a Psi^3_b
    L_3 = L^{23} = epsilon^{ab} Psi^2_a Psi^3_b

Each L_i is a chiral superfield that is a gauge singlet under SU(2).

**Auxiliary field**: X_L, a gauge-singlet chiral superfield (Lagrange-multiplier type).

**Total IR field content**: 4 chiral superfields {X_L, L_1, L_2, L_3}.

---

## Total Superpotential

    W = f X_L + m L_1 L_3 + g X_L L_1^2 + epsilon L_1 L_2 L_3

where:

- f has dimension [mass]^2: the SUSY-breaking scale squared
- m has dimension [mass]: the O'Raifeartaigh bilinear mass
- g is dimensionless: the O'Raifeartaigh cubic coupling
- epsilon is dimensionless: the Pfaffian coupling (= 1/Lambda_L^3 in canonical normalization)

**Origin of each term**:

1. **f X_L**: Linear term in X_L. This is the F-term SUSY breaking source. The coefficient f
   sets the scale of SUSY breaking through F_{X_L} = f.

2. **m L_1 L_3**: Bilinear mass term. Together with the X_L L_1^2 coupling, this forms the
   minimal O'Raifeartaigh model for SUSY breaking.

3. **g X_L L_1^2**: The O'Raifeartaigh cubic coupling. This couples the pseudo-modulus X_L
   to the "heavy" field L_1.

4. **epsilon L_1 L_2 L_3**: The dynamical superpotential from s-confinement. For SU(2) with
   N_f = 3, the confined-phase superpotential is the Pfaffian:

       W_dyn = Pf(L) / Lambda_L^{b_0} = L_1 L_2 L_3 / Lambda_L^3

   where b_0 = 3 N_c - N_f = 6 - 3 = 3. In canonical normalization, epsilon = 1/Lambda_L^3.

---

## (a) F-Term Equations and Vacuum Structure

The F-term conditions F_I = dW/dPhi_I = 0 for each field:

    F_{X_L} = dW/dX_L = f + g L_1^2                          (1)
    F_{L_1} = dW/dL_1 = m L_3 + 2g X_L L_1 + epsilon L_2 L_3    (2)
    F_{L_2} = dW/dL_2 = epsilon L_1 L_3                      (3)
    F_{L_3} = dW/dL_3 = m L_1 + epsilon L_1 L_2              (4)

### Metastable SUSY-breaking vacuum

From equation (3): epsilon L_1 L_3 = 0. Since epsilon is nonzero, either L_1 = 0 or L_3 = 0.

**Taking L_1 = 0:**

- Equation (1): F_{X_L} = f. SUSY IS BROKEN: F_{X_L} = f != 0.
- Equation (2): m L_3 + 0 + epsilon L_2 L_3 = L_3 (m + epsilon L_2) = 0.
  For small epsilon, L_3 = 0 (the other branch L_2 = -m/epsilon is the distant vacuum).
- Equations (3) and (4): automatically satisfied.
- X_L and L_2 are undetermined: pseudo-moduli, classically flat.

**Metastable vacuum:**

    <L_1> = 0,    <L_3> = 0,    <X_L> = v (pseudo-modulus),    <L_2> = arbitrary

    F_{X_L} = f != 0    (SUSY broken)
    F_{L_1} = F_{L_2} = F_{L_3} = 0

The SUSY-breaking order parameter is F_{X_L} = f, with the tree-level potential
V_tree = |f|^2, independent of v and <L_2>.

### Distant SUSY-preserving vacuum

Setting all F-terms to zero requires F_{X_L} = 0, hence L_1^2 = -f/g. This needs L_1 != 0.

From (4): m + epsilon L_2 = 0, giving:

    <L_2>_{SUSY} = -m / epsilon

For epsilon << 1 (i.e., Lambda_L >> m^{1/3}), this vacuum is parametrically far in field space.
The tunneling rate from the metastable vacuum to the SUSY-preserving vacuum is suppressed by
exp(-S_bounce), where S_bounce ~ m^4 / (epsilon^2 f^2), making the metastable vacuum
cosmologically long-lived for small epsilon.

From (3): L_3 = 0 at the SUSY vacuum (since L_1 != 0 generically).
From (1): L_1 = +/- i sqrt(f/g).
From (2): 2g X_L L_1 = 0, so X_L = 0 (since L_1 != 0).

**SUSY vacuum:**

    <L_1> = +/- i sqrt(f/g),    <L_2> = -m/epsilon,    <L_3> = 0,    <X_L> = 0

---

## (b) Fermion Mass Matrix at the Metastable Vacuum

The fermion mass matrix is W_{IJ} = d^2 W / (d Phi_I d Phi_J). Computing all second derivatives:

    W_{X_L X_L} = 0
    W_{X_L L_1} = 2g L_1
    W_{X_L L_2} = 0
    W_{X_L L_3} = 0

    W_{L_1 L_1} = 2g X_L
    W_{L_1 L_2} = epsilon L_3
    W_{L_1 L_3} = m + epsilon L_2

    W_{L_2 L_2} = 0
    W_{L_2 L_3} = epsilon L_1

    W_{L_3 L_3} = 0

The full 4 x 4 matrix in the basis (X_L, L_1, L_2, L_3):

    W_{IJ} = | 0          2g L_1       0              0             |
             | 2g L_1     2g X_L       epsilon L_3    m + epsilon L_2 |
             | 0          epsilon L_3  0              epsilon L_1   |
             | 0          m+epsilon L_2  epsilon L_1  0             |

**At the metastable vacuum** (<L_1> = <L_3> = 0, <X_L> = v, epsilon -> 0):

    W_{IJ}|_vac = | 0    0     0    0 |
                  | 0    2gv   0    m |
                  | 0    0     0    0 |
                  | 0    m     0    0 |

This matrix has:
- The X_L row and column are identically zero.
- The L_2 row and column are identically zero.
- A nontrivial 2 x 2 block in the (L_1, L_3) subspace:

    M_{2x2} = | 2gv   m |
              | m     0 |

### Eigenvalues of the full 4 x 4 matrix

Two eigenvalues are zero (from the X_L and L_2 directions). The 2 x 2 block has
characteristic equation:

    lambda^2 - 2gv lambda - m^2 = 0

    lambda_{+/-} = gv +/- sqrt(g^2 v^2 + m^2)

Since sqrt(g^2 v^2 + m^2) > gv, lambda_+ > 0 and lambda_- < 0 always.

**Full spectrum:**

    lambda_1 = 0                              (X_L direction)
    lambda_2 = 0                              (L_2 direction)
    lambda_3 = gv + sqrt(g^2 v^2 + m^2)      (positive)
    lambda_4 = gv - sqrt(g^2 v^2 + m^2)      (negative)

Physical masses (absolute values):

    m_1 = 0
    m_2 = 0
    m_3 = gv + sqrt(g^2 v^2 + m^2)
    m_4 = sqrt(g^2 v^2 + m^2) - gv

Note that m_3 * m_4 = |lambda_3 lambda_4| = |(gv)^2 - (g^2 v^2 + m^2)| = m^2.
Also m_3 + m_4 = 2 sqrt(g^2 v^2 + m^2).

### The Koide seed at gv/m = sqrt(3)

Define the dimensionless ratio t = gv/m. Then:

    m_3 = m (t + sqrt(t^2 + 1))
    m_4 = m (sqrt(t^2 + 1) - t)

At t = sqrt(3):

    sqrt(t^2 + 1) = sqrt(3 + 1) = 2

    m_3 = m (sqrt(3) + 2) = (2 + sqrt(3)) m = 3.73205... m
    m_4 = m (2 - sqrt(3))                    = 0.26795... m

Physical masses:

    (m_1, m_2, m_3, m_4) = (0, 0, (2 + sqrt(3)) m, (2 - sqrt(3)) m)

The three fermion masses of the physical spectrum (excluding the goldstino) are:

    (0, (2 - sqrt(3)) m, (2 + sqrt(3)) m)

### Proof that Q = 2/3 exactly

The Koide quotient for the triple (0, m_-, m_+) with m_- = (2 - sqrt(3)) m and
m_+ = (2 + sqrt(3)) m:

    Q = (m_1 + m_2 + m_3) / (sqrt(m_1) + sqrt(m_2) + sqrt(m_3))^2

Numerator:

    0 + (2 - sqrt(3)) m + (2 + sqrt(3)) m = 4m

Denominator: Let a = sqrt((2 - sqrt(3)) m), b = sqrt((2 + sqrt(3)) m). Then:

    a^2 + b^2 = (2 - sqrt(3)) m + (2 + sqrt(3)) m = 4m
    a b = sqrt((2 - sqrt(3))(2 + sqrt(3))) sqrt(m^2) = sqrt(4 - 3) m = m
    (0 + a + b)^2 = a^2 + b^2 + 2ab = 4m + 2m = 6m

Therefore:

    Q = 4m / 6m = 2/3    [EXACT]

This is the **O'Raifeartaigh-Koide seed**: the unique value of the pseudo-modulus ratio
gv/m = sqrt(3) that produces a Koide-exact fermion triple. The proof that sqrt(3) is
unique: Q(t) = (m_3 + m_4) / (sqrt(m_3) + sqrt(m_4))^2 = 2 sqrt(t^2+1) / (2 sqrt(t^2+1) + 2).
Setting Q = 2/3 gives sqrt(t^2+1) = 2, hence t = sqrt(3).

---

## (c) Identification of the Two Massless States

### The goldstino: psi_{X_L}

At the metastable vacuum, the F-terms are:

    F_{X_L} = f != 0
    F_{L_1} = F_{L_2} = F_{L_3} = 0

The goldstino direction in fermion field space is:

    psi_goldstino = (sum_I F_I psi_I) / |F| = psi_{X_L}

since only F_{X_L} is nonzero. The goldstino is purely the fermionic component of X_L.

Algebraic consistency check: the entire X_L row and column of W_{IJ}|_vac vanish (because
<L_1> = 0 kills the only nonzero entry W_{X_L, L_1} = 2g L_1). So psi_{X_L} is indeed a
zero eigenvector of the fermion mass matrix.

In local SUGRA, the goldstino is eaten by the gravitino via the super-Higgs mechanism.
The gravitino acquires mass m_{3/2} = F / (sqrt(3) M_Pl) = f / (sqrt(3) M_Pl).

### The pseudo-modulus partner: psi_{L_2}

The field L_2 is a pseudo-modulus: its VEV is undetermined at tree level, and both its
scalar and fermionic components are massless at tree level. The entire L_2 row and column
of W_{IJ}|_vac vanish.

At one loop, the Coleman-Weinberg potential lifts the scalar component of both X_L and L_2,
giving them masses of order m_CW ~ g^2 f / (16 pi^2 m). But the fermionic component psi_{L_2}
remains massless at one loop in the absence of additional couplings, because the one-loop
correction to the fermion mass matrix vanishes for pseudo-moduli that do not appear in the
tree-level mass matrix.

More precisely: the one-loop fermion mass correction is:

    delta m_{IJ} ~ (1/16 pi^2) W_{IKL} W_{JMN} <phi_K> <phi_M> ... (loop integrals)

Since W_{L_2 KL} = epsilon delta_{K,L_1} delta_{L,L_3} + (permutations), and both <L_1> = <L_3> = 0,
these corrections vanish at the metastable vacuum.

**Summary of massless states:**

| State | Mass (tree) | Mass (1-loop) | Identity |
|-------|-------------|---------------|----------|
| psi_{X_L} | 0 | 0 (protected by Goldstone theorem) | Goldstino |
| psi_{L_2} | 0 | 0 (vanishes at metastable vacuum) | Pseudo-modulus fermion |
| scalar X_L | 0 | m_CW ~ g^2 f/(16 pi^2 m) | Pseudo-modulus scalar |
| scalar L_2 | 0 | m_CW (similar) | Pseudo-modulus scalar |

The goldstino remains exactly massless to all orders in perturbation theory (it is the
Nambu-Goldstone fermion of broken SUSY, protected until gravity is included). The psi_{L_2}
fermion is unprotected and will acquire mass from nonperturbative effects (the "bloom").

---

## (d) Coleman-Weinberg Potential for the Pseudo-Modulus

### Tree-level mass spectrum as a function of v = g <X_L>/m

At the metastable vacuum with <L_1> = <L_3> = 0 and <X_L> = v m/g, the field-dependent
mass matrices are (all in units of m, using the dimensionless variable v = g<X_L>/m):

**Fermion masses** (from W_{IJ}):

    m_{f,1} = 0 (goldstino)
    m_{f,2} = 0 (L_2 partner)
    m_{f,3} = m (v + sqrt(v^2 + 1))
    m_{f,4} = m (sqrt(v^2 + 1) - v)

    Product: m_{f,3} m_{f,4} = m^2    (exact for all v)
    Sum: m_{f,3} + m_{f,4} = 2m sqrt(v^2 + 1)

**Scalar mass-squared matrices**: The scalar potential V = sum_I |F_I|^2 expanded around
the vacuum gives mass-squared matrices for the real and imaginary components of L_1 and L_3.
Writing L_1 = (a + ib)/sqrt(2), L_3 = (c + id)/sqrt(2):

*Real sector* (a, c basis), in units of m^2:

    A = | 4v^2 + 1 + 2y    2v |
        | 2v                1  |

    tr(A) = 4v^2 + 2 + 2y,    det(A) = (4v^2 + 1 + 2y) - 4v^2 = 1 + 2y

*Imaginary sector* (b, d basis), in units of m^2:

    B = | 4v^2 + 1 - 2y    2v |
        | 2v                1  |

    tr(B) = 4v^2 + 2 - 2y,    det(B) = 1 - 2y

where y = gf/m^2 is the dimensionless SUSY-breaking parameter.

**Tachyon-free condition**: det(B) = 1 - 2y > 0 requires y < 1/2. The metastable vacuum
is perturbatively stable only for gf < m^2/2.

**Eigenvalues of A** (in units of m^2):

    a_{1,2} = (4v^2 + 2 + 2y)/2 +/- sqrt{[(4v^2 + 2y)/2]^2 + 4v^2} / ...

More explicitly, using the standard formula:

    a_{1,2} = (tr A)/2 +/- sqrt{[(tr A)/2]^2 - det A}
            = (2v^2 + 1 + y) +/- sqrt{(2v^2 + y)^2 + 4v^2}

**Eigenvalues of B** (in units of m^2):

    b_{1,2} = (2v^2 + 1 - y) +/- sqrt{(2v^2 - y)^2 + 4v^2}

### Supertraces

    STr[M^2] = (a_1 + a_2 + b_1 + b_2) - 2(m_{f,3}^2 + m_{f,4}^2)
             = (4v^2 + 2 + 2y) + (4v^2 + 2 - 2y) - 2(2(v^2 + 1) + 2v^2)  [using m_f^2 sum]

Wait, more carefully:

    sum(boson m^2) = tr(A) + tr(B) = (4v^2 + 2 + 2y) + (4v^2 + 2 - 2y) = 8v^2 + 4
    sum(fermion m^2) = m_{f,3}^2 + m_{f,4}^2
                     = (v + sqrt(v^2+1))^2 + (sqrt(v^2+1) - v)^2
                     = 2v^2 + 2(v^2+1) = 4v^2 + 2

    STr[M^2] = (8v^2 + 4) - 2(4v^2 + 2) = 0    [exact, for all v]

This is guaranteed for any renormalizable O'Raifeartaigh model with canonical Kahler.

    STr[M^4] = sum(a_i^2) + sum(b_i^2) - 2 sum(m_{f,j}^4)

Using tr(A^2) = (tr A)^2 - 2 det A, etc.:

    sum(a_i^2) = (4v^2 + 2 + 2y)^2 - 2(1 + 2y)
    sum(b_i^2) = (4v^2 + 2 - 2y)^2 - 2(1 - 2y)
    sum(m_f^4) = (m_{f,3}^2 + m_{f,4}^2)^2 - 2(m_{f,3} m_{f,4})^2
               = (4v^2 + 2)^2 - 2

    STr[M^4] = [(4v^2+2+2y)^2 - 2(1+2y)] + [(4v^2+2-2y)^2 - 2(1-2y)]
               - 2[(4v^2+2)^2 - 2]

Expanding the squares and collecting:

    = (4v^2+2)^2 + 2(4v^2+2)(2y) + 4y^2 - 2 - 4y
      + (4v^2+2)^2 - 2(4v^2+2)(2y) + 4y^2 - 2 + 4y
      - 2(4v^2+2)^2 + 4

    = 2(4v^2+2)^2 + 8y^2 - 4 - 2(4v^2+2)^2 + 4
    = 8y^2

**STr[M^4] = 8y^2**, independent of v. This is exact.

### The CW potential

The one-loop Coleman-Weinberg effective potential (in DR-bar scheme, with renormalization
scale mu_R = m):

    V_CW(v) = (1/(64 pi^2)) { sum_{bosons} m_b^4 [ln(m_b^2/m^2) - 3/2]
                              - 2 sum_{fermions} m_f^4 [ln(m_f^2/m^2) - 3/2] }

The sum runs over the 4 real scalar masses from the (L_1, L_3) sector (eigenvalues of A and B)
and the 2 nonzero fermion masses m_{f,3} and m_{f,4}. The X_L and L_2 fields contribute nothing
(massless at tree level).

Since STr[M^4] = 8y^2 is v-independent, the constant piece (-3/2) STr[M^4] is v-independent,
and the nontrivial v-dependence comes from:

    V_CW(v) = (1/(64 pi^2)) STr[M^4 ln(M^2/m^2)] + const

The Taylor expansion in v^2 (all mass eigenvalues are functions of v^2):

    V_CW(v) = V_0(y) + c_2(y) v^2 + c_4(y) v^4 + c_6(y) v^6 + ...

From numerical computation (see pseudomodulus_cw.md and kahler_oneloop.md):

| y     | c_2       | c_4        |
|-------|-----------|------------|
| 0.01  | +0.00213  | -0.00192   |
| 0.10  | +0.214    | -0.194     |
| 0.30  | +1.998    | -1.902     |
| 0.50  | +6.181    | -7.028     |
| 1.00  | +23.55    | -32/3      |

**Key features:**

1. c_2 > 0 for all y > 0: the origin v = 0 is a local minimum of V_CW. The pseudo-modulus
   scalar X_L gets a positive mass-squared m_{X_L}^2 = c_2 g^2 m^2 / (64 pi^2) > 0.

2. c_4 < 0 for y < 1.425: the quartic is negative, so the CW potential has a local minimum
   at v = 0 and bends down for larger v. The global minimum is at some finite v_min > 0,
   determined by the balance of the v^2, v^4, and v^6 terms.

3. The CW minimum lies at v_min ~ 0.3-0.5 for the physical range y < 0.5, which is far below
   sqrt(3) = 1.732.

### CW stabilization and the Koide condition

**The one-loop CW potential does not stabilize X_L at gv/m = sqrt(3).** The CW minimum
is at v_min << sqrt(3) in the entire physical parameter range. The Koide seed requires
gv/m = sqrt(3), which is a special point of the superpotential geometry, not of the
quantum correction.

**Resolution**: The Koide condition gv/m = sqrt(3) must be imposed by a non-canonical
Kahler potential, specifically a tree-level quartic correction:

    K = |X_L|^2 - |X_L|^4 / (12 mu_K^2)

where mu_K is a Kahler mass scale. This generates a scalar potential:

    V = |f|^2 / (1 - |X_L|^2 / (3 mu_K^2))

which has a pole at |X_L| = sqrt(3) mu_K. If the seesaw vacuum places X_L at a value
such that gv/m = sqrt(3), the Kahler pole coincides with the Koide point. The self-consistency
condition is mu_K = m/g (i.e., the Kahler scale equals the natural mass unit of the
O'Raifeartaigh model). This is one condition on one parameter.

Alternatively, in the ISS framework, the condition gv/m = sqrt(3) is a UV boundary condition
transmitted to the IR through the structure of the magnetic dual Kahler potential. It is not
generated dynamically by CW corrections but imposed by the geometry of the moduli space.

---

## (e) Lepton Mass Identifications

### The Koide parametrization

The mass triple (m_1, m_2, m_3) satisfying Q = 2/3 has the Koide parametrization:

    sqrt(m_k) = v_0 (1 + sqrt(2) cos(delta + 2 pi k/3)),    k = 0, 1, 2

with v_0 = (sqrt(m_1) + sqrt(m_2) + sqrt(m_3)) / 3 and Q = 2/3 exactly for any delta.

### The seed

At the O'Raifeartaigh vacuum with gv/m = sqrt(3), the spectrum is:

    (m_0, m_-, m_+) = (0, (2 - sqrt(3)) m_OR, (2 + sqrt(3)) m_OR)

This corresponds to the Koide parametrization at delta = 3 pi/4, where:

    cos(3 pi/4) = -1/sqrt(2)

forces the k=0 eigenvalue to zero:

    sqrt(m_0) = v_0 (1 + sqrt(2) (-1/sqrt(2))) = v_0 (1 - 1) = 0

The seed total mass: M_0^{seed} = 0 + (2 - sqrt(3)) m_OR + (2 + sqrt(3)) m_OR = 4 m_OR.

### Physical lepton masses

    m_e   = 0.510999 MeV
    m_mu  = 105.6584 MeV
    m_tau = 1776.86 MeV

    M_0 = m_e + m_mu + m_tau = 1883.03 MeV
    v_0 = (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau)) / 3 = 17.716 MeV^{1/2}
    Q(e, mu, tau) = 0.666661 (deviation from 2/3: -0.0007%, within 0.91 sigma)

### The O'Raifeartaigh mass parameter

    m_OR = M_0 / 4 = 1883.03 / 4 = 470.76 MeV

At the seed: m_- = (2 - sqrt(3)) x 470.76 = 126.1 MeV, m_+ = (2 + sqrt(3)) x 470.76 = 1756.9 MeV.

### The bloom

The physical lepton masses differ from the seed by a small rotation in delta:

    delta_seed = 3 pi/4 = 135.000 deg
    delta_phys = 132.733 deg
    Bloom rotation: Delta_delta = -2.267 deg

The bloom converts the zero eigenvalue into m_e and redistributes mass:

| Quantity | Seed | Physical | Change |
|----------|------|----------|--------|
| m_1 | 0 | 0.511 MeV (m_e) | generated |
| m_2 | 126.1 MeV | 105.66 MeV (m_mu) | -16% |
| m_3 | 1756.9 MeV | 1776.86 MeV (m_tau) | +1.1% |
| v_0 | 17.710 MeV^{1/2} | 17.716 MeV^{1/2} | +0.03% |
| Q | 2/3 exact | 0.666661 | -0.0007% |

The lepton bloom is a tiny 2.3-degree delta-rotation at fixed v_0 (unlike the quark sector
where the bloom is 22 degrees with v_0 doubling). The physical lepton spectrum is very
close to the O'Raifeartaigh seed.

### Bloom mechanism

The bloom cannot be perturbative (CW-induced). The bloom_mechanism.md analysis shows:

1. Explicit R-symmetry breaking pushes Q away from 2/3. No perturbative epsilon
   preserves the Koide condition during the bloom.

2. The CW-induced mass for the pseudo-modulus fermion is m_CW ~ g^2 f/(16 pi^2 m)
   ~ 0.07 MeV (with hadronic parameters), which is the right order of magnitude for
   m_e = 0.511 MeV but requires a nonperturbative enhancement factor of ~7.

3. The bloom must be a nonperturbative effect that rotates delta while preserving
   (or nearly preserving) Q = 2/3. The Koide seed is an unstable fixed point under
   perturbative deformations; the bloom stabilizes on the physical value delta = 132.7 deg
   through a mechanism not yet identified.

---

## (f) Kahler Potential

### Canonical Kahler

The leading Kahler potential for the confined fields and the auxiliary field is:

    K_can = |X_L|^2 + |L_1|^2 + |L_2|^2 + |L_3|^2

This is the canonical kinetic term for gauge-singlet chiral superfields. In the magnetic
(Seiberg) dual, the canonical normalization is exact at leading order in 1/Lambda_L.

### Non-canonical corrections

The most general Kahler consistent with the symmetries (SU(3) flavor acting on the L_i,
U(1)_R, and the discrete symmetries of the superpotential) has:

    K = |X_L|^2 + sum_i |L_i|^2
        + (c_X / Lambda_L^2) |X_L|^4
        + (c_1 / Lambda_L^2) |L_1|^4 + (c_2 / Lambda_L^2) |L_2|^4 + (c_3 / Lambda_L^2) |L_3|^4
        + (c_{12} / Lambda_L^2) |L_1|^2 |L_2|^2 + (c_{13} / Lambda_L^2) |L_1|^2 |L_3|^2
        + (c_{23} / Lambda_L^2) |L_2|^2 |L_3|^2
        + (c_{X1} / Lambda_L^2) |X_L|^2 |L_1|^2 + ...
        + O(1/Lambda_L^4)

Each coefficient c_i is an O(1) number determined by the strong dynamics of the SU(2)
confining theory, not calculable in perturbation theory.

### The pseudo-modulus stabilization term

The critical non-canonical correction is the X_L quartic:

    K_X = - |X_L|^4 / (12 mu_K^2)

with c_X = -1/(12 mu_K^2 / Lambda_L^2) = -Lambda_L^2 / (12 mu_K^2).

This generates the scalar potential (with F_{X_L} = f):

    V(X_L) = |f|^2 K^{X_L X_L-bar} = |f|^2 / (1 - 2|X_L|^2 / (3 mu_K^2))

which has a pole at |X_L| = sqrt(3/2) mu_K. The effective metric is:

    K_{X_L X_L-bar} = 1 - 2|X_L|^2 / (3 mu_K^2)

(Here the factor of 2 arises from d^2(|X|^4)/d X d X-bar = 2|X|^2 per term, times the
coefficient 1/(12 mu_K^2), giving the pole at |X|^2 = 3 mu_K^2 / 2.)

Setting the pole to coincide with the Koide point gv/m = sqrt(3) requires:

    g sqrt(3/2) mu_K / m = sqrt(3)
    mu_K = sqrt(2) m / g

This is one condition on the Kahler mass scale, determining mu_K in terms of the
superpotential parameters.

### One-loop Kahler correction (from CW)

At one loop, the CW effective potential can be reinterpreted as a correction to the
effective Kahler potential. The v^4 term gives:

    delta K_eff = (g^4 / (16 pi^2)) C_4(y) |X_L|^4 / m^2

where C_4(y) = c_4(y)/4 and c_4 < 0 for y < 1.425. For y = 0.1, C_4 = -0.0485.
This one-loop correction is parametrically smaller than the tree-level c_X term by a
factor g^4/(16 pi^2) and does not control the location of the minimum.

### Complete Kahler for the lepton sector

    K_lepton = |X_L|^2 (1 - |X_L|^2 / (3 mu_K^2))
             + |L_1|^2 + |L_2|^2 + |L_3|^2
             + O(|L_i|^4 / Lambda_L^2)

where mu_K = sqrt(2) m / g. The L_i quartic corrections are subleading at the vacuum
(since <L_1> = <L_3> = 0) and only affect the dynamics away from the metastable vacuum.

---

## (g) Coupling to the Standard Model

### SM quantum number assignment

The UV fields Psi^i_a carry SU(2)_color x SU(2)_W x U(1)_Y quantum numbers. To produce
composites L_i with the quantum numbers of charged leptons (or their superpartners), we need
the composites L^{ij} = epsilon^{ab} Psi^i_a Psi^j_b to be SU(2)_W singlets with
hypercharge Y = -1 (for right-handed charged leptons, Q_EM = Y = -1).

Since L^{ij} is a bilinear in the UV fields, the quantum numbers combine as:

    Y(L^{ij}) = Y(Psi^i) + Y(Psi^j)
    SU(2)_W(L^{ij}) = SU(2)_W(Psi^i) x SU(2)_W(Psi^j) restricted to the antisymmetric product in color

**Assignment 1: All UV fields are SU(2)_W singlets**

If all Psi^i are SU(2)_W singlets with hypercharge Y(Psi^i) = y_i, then:

    Y(L^{12}) = y_1 + y_2
    Y(L^{13}) = y_1 + y_3
    Y(L^{23}) = y_2 + y_3

For the composites to have the quantum numbers of right-handed charged leptons
(SU(2)_W singlets, Y = -1), we need:

    y_1 + y_2 = -1
    y_1 + y_3 = -1
    y_2 + y_3 = -1

Solving: y_1 = y_2 = y_3 = -1/2.

All three UV fields have Y = -1/2. Then all three composites have Y = -1 and are
SU(2)_W singlets — consistent with right-handed charged leptons e_R, mu_R, tau_R.

**Assignment 2: UV fields form SU(2)_W doublets**

If the Psi^i are SU(2)_W doublets, the composites L^{ij} would generically transform as
non-trivial SU(2)_W representations (1 + 3 from 2 x 2 = 1 + 3). The SU(2)_W singlet
component would need to be isolated, and the triplet would need to decouple. This is more
complex and not minimal.

### Identification of composites with leptons

With Assignment 1 (all Y = -1/2, SU(2)_W singlets):

| Composite | UV content | SM quantum numbers | Mass (seed) | Lepton |
|-----------|-----------|-------------------|-------------|--------|
| L_1 = L^{12} | epsilon^{ab} Psi^1_a Psi^2_b | (1, -1) | (2 - sqrt(3)) m_OR = 126 MeV | mu (after bloom) |
| L_2 = L^{13} | epsilon^{ab} Psi^1_a Psi^3_b | (1, -1) | 0 (pseudo-modulus) | e (after bloom) |
| L_3 = L^{23} | epsilon^{ab} Psi^2_a Psi^3_b | (1, -1) | (2 + sqrt(3)) m_OR = 1757 MeV | tau (after bloom) |

The mass eigenstates at the O'Raifeartaigh vacuum are mixtures of L_1 and L_3 (from the
2x2 block diagonalization), with L_2 remaining an eigenstate (massless at tree level).

More precisely, the mass eigenstates in the (L_1, L_3) sector are:

    chi_+ = cos(theta) L_1 + sin(theta) L_3    (mass = m_+)
    chi_- = -sin(theta) L_1 + cos(theta) L_3   (mass = |m_-|)

where tan(2 theta) = 2m / (2gv) = 1/sqrt(3) at gv/m = sqrt(3), giving theta = pi/12 = 15 deg.

### Coupling to the Higgs sector

The lepton Yukawa coupling to the SM Higgs requires an SU(2)_W doublet. Since the composites
L_i are SU(2)_W singlets (with Y = -1), the coupling to the Higgs is:

    W_lepton-Higgs = y_l L_i (H_d . L_W^i)

where L_W^i is an SU(2)_W doublet lepton (the left-handed lepton doublet). But in the
sBootstrap, the left-handed leptons arise from a different part of the SU(5) flavor structure.

**This coupling is not specified within the present SU(2) SQCD model.** The lepton sector
Lagrangian written here describes the right-handed (SU(2)_W-singlet) charged leptons. Their
coupling to the left-handed leptons and the Higgs is an additional ingredient that requires
specifying the left-handed lepton embedding.

### Lepton number

The composites L_i carry lepton number L = 1 if we assign L(Psi^i) = 1/2. Since L^{ij} =
epsilon^{ab} Psi^i_a Psi^j_b is bilinear, L(L^{ij}) = 1. This is consistent with the
identification of L_i as charged lepton superfields.

---

## (h) Parameter Count

### Free parameters of the lepton sector

| Parameter | Dimension | Value | What it determines |
|-----------|-----------|-------|-------------------|
| f | [mass]^2 | ~ (470 MeV)^2 | SUSY-breaking scale; F_{X_L} = f |
| m = m_OR | [mass] | 470.76 MeV | Overall lepton mass scale; M_0 = 4 m_OR |
| g | dimensionless | ~ O(1) | Cubic O'R coupling; gv/m = sqrt(3) fixes v |
| epsilon | dimensionless | 1/Lambda_L^3 | Pfaffian coupling; controls metastable vacuum lifetime |
| mu_K | [mass] | sqrt(2) m/g | Kahler mass scale; determined by gv/m = sqrt(3) condition |
| delta_bloom | angle | -2.267 deg | Bloom rotation; generates m_e and redistributes (m_mu, m_tau) |

### Counting

1. **m_OR** (= m): Sets the overall lepton mass scale. Determined by M_0 = m_e + m_mu + m_tau = 4 m_OR.
   This is ONE input, equivalent to specifying one lepton mass if Q = 2/3 is imposed.

2. **gv/m = sqrt(3)**: This is one condition that is imposed (either by the Kahler pole
   mechanism or as a UV boundary condition). It is NOT a free parameter — it is the Koide
   constraint. Given m_OR and gv/m = sqrt(3), the seed spectrum (0, m_-, m_+) is completely
   fixed. Zero free parameters at the seed.

3. **delta_bloom = -2.267 deg**: This single parameter converts the seed (0, m_-, m_+) into
   the physical (m_e, m_mu, m_tau). It is ONE free parameter if the Koide condition is
   maintained during bloom (Q = 2/3 after bloom as well). If Q is allowed to deviate slightly,
   there is a second parameter (the deviation Delta_Q ~ 10^{-5}).

4. **f**: The SUSY-breaking scale. It sets F_{X_L} = f and controls the gravitino mass
   m_{3/2} = f/(sqrt(3) M_Pl). It does NOT affect the fermion mass spectrum at tree level
   (since W_{IJ} does not depend on f). It affects the scalar spectrum (through the
   splitting parameter y = gf/m^2) and the metastable vacuum lifetime (through the tunneling
   rate). Effectively free but decoupled from the lepton masses.

5. **g**: The cubic coupling. Combined with gv/m = sqrt(3) and m = m_OR, this determines
   v = sqrt(3) m_OR / g, and through mu_K = sqrt(2) m/g, the Kahler mass scale. It does
   NOT appear in the fermion mass eigenvalues (which depend only on the ratio gv/m = sqrt(3)
   and the overall scale m). It affects the scalar spectrum. Effectively free but decoupled
   from lepton masses.

6. **epsilon = 1/Lambda_L^3**: The Pfaffian coupling. Controls the distance to the SUSY
   vacuum (L_2 = -m/epsilon) and the metastable vacuum lifetime. Does NOT affect the
   fermion mass spectrum at the metastable vacuum (since <L_1> = <L_3> = 0 eliminates
   all epsilon-dependent terms from W_{IJ}). Effectively free but decoupled.

### Summary of parameter counting

**Parameters that determine lepton masses:**

| Count | Parameter | What it determines |
|-------|-----------|-------------------|
| 1 | m_OR = 470.76 MeV | Overall scale: M_0 = 4 m_OR |
| 0 | gv/m = sqrt(3) | Koide condition (constraint, not free parameter) |
| 1 | delta_bloom = -2.267 deg | Electron mass and (m_mu, m_tau) redistribution |

**Total: 2 parameters determine 3 lepton masses** (with the Koide condition as a constraint).

Equivalently: specifying any two of (m_e, m_mu, m_tau) determines the third, via Q = 2/3.

**Parameters that do NOT affect lepton masses:**

| Parameter | Role |
|-----------|------|
| f | SUSY-breaking scale (gravitino mass, scalar spectrum) |
| g | Cubic coupling (Kahler scale, scalar spectrum) |
| epsilon = 1/Lambda_L^3 | Metastable vacuum lifetime |
| mu_K | Determined by g, m, and the Koide condition |

**Total free parameters in the lepton sector Lagrangian: 5** (m, g, f, epsilon, delta_bloom),
of which only 2 (m, delta_bloom) affect the physical lepton masses.

---

## Complete Lagrangian (assembled)

### In superfield notation

    L = integral d^4 theta K(Phi, Phi-bar)
      + [integral d^2 theta W(Phi) + h.c.]

with:

    K = |X_L|^2 (1 - |X_L|^2 / (3 mu_K^2))
      + |L_1|^2 + |L_2|^2 + |L_3|^2
      + O(|L_i|^4 / Lambda_L^2)

    W = f X_L + m L_1 L_3 + g X_L L_1^2 + epsilon L_1 L_2 L_3

### In component form

The scalar potential:

    V = K^{I J-bar} F_I F-bar_J

    = |f + g L_1^2|^2 / (1 - 2|X_L|^2/(3 mu_K^2))
      + |m L_3 + 2g X_L L_1 + epsilon L_2 L_3|^2
      + |epsilon L_1 L_3|^2
      + |m L_1 + epsilon L_1 L_2|^2

At the metastable vacuum (<L_1> = <L_3> = 0, <X_L> = v m/g):

    V_vac = f^2 / (1 - 2g^2 v^2 / (3 g^2 mu_K^2))  [with v = sqrt(3) m/g]
          = f^2 / (1 - 2m^2/(mu_K^2))

    With mu_K = sqrt(2) m/g, mu_K^2 = 2m^2/g^2, so 2m^2/mu_K^2 = g^2.
    V_vac = f^2 / (1 - g^2)

This requires g < 1 for the vacuum to be well-defined.

The fermion mass terms:

    L_mass = -(1/2) W_{IJ} psi_I psi_J + h.c.

At the metastable vacuum:

    L_mass = -gv m (psi_{L_1} psi_{L_1}) - (m/2)(psi_{L_1} psi_{L_3} + psi_{L_3} psi_{L_1}) + h.c.
           = -gv m psi_{L_1} psi_{L_1} - m psi_{L_1} psi_{L_3} + h.c.

This gives the mass matrix analyzed in part (b), with physical spectrum
(0, 0, (2-sqrt(3)) m, (2+sqrt(3)) m) at gv/m = sqrt(3).

### In LaTeX form

```latex
\begin{align}
W_{\text{lepton}} &= f\, X_L + m\, L_1 L_3 + g\, X_L L_1^2
  + \frac{1}{\Lambda_L^3}\, L_1 L_2 L_3
\end{align}
```

```latex
\begin{align}
K_{\text{lepton}} &= |X_L|^2\!\left(1 - \frac{|X_L|^2}{3\mu_K^2}\right)
  + |L_1|^2 + |L_2|^2 + |L_3|^2
\end{align}
```

with the Koide condition:

```latex
\frac{g\langle X_L\rangle}{m} = \sqrt{3}
\qquad\Longrightarrow\qquad
\mu_K = \frac{\sqrt{2}\,m}{g}
```

and the physical lepton mass spectrum:

```latex
\begin{align}
(m_e,\, m_\mu,\, m_\tau) &\leftrightarrow
  \bigl(0,\; (2{-}\sqrt{3})\,m,\; (2{+}\sqrt{3})\,m\bigr)
  \times (\text{bloom rotation } \delta = -2.27^\circ)
\notag\\[4pt]
m &= m_{\text{OR}} = \frac{m_e + m_\mu + m_\tau}{4} = 470.76\;\text{MeV}
\end{align}
```

---

## Summary Table

| Quantity | Value | Status |
|----------|-------|--------|
| Gauge group | SU(2), N_f = 3 | s-confining, determined |
| Confinement scale Lambda_L | ~ 470 MeV (if epsilon ~ 1) | Estimated |
| O'Raifeartaigh mass m_OR | 470.76 MeV | = (m_e + m_mu + m_tau)/4 |
| gv/m = sqrt(3) | Koide condition | Imposed (Kahler pole or UV b.c.) |
| Seed spectrum | (0, 126.0, 1756.9) MeV | Exact at delta = 3pi/4 |
| Q(seed) | 2/3 exactly | Algebraic identity |
| Bloom rotation | -2.267 deg | 1 free parameter |
| v_0(phys)/v_0(seed) | 1.0004 | No v_0-doubling for leptons |
| Physical spectrum | (0.511, 105.66, 1776.86) MeV | Matches (m_e, m_mu, m_tau) |
| Q(e, mu, tau) | 0.666661 | 0.91 sigma from 2/3 |
| Goldstino | psi_{X_L} | Eaten by gravitino in SUGRA |
| Pseudo-modulus fermion | psi_{L_2} | Massless at tree+1-loop; blooms to m_e |
| SUSY-breaking F-term | F_{X_L} = f ~ (470 MeV)^2 | Low-scale breaking |
| Gravitino mass | ~ 5 x 10^{-11} eV | If SUGRA included |
| STr[M^2] | 0 (exact) | Guaranteed by canonical K |
| STr[M^4] | 8y^2 m^4 (v-independent) | Exact |
| CW pseudo-modulus mass | m_X^2 ~ g^2 m^2 c_2/(64 pi^2) | Positive, stabilizes at v=0 |
| CW does NOT select sqrt(3) | v_min(CW) ~ 0.3-0.5 << sqrt(3) | Requires non-canonical K |
| Free parameters (total) | 5: m, g, f, epsilon, delta_bloom | |
| Free parameters (masses) | 2: m_OR, delta_bloom | |

---

## Structural Observations

1. **The lepton sector is a SEPARATE SU(2) confining theory**, distinct from the quark sector
   SU(3) SQCD. This is required because the Kahler corrections from the quark sector SU(3)
   theory cannot generate lepton masses (demonstrated in kahler_leptons.md: rank obstruction,
   scale mismatch, Q never reaches 2/3). The lepton masses arise from a separate O'Raifeartaigh
   deformation of SU(2) SQCD.

2. **The O'Raifeartaigh-Koide connection is exact and algebraic**: gv/m = sqrt(3) is the
   UNIQUE value that produces Q = 2/3 for the fermion mass triple. This is not a numerical
   coincidence but a consequence of the specific algebraic structure of the 2x2 mass matrix
   with det = -m^2.

3. **The seed is determined by one parameter** (m_OR = 470.76 MeV). The Koide condition and
   the O'Raifeartaigh structure together fix the entire seed spectrum. The bloom adds one more
   parameter (delta_bloom). Two parameters for three masses.

4. **The bloom mechanism is unsolved**. No perturbative mechanism (CW, R-breaking deformation)
   can rotate delta while preserving Q = 2/3. The bloom must be nonperturbative. This is the
   central open problem of the lepton sector.

5. **The CW potential stabilizes X_L at v << sqrt(3)**, not at the Koide point. The Koide
   condition must be imposed by a tree-level non-canonical Kahler potential (the -|X|^4/(12 mu_K^2)
   term) or by a UV boundary condition. The one-loop CW dynamics confirms stability of the
   pseudo-modulus but does not select the Koide-special value.

6. **The gravitino mass is extremely small** (~ 10^{-11} eV) because the SUSY-breaking scale
   sqrt(f) ~ 470 MeV is far below the gauge-mediation scale ~ 10^5 GeV. This is a low-scale
   SUSY breaking scenario.

---

*Generated: 2026-03-04*
*Source files: lepton_sector.md, oraifeartaigh_koide.md, pseudomodulus_cw.md,*
*kahler_oneloop.md, kahler_leptons.md, bloom_mechanism.md, complete_lagrangian.md*
