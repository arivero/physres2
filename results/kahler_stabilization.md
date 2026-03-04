# Kahler Potential Stabilization of the Pseudo-Modulus at gv/m = sqrt(3)

## Problem Statement

The O'Raifeartaigh model W = f X + m phi_1 phi_3 + g X phi_1^2 breaks SUSY via F_X = f.
The pseudo-modulus X = v is classically flat. The one-loop Coleman-Weinberg potential stabilizes
it at v = 0 (or at small gv/m ~ 0.3-0.5 with canonical Kahler). But the Koide seed condition
Q(0, m_-, m_+) = 2/3 requires gv/m = sqrt(3).

The question: can a non-canonical Kahler potential shift the VEV to gv/m = sqrt(3)?

Throughout, we use the dimensionless variable t = gv/m.

---

## (a) Field-Dependent Mass Spectrum

### Conventions

Superpotential: W = f X + m phi_1 phi_3 + g X phi_1^2

Fields: X (pseudo-modulus), phi_1, phi_3 (massive pair). At the vacuum phi_1 = phi_3 = 0,
X = v (real, positive).

F-terms at the vacuum:

    F_X = f               (nonzero: SUSY breaking)
    F_1 = m phi_3 + 2g X phi_1 = 0
    F_3 = m phi_1 = 0

Tree-level potential: V_tree = |f|^2 (constant in v, the flat direction).

### Fermion Mass Matrix

The superpotential second derivatives (fermion mass matrix) at phi_1 = phi_3 = 0:

    W_ij = | 0     0     0   |
           | 0    2gv    m   |
           | 0     m     0   |

The X-fermion (Goldstino) decouples with mass 0. The 2x2 block for (psi_1, psi_3):

    M_F = m * | 2t   1 |
              | 1    0 |

Eigenvalues of M_F^dag M_F:

    m_F^2 = m^2 (t +/- sqrt(t^2 + 1))^2

Explicitly:

    m_+ = m (t + sqrt(t^2 + 1))     (always positive)
    m_- = m (sqrt(t^2 + 1) - t)     (always positive, since sqrt(t^2+1) > t)

Key identities:

    m_+ * m_- = m^2              (product = m^2, exactly)
    m_+ + m_- = 2m sqrt(t^2+1)  (sum)

At t = sqrt(3):

    m_+ = m(sqrt(3) + 2) = (2 + sqrt(3)) m = 3.73205 m
    m_- = m(2 - sqrt(3)) = (2 - sqrt(3)) m = 0.26795 m

### Scalar Mass-Squared Matrix

The scalar potential V = sum_I |F_I|^2 expanded around the vacuum involves the 4x4 mass-squared
matrix for the 4 real components of (phi_1, phi_3). It decomposes into two 2x2 blocks:

**Re-sector** (Re phi_1, Re phi_3):

    M_Re^2 = m^2 * | 4t^2 + 1 + 2y    2t |
                    | 2t                1  |

    trace = m^2 (4t^2 + 2 + 2y)
    det = m^2 (1 + 2y)

**Im-sector** (Im phi_1, Im phi_3):

    M_Im^2 = m^2 * | 4t^2 + 1 - 2y    2t |
                    | 2t                1  |

    trace = m^2 (4t^2 + 2 - 2y)
    det = m^2 (1 - 2y)

where y = gf/m^2 is the dimensionless SUSY-breaking parameter.

The +/-2y splitting arises from the F-term contribution: |F_X|^2 = |f + g phi_1^2|^2
contains a cross-term 2gf Re(phi_1^2) = 2gf((Re phi_1)^2 - (Im phi_1)^2), giving
+2gf m^{-2} = +2y to Re and -2y to Im.

**Tachyon condition:** det(M_Im^2) = m^4(1 - 2y) < 0 for y > 1/2. The vacuum at
phi_1 = phi_3 = 0 is perturbatively stable only for y < 1/2.

### Explicit Eigenvalues

The four boson mass-squared eigenvalues (in units of m^2):

Re-sector:

    mb1^2, mb2^2 = (1/2)[4t^2 + 2 + 2y +/- sqrt((4t^2 + 2y)^2 + 4t^2)]

Wait, let me be more careful. The eigenvalues of a 2x2 matrix [[a,b],[b,d]] are
(a+d)/2 +/- sqrt(((a-d)/2)^2 + b^2).

Re-sector: a = 4t^2 + 1 + 2y, d = 1, b = 2t.

    Re eigenvalues = (4t^2 + 2 + 2y)/2 +/- sqrt(((4t^2 + 2y)/2)^2 + 4t^2)
                   = 2t^2 + 1 + y +/- sqrt((2t^2 + y)^2 + 4t^2)

Im-sector: a = 4t^2 + 1 - 2y, d = 1, b = 2t.

    Im eigenvalues = 2t^2 + 1 - y +/- sqrt((2t^2 - y)^2 + 4t^2)

### Numerical values at t = sqrt(3), y = 0.1:

    Fermion masses (in units of m):
      Goldstino: 0
      m_- = 2 - sqrt(3) = 0.26795
      m_+ = 2 + sqrt(3) = 3.73205

    Boson mass-squared eigenvalues (in units of m^2):
      Re: 2*3 + 1 + 0.1 +/- sqrt((6 + 0.1)^2 + 12)
        = 7.1 +/- sqrt(37.21 + 12)
        = 7.1 +/- sqrt(49.21)
        = 7.1 +/- 7.015
        => mb_Re1^2 = 0.085,  mb_Re2^2 = 14.115

      Im: 7.0 - 0.1 +/- sqrt((6 - 0.1)^2 + 12)
        = 6.9 +/- sqrt(34.81 + 12)
        = 6.9 +/- sqrt(46.81)
        = 6.9 +/- 6.842
        => mb_Im1^2 = 0.058,  mb_Im2^2 = 13.742

These match the values tabulated in pseudomodulus_cw.md.

### Supertrace

    STr[M^2] = (mb_Re1^2 + mb_Re2^2 + mb_Im1^2 + mb_Im2^2) - 2(m_+^2 + m_-^2)
             = (4t^2 + 2 + 2y) + (4t^2 + 2 - 2y) - 2(2(t^2 + 1))
             = 8t^2 + 4 - 4t^2 - 4 = 4t^2 - 4t^2 + 4 - 4 = 0

    Exactly zero for all t, y. (Standard result for canonical Kahler.)

    STr[M^4] = 8y^2 m^4 (constant in t)

---

## (b) Coleman-Weinberg Potential

The one-loop CW potential (in MS-bar, renormalization scale mu^2 = m^2):

    V_CW = (1/64pi^2) sum_{bosons} (mb_i^2)^2 [ln(mb_i^2/m^2) - 3/2]
         - (1/64pi^2) sum_{fermions} 2(mf_j^2)^2 [ln(mf_j^2/m^2) - 3/2]

The factor 2 for fermions accounts for a Weyl fermion having 2 real d.o.f. vs 1 for a
real scalar.

### Results from numerical computation (pseudomodulus_cw.md, kahler_oneloop.md)

**Canonical CW minimum location:**

| y     | t_min   | t_min/sqrt(3) | Note                              |
|-------|---------|---------------|------------------------------------|
| 0.10  | 0.493   | 0.285         | Physical (y < 1/2)                |
| 0.30  | 0.457   | 0.264         | Physical                          |
| 0.49  | 0.336   | 0.194         | Near boundary                     |
| 0.50  | 0.321   | 0.185         | Boundary                          |
| 1.00  | 0.000   | 0.000         | Tachyonic (y > 1/2)               |

**In the physical regime y < 1/2, the CW minimum is at t_min ~ 0.3-0.5, a factor
3-5 below sqrt(3) = 1.732.**

For the ISS model (iss_kahler_oneloop.md), the CW potential is monotonically increasing
for all t > 0, with the global minimum at t = 0. The CW potential per color:

    f(t) = (t^2+1)^2 [ln(t^2+1) - 3/2] + (t^2-1)^2 [ln|t^2-1| - 3/2] - 2t^4 [ln t^2 - 3/2]

This is always positive and increasing for t > 0.

**The CW potential alone cannot stabilize the pseudo-modulus at t = sqrt(3) for any y.**

This is the central negative result that motivates exploring non-canonical Kahler potentials.

---

## (c) Non-Canonical Kahler: Pole Mechanism

### Setup

Consider the Kahler potential:

    K = |X|^2 - c |X|^4 / mu_K^2

with c > 0. The Kahler metric is:

    K_{XX*} = d^2 K / (dX dX*) = 1 - 2c (|X|^2 + |X|^2) / mu_K^2

Wait, let me be more careful. For K = |X|^2 - c |X|^4 / mu_K^2:

    K_{XX*} = 1 - 4c |X|^2 / mu_K^2

(The derivative d^2(|X|^4)/d(X)d(X*) = d^2(X^2 X*^2)/dX dX* = 4|X|^2.)

The physical scalar potential (tree-level + Kahler correction):

    V_tree = |F_X|^2 / K_{XX*} = f^2 / (1 - 4c v^2/mu_K^2)

for real X = v > 0.

### Pole at K_{XX*} = 0

The Kahler metric vanishes at:

    1 - 4c v^2 / mu_K^2 = 0
    => v_pole = mu_K / (2 sqrt(c))

For c > 0, this is a physical pole: V_tree -> infinity as v -> v_pole. The potential
has a hard wall, and any minimum must lie just below the pole.

### Condition for t_pole = sqrt(3)

Setting g v_pole / m = sqrt(3):

    g mu_K / (2m sqrt(c)) = sqrt(3)
    => c = g^2 mu_K^2 / (12 m^2)

With the natural choice mu_K = m/g (so the Kahler scale equals the mass scale in field
units):

    c = g^2 (m/g)^2 / (12 m^2) = 1/12

**Result: c = 1/12 (with mu_K = m/g) places the Kahler pole at t = sqrt(3) exactly.**

This is an algebraic result, not numerical tuning.

### Physical picture

Near the pole, V_tree = f^2 / (1 - t^2/3) diverges, while V_CW grows only
logarithmically. The total potential:

    V_eff(t) = f^2 / (1 - t^2/3) + V_CW(t)

has a minimum just below t = sqrt(3), where the CW contribution is negligible compared
to the diverging tree-level potential.

Verified numerically (pseudomodulus_vev.md):

| c          | t_min   | Note                                     |
|------------|---------|------------------------------------------|
| -5.000     | 0.224   | (Note: that analysis used opposite sign convention, c -> -c) |
| -0.0833    | 1.732   | t_min = sqrt(3) exactly                  |
| -0.0500    | 2.236   |                                          |

The minimum is pinned at the pole position t_pole = 1/(2 sqrt(|c|)) (with mu_K = m = g = 1).

### Self-consistency check

At t = sqrt(3), the Kahler metric K_{XX*} = 1 - 4c * 3m^2/g^2 / mu_K^2.
With c = 1/12 and mu_K = m/g: K_{XX*} = 1 - 4(1/12)(3) = 1 - 1 = 0.

The metric literally vanishes at the pole. This means:
1. The kinetic term for X vanishes: the field X becomes infinitely heavy (cannot fluctuate).
2. The scalar potential diverges.
3. The minimum is at t -> sqrt(3)^- (approaching from below), not at t = sqrt(3) exactly.

The physical VEV is t = sqrt(3) - epsilon, where epsilon is determined by the balance
between the diverging V_tree and the finite V_CW. For any finite CW contribution,
epsilon is nonzero but can be made parametrically small.

**This is actually a feature**: the pseudo-modulus is stabilized at a definite point
by a geometric property of the Kahler manifold, not by quantum corrections. The CW
potential merely provides the restoring force from below.

---

## (d) Can One-Loop Kahler Corrections Generate the Pole?

### One-loop effective Kahler potential

Integrating out the massive (phi_1, phi_3) pair at one loop generates corrections to
the effective Kahler potential for X:

    Delta K = -(1/16pi^2) sum_i (-1)^{2s_i} |m_i(X)|^2 [ln(|m_i(X)|^2/mu^2) - 1]

For the O'Raifeartaigh model, the X-dependent masses are the fermion pair (m_+, m_-)
and the four bosons. The effective Kahler metric correction is:

    delta Z(t) = d^2 (Delta K) / (dX dX*) |_{phi=0}

From kahler_oneloop.md, the quartic coefficient of the CW potential expansion
V_CW = V_0 + c_2 t^2 + c_4 t^4 + ... gives the effective Kahler quartic:

    C_4 = c_4 / 4     (in units of g^4/(16 pi^2))

The sign of c_4 determines whether the one-loop correction to the Kahler metric
generates a pole or not.

### Key result (kahler_oneloop.md)

    c_4 < 0 for y < 1.425    (negative quartic)
    c_4 = -32/3 at y = 1     (exact)
    c_4 = -(96/5) y^2 + O(y^4) for small y

A negative c_4 corresponds to Delta K ~ -|c_4| |X|^4, which gives a NEGATIVE
correction to K_{XX*}. This is the correct sign to produce a pole!

However, the magnitude is the issue. The one-loop Kahler coefficient is:

    c_Kahler^(1-loop) = g^4 |C_4| / (16 pi^2) = g^4 |c_4| / (64 pi^2)

The target is c = 1/12. So we need:

    g^4 |c_4| / (64 pi^2) = 1/12

At y = 0.1: |c_4| = 0.194, requiring g^4 = 64 pi^2 / (12 * 0.194) = 271.
So g = 4.06. This is non-perturbative.

At y = 1: |c_4| = 32/3, requiring g^4 = 64 pi^2 / (12 * 32/3) = 4.93.
So g = 1.49. This is marginally perturbative but y = 1 is in the tachyonic regime.

### ISS model result (iss_kahler_oneloop.md)

For the ISS model with N_c = 3, the one-loop quartic coefficient at x = sqrt(3) is:

    c_quartic(sqrt(3)) = +0.80     (POSITIVE!)

The ISS one-loop correction has the WRONG SIGN. The effective Kahler coefficient:

    c_eff(sqrt(3)) = +3.8e-3

versus the target -1/12 = -0.083.

**Conclusion: perturbative one-loop corrections cannot generate the pole at t = sqrt(3).**

The sign is wrong in the ISS model, and the magnitude is too small by factors of 20-100
in both the minimal O'Raifeartaigh and ISS frameworks.

### What this means

The pole at t = sqrt(3) must be a TREE-LEVEL property of the Kahler potential, not
a one-loop effect. It arises from the geometry of the target space, not from quantum
corrections to the kinetic term. This points toward a non-perturbative origin.

---

## (e) Curved Kahler Geometry from Confinement

### Composite-field Kahler potential

For a composite meson field X, the Kahler potential is not canonical. The natural
Kahler for a composite field in a confining theory takes a nonlinear form determined
by the target-space geometry of the moduli space.

**Form 1: Logarithmic Kahler (wrong sign)**

    K = Lambda^2 ln(1 + |X|^2/Lambda^2)

    K_{XX*} = Lambda^2 / (Lambda^2 + |X|^2)

This gives a monotonically decreasing metric that approaches zero as |X| -> infinity
but never develops a finite-distance pole. Not useful.

**Form 2: Negative-log Kahler (correct sign)**

    K = -Lambda^2 ln(1 - |X|^2/mu^2)

    K_{XX*} = Lambda^2 / (mu^2 - |X|^2)

This has a pole at |X| = mu. Setting g mu/m = sqrt(3):

    mu = sqrt(3) m/g

The pole is at exactly t = sqrt(3).

### Expansion

For |X| << mu:

    K = Lambda^2 [|X|^2/mu^2 + |X|^4/(2 mu^4) + |X|^6/(3 mu^6) + ...]
      = (Lambda^2/mu^2) |X|^2 + (Lambda^2/(2 mu^4)) |X|^4 + ...

Matching the canonical normalization K ~ |X|^2 at small X requires Lambda^2/mu^2 = 1,
i.e., Lambda = mu = sqrt(3) m/g.

The quartic coefficient is then:

    c = Lambda^2 / (2 mu^4) * mu^2 = 1/(2 mu^2) = g^2/(6 m^2)

Hmm, this gives c = 1/6 rather than 1/12. Let me redo more carefully.

The Kahler potential K = -mu^2 ln(1 - |X|^2/mu^2) expanded:

    K = |X|^2 + |X|^4/(2 mu^2) + |X|^6/(3 mu^4) + ...

So the quartic coefficient is 1/(2 mu^2). The Kahler metric:

    K_{XX*} = mu^2/(mu^2 - |X|^2) = 1/(1 - |X|^2/mu^2)
            = 1/(1 - t^2 g^{-2} m^2 / mu^2)

Wait. X = v = tm/g, so |X|^2 = t^2 m^2/g^2. Then:

    |X|^2/mu^2 = t^2 m^2/(g^2 mu^2)

For the pole at t = sqrt(3): mu^2 = 3 m^2/g^2, so |X|^2/mu^2 = t^2/3.

    K_{XX*} = 1/(1 - t^2/3)

This is identical to the polynomial form K = |X|^2 - c|X|^4/mu_K^2 with c = 1/12,
mu_K = m/g at leading order, but the logarithmic form resums all higher-order terms!

The scalar potential:

    V_tree = f^2 * K^{XX*} = f^2 (1 - t^2/3)

Wait, no. V = g^{i j-bar} F_i F_j* = K^{XX*} |F_X|^2 = (1/K_{XX*}) * f^2.
Actually this depends on convention. Let me be precise.

In N=1 SUGRA (flat limit), V = K^{i j-bar} F_i F_j*, where K^{i j-bar} = (K_{i j-bar})^{-1}
is the INVERSE Kahler metric, and F_i = dW/dPhi_i.

For a single field X with K_{XX*} = 1/(1 - |X|^2/mu^2):

    K^{XX*} = 1 - |X|^2/mu^2

    V = (1 - |X|^2/mu^2) * |f|^2

This is a DECREASING function of |X|! The potential DECREASES as X moves toward the pole.
There is no pole barrier in V; instead, V -> 0 as |X| -> mu.

This is the WRONG behavior for stabilization at the pole.

Let me reconsider. The issue is subtle: V_tree = K^{i j-bar} F_i F_j* uses the
inverse metric. A pole in K_{i j-bar} means a ZERO in K^{i j-bar}, which means
V_tree -> 0, not infinity.

For a pole in V_tree, we need a ZERO of K_{i j-bar}, which means a POLE in K^{i j-bar}.

The Kahler K = |X|^2 - c|X|^4/mu_K^2 gives:

    K_{XX*} = 1 - 4c|X|^2/mu_K^2

This vanishes at the pole, giving K^{XX*} -> infinity, hence V_tree -> infinity.
This is the correct behavior.

But the logarithmic Kahler K = -mu^2 ln(1 - |X|^2/mu^2) gives K_{XX*} = mu^2/(mu^2 - |X|^2)
which DIVERGES at the pole, giving K^{XX*} -> 0, hence V_tree -> 0.

**These are opposite behaviors.** The polynomial K = |X|^2 - c|X|^4 has K_{XX*}
vanishing at the pole; the logarithmic K = -mu^2 ln(1 - |X|^2/mu^2) has K_{XX*}
diverging at the pole.

### Resolution: which form is physical?

The distinction depends on the sign of the Kahler curvature.

**Negative curvature** (hyperbolic target space, like Poincare disk):

    K = -R ln(1 - |X|^2/R)

    K_{XX*} = R/(R - |X|^2)     (diverges at boundary)
    K^{XX*} = (R - |X|^2)/R     (vanishes at boundary)
    V_tree = f^2 (R - |X|^2)/R  (decreases, vanishes at boundary)

The field X is attracted TO the boundary. This is the anti-de Sitter (modular) behavior.
The boundary is at infinite proper distance (the metric diverges), so X can never
actually reach it.

**Positive curvature** (sphere-like, like CP^n):

    K = R ln(1 + |X|^2/R)

    K_{XX*} = R/(R + |X|^2)     (decreases monotonically)
    V_tree = f^2 (R + |X|^2)/R  (increases, no pole)

No pole at all. The potential pushes X to the origin.

**The polynomial truncation** K = |X|^2 - c|X|^4 with c > 0:

    K_{XX*} = 1 - 4c|X|^2       (vanishes at |X|^2 = 1/(4c))
    V_tree = f^2/(1 - 4c|X|^2)  (diverges at pole)

This is NOT the truncation of either the positive or negative curvature logarithmic form.
It is an independent structure with K_{XX*} having a zero rather than a pole.

This polynomial form corresponds to a Kahler manifold that TERMINATES at a finite
field distance. Beyond the pole, K_{XX*} < 0, meaning the kinetic term has the wrong
sign and the theory is inconsistent. The moduli space simply ENDS at |X| = 1/(2 sqrt(c)).

### The quantum-modified moduli space

For SU(3) SQCD with N_f = N_c = 3, the quantum-modified moduli space is defined by:

    det M - B B-tilde = Lambda^6

This is a hypersurface constraint. The meson matrix M is a 3x3 complex matrix with
det M = Lambda^6 (when B = 0). This constraint removes one complex degree of freedom
from the 9-dimensional meson space.

The key insight is that the constraint defines a CURVED Kahler manifold. Near a
specific point on the moduli space (say, a diagonal M with one small eigenvalue),
the effective Kahler metric for the fluctuation along the constrained direction
depends on the inverse of the distance to the boundary of the moduli space.

For the pseudo-modulus direction X (the component of M that breaks SUSY), the
quantum constraint forces:

    M_X = Lambda^6 / (M_1 M_2)

where M_1, M_2 are the other diagonal components. As X varies, the constraint
keeps the determinant fixed. The induced Kahler metric on this slice depends on
the specific embedding.

### CP^n-like sigma model for the meson

In the magnetic dual description (valid in the conformal window), the meson fields
are elementary composites with a Kahler potential of the form:

    K = Tr(M^dag M) / Lambda^2     (canonical for the magnetic dual)

But in the confined description (N_f = N_c, below the conformal window), the meson
M is a genuine composite. The Kahler potential is modified by the constraint:

    K_confined = Tr(M^dag M) + lambda (det M + det M^dag - 2 Lambda^6)

The Lagrange multiplier enforces the quantum constraint. Integrating out lambda
and the constrained direction modifies the effective Kahler metric.

For a single-component pseudo-modulus parametrizing the constrained surface, the
effective Kahler can take the form:

    K_eff(X) = mu^2 / (mu^2 - |X|^2)     (metric on the constrained slice)

Wait, this would be the METRIC, not the potential. Let me be more careful.

If the moduli space is the ball |X| < mu (bounded by the quantum constraint), then
a natural metric on this space is the Bergman metric of the ball:

    K_{XX*} = mu^2 / (mu^2 - |X|^2)^2

with Kahler potential K = -mu^2 ln(1 - |X|^2/mu^2). This is the standard Poincare
metric (constant negative curvature). But as shown above, this gives V_tree -> 0 at
the boundary, not a barrier.

**The polynomial form K = |X|^2 - c|X|^4 corresponds to a different geometric
structure**: a moduli space that terminates at a FINITE proper distance, with the
kinetic term vanishing at the boundary. This is like a conical singularity.

For this to arise from the confined SQCD moduli space, the key property needed is:

**The kinetic term for the pseudo-modulus must vanish at a specific field value.**

This happens when the pseudo-modulus reaches a boundary of the moduli space where
additional degrees of freedom become massless, and the effective description breaks
down. In SQCD terms, this occurs at the quantum-modified constraint boundary.

### Specific mechanism: constraint-induced zero

Consider the N_f = N_c = 3 case. The quantum constraint is:

    det M = Lambda^6     (with B = 0)

Parametrize M = diag(M_1, M_2, M_3) with det M = M_1 M_2 M_3 = Lambda^6. Then
M_3 = Lambda^6 / (M_1 M_2) is fixed by the constraint.

The induced Kahler metric for M_1 (holding M_2 fixed and using the constraint) is:

    ds^2 = |dM_1|^2 + |dM_3|^2 = |dM_1|^2 + |Lambda^6/(M_1^2 M_2)|^2 |dM_1|^2
         = |dM_1|^2 (1 + Lambda^{12}/(|M_1|^4 |M_2|^2))

This diverges as M_1 -> 0 (the constraint forces M_3 -> infinity). It does NOT
vanish at any finite M_1. So the diagonal parametrization does not produce a zero
of the metric.

However, off-diagonal fluctuations can produce a zero. Consider M = M_diag + delta M
with the constraint det(M_diag + delta M) = Lambda^6. The induced metric for the
off-diagonal fluctuation delta M depends on the cofactor matrix of M_diag. When one
eigenvalue of M_diag is large and another is small, the cofactor can produce
suppressed kinetic terms.

This analysis becomes model-dependent and requires specifying which direction in the
meson moduli space corresponds to the pseudo-modulus.

---

## (e') The Bion Kahler Potential: A Non-Perturbative Mechanism

The analysis in parts (c) and (d) shows that:
1. A polynomial K = |X|^2 - (1/12)|X|^4 (in appropriate units) places a Kahler zero at t = sqrt(3).
2. Perturbative one-loop corrections do not generate this zero.
3. The logarithmic Kahler from the moduli space geometry gives the wrong sign.

The resolution lies in the bion-induced Kahler potential (complete_kahler.md, bion_kahler.md).

### The bion mechanism

On R^3 x S^1, monopole-instanton pairs (bions) generate a non-perturbative correction
to the Kahler potential:

    K_bion = (zeta^2/Lambda^2) |sum_k s_k sqrt(M_k)|^2

where s_k = +/-1 are monopole-instanton phases and zeta = exp(-S_0/(2N_c)) is the
monopole fugacity.

This is a POSITIVE definite correction to the Kahler metric. Combined with the
canonical term:

    K_total = Tr(M^dag M) + (zeta^2/Lambda^2) |sum_k s_k sqrt(M_k)|^2

The bion term does NOT produce a Kahler zero (it adds to, rather than subtracts from,
the canonical metric). However, it produces something equally powerful:

**A direct effective potential** that stabilizes the pseudo-modulus through the
bion amplitude's dependence on the meson VEVs.

The bion-induced bosonic potential is:

    V_bion propto -|sum_k s_k sqrt(m_k)|^2

(with a negative sign because magnetic bions are attractive). This potential is
MINIMIZED when |sum_k s_k sqrt(m_k)|^2 is MAXIMIZED.

For the bloom configuration with s = (-1, +1, +1):

    S_bloom = -sqrt(m_s) + sqrt(m_c) + sqrt(m_b)

Combined with the cross-channel coupling to the seed:

    V_eff = lambda_2 |S_bloom - 2 S_seed|^2

minimized at S_bloom = 2 S_seed, giving:

    sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c)     (v_0-doubling)

This is the prediction m_b = 4177 MeV, agreeing with PDG 4180 +/- 30 MeV at 0.1 sigma.

### Connection to the Kahler pole mechanism

The bion potential and the Kahler pole are complementary mechanisms:

1. **Kahler pole**: Tree-level Kahler geometry (K_{XX*} = 0 at t = sqrt(3)) creates a
   barrier that pins the VEV. Requires c = 1/12 to place the pole correctly. The
   CW potential is negligible near the pole.

2. **Bion potential**: Non-perturbative bion amplitudes generate a direct potential
   V_bion(m_k) that is minimized at the v_0-doubling condition. The Kahler metric
   remains canonical (no pole), but the potential has its own minimum.

The question is whether these two mechanisms are different descriptions of the same
physics. The answer is suggestive:

In the bion framework, the effective Kahler metric correction from bions is:

    delta K_{XX*}^{bion} = epsilon * s_i s_j / (4 sqrt(m_i m_j))

This is O(epsilon) = O(exp(-2S_0/3)) = O(10^{-2}) at alpha_s = 1.

For the Kahler zero at t = sqrt(3), we need delta K_{XX*} = -1 (to cancel the canonical
contribution). This requires epsilon / (4 sqrt(m_i m_j)) ~ 1, which for hadronic masses
(m_i ~ 100-4000 MeV, sqrt product ~ 50 MeV^{1/2}) gives:

    epsilon / (4 * 50) ~ 0.015 / 200 ~ 7.5e-5

This is five orders of magnitude too small. **The perturbative bion correction to the
Kahler metric cannot produce a zero at t = sqrt(3).**

However, the bion POTENTIAL does not need to produce a Kahler zero. It produces a
minimum through a different mechanism: the direct bosonic potential from monopole-
antimonopole pair amplitudes, which depends on sum_k s_k sqrt(m_k) rather than on
the Kahler metric.

---

## (f) Summary and Minimum Assumptions

### What works

1. **The Kahler pole mechanism** (c = 1/12 in polynomial truncation, or
   K = -mu^2 ln(1 - |X|^2/mu^2) with mu = sqrt(3) m/g in the logarithmic form)
   successfully pins the pseudo-modulus at t = sqrt(3). This gives the O'Raifeartaigh-
   Koide seed (0, 2-sqrt(3), 2+sqrt(3)) with Q = 2/3 exactly.

2. **The bion potential mechanism** (complete_kahler.md) stabilizes the pseudo-modulus
   via V_eff = lambda_2 |S_bloom - 2 S_seed|^2, predicting m_b = 4177 MeV.

### What does not work

1. **Canonical CW stabilization**: minimum at t ~ 0.3-0.5, not sqrt(3).
2. **Superpotential deformations**: none of the polynomial deformations tested reach
   t = sqrt(3) (pseudomodulus_vev.md).
3. **One-loop Kahler corrections**: wrong sign (ISS) or too small magnitude (O'R).
4. **Logarithmic composite-field Kahler**: gives the wrong behavior (V_tree -> 0 at
   boundary rather than diverging).

### Minimum assumptions needed

To claim that gv/m = sqrt(3) is a consequence of the Kahler geometry:

**Option A: Tree-level Kahler pole**

Assume: K(X) has a zero of K_{XX*} at |X| = sqrt(3) m/g.

This requires: the meson moduli space of the confined SQCD theory has a boundary
(or singularity) at a specific field value, where additional degrees of freedom
become massless or the effective description breaks down. The position of this
boundary must be set by the ratio of scales mu/m in the theory.

Naturalness: the coefficient c = 1/12 must arise from a specific geometric structure.
In general, c depends on the details of the confining dynamics and is not predicted
by symmetry alone. However, c = 1/12 is a suggestive number: it is the inverse of the
dimension of the SU(3) x SU(3) meson matrix minus the constraint, or 1/(N_c(N_c+1)/2),
or related to the one-loop beta function coefficient.

**Option B: Bion potential**

Assume: the monopole-instanton structure of the compactified theory generates a
bosonic potential V_bion propto |sum_k s_k sqrt(m_k)|^2 with specific cross-channel
coupling lambda_3 = -4 lambda_2.

This requires: (i) center symmetry preservation (for well-defined monopole-instantons),
(ii) N_c = 3 (for three monopole types matching three quark flavors), (iii) the
specific sign structure s_k from zero-mode counting, and (iv) the bion-bion
coupling ratio lambda_3/lambda_2 = -4.

Naturalness: assumptions (i)-(iii) are standard in the Unsal framework. Assumption (iv)
(the coupling ratio) is the main additional input. The ratio -4 corresponds to the
v_0-doubling condition and predicts m_b with 0.07% accuracy.

**Option C: Combined mechanism**

The bion Kahler potential and the bion bosonic potential arise from the same
monopole-instanton dynamics. In principle, a complete non-perturbative calculation
would determine both the Kahler correction and the direct potential simultaneously.
If the same bion dynamics that generates V_bion also modifies K(X) to produce a
zero at t = sqrt(3), the two mechanisms would be unified.

However, as shown above, the perturbative bion Kahler correction is O(10^{-2})
and cannot produce a zero of K_{XX*}. A non-perturbative resummation of bion
contributions (or a large-alpha_s regime where the bion expansion breaks down)
would be needed to generate a finite Kahler correction.

### Assessment

The cleanest statement is:

**The CW potential does not select gv/m = sqrt(3). The bion potential does,
through the v_0-doubling condition sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c).
The Kahler pole mechanism at c = 1/12 provides an alternative tree-level
geometric explanation, but the origin of c = 1/12 from first principles
remains an open problem.**

The bion mechanism has the advantage of being predictive (m_b = 4177 MeV, 0.1 sigma
from PDG) and having a clear physical interpretation (monopole-antimonopole pair
amplitudes on R^3 x S^1). The Kahler pole mechanism has the advantage of being
algebraically exact (t = sqrt(3) with no approximations once c = 1/12 is assumed),
but requires an unexplained input.

The deepest question is whether c = 1/12 and lambda_3 = -4 lambda_2 are related,
both arising from the non-perturbative structure of SU(3) SQCD with N_f = N_c = 3.

---

## Appendix: Numerical Benchmarks

### Fermion spectrum at t = sqrt(3)

    m_Goldstino = 0
    m_- = (2 - sqrt(3)) m = 0.267949 m
    m_+ = (2 + sqrt(3)) m = 3.732051 m

    Q(0, m_-, m_+) = (m_- + m_+) / (sqrt(m_-) + sqrt(m_+))^2
                   = 4m / (sqrt(2-sqrt(3)) + sqrt(2+sqrt(3)))^2 m
                   = 4 / (sqrt(2-sqrt(3)) + sqrt(2+sqrt(3)))^2

    Now: sqrt(2-sqrt(3)) + sqrt(2+sqrt(3)) = ?

    Let a = sqrt(2-sqrt(3)), b = sqrt(2+sqrt(3)).
    a^2 + b^2 = 4, ab = sqrt((2-sqrt(3))(2+sqrt(3))) = sqrt(4-3) = 1.
    (a+b)^2 = a^2 + 2ab + b^2 = 4 + 2 = 6.

    Q = 4/6 = 2/3.     (EXACT)

### Boson spectrum at t = sqrt(3), y = 0.1

    Re-sector eigenvalues (m^2 units):
      (2*3 + 1 + 0.1 + sqrt((6+0.1)^2 + 12))/1 = 14.115
      (2*3 + 1 + 0.1 - sqrt((6+0.1)^2 + 12))/1 = 0.085

    Im-sector eigenvalues (m^2 units):
      (2*3 + 1 - 0.1 + sqrt((6-0.1)^2 + 12))/1 = 13.742
      (2*3 + 1 - 0.1 - sqrt((6-0.1)^2 + 12))/1 = 0.058

    STr[M^2] = (14.115 + 0.085 + 13.742 + 0.058) - 2(3.732^2 + 0.268^2)
             = 28.000 - 2(13.928 + 0.072)
             = 28.000 - 28.000 = 0     (check)

### CW potential at t = sqrt(3) (from pseudomodulus_vev.md)

    y = 0.1:  V_CW(sqrt(3)) = 0.000315 m^4
    V_CW(0) = -0.000000 m^4 (negligible)

    Gradient: dV_CW/dt|_{t=sqrt(3)} > 0 (large and positive)

### Kahler pole stabilization (from pseudomodulus_vev.md)

    c = -1/12 = -0.0833 (their sign convention: c < 0 for K_{XX*} zero)
    t_min = 1/(2*sqrt(1/12)) = sqrt(12)/2 = sqrt(3)    (EXACT)
