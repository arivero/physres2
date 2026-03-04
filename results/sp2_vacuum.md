# Sp(2) SQCD with N_f = 3: Complete Vacuum Analysis

## Convention

Sp(2) denotes the rank-2 symplectic group, equivalent to Sp(4) in the matrix convention (4x4 antisymmetric symplectic form). The fundamental representation is 4-dimensional. With N_f = 3 flavors of quarks Q^i (i = 1,...,3) in the fundamental, the meson matrix M^{ij} = Q^i J Q^j carries 2N_f = 6 flavor indices (from the pseudo-real doubling). The antisymmetric 6x6 matrix M has 15 independent components.

This is the s-confining case: N_f = N_c + 1 for Sp(N_c) with N_c = 2.


## 1. Confined Spectrum

**Mesons:** M^{ij} = Q^i_a J^{ab} Q^j_b, contracted over color (a,b = 1,...,4). Since the symplectic form J is antisymmetric, M^{ij} = -M^{ji}.

**Independent components:** C(6,2) = **15**.

**Quantum-modified constraint (Pfaffian):**

    Pf(M) = Lambda^{2(N_c + 1)} = Lambda^6

where Pf(M) is the Pfaffian of the antisymmetric 6x6 matrix. For any 2n x 2n antisymmetric matrix, Pf(M)^2 = det(M). Here 2n = 6, so Pf(M) is a degree-3 polynomial in the M^{ij}. This is the Sp analogue of the quantum-modified constraint det(M) = Lambda^{2N_f} for SU(N_c).

**Confined superpotential:**

    W_conf = X (Pf(M) - Lambda^6)

where X is a Lagrange multiplier superfield enforcing the quantum constraint. This is the exact low-energy description of the s-confining Sp(2) theory with N_f = 3 (Intriligator-Pouliot 1995).

**One-loop beta coefficient:** b_0 = 3(N_c + 1) - N_f = 9 - 3 = 6 (asymptotically free).


## 2. O'Raifeartaigh Vacuum

### Naive deformation (SUSY preserved)

Adding W_tree = g Tr(M A) + m Tr(M) with antisymmetric background A = diag(a_1 J_2, a_2 J_2, a_3 J_2), the total superpotential is:

    W = X(Pf(M) - Lambda^6) + g Tr(M A) + m Tr(M J)

F-term equations in the Pfaffian-eigenvalue basis:

    F_X = lambda_1 lambda_2 lambda_3 - Lambda^6 = 0
    F_{lambda_k} = X Lambda^6 / lambda_k + g a_k + m = 0

From the second equation: lambda_k = -X Lambda^6 / (g a_k + m). Substituting into the first:

    X^3 = -Prod_k(g a_k + m) / Lambda^{12}

This cubic always has a real solution. **SUSY is preserved.** The naive construction fails because the Pfaffian structure is too accommodating -- unlike the ISS mechanism, there is no rank condition (no magnetic quarks in the s-confining theory).

### Working construction (SUSY broken)

Embed a standard O'Raifeartaigh model into the meson sector:

    W = f X + m M_{12} M_{34} + g X M_{12}^2 + (eps/Lambda^3) Pf(M)

At eps = 0 (pure O'Raifeartaigh):

    F_X    = f + g M_{12}^2  => wants M_{12} != 0
    F_{M34} = m M_{12}        => wants M_{12} = 0
    CONFLICT => SUSY BROKEN

The vacuum is at M_{12} = M_{34} = 0, X = v (pseudo-modulus, flat at tree level). The vacuum energy is V = |f|^2.

At eps != 0: A distant SUSY-restoring vacuum appears at M ~ Lambda^3/eps, but the metastable SUSY-breaking vacuum at the origin persists. The tunneling rate Gamma ~ exp(-(Lambda^3/eps)^4) is exponentially suppressed for small eps. The Pfaffian deformation is cubic in the meson fields, so all its second derivatives vanish at M = 0 -- it does NOT modify the fermion mass matrix at the metastable vacuum.

**SUSY-breaking conditions:**
1. f, g, m all nonzero
2. F_X and F_{M34} are incompatible (quadratic vs linear in M_{12})
3. Stability: gf/m^2 < 1/2 (no tachyons)


## 3. Pfaffian Eigenvalues

The antisymmetric 6x6 matrix M is brought to Darboux (Williamson) normal form by a congruence U^T M U:

    M_canonical = diag(lambda_1 J_2, lambda_2 J_2, lambda_3 J_2)

with J_2 = [[0,1],[-1,0]] and lambda_k >= 0 the **Pfaffian eigenvalues**. These are related to the ordinary (purely imaginary) eigenvalues by: eig(M) = {+/- i lambda_k, k = 1,2,3}.

**Pfaffian:** Pf(M) = lambda_1 lambda_2 lambda_3

**Quantum constraint:** lambda_1 lambda_2 lambda_3 = Lambda^6

**At the O'Raifeartaigh vacuum:**
- The naive construction gives lambda_k = -X Lambda^6/(g a_k + m)
- With a_k = 0 (no flavor breaking): all equal, lambda_k = Lambda^2
- With nonzero a_k: eigenvalues split according to the background


## 4. Koide Condition at gv/m = sqrt(3)

At the metastable O'Raifeartaigh vacuum, the fermion mass matrix in the (M_{12}, M_{34}) sector is:

    W_F = [[2gv,  m],
           [m,    0]]

with eigenvalues m_{pm} = m(t +/- sqrt(t^2 + 1)), t = gv/m.

Full fermion spectrum: **(0 [Goldstino], 0 [flat directions], m_-, m_+)**

The Koide quality factor for the triple (0, m_-, m_+):

    Q = (m_- + m_+) / (sqrt(m_-) + sqrt(m_+))^2
      = sqrt(t^2+1) / (sqrt(t^2+1) + 1)

**Q = 2/3 uniquely at t = sqrt(3):**

Setting Q = 2/3: sqrt(t^2+1) = 2, so t^2 = 3, t = sqrt(3).

At t = sqrt(3):

    m_+ = (2 + sqrt(3)) m = 3.7321 m
    m_- = (2 - sqrt(3)) m = 0.2679 m
    m_+ m_- = m^2

**Proof that Q = 2/3 exactly:**
- Sum(m_i) = 0 + (2 - sqrt(3)) + (2 + sqrt(3)) = 4 (in units of m)
- sqrt(m_-) sqrt(m_+) = sqrt(m_- m_+) = sqrt(1) = 1
- (Sum sqrt(m_i))^2 = (sqrt(2-sqrt(3)) + sqrt(2+sqrt(3)))^2 = 4 + 2 = 6
- Q = 4/6 = 2/3. QED.

**Numerically:** |Q - 2/3| < 10^{-15} (machine epsilon).

The Koide condition gv/m = sqrt(3) is a single algebraic constraint linking the Yukawa coupling g, the pseudo-modulus VEV v, and the mass parameter m.


## 5. Bloom Rotation

### 5a. Q-preserving rotation

The Koide parametrization z_k = 1 + sqrt(2) cos(delta + 2 pi k/3) with m_k = M_0 z_k^2 satisfies **Q = 2/3 identically for all delta** (when z_k > 0 or signed roots are used consistently).

Proof: Sum(z_k) = 3 (cosines cancel), Sum(z_k^2) = 6 (from cos^2 sum = 3/2), so Q = 6/9 = 2/3.

**Result:** ANY pure delta-rotation at fixed M_0 preserves Q = 2/3 exactly. The bloom is a Q-preserving delta-rotation.

### 5b. Delta for lepton hierarchy

The seed is at delta_0 = 3 pi/4 = 135 deg (where z_0 = 0, giving mass = 0 for the first eigenvalue).

Fitting (m_e, m_mu, m_tau) = (0.511, 105.658, 1776.86) MeV to the Koide parametrization:

| Parameter | Value |
|---|---|
| M_0 | 313.84 MeV |
| delta | 132.73 deg |
| Bloom (delta - delta_0) | -2.27 deg |

The physical lepton hierarchy is achieved by rotating delta by only **2.27 degrees** from the seed. The small rotation reflects the extreme hierarchy m_e << m_mu << m_tau.

Bloom parameters for all three Koide triples:

| Triple | M_0 (MeV) | delta (deg) | Bloom (deg) | Q |
|---|---|---|---|---|
| (e, mu, tau) | 313.84 | 132.73 | -2.27 | 0.66666 |
| (-s, c, b) | 912.56 | 157.21 | +22.21 | 0.4585 |
| (c, b, t) | 29576.44 | 123.93 | -11.07 | 0.6695 |

Note: Q values differ from 2/3 for the quark triples because the physical masses do not sit exactly on the Koide manifold. The lepton triple is within 0.001% of Q = 2/3.

### 5c. M_0 for the lepton sector

    M_0 = (sum sqrt(m_k) / 3)^2 = 313.84 MeV

The seed spectrum at delta = 3 pi/4:

    (0, (2-sqrt(3)) M_0, (2+sqrt(3)) M_0) = (0, 84.09, 1171.27) MeV


## 6. Bion Corrections: Sp(2) vs SU(3)

### Root structure comparison

**SU(3) (A_2):**
- Extended Dynkin diagram: triangle (all 3 pairs of nodes adjacent)
- Center symmetry: Z_3
- Magnetic bion types: 3 (related by Z_3)
- Bion potential in delta-space: V ~ cos(3 delta), with 3 equivalent minima per period

**Sp(2) = USp(4) (C_2):**
- Simple roots: alpha_1 = e_1 - e_2 (short), alpha_2 = 2 e_2 (long)
- Affine root: alpha_0 = -2 e_1
- Extended Dynkin diagram: **chain** (alpha_0 =>= alpha_1 =<= alpha_2), NOT a cycle
- Center symmetry: Z_2
- Magnetic bion types: **2** (only adjacent pairs: B_{01} and B_{12}; alpha_0 and alpha_2 are NOT adjacent)
- Bion potential in delta-space: V ~ cos(2 delta), with 2 equivalent minima per period
- Dual Coxeter number: h^v = 3 (comarks all equal to 1)

### v_0-doubling ratio

The v_0-doubling is defined through the quark triple (-s, c, b):

    v_0(seed) = (sqrt(m_s) + sqrt(m_c)) / 3 = 15.10 MeV^{1/2}
    v_0(full) = (-sqrt(m_s) + sqrt(m_c) + sqrt(m_b)) / 3 = 30.21 MeV^{1/2}
    Ratio = 2.0005

**This ratio is a kinematic result, independent of the gauge group.**

The factor of 2 arises from the Koide algebra: when the bloom converts (s, c, 0) to (-s, c, b) by flipping the sign of sqrt(m_s) and adding sqrt(m_b), the sum doubles when sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c). This is a relation among the physical masses, not a property of the gauge dynamics.

**Prediction:** sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c) gives m_b = 4177 MeV (PDG: 4180 +/- 30, deviation 0.07%).

**Result: v_0(full)/v_0(seed) = 2.0005 for BOTH SU(3) and Sp(2).**

The gauge group affects the dynamical mechanism of the bloom (which potential V(delta) selects the physical angle) but not the kinematic v_0 ratio:

- SU(3): cos(3 delta) potential from Z_3 center, 3 minima per period. Seed at delta = 135 deg is 75 deg from the nearest minimum at 120 deg.
- Sp(2): cos(2 delta) potential from Z_2 center, 2 minima per period. Seed at delta = 135 deg is **45 deg** from the nearest minimum at 180 deg.

The Sp(2) bloom rotation toward the nearest bion minimum is shorter (45 deg vs 75 deg), suggesting the Sp(2) dynamics may be more "natural" for the bloom, but the endpoint mass spectrum and v_0 ratio are determined by the physical masses, not the gauge group.

For comparison, the lepton sector has v_0(full)/v_0(seed) = 1.014 (not 2), because the sign flip is absent and sqrt(m_e) is negligible compared to the other roots.


## Summary of Main Results

| Part | Result |
|---|---|
| 1. Mesons | 15 independent components (antisymmetric 6x6) |
| 1. Constraint | Pf(M) = Lambda^6 |
| 1. Superpotential | W = X(Pf(M) - Lambda^6) |
| 2. SUSY breaking | W = fX + m M_{12} M_{34} + g X M_{12}^2 (O'R mechanism) |
| 2. F-term | F_X = f != 0 at metastable vacuum |
| 3. Pfaffian evals | lambda_1 lambda_2 lambda_3 = Lambda^6; split by background |
| 4. Koide Q | Q = 2/3 exactly at gv/m = sqrt(3) |
| 4. Spectrum | (0, (2-sqrt(3))m, (2+sqrt(3))m) |
| 5a. Q-preserving | All delta-rotations at fixed M_0 preserve Q = 2/3 |
| 5b. Lepton bloom | delta = 132.73 deg, bloom = -2.27 deg from seed |
| 5c. Lepton M_0 | 313.84 MeV |
| 6. Sp(2) bions | 2 types (chain Dynkin), Z_2 center, cos(2 delta) potential |
| 6. v_0 doubling | 2.0005 (same as SU(3); kinematic, not gauge-dependent) |
