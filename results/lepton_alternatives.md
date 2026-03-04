# Alternative Mechanisms for Lepton Mass Generation in the SQCD Framework

## Preamble: What Has Failed

The direct approach -- identifying charged leptons with light mesino eigenstates of the
W_IJ fermion mass matrix at the Seiberg vacuum -- has been conclusively ruled out:

- **Perturbative mesino masses**: The 6x6 W_IJ gives three light eigenvalues at
  (7.7 x 10^{-4}, 7.8 x 10^{-4}, 1.04 x 10^{-2}) MeV with Q = 0.481 (not 2/3).
  The mass ratios are completely unlike (m_e : m_mu : m_tau).

- **Kahler corrections**: Bion-induced delta_K_ij ~ epsilon * sqrt(m_i m_j)/Lambda^2 ~ 10^{-7},
  negligible. A complete scan over flavor-dependent Kahler parameters finds Q in [0.35, 0.50],
  never reaching 2/3. This is a rank invariance obstruction: K^{-1/2} W K^{-1/2} cannot
  change the rank of W_IJ.

- **ISS Coleman-Weinberg**: CW masses for broken mesons are nearly degenerate
  (m_CW ~ 27 MeV for all flavors). Q(CW) = 1/3 (degenerate limit). The affine shift
  f_a = h mu^2 - m_a does not preserve Koide.

The lepton mass mechanism must come from somewhere else.

---

## Mechanism 1: Separate Confining Sector for Leptons

### Assessment

**Principle**: Leptons arise from a SEPARATE confining gauge group G_L, distinct from
the SU(3) SQCD that produces mesons/baryons. This sector has its own confinement scale
Lambda_L and its own O'Raifeartaigh dynamics that generate Q = 2/3 through the
gv/m = sqrt(3) mechanism.

**Does it work in principle?** Yes. The O'Raifeartaigh-Koide mechanism is algebraically
clean: any O'Raifeartaigh model W = f Phi_0 + m Phi_a Phi_b + g Phi_0 Phi_a^2 produces
the seed triple (0, 2-sqrt(3), 2+sqrt(3)) x m when gv/m = sqrt(3). This is a
one-parameter constraint and is robust.

**Group and matter content**: The minimal choice is SU(2)_L with N_f = 2 flavors,
giving a confining theory with mesons M^i_j (2x2 matrix), baryons, and a quantum
constraint det M - BB_tilde = Lambda_L^4. The meson matrix has 4 components; with the
constraint and Lagrange multiplier X, the superpotential has the same O'Raifeartaigh
structure. However, N_f = N_c = 2 gives only a 2x2 meson matrix (3 eigenvalues including
the singlet). For three charged leptons, you need at least three independent mass eigenvalues.

**Alternative**: Sp(2) = SU(2) with N_f = 3 antisymmetric flavors. The mesons form a
6x6 antisymmetric matrix (15 components in the symmetric product). The quantum constraint
is Pf(M) = Lambda^4. This gives 3 independent eigenvalues from the 3 Pfaffian
eigenvalue pairs. The O'Raifeartaigh structure generalizes naturally with three mass
parameters (m_1, m_2, m_3).

**Does it preserve Q = 2/3?** Yes, if the Sp(2) O'Raifeartaigh vacuum satisfies the
analogue of gv/m = sqrt(3) for each eigenvalue. The key point is that the Koide seed
is a property of a 3-field O'Raifeartaigh model with one coupling condition. The Sp(2)
model with N_f = 3 naturally provides three mass eigenvalues.

**Free parameters**: Lambda_L (sets overall lepton mass scale M_0 = 313.8 MeV),
delta (the bloom phase, observed at 132.73 deg), and the coupling condition gv/m = sqrt(3).
The hierarchy m_tau/m_e ~ 3475 comes from the bloom rotation (2.27 deg from seed).

**Analogue of meson = quark bilinear**: In the Sp(2) model, "leptonic mesons"
L^{ij} ~ psi^i psi^j (antisymmetric) are composites of the confining preons psi^i.
The diagonal VEVs L^{ii} = Lambda_L^2 / m_i^{preon} give the lepton mass eigenvalues
through the seesaw. The preon masses m_i^{preon} are the UV parameters; the Koide
condition constrains their ratios.

**Verdict**: Feasible in principle. Natural group-theoretic structure. But introduces
an entirely new confining sector with no direct connection to the quark SQCD. The
Koide relation for leptons and quarks would then be coincidental (both arising from
separate O'Raifeartaigh mechanisms) unless a UV unification links the two sectors.

---

## Mechanism 2: Non-Renormalizable Operators

### Assessment

**Principle**: Add higher-dimension operators coupling the SQCD meson fields to lepton
superfields:

    W_NR = (c_ij / M_*) M^i_j L_bar_i e_j

where L, e are lepton superfields and M_* is a UV cutoff (e.g., M_Pl or a GUT scale).

At the Seiberg vacuum, the meson VEVs are <M^i_i> = C/m_i (the seesaw). This gives
lepton Yukawa couplings:

    y_{ij}^{lep} = c_ij * C / (m_i * M_*)

**Does it generate the hierarchy?** The meson VEVs span:

    <M_u> = Lambda^2/m_u = 41667 MeV
    <M_d> = Lambda^2/m_d = 19272 MeV
    <M_s> = Lambda^2/m_s = 964 MeV

The ratio M_u/M_s = m_s/m_u = 43.2. This is much smaller than m_tau/m_e = 3475.
With a single meson insertion, the hierarchy is insufficient by two orders of magnitude.
Even with the Wilson coefficient c_ij, you would need c_tau/c_e ~ 80 to compensate.
This is not naturally explained.

With TWO meson insertions (dimension-7 operator):

    W_NR = (c_ij / M_*^3) M^i_k M^k_j L_bar_i e_j

the effective coupling goes as (Lambda^2/m_i)^2 / M_*^3, and the ratio becomes
(m_s/m_u)^2 ~ 1870, closer to the lepton hierarchy. But this still requires
fine-tuning of c_ij.

**Does it produce Q = 2/3?** Not generically. The lepton masses become:

    m_l^i = v * c_ii * C / (m_i * M_*)

For diagonal c_ij, the lepton masses are proportional to 1/m_i (the seesaw masses).
Q(1/m_u, 1/m_d, 1/m_s) = 0.443, far from 2/3. Q(1/m_d, 1/m_s, 1/m_b) = 0.665,
close to 2/3, but this uses down-type quarks (d, s, b), not the (u, d, s) of the SQCD.

For Q = 2/3 to emerge from c_ij, one needs

    c_ij such that (c_11/m_u, c_22/m_d, c_33/m_s) lies on the Koide manifold.

This is three free parameters constrained by one equation. It works but explains nothing.

**Free parameters**: c_ij (9 complex parameters for a general matrix), M_*. Over-parametrized.

**Verdict**: Technically possible but over-parametrized. The Koide relation becomes
an accidental tuning of Wilson coefficients rather than a dynamical output. Does not
explain why Q = 2/3.

---

## Mechanism 3: Radiative Lepton Masses

### Assessment

**Principle**: Leptons are massless at tree level and acquire masses through loop diagrams
involving the SQCD mesons/mesinos. The loop suppression factor alpha/(4 pi) ~ 10^{-3}
could naturally explain the smallness of lepton masses relative to the quark condensate
scale.

**Loop topology**: The simplest diagram has a lepton line emitting a virtual Higgs, which
couples to a meson loop. The meson propagator involves the VEV <M^i_i> = Lambda^2/m_i.
The one-loop contribution is:

    delta_m_l ~ (y_l * y_M / (16 pi^2)) * m_mesino * log(Lambda^2/m_mesino^2)

where y_l is the lepton-Higgs-meson coupling and y_M is the meson self-coupling. The
mesino masses are the eigenvalues of W_IJ, which are of order Lambda^5/(m_i m_j) ~ MeV to GeV.

**Flavor structure**: The loop factor is generically flavor-dependent through the meson
propagator. If the dominant contribution comes from the diagonal mesons M^i_i, the
lepton mass matrix is diagonal with entries proportional to f(m_i). The flavor dependence
enters through the meson mass m_mesino(i) ~ Lambda^4/(m_i m_j), which is inversely
proportional to quark masses. This produces a seesaw-like hierarchy.

**Does it preserve Q = 2/3?** The loop function f(m_i) would need to map the quark mass
spectrum to a Koide triple. Since the loop involves logarithms and power-law functions
of the meson masses, the Koide condition would require a specific functional form of
the loop integrand. There is no generic reason for this to happen.

However, if the dominant loop topology is a SUSY-breaking insertion on the meson line
(from the F-term F_X), the mass formula becomes:

    delta_m_l^i ~ (alpha/(4 pi)) * F_X / m_mesino^i

With m_mesino^i ~ 1/m_q^i (seesaw), this gives delta_m_l^i ~ m_q^i, reproducing the
quark hierarchy in the lepton sector. But Q(m_u, m_d, m_s) = 0.567, not 2/3.

**Free parameters**: The coupling constants at the loop vertices (at least 2), the
SUSY-breaking insertion F_X, and the UV cutoff.

**Verdict**: The loop suppression factor is qualitatively right for the lepton mass
scale. But the flavor structure generically maps the quark spectrum or its seesaw
to the leptons, neither of which gives Q = 2/3 for (u, d, s) flavors. Would need
specific loop topology engineered to produce Koide, which is unnatural.

---

## Mechanism 4: Koide as UV Boundary Condition

### Assessment

**Principle**: The Koide relation Q = 2/3 is imposed at a UV scale (the confinement
scale Lambda or a higher scale) as a boundary condition on the lepton Yukawa matrix.
The question is whether Q = 2/3 survives RG running to the lepton pole masses.

**RG running of Q**: Under one-loop SUSY RG (gauge + Yukawa), the mass ratios
evolve as:

    d(m_i)/dt = (gamma_i + gamma_H) m_i

where gamma_i are anomalous dimensions. In the MSSM, gauge contributions are
flavor-universal and cancel in Q. The Yukawa contributions are proportional to
y_i^2 ~ m_i^2/v^2, so heavier leptons run faster.

The key result (established in the project): Q flows toward 1/3 under one-loop
SUSY RG by the Cauchy-Schwarz inequality. Specifically:

    dQ/dt = -(2 Q / (sum sqrt(m))^2) * sum_i gamma_i^Y * (sqrt(m_i) - v0)^2

where gamma_i^Y is the Yukawa anomalous dimension. Since (sqrt(m_i) - v0)^2 >= 0,
Q always decreases (Q -> 1/3 = degenerate limit).

**Running from Lambda to m_tau**: For the charged leptons with MSSM-like running:

    delta_Q ~ -(y_tau^2 / (16 pi^2)) * ln(Lambda/m_tau) * Q * (spread factor)

With y_tau ~ 0.01, Lambda ~ 300 MeV (or higher), ln(Lambda/m_tau) ~ 5:

    delta_Q ~ -(10^{-4}) * 5 * (2/3) * O(1) ~ -3 x 10^{-4}

This is a 0.05% shift, comparable to the experimental deviation of Q from 2/3
(which is 0.001%). So the shift is about 50x larger than the observed deviation.

**Does Q = 2/3 survive?** For purely leptonic running from Lambda = 300 MeV,
the shift is small but non-negligible. From a GUT scale (Lambda ~ 10^{16} GeV),
ln(Lambda/m_tau) ~ 40, and delta_Q ~ -0.003, a 0.4% shift. This is 400x larger
than the observed 0.001% deviation. Q = 2/3 at the GUT scale would NOT survive
to Q = 2/3 at low energy.

However, if the boundary condition is set at a LOW scale (e.g., Lambda ~ few hundred MeV,
the confinement scale), the running is minimal. This is consistent with the SQCD framework
where confinement produces the Koide seed at Lambda ~ 300 MeV.

**Conditions for survival**: The boundary scale must be O(GeV) or lower, so that
the RG running is insufficient to destroy Q = 2/3. In the SQCD framework, Lambda ~ 300 MeV
is natural. Below that scale, only QED running applies, which shifts individual masses
by O(alpha/pi) ~ 0.2% but preserves Q to O(alpha^2/pi^2) ~ 10^{-5} (because
QED running is nearly flavor-universal for charged leptons).

**Free parameters**: The boundary scale Lambda_BC and the initial values of the
lepton Yukawa ratios at Lambda_BC. The Koide condition is one constraint, leaving
two free parameters (e.g., M_0 and delta).

**Verdict**: This is the most natural interpretation if the confinement scale is
low (O(100 MeV) to O(GeV)). Q = 2/3 set at the confinement scale survives to
low energy because QED running is nearly universal. However, this mechanism does not
EXPLAIN why Q = 2/3; it merely states that whatever dynamics produces Q = 2/3 at
the confinement scale, the condition is preserved by low-energy running. The
explanation is deferred to the UV dynamics (which must be one of the other mechanisms).

---

## Mechanism 5: Lepton-Meson Duality

### Assessment

**Principle**: In SU(5) flavor, mesons are in 5 x 5_bar = 1 + 24. Leptons are in
the 5_bar (or 10) of SU(5). Is there a group-theoretic map from the meson eigenvalue
spectrum to lepton masses?

**The 24 decomposition under SU(3)_flavor x U(1)**: The adjoint 24 of SU(5) decomposes
as 24 -> (8,0) + (3, 5/3) + (3_bar, -5/3) + (1, 0) under SU(3) x SU(2) x U(1).
The singlet component is the traceless part of the diagonal. This decomposition does
not directly relate meson eigenvalues to lepton quantum numbers.

**Casimir approach**: The quadratic Casimir of SU(3) is C_2 = (p^2 + pq + q^2 + 3p + 3q)/3
for representation (p,q). For the fundamental 3: C_2 = 4/3. For the adjoint 8: C_2 = 3.
These are group-theory constants, not mass ratios.

A more interesting approach: the EIGENVALUES of the meson matrix M^i_j at the vacuum
are Lambda^2/m_i. The sum of eigenvalues is Tr(M) = Lambda^2 sum(1/m_i). The
characteristic polynomial of M has roots Lambda^2/m_u, Lambda^2/m_d, Lambda^2/m_s.
The elementary symmetric polynomials are e_1 = Tr(M), e_2 = sum_{i<j} M_i M_j,
e_3 = det(M) = Lambda^6.

The Koide quotient of the seesaw eigenvalues is Q(1/m_u, 1/m_d, 1/m_s) = 0.443.
For Q(1/m_d, 1/m_s, 1/m_b) = 0.665, close to 2/3. But this uses the DOWN-type
quarks, not the light quarks (u,d,s) of the N_f = 3 SQCD.

**Free parameters**: None beyond the quark masses (already specified). But the
connection is numerical, not group-theoretic.

**Verdict**: No clean group-theoretic map from meson eigenvalues to lepton masses
has been identified. The dual Koide Q(1/m_d, 1/m_s, 1/m_b) = 0.665 is suggestive but
uses different quarks than the SQCD with N_f = 3 light flavors. If one extends to
N_f = 5 (including b), the Seiberg seesaw produces meson VEVs M_i ~ 1/m_i whose
Koide quotient for the down-type subset (d,s,b) is near 2/3. But this is the
dual Koide observation, not a lepton mass mechanism.

---

## Mechanism 6: Seiberg Duality Cascade

### Assessment

**Principle**: The lepton sector arises from a DIFFERENT SQCD with its own Seiberg
duality. For example, if leptons come from SU(2) SQCD with N_f = 4 (ISS window:
N_c < N_f < 3N_c/2, i.e., 2 < 4 < 3), the magnetic dual has SU(N_f - N_c) = SU(2)
magnetic colors and its own meson fields.

Alternatively, the quarks and leptons both arise from a SINGLE UV SQCD but at
different stages of a duality cascade, where the effective N_c and N_f change at
each step.

**ISS for SU(2), N_f = 4**: The magnetic dual has SU(2) gauge group, 4 flavors,
16 meson singlets. The ISS vacuum has rank condition: 2 flavors get VEVs, 2 have
F != 0. The CW potential lifts the broken directions. The meson spectrum has
O(4) symmetry broken to O(2) x O(2), giving 2 pseudo-Goldstone pairs and
heavy mesons. This produces two mass scales, not three.

For three lepton generations, one needs at least 3 independent mass eigenvalues.
SU(2) with N_f = 4 gives 2 independent CW masses (from the 2 F-term directions).
Not enough.

**SU(2), N_f = 5**: This is outside the ISS window (N_f = 5 > 3N_c/2 = 3), so
SUSY is not broken at tree level in the ISS sense. The magnetic dual is SU(3) with
N_f = 5, which is the SAME theory as the quark sector. No new physics.

**SU(3), N_f = 4**: ISS window (3 < 4 < 4.5). Magnetic SU(1) = U(1), 4 flavors.
This was already analyzed (iss_cw_koide.md): CW masses are nearly degenerate
(Q ~ 1/3). Does not produce Koide.

**Duality cascade**: In a Klebanov-Strassler type cascade, the effective (N_c, N_f)
changes at each step: (N_c, N_f) -> (N_f - N_c, N_f) -> ... The cascade
terminates when N_c = 0 or N_f < N_c. Starting from (3, 5):
Step 1: (3,5) -> (2,5) (outside ISS for SU(2))
Step 2: (2,5) -> (3,5) (back to original)

The cascade is periodic, not terminating. For (3,4):
Step 1: (3,4) -> (1,4) (U(1) with 4 flavors, confines)
Terminates after one step. No interesting cascade.

**Does it produce Q = 2/3?** The ISS CW mechanism generically produces Q ~ 1/3
because the affine shift F_a = h mu^2 - m_a is degenerate when all m_a << h mu^2.
A cascade does not change this.

**Free parameters**: The (N_c, N_f) of the lepton sector, Lambda_L, and the
lepton-sector quark masses.

**Verdict**: Cascades from realistic SQCD theories either cycle or terminate too
quickly to produce new structure. The ISS CW mechanism generically gives Q ~ 1/3,
not 2/3. This approach does not solve the lepton Koide problem.

---

## Identification of the Two Most Promising Approaches

### Ranking

| Mechanism | Q = 2/3? | Natural hierarchy? | Predictive? | Overall |
|-----------|----------|-------------------|-------------|---------|
| 1. Separate sector | Yes (by construction) | Yes (bloom) | Moderate | **Best** |
| 2. Non-renormalizable | Requires tuning | Insufficient | Over-parametrized | Poor |
| 3. Radiative | No generic reason | Qualitatively right | Low | Moderate |
| 4. UV boundary | Preserved if low-scale | Deferred | High (if UV known) | **Second best** |
| 5. Lepton-meson duality | Wrong quarks | No | None | Poor |
| 6. Duality cascade | No (Q ~ 1/3) | No | None | Poor |

**Most promising**: Mechanism 1 (Separate confining sector) and Mechanism 4 (UV boundary
condition), with Mechanism 4 being understood as the condition that Mechanism 1 must
satisfy.

---

## Detailed Analysis: Mechanism 1 + 4 Combined

### The Proposal

Leptons arise from a separate Sp(2) SQCD sector with N_f = 3 antisymmetric flavors.
At the confinement scale Lambda_L ~ O(100 MeV), the O'Raifeartaigh-Koide mechanism
sets Q = 2/3 as a SEED (with one zero mass). The bloom mechanism then rotates delta
from the seed value 3 pi/4 to the physical value 132.73 deg, generating all three
nonzero lepton masses while preserving Q = 2/3 identically (since Q is a geometric
property of the Koide manifold, not a dynamical condition).

Below the confinement scale, QED running is nearly flavor-universal and preserves
Q = 2/3 to O(alpha^2/pi^2) ~ 10^{-5}, consistent with the observed deviation of
9 ppm.

### Sp(2) Model Specification

**Gauge group**: Sp(2) = SU(2) (isomorphic for rank 1)

**Matter content**: 3 chiral superfields Psi^i (i = 1,2,3) in the fundamental 2 of Sp(2).

**Composite operators**: The antisymmetric meson matrix

    L^{ij} = epsilon^{ab} Psi^i_a Psi^j_b / Lambda_L

with L^{ij} = -L^{ji}. This has 3 independent components: L^{12}, L^{13}, L^{23}.

**Quantum constraint** (for Sp(2) with N_f = 3):

    Pf(L) = Lambda_L^{2(N_f - N_c - 1)} = Lambda_L^2

where Pf is the Pfaffian of the antisymmetric 6x6 matrix (extended by the symplectic
form). More precisely, for Sp(2) with N_f = 3, the instanton-generated superpotential
replaces the constraint:

    W_dyn = 1/Lambda_L^3 * Pf(L)

This is the Affleck-Dine-Seiberg (ADS) superpotential for Sp(N_c) with N_f = N_c + 1 = 3.

**Mass deformation**:

    W = sum_{i<j} mu_{ij} L^{ij} + W_dyn + W_{O'R}

where mu_{ij} are preon mass parameters and W_{O'R} is the O'Raifeartaigh coupling.

### O'Raifeartaigh Structure

The superpotential with a singlet X (Lagrange multiplier/NMSSM-like):

    W = sum_{i<j} mu_{ij} L^{ij} + X * (Pf(L) / Lambda_L^3 - v_L^2) + g X L^{12} L^{12}

At the vacuum with L^{ij} != 0, the fermion mass matrix has the same structure as
the standard O'Raifeartaigh model. The condition gv_X/mu = sqrt(3) produces the
Koide seed (0, 2-sqrt(3), 2+sqrt(3)) x mu.

### Identification of Leptons

The three independent Pfaffian eigenvalues of L^{ij} map to the three charged lepton
masses:

    m_e <-> L^{23} (lightest, from bloom)
    m_mu <-> L^{13} (intermediate)
    m_tau <-> L^{12} (heaviest)

The lepton-Higgs coupling is mediated by the operator:

    W_Yuk = (h_L / M_*) L^{ij} H_d e_R^k epsilon_{ijk}

where e_R are the right-handed lepton superfields and M_* is the mediation scale.
The effective Yukawa couplings are:

    y_l^k = h_L * <L^{ij}> / M_* * epsilon_{ijk}

### Numerical Check: Koide from the Seed

Set mu = 1 (arbitrary units), gv/mu = sqrt(3).

Seed eigenvalues: (0, (2-sqrt(3)), (2+sqrt(3))) = (0, 0.2679, 3.7321).

    Sum m = 4.0
    Sum sqrt(m) = sqrt(0.2679) + sqrt(3.7321) = 0.5176 + 1.9319 = 2.4495
    (Sum sqrt(m))^2 = 5.9999
    Q = 4.0 / 6.0 = 2/3 exactly.

This matches the O'Raifeartaigh-Koide result (oraifeartaigh_koide.md).

### Bloom to Physical Masses

On the Koide manifold, z_k = 1 + sqrt(2) cos(delta + 2 pi k/3), and
m_k = M_0 * z_k^2. The physical lepton delta = 132.73 deg, M_0 = 313.84 MeV.

The bloom rotation delta_seed -> delta_phys = 135.00 deg -> 132.73 deg is a
2.27 deg shift. The bion + instanton potential V(delta) governs this:

    V(delta) = A * (9/2)[5 - 4 cos(delta_0 - delta)] + B cos(9 delta + phi)

The bion restoring force dominates (delta sits 2.27 deg from seed, 12.73 deg from
the nearest instanton minimum at 120 deg).

### The v_0-Doubling Test

For leptons, v_0(full)/v_0(seed) = 1.014, NOT 2.0. The v_0-doubling relation
does NOT hold for leptons. This is a critical distinction from quarks (where the
ratio is 2.0005).

**Interpretation**: The separate leptonic sector has different bion dynamics.
In the quark sector, the bion Kahler generates the effective potential
V_mon ~ |S_bloom - 2 S_seed|^2 whose factor of 2 produces v_0-doubling.
In the lepton sector, the analogous bion has a different coefficient (near 1,
not 2), reflecting the different gauge group (Sp(2) vs SU(3)) and different
instanton structure.

Specifically, for Sp(N_c) the monopole-instanton action is S_0 = 2 pi / (alpha_s N_c),
with N_c = 1 for Sp(2). The bion fugacity is exp(-2 S_0/N_c) = exp(-4 pi/alpha_s),
which is more suppressed than in SU(3) (exp(-2 pi/(3 alpha_s))). The weaker bion
restoring force means the bloom is a smaller angular rotation, consistent with
the observed 2.27 deg (vs ~ 22 deg for the quark bloom from seed (0,s,c) to (-s,c,b)).

### Mass Scale Determination

The overall scale M_0 = 313.84 MeV is set by Lambda_L. In the O'Raifeartaigh model:

    M_0 = mu^2 / (something involving Lambda_L and couplings)

For the lepton masses to be correct, we need M_0 v_EW / M_* ~ m_tau, so:

    M_* ~ M_0 v_EW / m_tau = 313.84 * 246220 / 1776.86 = 43,480 MeV ~ 43.5 GeV

This is a low mediation scale, below the weak scale. This could be problematic
unless the mediator is part of the Higgs sector itself.

Alternatively, if the lepton masses are directly the mesino masses (not mediated by a
Higgs), then M_0 = 313.84 MeV is the mass scale and the lepton masses are:

    m_e = 313.84 * z_e^2 = 0.511 MeV
    m_mu = 313.84 * z_mu^2 = 105.66 MeV
    m_tau = 313.84 * z_tau^2 = 1776.86 MeV

where z_k are the Koide parametrization values. This works BY DEFINITION of the
parametrization. The question is what dynamics sets M_0 = 313.84 MeV.

If Lambda_L ~ 300 MeV (same order as Lambda_QCD ~ 300 MeV), then M_0 ~ Lambda_L
naturally. The near-coincidence of the lepton M_0 with Lambda_QCD is suggestive
of a common dynamical origin.

### Connection Between Quark and Lepton Sectors

The quark sector uses SU(3) SQCD with N_f = 3 (or 5). The lepton sector uses
Sp(2) with N_f = 3. Both produce O'Raifeartaigh-Koide seeds. The coupling between
the sectors occurs through the Higgs fields H_u, H_d, which are shared:

    W = W_quark(M, B, B_tilde, X, H_u, H_d) + W_lepton(L, X_L, H_d, e_R) + W_Higgs

The Higgs sector mediates EWSB to both sectors. The Koide relation in each sector
is independent (each has its own O'Raifeartaigh mechanism), but the EW scale v = 246 GeV
is shared.

### Critical Question: Why Are There Two Sectors?

In the SU(5) flavor framework, quarks are in the 5 and leptons are in the 5_bar.
The separate confining sector for leptons requires explaining why the 5_bar leptons
confine under a DIFFERENT gauge group than the 5 quarks. This is not automatic in
SU(5) GUT unification.

A possible resolution: the SU(5) flavor is a GLOBAL symmetry, not a gauge symmetry.
The quarks are charged under SU(3)_color and the leptons are charged under a
DIFFERENT gauge group (Sp(2)_L or similar). The SU(5)_flavor is a classification
symmetry that organizes both sectors, not a dynamical gauge symmetry.

This is consistent with the sBootstrap framework where SU(5)_flavor is the flavor
symmetry of 5 light quarks (u, c, d, s, b), not a GUT gauge group.

### Numerical Predictions

If the lepton sector is Sp(2) with the O'Raifeartaigh-Koide mechanism:

1. **Q = 2/3 exactly** at the confinement scale, preserved to 10^{-5} by QED running.
   Prediction: Q(physical) = 2/3 - O(10^{-5}). Observed: Q = 0.666661 (9 ppm from 2/3).
   Consistent.

2. **delta = 132.73 deg** is a dynamical output of the bion + instanton potential.
   The closeness of delta mod (2pi/3) to 2/9 (33 ppm) suggests a rational structure
   in the potential.

3. **M_0 = 313.84 MeV** is set by the Sp(2) confinement scale Lambda_L.
   Prediction: Lambda_L = O(300 MeV). This is the same order as Lambda_QCD.

4. **Free parameters**: Lambda_L (1 parameter), delta (1 parameter, dynamically determined).
   All three lepton masses are predicted from these two inputs. Alternatively:
   one input mass (e.g., m_tau) plus the Koide condition determines the other two.
   This leaves ONE free parameter for 3 masses, i.e., 2 predictions.

---

## Summary

| Approach | Explains Q = 2/3? | Hierarchy? | Parameters | Status |
|----------|-------------------|------------|------------|--------|
| Separate Sp(2) sector | Yes (O'R-Koide) | Yes (bloom) | 2 (Lambda_L, delta) | Most promising |
| UV boundary at Lambda ~ 300 MeV | Preserved (QED universal) | Deferred | 2 (M_0, delta) | Complementary |
| Non-renormalizable ops | Tuning only | Poor | 9+ | Ruled out as explanation |
| Radiative | No | Qualitative only | 4+ | Unlikely |
| Lepton-meson duality | Wrong quarks | No | 0 | Ruled out |
| Duality cascade | No (Q ~ 1/3) | No | 3+ | Ruled out |

**Recommended path forward**: Develop the Sp(2) O'Raifeartaigh-Koide model in detail.
The key computations needed are:

1. **Sp(2) N_f = 3 vacuum structure**: Verify that the ADS superpotential + mass
   deformation produces an O'Raifeartaigh vacuum with gv/mu = sqrt(3).

2. **Bion potential for Sp(2)**: Compute the monopole-instanton structure on R^3 x S^1
   and verify that it produces a delta-dependent potential with a minimum near
   delta = 132.73 deg.

3. **Mediation to SM leptons**: Specify the coupling between the Sp(2) composites
   and the SM lepton fields. Determine whether the mediator scale M_* is consistent
   with EWSB and flavor constraints.

4. **FCNC constraints**: Check whether the Sp(2) sector generates dangerous
   lepton-flavor-violating operators (mu -> e gamma, etc.) at acceptable levels.
