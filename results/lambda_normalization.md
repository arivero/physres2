# Lambda Normalization and the Higgs Mass in Seiberg-Confined SQCD

## Part (a): Dimensions of lambda_0

The superpotential has dimension [W] = mass^3. The term in question is:

    W superset lambda_0 X H_u . H_d

We know:
- [X] = mass^{-3} (Lagrange multiplier enforcing a dimension-6 constraint)
- [H_u . H_d] = mass^2 (product of two canonical scalars)

Therefore:

    [lambda_0] = [W] / ([X][H_u . H_d]) = mass^3 / (mass^{-3} . mass^2) = mass^4

Writing lambda_0 = lambda_hat . Lambda^n with [lambda_hat] = mass^0:

    **n = 4, so lambda_0 = lambda_hat . Lambda^4**

This makes physical sense: the coupling carries four powers of the confinement scale because X is a highly composite operator (dimension -3), and the coupling must compensate to make the superpotential dimension 3.


## Part (b): Canonical field X_c

We want X_c = X . Lambda^p with [X_c] = mass (canonical scalar dimension).

Since [X] = mass^{-3}:

    [X_c] = [X] . [Lambda^p] = mass^{-3} . mass^p = mass^{p-3}

Setting this equal to mass^1:

    **p = 4, so X_c = X . Lambda^4**

Equivalently, X = X_c / Lambda^4. Substituting into the superpotential:

    W = sum_i m_i M^i_i + (X_c / Lambda^4)(det M - BB_tilde - Lambda^6) + lambda_0 (X_c / Lambda^4) H_u . H_d + m_B BB_tilde

The constraint term becomes:

    (X_c / Lambda^4)(det M - BB_tilde - Lambda^6)

And the Higgs coupling term becomes:

    lambda_0 / Lambda^4 . X_c H_u . H_d = (lambda_hat . Lambda^4 / Lambda^4) X_c H_u . H_d = lambda_hat . X_c H_u . H_d

**The coupling of the canonical singlet X_c to H_u H_d is simply lambda_hat (dimensionless).**


## Part (c): The NMSSM Higgs mass

In the NMSSM, the superpotential contains:

    W superset lambda_NMSSM S H_u . H_d

where S is a canonical singlet ([S] = mass) and lambda_NMSSM is dimensionless.

From Part (b), when we identify S = X_c, the superpotential term is:

    lambda_hat . X_c H_u . H_d

Comparing directly:

    **lambda_NMSSM = lambda_hat = lambda_0 / Lambda^4**

The tree-level Higgs mass at tan beta = 1 is:

    m_h^2 = lambda_NMSSM^2 v^2 = lambda_hat^2 v^2


## Part (d): Numerical evaluation

**Step 1: lambda_NMSSM from the Higgs mass.**

    m_h = lambda_NMSSM . v . |sin(2 beta)|

At tan beta = 1, sin(2 beta) = 1, so:

    lambda_NMSSM = m_h / v = 125400 / 246220 = 0.5093

**Step 2: lambda_hat.**

    lambda_hat = lambda_NMSSM = 0.5093

This is just a restatement: lambda_hat IS lambda_NMSSM.

**Step 3: Consistency check.**

The constraint for positive meson VEVs requires lambda_hat < 2.97 x 10^{-6}.

We need lambda_hat = 0.5093.

    **0.5093 >> 2.97 x 10^{-6} by a factor of 1.7 x 10^5**

This is a spectacular failure: the Higgs mass demands a coupling that is five orders of magnitude larger than what the meson vacuum can tolerate.


## Part (e): The conflict

The required lambda_hat ~ 0.5 versus the constraint lambda_hat < 3 x 10^{-6} means the naive identification X_c = S (NMSSM singlet) cannot simultaneously give the correct Higgs mass and maintain a healthy meson vacuum. Let us assess the three options.

### Option 1: X is not the NMSSM singlet

This is the cleanest resolution. One introduces a separate elementary singlet S with a standard NMSSM coupling lambda_NMSSM S H_u H_d, while X remains a non-dynamical Lagrange multiplier enforcing the Seiberg constraint. The Higgs mass comes from the NMSSM mechanism with lambda_NMSSM ~ 0.5, and X has either no coupling to the Higgs sector or a tiny coupling lambda_hat ~ 10^{-6} consistent with positive meson VEVs.

Cost: one additional elementary field S and its coupling lambda_NMSSM as a free parameter. The Seiberg-confined sector decouples from EWSB at tree level.

Assessment: m_h = 125 GeV is straightforwardly obtained. This is the standard NMSSM with an added confined sector.

### Option 2: X acquires a kinetic term

Radiative corrections or Kahler potential terms (e.g., K superset c X^dagger X with a coefficient determined by strong dynamics) could give X a propagating degree of freedom. If the kinetic normalization is Z_X, the physical canonical field is X_c = sqrt(Z_X) . X, and the effective coupling becomes:

    lambda_eff = lambda_0 / sqrt(Z_X)

For this to give lambda_eff ~ 0.5 with lambda_0 = lambda_hat . Lambda^4 and lambda_hat << 10^{-6}, one needs Z_X to be enormously large, which is unnatural.

Alternatively, if Z_X is of order Lambda^8 (which is what dimensional analysis for a composite kinetic term would suggest), then X_c ~ X . Lambda^4 as before, and we recover the same problem.

Cost: adds a Kahler potential parameter and makes the constraint only approximate (F_X != 0). The meson vacuum shifts, potentially destabilizing the theory.

Assessment: does not resolve the conflict unless Z_X takes an unnaturally large value.

### Option 3: Higgs mass from radiative corrections (MSSM-like)

In the MSSM without any singlet coupling, the Higgs mass at tree level at tan beta = 1 is:

    m_h^{tree} = m_Z |cos(2 beta)| = 0

The full 125 GeV must come from radiative corrections, primarily the top-stop loop:

    Delta m_h^2 ~ (3 y_t^4 v^2) / (8 pi^2) ln(m_stop^2 / m_t^2)

This requires m_stop ~ 1-10 TeV depending on mixing, with no direct connection to the Seiberg sector.

Cost: heavy stops (fine-tuning), and the Seiberg sector contributes nothing to the Higgs mass. The coupling lambda_0 X H_u H_d can then be dropped entirely (lambda_hat = 0) or kept tiny, preserving the meson vacuum.

Assessment: m_h = 125 GeV is obtained via the standard MSSM mechanism. This decouples the confined SQCD sector from EWSB entirely. Predictivity from the Seiberg sector regarding the Higgs is lost.

### Summary

The conflict is genuine and severe. The most economic resolution is **Option 1**: X is not the NMSSM singlet. The Seiberg-confined sector and the electroweak Higgs sector require separate singlets operating at vastly different scales (Lambda ~ 300 MeV versus v ~ 246 GeV). Attempting to unify them into a single field creates an irreconcilable tension of five orders of magnitude.


## Part (f): F-term with and without kinetic term

### Case 1: X has no kinetic term (Lagrange multiplier)

The equation of motion for X (varying the action with respect to X) gives:

    dW/dX = det M - BB_tilde - Lambda^6 + lambda_0 H_u . H_d = 0

This is an algebraic constraint, not a dynamical equation. It enforces:

    det M = Lambda^6 + BB_tilde - lambda_0 H_u . H_d

Now, F_X is defined as:

    F_X^* = -dW/dX

But since X has no kinetic term, the scalar potential does not contain |F_X|^2. The standard formula V = sum |F_i|^2 sums only over fields with kinetic terms. Since X is not dynamical:

    **F_X does not contribute to the potential. The constraint is enforced exactly, and there is no F_X^2 energy cost.**

More precisely: for a non-dynamical field, the "F-term" is not an independent degree of freedom. The equation dW/dX = 0 is a constraint, not an energy contribution. We can say F_X = 0 in the sense that the constraint ensures dW/dX = 0.

If lambda_0 = 0 (or is negligible), the constraint reduces to:

    det M = Lambda^6 + BB_tilde

Setting B = B_tilde = 0, the meson vacuum satisfies det M = Lambda^6 with meson VEVs M^i_i = Lambda^2 (for diagonal M).

If lambda_0 != 0, the constraint becomes det M = Lambda^6 - lambda_0 <H_u . H_d> (setting B = 0). With <H_u . H_d> = v^2/2, this shifts the meson VEVs. For the meson VEVs to remain positive, we need:

    lambda_0 v^2/2 << Lambda^6, i.e., lambda_hat << 2 Lambda^2 / v^2 ~ 3 x 10^{-6}

which is precisely the bound stated in the problem.

### Case 2: X has a canonical Kahler potential K superset X^dagger X

Now X is a dynamical field with a propagating degree of freedom. The scalar potential includes:

    V superset |F_X|^2 = |dW/dX|^2 = |det M - BB_tilde - Lambda^6 + lambda_0 H_u . H_d|^2

The constraint is no longer enforced exactly. Instead, the potential has a minimum where F_X is minimized but generically nonzero:

    F_X = det M - Lambda^6 + lambda_0 H_u . H_d (setting B = 0)

At the electroweak vacuum with <H_u . H_d> = v^2/2:

    F_X = det M - Lambda^6 + lambda_0 v^2/2

For F_X = 0, we would need det M = Lambda^6 - lambda_0 v^2/2, but the other F-term conditions (from dW/dM^i_j) also constrain the meson VEVs. In general, all F-terms cannot vanish simultaneously, and:

    **F_X != 0 generically. SUSY is broken by the meson-Higgs coupling.**

The physical implications are:
- The meson VEVs deviate from det M = Lambda^6.
- F_X != 0 breaks SUSY, generating soft masses proportional to F_X / M_Pl (if gravity-mediated) or through gauge mediation if the mesons are charged.
- The |F_X|^2 term in the potential creates a runaway direction unless stabilized by other terms (e.g., the meson mass terms m_i M^i_i).
- With lambda_hat ~ 0.5, the F_X contribution would be of order Lambda^6, completely destabilizing the meson vacuum.

This reinforces the conclusion from Part (e): the Seiberg Lagrange multiplier X should not be promoted to a dynamical NMSSM singlet.
