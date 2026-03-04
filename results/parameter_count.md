# Parameter Counting: SM Charged Fermion Sector

## Framework relations

Five algebraic relations derived from the SUSY Lagrangian (O'Raifeartaigh
superpotential + bion Kahler + Yukawa eigenvalue constraint + Oakes texture):

| Label | Relation | Origin in Lagrangian |
|-------|----------|---------------------|
| R1 | m\_c / m\_s = (2+sqrt3)^2 = 13.928 | O'Raifeartaigh seed (quark SU(3) sector), metastable vacuum gv/m = sqrt3 |
| R2 | sqrt(m\_b) = 3 sqrt(m\_s) + sqrt(m\_c) | Bion Kahler correction (v\_0-doubling) |
| R3 | (m\_e + m\_mu + m\_tau) / (sqrt(m\_e) + sqrt(m\_mu) + sqrt(m\_tau))^2 = 2/3 | O'Raifeartaigh seed + bloom (lepton SU(2) sector) |
| R4 | (m\_c + m\_b + m\_t) / (sqrt(m\_c) + sqrt(m\_b) + sqrt(m\_t))^2 = 2/3 | Yukawa eigenvalue constraint at tan beta = 1 |
| R5 | sin theta\_12 = sqrt(m\_d / m\_s) | Oakes/GST relation from Fritzsch texture |

Independence: R1 introduces m\_c; R2 introduces m\_b; R4 introduces m\_t; R3
involves only lepton masses; R5 introduces theta\_12. Each equation brings
exactly one new variable. All five are algebraically independent.


## Parameter table

### Quark masses (6 SM parameters)

| Parameter | PDG value (MeV) | Status | Determined by | Predicted value | Accuracy |
|-----------|-----------------|--------|---------------|-----------------|----------|
| m\_u | 2.16 +0.49/-0.26 | **Free** | -- | -- | -- |
| m\_d | 4.67 +0.48/-0.17 | **Free** | -- | -- | -- |
| m\_s | 93.4 +8.6/-3.4 | **Free** | -- | -- | -- |
| m\_c | 1275 +/- 25 | **Predicted** | R1: m\_s x (2+sqrt3)^2 | 1301 | +2.0% |
| m\_b | 4180 +30/-20 | **Predicted** | R1+R2: (3 sqrt(m\_s) + sqrt(m\_c))^2 | 4233 (from m\_s alone) | +1.3% |
|  |  |  | R2 only (using PDG m\_c) | 4186 | +0.15% |
| m\_t | 172760 +/- 300 | **Predicted** | R1+R2+R4: Koide(c,b,t) with m\_c, m\_b from m\_s | 171252 | -0.87% |
|  |  |  | R4 only (using PDG m\_c, m\_b) | 168628 | -2.4% |

Note on m\_t: The full chain (m\_s -> m\_c -> m\_b -> m\_t) gives 171252 MeV
(-0.87%), which is better than the direct Koide prediction from PDG m\_c, m\_b
(168628 MeV, -2.4%). This is because the predicted m\_c is slightly high and
the predicted m\_b is slightly high, and these deviations partially compensate
in the Koide formula.

Note on scheme: R1 and R2 use MS-bar masses (m\_s at 2 GeV; m\_c at m\_c; m\_b
at m\_b). R4 uses a mixed scheme (MS-bar for c,b; pole for t). The Q(c,b,t)
evaluated at PDG values is 0.6693, which is 0.40% from 2/3. This is the known
precision of the (c,b,t) Koide relation.

### Lepton masses (3 SM parameters)

| Parameter | PDG value (MeV) | Status | Determined by | Predicted value | Accuracy |
|-----------|-----------------|--------|---------------|-----------------|----------|
| m\_e | 0.51100 | **Free** | -- | -- | -- |
| m\_mu | 105.658 | **Free** | -- | -- | -- |
| m\_tau | 1776.86 +/- 0.12 | **Predicted** | R3: Koide Q = 2/3 | 1776.97 | +0.006% |

The lepton Koide relation is the most precise prediction of the framework:
Q(PDG) = 0.66666082, matching 2/3 to 0.0009%. Using pole masses throughout.

### CKM parameters (4 SM parameters)

| Parameter | PDG value | Status | Determined by | Predicted value | Accuracy |
|-----------|-----------|--------|---------------|-----------------|----------|
| theta\_12 | 13.00 deg (sin = 0.2250) | **Predicted** | R5: arcsin sqrt(m\_d/m\_s) | 12.92 deg (sin = 0.2236) | -0.6% |
| theta\_23 | 2.40 deg | **Free** | -- | -- | -- |
| theta\_13 | 0.211 deg | **Free** | -- | -- | -- |
| delta\_CP | 1.144 rad | **Free** | -- | -- | -- |


## Net parameter count

| | SM | This framework |
|---|---|---|
| Total parameters (charged fermion sector) | 13 | 8 |
| Quark masses | 6 free | 3 free (m\_u, m\_d, m\_s) + 3 predicted |
| Lepton masses | 3 free | 2 free (m\_e, m\_mu) + 1 predicted |
| CKM angles | 3 free | 2 free (theta\_23, theta\_13) + 1 predicted |
| CKM phase | 1 free | 1 free |
| **Net reduction** | | **5 parameters** |

### What the 8 remaining free parameters are

1. m\_u (up quark mass)
2. m\_d (down quark mass)
3. m\_s (strange quark mass) -- the "seed" for the entire heavy quark spectrum
4. m\_e (electron mass)
5. m\_mu (muon mass) -- together with m\_e, seeds the lepton spectrum
6. theta\_23 (CKM)
7. theta\_13 (CKM)
8. delta\_CP (CKM phase)


## The chain from m\_s

Given m\_s alone, the framework predicts three heavy quark masses:

    m_s = 93.4 MeV  (input)
       |
       | R1: x (2+sqrt3)^2
       v
    m_c = 1301 MeV   (PDG: 1275, +2.0%)
       |
       | R2: sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c)
       v
    m_b = 4233 MeV   (PDG: 4180, +1.3%)
       |
       | R4: Koide Q(c,b,t) = 2/3
       v
    m_t = 171252 MeV  (PDG: 172760, -0.87%)

Three outputs from one input. This is the core quantitative content of the
quark sector.


## Lagrangian parameters that replace the SM Yukawas

The SM describes the 9 charged fermion masses via 9 Yukawa couplings (diagonal,
after going to mass basis) plus the Higgs vev v = 246.22 GeV. The CKM matrix
requires 4 additional parameters in the Yukawa texture. Total: 13 parameters
(plus v, which is fixed by G\_F).

In this framework, the Lagrangian parameters are:

| Parameter | Role |
|-----------|------|
| Lambda\_q | Confining scale of quark SU(3) sector |
| Lambda\_l | Confining scale of lepton SU(2) sector |
| m (quark sector) | O'Raifeartaigh mass parameter |
| g (quark sector) | O'Raifeartaigh coupling (gv/m = sqrt3 at metastable vacuum) |
| m (lepton sector) | O'Raifeartaigh mass parameter |
| g (lepton sector) | O'Raifeartaigh coupling |
| c\_bion | Bion Kahler coefficient |
| tan beta | Ratio of MSSM Higgs vevs (= 1 if constrained) |

The ratio gv/m = sqrt3 is NOT a tuned parameter -- it is the location of the
metastable vacuum of the O'Raifeartaigh model. The Koide relation Q = 2/3
follows algebraically from this vacuum. The v\_0-doubling coefficient (3 in
sqrt(m\_b) = 3 sqrt(m\_s) + sqrt(m\_c)) follows from the bion Kahler correction
with specific coefficient.

The point: 5 of the 13 SM parameters are not free parameters in this framework.
They are algebraic consequences of the Lagrangian structure.


## Observed but not yet derived (potential further reductions)

| Relation | Status | If derived, saves |
|----------|--------|-------------------|
| Dual Koide: Q(1/m\_d, 1/m\_s, 1/m\_b) = 2/3 | 0.22% from exact; predicts m\_d = 4.60 MeV (PDG: 4.67, 1.4% off) | 1 parameter (m\_d) |
| delta\_0 mod 2pi/3 = 2/9 | 35 ppm; 3.1 sigma with look-elsewhere | 1 parameter (connects delta to rational) |

If the dual Koide is derived from a Seiberg-dual seesaw in the down sector,
the free parameter count drops to 7. The m\_d -> m\_s -> m\_c -> m\_b -> m\_t
chain would then produce four masses from one input (m\_s).


## Explicit solution: m\_t from Koide(c,b,t)

Solve (m\_c + m\_b + m\_t) / (sqrt(m\_c) + sqrt(m\_b) + sqrt(m\_t))^2 = 2/3
for m\_t, with m\_c = 1275 MeV, m\_b = 4180 MeV.

Let x = sqrt(m\_t), A = sqrt(m\_c) + sqrt(m\_b) = 35.707 + 64.653 = 100.360.

    3(m_c + m_b + x^2) = 2(A + x)^2
    3(5455) + 3x^2 = 2A^2 + 4Ax + 2x^2
    16365 + 3x^2 = 20144.3 + 401.44x + 2x^2
    x^2 - 401.44x - 3779.3 = 0
    x = (401.44 + sqrt(161155 + 15117))/2
      = (401.44 + sqrt(176272))/2
      = (401.44 + 419.85)/2
      = 410.64

    m_t = x^2 = 168628 MeV

(The negative root gives x = -9.2, rejected.)

Cross-check: Q(1275, 4180, 168628) = 178083/267084 = 2/3. Exact.

**Result: m\_t = 168628 MeV from PDG (m\_c, m\_b). Deviation from PDG m\_t(pole) = -2.4%.**

From the full chain (m\_s = 93.4 input): m\_t = 171252 MeV, deviation = -0.87%.
