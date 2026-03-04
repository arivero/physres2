# Bloom Mechanism: Seed to Full Triple

## Parametrization
The Koide parametrization z_k = v0 + r*cos(2*pi*k/3 + delta) gives:
- Seed: delta = 3*pi/4 (where cos(delta) = -1/sqrt(2) forces z_0 = 0)
- r = sqrt(2)*v0 when Q = 2/3

## Key Results

### 1. Q-Preserving Bloom
A Koide-preserving bloom is a pure delta-rotation at fixed v0 and r = sqrt(2)*v0.
Total mass M = 6*v0^2 is conserved. Mass redistributes among eigenvalues.

### 2. Sign Flip Is Instantaneous
The seed sits exactly at the boundary z_0 = 0. Any infinitesimal delta-rotation
generates nonzero z_0 with definite sign.

### 3. R-Breaking Instability
Explicit R-symmetry breaking (delta_W ~ epsilon*Phi_0*Phi_2) lifts goldstino mass
as m_0 ~ epsilon^2/m, but pushes Q away from 2/3. No nonzero epsilon preserves
the Koide condition. The seed is an UNSTABLE FIXED POINT.

### 4. CW Mass Too Small
Coleman-Weinberg mass m_CW ~ g^2*f/(16*pi^2*m) ~ 0.073 MeV with hadronic
parameters. Physical m_b = 4180 MeV is 5 orders of magnitude larger.
Bloom must be nonperturbative.

### 5. v0 Doubling (Quark Sector Only)
v0(full)/v0(seed) = 2.0005 for quarks. Predicts m_b = 4177 MeV (0.1 sigma).
Does NOT hold for mesons (ratio 1.87) or leptons (ratio 1.01).

### 6. ISS Transmits Koide
The ISS mechanism (N_c=3, N_f=5) produces pseudo-moduli with CW masses
proportional to quark masses. It TRANSMITS the Koide condition from the UV
quark mass spectrum to the IR, but does not generate it.

## Bloom Parameters (Quark Sector)
| Parameter | Seed (0,s,c) | Full (-s,c,b) | Change |
|-----------|-------------|---------------|--------|
| v0 | 15.10 | 30.21 | x2.00 |
| r | 21.28 | 43.25 | x2.03 |
| delta | 135.2 deg | 157.2 deg | +22 deg |
| Q | 0.664 | 0.675 | +1.6% |
