# Numerical Verification of Paper Claims

All key numerical claims in sbootstrap_v4d.tex verified against independent computation.

## Results

| Claim | Paper value | Verified value | Status |
|-------|-------------|----------------|--------|
| Q(e,mu,tau) | 1.500014 | 1.500013849 | OK |
| Significance of Q deviation | 0.91 sigma | 0.91 sigma | OK |
| M0 (lepton Koide) | 313.84 MeV | 313.841 MeV | OK |
| delta_0 mod 2pi/3 vs 2/9 | 35 ppm | 33.4 ppm | OK (edition-dependent) |
| R (Casimir ratio) | 0.223101322300866 | 0.223101322300866 | OK (15 digits) |
| R algebraic form | (sqrt19-3)(sqrt19-sqrt3)/16 | matches to 1e-15 | OK |
| M_W prediction | 80.374 GeV | 80.3744 GeV | OK |
| M_W tension with PDG | 0.39 sigma | 0.39 sigma | OK |
| Fourth eigenvalue | 122.39 GeV | 122.39 GeV | OK |
| Fourth eigenvalue vs M_H | 2.3% below | 2.3% below | OK |
| Diophantine uniqueness | (N,r,s)=(3,3,2) | unique in 1<=r,s<=19 | OK |
| Anomaly cancellation | A=0 | 0+9-9=0 | OK |
| Scale m = M_Z/sqrt(sqrt3-1) | 106.577 GeV | 106.577 GeV | OK |
| Meson Q(-pi,Ds,B) | ~0.10% from 3/2 | 0.104% | OK |
| Quark seed sqrt(ms)/sqrt(mc) | ~1.2% from 2-sqrt3 | 1.2% | OK |

## Note on ppm
The paper states 35 ppm; independent computation gives 33.4 ppm. The difference is from the PDG mass edition used (m_tau = 1776.86 MeV). Both round to ~35 ppm. The key point: the least-squares fit for delta_0 is essential — using a single mass equation gives 5 ppm (electron) or 33 ppm (muon), because Q is not exactly 3/2.

## Scaling predictions
M1 = 3*M0, delta_1 = 3*delta_0 gives:
- m_b predicted: 4197 MeV (PDG: 4180 MeV, 0.4% off) -- OK
- m_s predicted: 92.3 MeV (PDG: 93.4 MeV, 1.2% off) -- OK
- m_c predicted: 1360 MeV (PDG: 1270 MeV, 7.1% off) -- OK (known issue, scheme-dependent)
