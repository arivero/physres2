# O'Raifeartaigh Model Produces Exact Koide Seed

## Result
The minimal O'Raifeartaigh superpotential W = f Phi_0 + m Phi_a Phi_b + g Phi_0 Phi_a^2 produces a Koide seed triple (0, m_-, m_+) with Q = 2/3 EXACTLY when gv/m = sqrt(3), where v is the pseudo-modulus VEV.

## Proof
At the O'R vacuum (phi_a = 0, phi_0 = v), the fermionic mass matrix is:
W_ij = [[0, 0, 0], [0, 2gv, m], [0, m, 0]]

Eigenvalues: 0 (goldstino), m_pm = sqrt(g^2 v^2 + m^2) pm gv.

Setting sqrt(m_-)/sqrt(m_+) = 2 - sqrt(3) (the Koide seed ratio) requires:
sqrt(t^2+1) - t = 2 - sqrt(3), where t = gv/m.

Solving: t = (1 - (2-sqrt(3))^2) / (2(2-sqrt(3))) = (2sqrt(3)-3)/(2-sqrt(3)).
Rationalizing with (2+sqrt(3)): t = sqrt(3).

Then t^2+1 = 4, so m_pm = (2 pm sqrt(3)) m.

The Koide ratio: Q = (0 + (2-sqrt(3))m + (2+sqrt(3))m) / (0 + sqrt((2-sqrt(3))m) + sqrt((2+sqrt(3))m))^2
= 4m / (sqrt(2-sqrt(3)) + sqrt(2+sqrt(3)))^2 m = 4/6 = 2/3.

The "6" comes from: (a+b)^2 = a^2 + b^2 + 2ab = (2-sqrt(3)) + (2+sqrt(3)) + 2*1 = 6,
using a*b = sqrt((2-sqrt(3))(2+sqrt(3))) = sqrt(4-3) = 1.

## Significance
- The Koide condition on the seed is ONE algebraic constraint on the superpotential: gv = sqrt(3) m
- The hierarchy m_s/m_c = (2-sqrt(3))^2 ~ 0.07 emerges from the MINIMAL model, not from "additional dynamics"
- The constraint links the cubic coupling g, the pseudo-modulus VEV v, and the bilinear mass m
- This is the microscopic content of the energy-balance condition applied to the seed

## Numerical verification
- t = gv/m = 1.7320508076 = sqrt(3) to machine epsilon
- m_- = 0.267949 m, m_+ = 3.732051 m
- Q = 0.6666666667 = 2/3 exactly
- sqrt(m_-)/sqrt(m_+) = 0.2679491924 = 2-sqrt(3) exactly
