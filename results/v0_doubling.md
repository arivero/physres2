# v0-Doubling in the Quark Bloom

## Result
The vacuum charge v0 = (1/3) sum(z_k) doubles during the quark seed-to-full transition:
- v0(seed) = (sqrt(m_s) + sqrt(m_c))/3 = 15.10 MeV^{1/2}
- v0(full) = (-sqrt(m_s) + sqrt(m_c) + sqrt(m_b))/3 = 30.21 MeV^{1/2}
- Ratio: 2.0005

## Prediction
If v0 doubles exactly: sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c)
- Predicted m_b = 4177 MeV (PDG: 4180 +/- 30 MeV, 0.1 sigma)
- Precision: 0.07% (more precise than overlap prediction at 0.5%)

## Look-elsewhere
Among 5 random log-uniform masses in [2, 5000] MeV, searching all
relations sqrt(m_i) = a*sqrt(m_j) + b*sqrt(m_k) with |a|,|b| <= 5:
- 24.9% produce a relation as good as 0.07%
- The v0-doubling is NOT statistically significant in isolation
- Significance comes from the structural interpretation within sBootstrap

## Other integer relations among quark sqrt-masses (within 1%)
1. sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c): 0.035%  [THE v0-doubling]
2. sqrt(m_c) = -3*sqrt(m_s) + sqrt(m_b): 0.064%  [same relation rearranged]
3. sqrt(m_b) = -3*sqrt(m_d) + 2*sqrt(m_c): 0.21%
4. sqrt(m_c) = -2*sqrt(m_u) + 4*sqrt(m_s): 0.23%
5. sqrt(m_c) = 3*sqrt(m_d) + 3*sqrt(m_s): 0.45%

## Meson sector
The v0-doubling does NOT hold for mesons:
v0(full)/v0(seed) = 1.87 (not 2)

## Physical interpretation
The bloom doubles v0, which means the total signed charge sum(z_k)
also doubles. This is NOT a consequence of the Koide condition (which
constrains ratios, not sums). It is an additional, independent constraint
on the mass spectrum.
