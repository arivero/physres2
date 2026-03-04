# Meson Koide Look-Elsewhere Analysis

## Setup
- 11 pseudoscalar mesons (pi, K, eta, eta', D, D_s, B, B_s, B_c, eta_c, eta_b)
- C(11,3) = 165 triples x 4 sign patterns = 660 combinations
- Q = (sum m) / (sum s_k sqrt(m_k))^2 -- standard Koide convention (Q = 2/3)

## Results

### Top hits (sorted by deviation from 2/3)
| Triple | Signs | Q | Deviation |
|--------|-------|---|-----------|
| (-pi, D_s, B) | -++ | 0.667359 | 0.10% |
| (-pi, D_s, B_s) | -++ | 0.667632 | 0.14% |
| (-pi, D_s, eta_c) | -++ | 0.669988 | 0.50% |
| (-pi, D_s, B_c) | -++ | 0.671068 | 0.66% |
| (-pi, D, B) | -++ | 0.672774 | 0.92% |
| (-pi, D, B_s) | -++ | 0.673073 | 0.96% |

Within 0.5%: 3 triples. Within 0.2%: 2 triples.

### Pattern
ALL top-6 hits have structure: (-pi, D_(s), B_(s,c)) -- negative pion + D-meson + B-meson.
This is the flavor structure predicted by the sBootstrap.

### Monte Carlo (10M random triples from [100, 10000] MeV)
- p(within 0.10% of 2/3) = 0.021% per trial
- p(within 0.50% of 2/3) = 0.10% per trial

### Look-elsewhere corrected significance
- N_trials = 660
- p_single = 0.021%
- p_corrected = 0.130 (1.12 sigma)

Not statistically significant after look-elsewhere correction.

### Caveats
1. The physically motivated subset (negative pion + D + B mesons) gives
   N_trials ~ 16, which would give p_corrected ~ 0.003 (2.7 sigma).
   But this restriction is post-hoc.
2. The (-pi, D_s, B_s) near-miss at 0.14% is physically distinct from
   (-pi, D_s, B) -- it substitutes b -> bs -- and independently close.
