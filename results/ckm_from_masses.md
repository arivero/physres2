# CKM from Mass Predictions + Fritzsch Texture

## Mass Predictions (from m_s = 93.4 MeV input)
| Quark | Predicted (MeV) | PDG (MeV) | Pull |
|-------|----------------|-----------|------|
| m_c | 1301 | 1270±20 | +1.5σ |
| m_b | 4233 | 4180±30 | +1.8σ |

Using v₀-doubling for m_b: 4177 MeV (+0.07%, −0.1σ) is more precise.

## CKM Results
| Element | Predicted | PDG | Status |
|---------|-----------|-----|--------|
| |V_us| | 0.2236 | 0.2243±0.0008 | **−0.9σ** ✓ |
| |V_cb| | ≥0.059 | 0.0422±0.0008 | **FAILS** (40% too large) |
| |V_ub| | ~0.033 | 0.0039±0.0004 | **FAILS** (8.4× too large) |

## Diagnosis
V_cb failure is STRUCTURAL in the 6-texture-zero Fritzsch ansatz:
√(m_s/m_b) = 0.149, √(m_c/m_t) = 0.088, difference = 0.059
No CP phase can bring |V_cb| below 0.059 (cos ψ > 1 needed)

Known since late 1980s. Requires modified texture (5 zeros, Roberts-Ross)
or RG running from high scale.

## Parameter Count
SM: 10 free (6 masses + 3 angles + 1 phase)
Theory: 4 inputs (m_s, m_u, m_d, m_t) → predicts m_c, m_b, V_us
Reduction: 10 → 4 (but CKM beyond V_us requires better texture)
