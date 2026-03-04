# Chiral Spectrum of T⁶/Z₃ Type I Orientifold

## Setup
- T⁶/Z₃ with SU(3) root lattice, twist (1/3,1/3,1/3)
- γ_θ = diag(ω·1₅, ω²·1₅, 1₃, 1₃)
- Ω swaps 5-blocks (A↔B) and 3-blocks (C↔D)
- Gauge group: SU(5) × SU(3) × U(1)²

## D9-Brane Chiral Spectrum

| Mult | SU(5) | SU(3) | U(1)_A | U(1)_C | Source |
|------|-------|-------|--------|--------|--------|
| 3    | 10    | 1     | +2     | 0      | (A,B) sector |
| 3    | 5̄    | 3     | -1     | +1     | (C,A)/(D,B) sector |
| 3    | 5     | 3̄    | +1     | -1     | (B,C)/(A,D) sector |

## D5-Brane Contribution (from anomaly cancellation)

| Mult | SU(5) | SU(3) | Source |
|------|-------|-------|--------|
| 3    | 5̄    | 1     | D5-branes at fixed points |

## Anomaly Check
- SU(5)³: 3(1) + 9(-1) + 9(1) = 3 from D9 → needs D5 contribution
- With D5: 3 + 3(-1) = 0 ✓
- SU(3)³: 9 - 9 = 0 ✓

## Generation count: 3
From the 3 complex dimensions of T⁶, each contributing one chiral multiplet.

## CRITICAL TENSION: 10 vs 15
The (A,B) sector gives ANTISYMMETRIC 10 of SU(5), not symmetric 15.
The orientifold sign argument (odd n=5 forces symmetric) applies to the
(A,A) sector, but (A,A) has no chiral matter (Z₃ phase = 0).
The chiral matter comes from (A,B), where the sign is determined by the
SO(32) parent's antisymmetric form, giving 10.

This is the standard SU(5) GUT spectrum: 10 + 5̄, not the sBootstrap's 15 + 5̄.
