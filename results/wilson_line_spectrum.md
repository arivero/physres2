# Wilson Line Breaking SU(5)_f → SU(3)'×SU(2)×U(1)_Y

## Result: Quantum number mismatch

Breaking SU(5)_f by Wilson line produces states with correct SU(2)×U(1)_Y quantum numbers
BUT wrong SU(3)_c assignments. Only e_R^c matches (it's singlet under both).

## The problem
SM requires correlated color-electroweak quantum numbers (Q_L is simultaneously 3_c AND 2_L).
In SU(5)_f × SU(3)_c, states with right (SU(2),Y) are color singlets (from 10,1),
while color-charged states (from 5̄,3) have wrong Y values.

## Key table
| SM field | Needed (3_c, 2, Y) | Found? | Problem |
|----------|-------------------|--------|---------|
| Q_L | (3, 2, 1/6) | NO | State B is (1_c, 2, 1/6) — right EW, wrong color |
| u_R^c | (3̄, 1, -2/3) | NO | State A is (1_c, 1, -2/3) — color singlet |
| d_R^c | (3̄, 1, 1/3) | NO | State D is (3_c, 1, 1/3) — color triplet not antitriplet |
| L | (1, 2, -1/2) | NO | State E is (3_c, 2, -1/2) — has color |
| e_R^c | (1, 1, 1) | YES | State C matches exactly |

## Resolution needed
Wilson line must act on full SO(32), not just SU(5) factor.
SU(5)_f is gauged in string construction (Chan-Paton), not global.
The SM embedding requires understanding the full SO(32) breaking chain.

## Note
This does NOT invalidate the sBootstrap — it means the SM spectrum identification
requires the full string-theoretic breaking, not a naive factored Wilson line.
