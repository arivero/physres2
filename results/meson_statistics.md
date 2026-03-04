# Meson Koide Triple — Look-Elsewhere Statistics

## Target triple
(-pi+, Ds, B+): Q = 1.49844312, |Q-3/2| = 1.557e-3, deviation = 0.104%

## Scan of 15 common pseudoscalar mesons
With all sign choices (7 per triplet, excluding all-positive), 3185 total trials.

### Top 10 closest to Q = 3/2:
| Rank | |Q-3/2| | Q | Triple |
|------|---------|---|--------|
| 1 | 6.37e-4 | 1.50064 | (+pi0, -Ds, -eta_c) |
| 2 | 6.37e-4 | 1.50064 | (-pi0, +Ds, +eta_c) |
| 3 | 1.56e-3 | 1.49844 | (+pi+, -Ds, -B+) |
| 4 | 1.56e-3 | 1.49844 | (-pi+, +Ds, +B+) |
| 5 | 1.56e-3 | 1.49844 | (+pi+, -Ds, -B0) |
| 6 | 1.56e-3 | 1.49844 | (-pi+, +Ds, +B0) |
| 7 | 2.17e-3 | 1.49783 | (+pi+, -Ds, -Bs) |
| 8 | 2.17e-3 | 1.49783 | (-pi+, +Ds, +Bs) |
| 9 | 2.66e-3 | 1.50266 | (+pi0, -Ds, -J/psi) |
| 10 | 2.66e-3 | 1.50266 | (-pi0, +Ds, +J/psi) |

## Key finding: (pi0, Ds, eta_c) is BETTER
The (pi0, -Ds, -eta_c) triple has |Q-3/2| = 6.37e-4, which is 2.4x closer than (-pi, Ds, B).
This is a new observation not in the paper.

## Look-elsewhere Monte Carlo
50,000 trials with 15 random masses (log-uniform in [100, 10000] MeV), scanning all C(15,3) * 7 = 3185 sign-choice triplets per trial.

Result: 44.7% of random mass sets produce a triplet as good as (-pi, Ds, B).
Significance: ~0.1 sigma. NOT statistically significant after look-elsewhere.

## Interpretation
The 0.10% precision of the (-pi, Ds, B) triple is not statistically surprising given the number of trials (3185 sign-choice triplets from 15 mesons). The physical significance, if any, must come from the SPECIFIC particles involved (pi, Ds, B are the lightest pseudoscalars with those flavor quantum numbers), not from the closeness to 3/2 alone.

The (pi0, Ds, eta_c) triple deserves investigation — it's closer to 3/2 and involves physically meaningful particles (the lightest unflavored, charmed-strange, and charmonium pseudoscalars).
