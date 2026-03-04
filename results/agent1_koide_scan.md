# Agent 1: Mass Triplet Formula — Exhaustive Scan Results

## Method

Scanned all C(9,3) = 84 triplets × 8 sign combinations = 672 evaluations.
Self-check: direct formula vs angle formula, max discrepancy = 6.66 × 10⁻¹⁶ (machine epsilon).

## Charged Lepton Triplet (e, μ, τ) — High Precision

| Quantity | Value |
|----------|-------|
| Q | 1.500013849 |
| |Q - 3/2| | 1.384930375 × 10⁻⁵ |
| θ | 44.999735° |
| |θ - 45°| | 2.645 × 10⁻⁴ degrees |

## Ranked Table: All triplets with |Q - 3/2| < 0.05

| Rank | Triplet | Signs | Q | |Q - 3/2| | θ (deg) |
|------|---------|-------|---|----------|---------|
| 1 | (e, μ, τ) | +++ | 1.5000138493 | 1.385e-05 | 44.9997 |
| 2 | (e, μ, τ) | --- | 1.5000138493 | 1.385e-05 | 135.000 |
| 3 | (c, b, t) | +++ | 1.4937748943 | 6.225e-03 | 45.119 |
| 4 | (c, b, t) | --- | 1.4937748943 | 6.225e-03 | 134.881 |
| 5 | (e, μ, c) | -++ | 1.4846663848 | 1.533e-02 | 45.293 |
| 6 | (e, μ, c) | +-- | 1.4846663848 | 1.533e-02 | 134.707 |
| 7 | (τ, u, s) | +++ | 1.5164821100 | 1.648e-02 | 44.685 |
| 8 | (τ, u, s) | --- | 1.5164821100 | 1.648e-02 | 135.315 |
| 9 | (s, c, b) | -++ | 1.4815819582 | 1.842e-02 | 45.352 |
| 10 | (s, c, b) | +-- | 1.4815819582 | 1.842e-02 | 134.648 |
| 11 | (τ, s, c) | +-+ | 1.4779301658 | 2.207e-02 | 45.422 |
| 12 | (τ, s, c) | -+- | 1.4779301658 | 2.207e-02 | 134.578 |
| 13 | (e, τ, s) | +++ | 1.4751211281 | 2.488e-02 | 45.475 |

Note: entries come in ±± pairs (all-positive and all-negative signs give the same Q).

## Summary Statistics

- Total combinations scanned: 672
- |Q - 3/2| < 0.001: 2 (both the lepton triplet)
- |Q - 3/2| < 0.01: 4 (leptons + (c,b,t))
- |Q - 3/2| < 0.05: 26

## Key Finding

The charged lepton triplet (e, μ, τ) is the only triplet within 0.001 of Q = 3/2.
The quark triplet (c, b, t) is second but 450× worse.
The sign-flipped triplet (−√s, √c, √b) appears as rank 9 with |Q - 3/2| ≈ 0.018.
