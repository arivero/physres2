# Z₃ Orbifold Projection on the 15 of SU(5)

## Question
Can the Z₃ orbifold project out the (1,3) component of the 15 while keeping (6,1) and (3,2)?

## Answer: NO — arithmetic obstruction

For γ = diag(ω^p, ω^p, ω^p, ω^q, ω^q), the Z₃ charges are:
- (6,1): ω^{2p}
- (3,2): ω^{p+q}
- (1,3): ω^{2q}

The constraint 2p + 2q ≡ 2(p+q) mod 3 means: if (6,1) and (3,2) survive (charges = 0 mod 3), then (1,3) must also survive.

## Exhaustive table (p,q) mod 3

| (p,q) | (6,1) | (3,2) | (1,3) |
|-------|-------|-------|-------|
| (0,0) |  ✓    |  ✓    |  ✓    |
| (0,1) |  ✓    |  ✗    |  ✗    |
| (0,2) |  ✓    |  ✗    |  ✗    |
| (1,0) |  ✗    |  ✗    |  ✓    |
| (1,1) |  ✗    |  ✗    |  ✗    |
| (1,2) |  ✗    |  ✓    |  ✗    |
| (2,0) |  ✗    |  ✗    |  ✓    |
| (2,1) |  ✗    |  ✓    |  ✗    |
| (2,2) |  ✗    |  ✗    |  ✗    |

## Type I construction
γ_θ = ω·1₅ → p=q=1 → ALL components projected out.
SU(5)→SU(3)×SU(2) breaking must come from Wilson lines or Higgsing, not the orbifold.
