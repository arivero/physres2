# Assembly: Composite Scalar Counting and Generation Uniqueness

## Section Draft

We consider a model in which the scalar sector is composed of bound states formed from r constituents of down-type charge (−1/3) and s constituents of up-type charge (+2/3), distributed across N generations. The requirement that the composite scalars reproduce the quantum numbers needed for Yukawa couplings to the known fermions imposes a system of Diophantine constraints on the triple (N, r, s).

Specifically, charge +1/3 composites formed from one up-type and one down-type constituent must account for 2N real degrees of freedom, giving

  rs = 2N.                              (1)

Charge −2/3 composites, formed as symmetric pairs of down-type constituents, yield

  r(r+1)/2 = 2N.                        (2)

Finally, neutral composites — including those providing a right-handed neutrino — require

  r² + s² − 1 = 4N.                     (3)

**Theorem.** The system (1)–(3) admits a unique solution over the positive integers: (N, r, s) = (3, 3, 2).

*Proof.* Eliminating N between (1) and (2) gives s = (r+1)/2. Substituting into (3) yields (r−3)(r+1) = 0. Since r > 0, the unique solution is r = 3, s = 2, N = 3. □

The five constituents r + s = 5 naturally furnish a fundamental of SU(5). The charge census for 24 ⊕ 15 ⊕ 15̄ (54 states) with Q = αY₁ + βT₃ restricted to {0, ±1/3, ±2/3, ±1, ±4/3} gives:

| |Q| | 0 | 1/3 | 2/3 | 1 | 4/3 |
|-----|---|-----|-----|---|-----|
| mult | 16 | 2×4 | 2×10 | 2×4 | 2×1 |

The ±4/3 states live exclusively in (1,3)₁ ⊂ 15 and (1,3)₋₁ ⊂ 15̄.

---

## GAPS IDENTIFIED (for Round 2)

1. **Physical derivation of Diophantine equations**: Why rs (not 2rs)? Why symmetric pairing in eq 2? Where does the "−1" in eq 3 come from?

2. **Factor-of-two conventions**: 2N counts real d.o.f. but SU(5) reps count complex d.o.f. Bridge needed.

3. **Explicit mapping**: Which SU(5) components correspond to which equation? Does (6,1)₋₂/₃ or (3̄,1)₋₂/₃ satisfy eq 2?

4. **Why 24 ⊕ 15 ⊕ 15̄?**: Why not 24 ⊕ 10 ⊕ 10̄? Is this forced by the symmetric product (eq 2)?

5. **Asymptotic freedom bound**: n_f = 2(r+s)−1 < 33 not derived. What gauge group? Standard QCD gives n_f < 16.5.

6. **Second solution (7,4,28)**: n_f = 11, would give SU(11). Physical significance?

7. **Charge operator normalization**: (α,β) = (±1, ±1/3) vs standard Gell-Mann–Nishijima Q = Y/2 + T₃. Convention clarification needed.

8. **Charge census verification**: 16 neutral states — how many from 24 vs from 15 ⊕ 15̄?

9. **Anomaly cancellation**: Does 24 ⊕ 15 ⊕ 15̄ introduce gauge anomalies?

10. **Neutral state discrepancy**: Eq 3 gives 12 neutrals, but census gives 16. The extra 4 presumably from the 24 adjoint (gauge/Cartan). Needs explicit accounting.
