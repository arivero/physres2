#!/usr/bin/env python3
"""
Orientifold projection: antisymmetric 10 vs symmetric 15

The Z_3 orbifold on T^6 with Type I (orientifold) gives chiral matter
from the D9-brane sector. The representation (symmetric or antisymmetric
of SU(5)_f) depends on the orientifold projection type.

This script examines:
1. Why the standard construction gives antisymmetric 10
2. What modifications could give symmetric 15
3. Tadpole cancellation constraints
"""
import numpy as np
from itertools import product

print("=" * 70)
print("ORIENTIFOLD PROJECTION: 10 vs 15 OF SU(5)")
print("=" * 70)

# =====================================================================
# 1. Standard Z_3 orientifold: tadpole constraints
# =====================================================================
print("\n--- 1. Tadpole cancellation for Z_3 on T^6 ---")
print()
print("D9-brane sector (no D5-branes):")
print("  n_0 + 2*n_1 = 16  (untwisted RR tadpole: 32 CP labels / 2 for orientifold)")
print("  n_0 - n_1 = 1     (twisted RR tadpole: Tr(gamma_theta) = 1)")
print()

# Solve
# n_0 + 2*n_1 = 16
# n_0 - n_1 = 1  =>  n_0 = 1 + n_1
# => 1 + 3*n_1 = 16 => n_1 = 5
n1 = 5
n0 = 1 + n1
print(f"  Unique solution: n_0 = {n0}, n_1 = n_2 = {n1}")
print(f"  Total: {n0} + 2*{n1} = {n0 + 2*n1} = 16  ✓")
print(f"  Tr(gamma) = {n0} - {n1} = {n0 - n1} = 1  ✓")

# =====================================================================
# 2. Orientifold action on matter representations
# =====================================================================
print("\n--- 2. Orientifold projection on SU(5) matter ---")
print()
print("The n_1 = 5 branes with eigenvalue ω carry an SU(5) gauge factor.")
print("The orientifold Ω maps ω-branes ↔ ω²-branes (conjugate pair).")
print()
print("For strings with BOTH ends on the ω-block (before Ω projection):")
print("  5 ⊗ 5 = 10_A + 15_S  (antisymmetric + symmetric)")
print()
print("The Ω projection keeps one of these, determined by the sign ε:")
print("  ε = +1 (O⁺, SO-type): keeps ANTISYMMETRIC → 10")
print("  ε = -1 (O⁻, Sp-type): keeps SYMMETRIC → 15")
print()
print("Which sign is forced?")
print()
print("For O⁺ planes (standard Type I): the projection matrix γ_Ω")
print("  on n-dimensional blocks satisfies γ_Ω^T = +γ_Ω (symmetric).")
print("  This works for ANY n.")
print()
print("For O⁻ planes (Sp-type): γ_Ω^T = -γ_Ω (antisymmetric).")
print("  An antisymmetric n×n matrix has rank ≤ n-1 for odd n,")
print("  so it is SINGULAR for odd n. Sp(n) requires even n.")
print()
print(f"The SU(5) block has n_1 = {n1} (ODD).")
print("  → O⁺ is the ONLY consistent choice → antisymmetric 10.")
print("  → Sp-type projection is ALGEBRAICALLY IMPOSSIBLE for n=5.")

# =====================================================================
# 3. Can D5-branes change the game?
# =====================================================================
print("\n--- 3. Modified tadpole with D5-branes ---")
print()
print("With N_5 D5-branes at Z_3 fixed points, the twisted tadpole becomes:")
print("  Tr(γ_{θ,9}) + Tr(γ_{θ,5}) = c  (c depends on O-plane charge)")
print()
print("If the D5 contribution shifts the RHS, n_1 could change.")
print("We need n_1 EVEN for Sp-type. Scanning possible solutions:")
print()

# With D5 branes, twisted tadpole: n0_9 - n1_9 + n0_5 - n1_5 = c
# Untwisted: n0_9 + 2*n1_9 = 16 (D9 count fixed by anomaly)
# D5 branes don't affect D9 untwisted tadpole
# But twisted tadpole changes to: n0_9 - n1_9 = c - (n0_5 - n1_5)

for c_twist in [0, 1, 2, 4]:
    print(f"  Twisted tadpole RHS = {c_twist}:")
    # n0 + 2*n1 = 16, n0 - n1 = c_twist => n0 = c_twist + n1
    # c_twist + 3*n1 = 16 => n1 = (16 - c_twist)/3
    n1_try = (16 - c_twist) / 3
    if n1_try == int(n1_try) and n1_try > 0:
        n1_try = int(n1_try)
        n0_try = c_twist + n1_try
        parity = "EVEN" if n1_try % 2 == 0 else "ODD"
        sp_ok = "✓ Sp possible!" if n1_try % 2 == 0 else "✗ Sp blocked"
        gauge = f"SU({n1_try})×SU({n0_try//2})" if n0_try % 2 == 0 else f"SU({n1_try})×..."
        print(f"    n_1 = {n1_try} ({parity}), n_0 = {n0_try}  {sp_ok}")
        if n1_try % 2 == 0:
            print(f"    → SU({n1_try}) flavor × (n_0 = {n0_try} color sector)")
            print(f"    → But SU({n1_try}) ≠ SU(5)! Need different flavor group.")
    else:
        print(f"    n_1 = {n1_try:.2f} — not integer, no solution")
    print()

print("Summary: integer solutions with EVEN n_1 require")
print("  c = 4 → n_1 = 4 → SU(4), not SU(5)")
print("  c = -2 → n_1 = 6 → SU(6), not SU(5)")
print("No solution gives SU(5) with even-dimensional blocks.")

# =====================================================================
# 4. Discrete torsion
# =====================================================================
print("\n--- 4. Discrete torsion ---")
print()
print("For Z_N orientifolds, a discrete torsion phase ε_k = ±1")
print("in the k-th twisted sector can flip the projection:")
print("  ε_1 = +1: standard → antisymmetric (10)")
print("  ε_1 = -1: non-standard → symmetric (15)")
print()
print("However, for Z_3:")
print("  - ε_1 and ε_2 are related by ε_2 = ε_1* (consistency)")
print("  - The choice ε_1 = -1 modifies the twisted tadpole condition")
print("  - New twisted tadpole: Tr(γ_θ) = -1 (sign flip)")
print("    → n_0 - n_1 = -1 → n_0 = n_1 - 1")
print("    → n_1 - 1 + 2*n_1 = 16 → 3*n_1 = 17")
print()
n1_dt = 17/3
print(f"  n_1 = 17/3 = {n1_dt:.4f} — NOT INTEGER!")
print()
print("Discrete torsion ε = -1 has NO integer solution for Z_3.")
print("This option is INCONSISTENT.")
print()
print("(For Z_2 × Z_2 orientifolds, discrete torsion works because")
print("the tadpole structure is different. Z_3 is too rigid.)")

# =====================================================================
# 5. Magnetized branes / branes at angles
# =====================================================================
print("\n--- 5. Magnetized branes (intersecting brane models) ---")
print()
print("In magnetized D-brane models, the chiral spectrum comes from")
print("the intersection numbers of brane stacks wrapped on cycles.")
print()
print("For two stacks a, b with intersection number I_{ab}:")
print("  I_{ab} > 0 → |I_{ab}| copies of (□_a, □̄_b)")
print("  I_{ab} < 0 → |I_{ab}| copies of (□̄_a, □_b)")
print()
print("For the SYMMETRIC representation:")
print("  The intersection of stack a with its ORIENTIFOLD IMAGE a'")
print("  gives matter in the two-index representation of U(n_a).")
print("  The type depends on the O-plane:")
print("    I_{aa'} with O⁺: antisymmetric □□_A (10 of SU(5))")
print("    I_{aa'} with O⁻: symmetric □□_S (15 of SU(5))")
print()
print("In magnetized brane models on T^6:")
print("  - Multiple types of O-planes can coexist (O9, O5, etc.)")
print("  - The net O-plane charge determines the projection type")
print("  - With appropriate magnetic fluxes, the EFFECTIVE projection")
print("    can be Sp-type even with overall O⁺ charge")
print()
print("This is the most promising route to the symmetric 15.")
print("The Z_3 orbifold point is a SPECIAL CASE of magnetized branes")
print("(zero flux). Moving away from the orbifold point into the")
print("moduli space of magnetized branes may allow the 15.")

# =====================================================================
# 6. Gimon-Polchinski variant: O5⁻ at fixed points
# =====================================================================
print("\n--- 6. O5⁻ planes at Z_3 fixed points ---")
print()
print("T^6/Z_3 has 27 fixed points (3^3 on each T^2 factor).")
print("At each fixed point, one can place O5⁺ or O5⁻ planes.")
print()
print("O5⁺: contributes +1 to twisted tadpole (per fixed point)")
print("O5⁻: contributes -1 to twisted tadpole")
print()
print("With k O5⁻ planes at fixed points (and 27-k O5⁺):")
print("  Net twisted tadpole shift: (27-k)·(+1) + k·(-1) = 27 - 2k")
print()
print("Modified D9 tadpole: Tr(γ_θ) = 1 - (27-2k)/27 correction")
print("(The exact coefficient depends on the model details.)")
print()
print("If the O5⁻ planes dominate at the SU(5) locus,")
print("the LOCAL projection can be Sp-type → symmetric 15,")
print("even though the GLOBAL projection is SO-type → SO(32).")
print()
print("This is the mechanism proposed by Antoniadis, Bachas, Dudas (1999)")
print("for getting symmetric representations in Type I models.")

# =====================================================================
# 7. Summary of options
# =====================================================================
print("\n" + "=" * 70)
print("SUMMARY: ROUTES FROM 10 TO 15")
print("=" * 70)
print()
print("Option                    | Status         | Problem")
print("-" * 70)
print("Standard Z_3 orientifold  | Gives 10       | n_1=5 (odd) → no Sp")
print("D5-branes shift tadpole   | n_1=4 or 6     | Loses SU(5)")
print("Discrete torsion (ε=-1)   | No solution    | 17/3 not integer")
print("Magnetized branes         | POSSIBLE        | Leaves orbifold point")
print("Local O5⁻ at fixed pts   | POSSIBLE        | Requires explicit model")
print()
print("The two viable options both involve going BEYOND the strict")
print("Z_3 orbifold point in moduli space:")
print()
print("  A) Magnetized branes: replace orbifold by magnetic fluxes.")
print("     The Z_3 structure is recovered as a special flux choice.")
print("     Symmetric reps arise from branes wrapping cycles that")
print("     intersect O⁻ planes. Three generations from I_{aa'} = 3.")
print()
print("  B) Local O5⁻ planes: keep the orbifold but place O5⁻ at")
print("     the fixed points where the SU(5) matter lives. This")
print("     locally flips the projection from antisymmetric to symmetric.")
print("     Global consistency (tadpole cancellation) must be checked.")
print()
print("In both cases, the key insight is that the global Type I")
print("gauge group remains SO(32), but the LOCAL projection at the")
print("orbifold fixed points can differ from the global one.")
print()
print("The sBootstrap's requirement of the symmetric 15 thus")
print("CONSTRAINS the string moduli: it selects a specific region")
print("of the Type I moduli space where O⁻ effects dominate locally.")
