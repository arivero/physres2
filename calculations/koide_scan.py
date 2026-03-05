"""
Exhaustive Koide Q scan over all quark triplets with sign and inverse variations.

Koide Q = sum(y_k^2) / (sum(y_k))^2,  target Q = 2/3

For each triplet (a,b,c) from {u,d,s,c,b,t}:
  - each quark can be normal (y = sqrt(m)) or inverse (y = 1/sqrt(m))
  - each y_k can have sign +1 or -1
  Note: (s1,s2,s3) and (-s1,-s2,-s3) give same Q (overall sign flip)
  → 4 distinct sign combos per mass-type choice
  → 2^3 mass-type choices × 4 sign combos × C(6,3) triplets = 640 total cases
"""

import itertools
import math

# PDG 2024 MSbar masses in MeV
quarks = {
    'u': 2.16,
    'd': 4.67,
    's': 93.4,
    'c': 1270.0,
    'b': 4183.0,
    't': 162500.0,   # MS-bar at m_t scale
}
names = list(quarks.keys())
masses = [quarks[n] for n in names]

TARGET = 2.0/3.0

results = []

for idx in itertools.combinations(range(6), 3):
    qa, qb, qc = names[idx[0]], names[idx[1]], names[idx[2]]
    ma, mb, mc = masses[idx[0]], masses[idx[1]], masses[idx[2]]

    # iterate over inverse choices: 0=normal (sqrt(m)), 1=inverse (1/sqrt(m))
    for inv in itertools.product([0, 1], repeat=3):
        ra = math.sqrt(ma) if inv[0]==0 else 1.0/math.sqrt(ma)
        rb = math.sqrt(mb) if inv[1]==0 else 1.0/math.sqrt(mb)
        rc = math.sqrt(mc) if inv[2]==0 else 1.0/math.sqrt(mc)

        # label for mass type
        la = ('1/'+qa if inv[0] else qa)
        lb = ('1/'+qb if inv[1] else qb)
        lc = ('1/'+qc if inv[2] else qc)

        # 4 distinct sign combos (fix first sign = +1 to remove degeneracy)
        for signs in [(+1,+1,+1),(+1,+1,-1),(+1,-1,+1),(+1,-1,-1)]:
            sa, sb, sc = signs
            ya = sa * ra
            yb = sb * rb
            yc = sc * rc

            denom = ya + yb + yc
            if abs(denom) < 1e-30:
                continue

            numer = ya**2 + yb**2 + yc**2
            Q = numer / denom**2

            # Q must be >= 1/3 by Cauchy-Schwarz (when denom != 0)
            # Negative Q can occur if denom^2 is computed wrong -- shouldn't happen
            # but skip unphysical cases
            if Q < 0.3 or Q > 10.0:
                continue

            # label signs
            def slabel(s, lbl):
                return ('-' + lbl) if s < 0 else lbl

            label = f"({slabel(sa,la)}, {slabel(sb,lb)}, {slabel(sc,lc)})"

            # count how many are inverted (for categorization)
            n_inv = sum(inv)
            # count how many are negative
            n_neg = sum(1 for s in signs if s < 0)

            # Dimensional validity: all-direct (n_inv=0) or all-inverse (n_inv=3)
            # are dimensionally consistent (all MeV^{+1/2} or all MeV^{-1/2}).
            # Mixed (n_inv=1 or 2) adds MeV^{+1/2} and MeV^{-1/2} -- ill-defined.
            dim_ok = (n_inv == 0 or n_inv == 3)

            deviation = abs(Q - TARGET)
            results.append({
                'label': label,
                'Q': Q,
                'dev': deviation,
                'n_inv': n_inv,
                'n_neg': n_neg,
                'dim_ok': dim_ok,
            })

# Sort by closeness to 2/3
results.sort(key=lambda r: r['dev'])

# Write output
lines = []
lines.append("Koide Q scan — all quark triplets with sign and inverse variations")
lines.append("Target Q = 2/3 = 0.666667")
lines.append("Quark masses (PDG 2024 MSbar, MeV): u=2.16, d=4.67, s=93.4, c=1270, b=4183, t=162500")
n_dim_ok = sum(1 for r in results if r['dim_ok'])
n_dim_bad = sum(1 for r in results if not r['dim_ok'])
lines.append("")
lines.append(f"Total cases enumerated: {len(results)}")
lines.append(f"  Dimensionally valid   (n_inv=0 or 3, all MeV^{{+/-1/2}}): {n_dim_ok}")
lines.append(f"  Dimensionally INVALID (n_inv=1 or 2, mixes MeV^{{+1/2}} and MeV^{{-1/2}}): {n_dim_bad} [marked *]")
lines.append("")
lines.append("NOTE: Mixed triples (n_inv=1 or 2) sum sqrt(m) and 1/sqrt(m) with")
lines.append("      incompatible dimensions -- Q has no physical meaning for these.")
lines.append("      UNIT DEPENDENCE TEST: Q(s,c,1/b) = 0.6639 in MeV but 0.4340 in GeV.")
lines.append("      The 'near-2/3' values for mixed triples are a MeV-unit artifact.")
lines.append("      All-direct and all-inverse Q values are unit-independent (scale-invariant).")
lines.append("")
lines.append("Legend:")
lines.append("  1/x  = inverse mass used for quark x")
lines.append("  -x   = negative sqrt (sign flip in Koide root)")
lines.append("  n_inv = number of inverted quarks, n_neg = number with negative sqrt")
lines.append("  *    = dimensionally invalid (mixed direct and inverse)")
lines.append("")
lines.append(f"{'Rank':<5} {'Triplet':<32} {'Q':>10} {'|Q-2/3|':>12} {'n_inv':>6} {'n_neg':>6} {'dim':>5}")
lines.append("-" * 85)

for i, r in enumerate(results[:200], 1):
    dim_flag = "  ok" if r['dim_ok'] else "  * "
    lines.append(f"{i:<5} {r['label']:<32} {r['Q']:>10.6f} {r['dev']:>12.6f} {r['n_inv']:>6} {r['n_neg']:>6} {dim_flag}")

# Also write a summary of the best hits grouped by category
lines.append("")
lines.append("=" * 80)
lines.append("TOP HITS BY CATEGORY (best per category, |Q-2/3| < 0.05)")
lines.append("=" * 80)

cats = {}
for r in results:
    key = (r['n_inv'], r['n_neg'])
    if key not in cats:
        cats[key] = []
    cats[key].append(r)

for key in sorted(cats.keys()):
    n_inv, n_neg = key
    best = [r for r in cats[key] if r['dev'] < 0.05]
    if not best:
        continue
    inv_str = "all-normal" if n_inv==0 else f"{n_inv}-inverse"
    neg_str = "all-positive" if n_neg==0 else f"{n_neg}-negative"
    lines.append(f"\n--- n_inv={n_inv} ({inv_str}), n_neg={n_neg} ({neg_str}) ---")
    lines.append(f"  {'Triplet':<32} {'Q':>10} {'|Q-2/3|':>12}")
    for r in best[:20]:
        lines.append(f"  {r['label']:<32} {r['Q']:>10.6f} {r['dev']:>12.6f}")

output = "\n".join(lines)
print(output)

with open("koide_scan_results.txt", "w") as f:
    f.write(output)

lines.append("")
lines.append("=" * 85)
lines.append("DIMENSIONALLY VALID CASES ONLY (n_inv=0 all-direct, or n_inv=3 all-inverse)")
lines.append("Ranked by |Q-2/3|")
lines.append("=" * 85)
valid = [r for r in results if r['dim_ok']]
lines.append(f"\n{'Rank':<5} {'Triplet':<32} {'Q':>10} {'|Q-2/3|':>12} {'type':>12}")
lines.append("-" * 75)
for i, r in enumerate(valid[:50], 1):
    t = "all-direct" if r['n_inv']==0 else "all-inverse"
    lines.append(f"{i:<5} {r['label']:<32} {r['Q']:>10.6f} {r['dev']:>12.6f} {t:>12}")

output = "\n".join(lines)

with open("koide_scan_results.txt", "w") as f:
    f.write(output)

print(f"Saved to koide_scan_results.txt")
print(f"\nTop 10 (all cases, * = dimensionally invalid):")
for r in results[:10]:
    flag = "" if r['dim_ok'] else " *"
    print(f"  {r['label']:<32}  Q={r['Q']:.8f}  dev={r['dev']:.2e}{flag}")
print(f"\nTop 10 dimensionally valid only:")
for r in valid[:10]:
    t = "all-direct" if r['n_inv']==0 else "all-inverse"
    print(f"  {r['label']:<32}  Q={r['Q']:.8f}  dev={r['dev']:.2e}  ({t})")
