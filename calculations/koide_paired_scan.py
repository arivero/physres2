"""
Koide Q scan structured by quark pair decomposition.

The 6 quarks pair as (t,s) × (u,b) × (c,d).
Each Koide triple picks one from each pair.
Extra options: u → 0 (massless limit), t → ∞ (decoupled limit).
This gives 3 × 3 × 2 = 18 base mass triples.

For each: signs (4 combos) × seesaw (direct/inverse per entry,
constrained by 0 and ∞). Dimensional validity: all finite nonzero
entries must be same type (all-direct or all-inverse).

Output:
  1) Global top 12
  2) Best 2 per each of the 8 original (t/s × u/b × c/d) triples
"""

import itertools
import math

# PDG 2024 MSbar masses in MeV
mass_vals = {
    'u': 2.16, 'd': 4.67, 's': 93.4, 'c': 1270.0, 'b': 4183.0, 't': 162500.0,
    '0': 0.0,
}

pair1 = ['t', 's', 'inf']   # (t,s) pair + t→∞
pair2 = ['u', 'b', '0']     # (u,b) pair + u→0
pair3 = ['c', 'd']           # (c,d) pair

TARGET = 2.0 / 3.0
results = []

for q1 in pair1:
    for q2 in pair2:
        for q3 in pair3:
            # Which of the 8 original triples does this belong to?
            orig_q1 = 't' if q1 in ('t', 'inf') else 's'
            orig_q2 = 'u' if q2 in ('u', '0') else 'b'
            orig_q3 = q3
            orig_key = f'{orig_q1}{orig_q2}{orig_q3}'

            m1 = float('inf') if q1 == 'inf' else mass_vals[q1]
            m2 = mass_vals[q2]
            m3 = mass_vals[q3]

            for inv1, inv2, inv3 in itertools.product([False, True], repeat=3):
                # Constraints: m=0 can't be inverted; m=inf can't be direct
                if m1 == 0 and inv1: continue
                if m2 == 0 and inv2: continue
                if m3 == 0 and inv3: continue
                if m1 == float('inf') and not inv1: continue
                if m2 == float('inf') and not inv2: continue
                if m3 == float('inf') and not inv3: continue

                def get_y(m, inv):
                    if m == 0: return 0.0
                    if m == float('inf'): return 0.0  # 1/sqrt(inf) = 0
                    return 1.0 / math.sqrt(m) if inv else math.sqrt(m)

                r1, r2, r3 = get_y(m1, inv1), get_y(m2, inv2), get_y(m3, inv3)

                # Dimensional validity: all finite nonzero entries same type
                finite_invs = []
                for m, inv in [(m1, inv1), (m2, inv2), (m3, inv3)]:
                    if m not in (0, float('inf')):
                        finite_invs.append(inv)
                dim_ok = len(set(finite_invs)) <= 1

                def mk_label(q, inv):
                    if q == 'inf': return '(t->inf)'
                    if q == '0': return '(u->0)'
                    return f'1/{q}' if inv else q

                l1, l2, l3 = mk_label(q1, inv1), mk_label(q2, inv2), mk_label(q3, inv3)

                for signs in [(1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1)]:
                    ya, yb, yc = signs[0]*r1, signs[1]*r2, signs[2]*r3
                    denom = ya + yb + yc
                    if abs(denom) < 1e-30:
                        continue
                    Q = (ya**2 + yb**2 + yc**2) / denom**2
                    if Q < 0.3 or Q > 5.0:
                        continue

                    def sl(s, lbl):
                        return f'-{lbl}' if s < 0 else lbl

                    trip = f'({sl(signs[0],l1)}, {sl(signs[1],l2)}, {sl(signs[2],l3)})'
                    n_neg = sum(1 for s in signs if s < 0)

                    results.append({
                        'label': trip,
                        'Q': Q,
                        'dev': abs(Q - TARGET),
                        'orig': orig_key,
                        'dim_ok': dim_ok,
                        'n_neg': n_neg,
                        'has_limit': q1 == 'inf' or q2 == '0',
                    })

results.sort(key=lambda r: r['dev'])

# ── Output ──
lines = []
lines.append('Koide Q scan — paired structure (t,s)x(u,b)x(c,d)')
lines.append('With u->0 and t->inf limits')
lines.append('Quark masses (PDG 2024 MSbar, MeV): u=2.16, d=4.67, s=93.4, c=1270, b=4183, t=162500')
lines.append(f'Target Q = 2/3 = {TARGET:.6f}')
lines.append(f'Total cases: {len(results)}')
lines.append('')
lines.append('Dimensional validity: "ok" if all finite nonzero entries are same type')
lines.append('                      "*"  if mixed direct/inverse (dimensionally ill-defined)')
lines.append('')

# 1) Global top 12
lines.append('=' * 90)
lines.append('GLOBAL TOP 12')
lines.append('=' * 90)
hdr = f'{"Rk":<4} {"Triplet":<45} {"Q":>9} {"|Q-2/3|":>10} {"orig":>5} {"dim":>4} {"lim":>4}'
lines.append(hdr)
lines.append('-' * 90)
for i, r in enumerate(results[:12], 1):
    d = 'ok' if r['dim_ok'] else ' *'
    lim = 'lim' if r['has_limit'] else ''
    lines.append(f'{i:<4} {r["label"]:<45} {r["Q"]:>9.6f} {r["dev"]:>10.6f} {r["orig"]:>5} {d:>4} {lim:>4}')

# 2) Best 2 per original triple
lines.append('')
lines.append('=' * 90)
lines.append('BEST 2 PER ORIGINAL TRIPLE (from 8 base triples: t/s x u/b x c/d)')
lines.append('=' * 90)

orig_keys = ['tuc', 'tud', 'tbc', 'tbd', 'suc', 'sud', 'sbc', 'sbd']
for key in orig_keys:
    subset = [r for r in results if r['orig'] == key]
    lines.append(f'\n--- {key.upper()} ({key[0]}/{key[0]}*) x ({key[1]}/{key[1]}*) x ({key[2]}) ---')
    lines.append(f'  {"Triplet":<45} {"Q":>9} {"|Q-2/3|":>10} {"dim":>4} {"lim":>4}')
    for r in subset[:2]:
        d = 'ok' if r['dim_ok'] else ' *'
        lim = 'lim' if r['has_limit'] else ''
        lines.append(f'  {r["label"]:<45} {r["Q"]:>9.6f} {r["dev"]:>10.6f} {d:>4} {lim:>4}')

# 3) Dimensionally valid only, full list
lines.append('')
lines.append('=' * 90)
lines.append('ALL DIMENSIONALLY VALID CASES, RANKED')
lines.append('=' * 90)
valid = [r for r in results if r['dim_ok']]
lines.append(f'{"Rk":<4} {"Triplet":<45} {"Q":>9} {"|Q-2/3|":>10} {"orig":>5} {"lim":>4}')
lines.append('-' * 80)
for i, r in enumerate(valid[:40], 1):
    lim = 'lim' if r['has_limit'] else ''
    lines.append(f'{i:<4} {r["label"]:<45} {r["Q"]:>9.6f} {r["dev"]:>10.6f} {r["orig"]:>5} {lim:>4}')

output = '\n'.join(lines)
print(output)

with open('calculations/koide_paired_scan.txt', 'w') as f:
    f.write(output)

print(f'\nSaved to calculations/koide_paired_scan.txt')
