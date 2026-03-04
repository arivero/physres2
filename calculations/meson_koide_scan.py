#!/usr/bin/env python3
"""
Systematic scan of Koide ratio Q = 2/3 across ALL pseudoscalar meson triples
with ALL sign assignments.

Q = (m1 + m2 + m3) / (s1*sqrt(m1) + s2*sqrt(m2) + s3*sqrt(m3))^2

PDG 2024 pseudoscalar meson masses in MeV.
"""

import itertools
import math
import random
import os

# PDG 2024 pseudoscalar meson masses (MeV)
MESONS = {
    "pi":    139.570,
    "K":     493.677,
    "eta":   547.862,
    "eta'":  957.78,
    "D":    1869.66,
    "D_s":  1968.35,
    "B":    5279.34,
    "B_s":  5366.92,
    "B_c":  6274.47,
    "eta_c": 2983.9,
    "eta_b": 9398.7,
}

TARGET_Q = 2.0 / 3.0

def koide_Q(m1, m2, m3, s1, s2, s3):
    """Compute Koide ratio: Q = (m1+m2+m3) / (s1*sqrt(m1)+s2*sqrt(m2)+s3*sqrt(m3))^2"""
    numerator = m1 + m2 + m3
    denominator = (s1 * math.sqrt(m1) + s2 * math.sqrt(m2) + s3 * math.sqrt(m3)) ** 2
    if denominator < 1e-30:
        return float('inf')
    return numerator / denominator

def sign_label(s):
    return "+" if s > 0 else "-"

def triple_label(s1, s2, s3, n1, n2, n3):
    return f"({sign_label(s1)}{n1}, {sign_label(s2)}{n2}, {sign_label(s3)}{n3})"

def main():
    names = list(MESONS.keys())
    masses_list = [MESONS[n] for n in names]
    n = len(names)

    # Sign assignments: fix s1=+1 (since global flip gives same Q), vary s2,s3
    # This gives 4 independent combos per triple
    sign_combos = [(+1, +1, +1), (+1, +1, -1), (+1, -1, +1), (+1, -1, -1)]

    results = []
    triple_count = 0

    for combo in itertools.combinations(range(n), 3):
        triple_count += 1
        i, j, k = combo
        n1, n2, n3 = names[i], names[j], names[k]
        m1, m2, m3 = masses_list[i], masses_list[j], masses_list[k]

        for s1, s2, s3 in sign_combos:
            Q = koide_Q(m1, m2, m3, s1, s2, s3)
            dev = abs(Q - TARGET_Q)
            label = triple_label(s1, s2, s3, n1, n2, n3)
            results.append({
                'label': label,
                'names': (n1, n2, n3),
                'masses': (m1, m2, m3),
                'signs': (s1, s2, s3),
                'Q': Q,
                'dev': dev,
                'dev_pct': dev / TARGET_Q * 100,
            })

    results.sort(key=lambda r: r['dev'])
    total = len(results)

    # --- Threshold counts ---
    count_001 = sum(1 for r in results if r['dev'] < 0.001)
    count_0007 = sum(1 for r in results if r['dev'] < 0.0007)
    count_005 = sum(1 for r in results if r['dev'] < 0.005)

    # --- Find (-pi, D_s, B) triple ---
    # With our convention s1=+1, the triple (pi, D_s, B) has pi first.
    # (-pi, +D_s, +B) is equivalent to (+pi, -D_s, -B) by global flip.
    # So we look for names=(pi, D_s, B) with signs=(+1, -1, -1).
    target_triple = None
    target_rank = None
    for idx, r in enumerate(results):
        if r['names'] == ("pi", "D_s", "B") and r['signs'] == (1, -1, -1):
            target_triple = r
            target_rank = idx + 1
            break

    # All sign combos for (pi, D_s, B)
    pi_ds_b_all = [r for r in results if r['names'] == ("pi", "D_s", "B")]

    # --- Monte Carlo ---
    random.seed(42)
    mc_trials = 10000
    mc_hits_001 = 0
    mc_hits_0007 = 0

    for _ in range(mc_trials):
        rm1 = random.uniform(100, 10000)
        rm2 = random.uniform(100, 10000)
        rm3 = random.uniform(100, 10000)
        best_dev = float('inf')
        for s1, s2, s3 in sign_combos:
            Q = koide_Q(rm1, rm2, rm3, s1, s2, s3)
            dev = abs(Q - TARGET_Q)
            if dev < best_dev:
                best_dev = dev
        if best_dev < 0.001:
            mc_hits_001 += 1
        if best_dev < 0.0007:
            mc_hits_0007 += 1

    mc_rate_001 = mc_hits_001 / mc_trials
    mc_rate_0007 = mc_hits_0007 / mc_trials

    # Expected counts
    expected_001 = mc_rate_001 * triple_count
    expected_0007 = mc_rate_0007 * triple_count

    # --- Console output ---
    print(f"Total mesons: {n}")
    print(f"Total triples C({n},3): {triple_count}")
    print(f"Sign assignments per triple (independent): 4")
    print(f"Total configurations scanned: {total}")
    print()

    print("--- TOP 20 ---")
    for idx, r in enumerate(results[:20]):
        print(f"  {idx+1:3d}. {r['label']:35s}  Q={r['Q']:.8f}  |dQ|={r['dev']:.8f}  ({r['dev_pct']:.4f}%)")
    print()

    print("--- THRESHOLD COUNTS ---")
    print(f"  |Q - 2/3| < 0.005:  {count_005} configs")
    print(f"  |Q - 2/3| < 0.001:  {count_001} configs")
    print(f"  |Q - 2/3| < 0.0007: {count_0007} configs")
    print()

    print("--- (-pi, D_s, B) TRIPLE ---")
    if target_triple:
        print(f"  Label: {target_triple['label']}")
        print(f"  Q = {target_triple['Q']:.10f}")
        print(f"  |Q - 2/3| = {target_triple['dev']:.10f}")
        print(f"  Deviation: {target_triple['dev_pct']:.4f}%")
        print(f"  Rank: {target_rank} out of {total}")
    else:
        print("  NOT FOUND - listing all (pi, D_s, B) combos:")
    print()
    print("  All (pi, D_s, B) sign combos:")
    for r in pi_ds_b_all:
        rank = results.index(r) + 1
        print(f"    {r['label']:35s}  Q={r['Q']:.8f}  |dQ|={r['dev']:.8f}  rank={rank}")
    print()

    print("--- MONTE CARLO ---")
    print(f"  Trials: {mc_trials}")
    print(f"  Background rate (best of 4 signs) |dQ|<0.001: {mc_hits_001}/{mc_trials} = {mc_rate_001*100:.2f}%")
    print(f"  Background rate (best of 4 signs) |dQ|<0.0007: {mc_hits_0007}/{mc_trials} = {mc_rate_0007*100:.2f}%")
    print(f"  Expected triples |dQ|<0.001: {expected_001:.1f}, observed: {count_001}")
    print(f"  Expected triples |dQ|<0.0007: {expected_0007:.1f}, observed: {count_0007}")
    print()

    # --- Build markdown report ---
    L = []
    L.append("# Koide Ratio Scan: All Pseudoscalar Meson Triples")
    L.append("")
    L.append("## Definition")
    L.append("")
    L.append("Koide ratio with signed square roots:")
    L.append("")
    L.append("    Q = (m1 + m2 + m3) / (s1*sqrt(m1) + s2*sqrt(m2) + s3*sqrt(m3))^2")
    L.append("")
    L.append("where si = +/-1. Target: Q = 2/3.")
    L.append("")
    L.append("Since Q is invariant under simultaneous flip of all signs, there are")
    L.append("4 independent sign assignments per triple (not 8).")
    L.append("")

    L.append("## Pseudoscalar Meson Masses (PDG 2024)")
    L.append("")
    L.append("| Meson | Mass (MeV) |")
    L.append("|-------|------------|")
    for name in names:
        L.append(f"| {name} | {MESONS[name]} |")
    L.append("")

    L.append("## Combinatorics")
    L.append("")
    L.append(f"- Mesons: {n}")
    L.append(f"- Triples: C({n},3) = {triple_count}")
    L.append(f"- Independent sign assignments per triple: 4")
    L.append(f"- Total configurations: {total}")
    L.append("")

    # Top 20
    L.append("## Top 20 Closest to Q = 2/3")
    L.append("")
    L.append("| Rank | Triple | Q | |Q - 2/3| | Deviation % |")
    L.append("|------|--------|---|----------|-------------|")
    for idx, r in enumerate(results[:20]):
        L.append(f"| {idx+1} | {r['label']} | {r['Q']:.8f} | {r['dev']:.8f} | {r['dev_pct']:.4f}% |")
    L.append("")

    # Extended top 50
    L.append("## Extended: Top 50")
    L.append("")
    L.append("| Rank | Triple | Q | |Q - 2/3| | Dev % |")
    L.append("|------|--------|---|----------|-------|")
    for idx, r in enumerate(results[:50]):
        L.append(f"| {idx+1} | {r['label']} | {r['Q']:.8f} | {r['dev']:.8f} | {r['dev_pct']:.4f}% |")
    L.append("")

    # Threshold counts
    L.append("## Threshold Counts")
    L.append("")
    L.append(f"| Threshold | Count | Fraction of {total} |")
    L.append("|-----------|-------|---------------------|")
    L.append(f"| |Q - 2/3| < 0.005 (0.75%) | {count_005} | {count_005/total:.4f} |")
    L.append(f"| |Q - 2/3| < 0.001 (0.15%) | {count_001} | {count_001/total:.4f} |")
    L.append(f"| |Q - 2/3| < 0.0007 (0.1%) | {count_0007} | {count_0007/total:.4f} |")
    L.append("")

    # (-pi, D_s, B) section
    L.append("## The (-pi, D_s, B) Triple")
    L.append("")
    if target_triple:
        L.append(f"- **Signed triple**: {target_triple['label']}")
        L.append(f"- **Q value**: {target_triple['Q']:.10f}")
        L.append(f"- **|Q - 2/3|**: {target_triple['dev']:.10f}")
        L.append(f"- **Deviation**: {target_triple['dev_pct']:.4f}%")
        L.append(f"- **Rank**: {target_rank} out of {total} total configurations")
        L.append(f"- **Rank among triples**: top {target_rank} (counting each sign assignment separately)")
    else:
        L.append("Not found with expected sign convention. See all combos below.")
    L.append("")

    L.append("### All sign assignments for (pi, D_s, B)")
    L.append("")
    L.append("| Signs | Q | |Q - 2/3| | Dev % | Rank |")
    L.append("|-------|---|----------|-------|------|")
    for r in sorted(pi_ds_b_all, key=lambda x: x['dev']):
        rank = results.index(r) + 1
        L.append(f"| {r['label']} | {r['Q']:.8f} | {r['dev']:.8f} | {r['dev_pct']:.4f}% | {rank} |")
    L.append("")

    # Monte Carlo
    L.append("## Monte Carlo Background Estimate")
    L.append("")
    L.append(f"- **Trials**: {mc_trials} random triples")
    L.append(f"- **Mass range**: uniform in [100, 10000] MeV")
    L.append(f"- **Method**: for each random triple, test all 4 sign assignments, keep best")
    L.append("")
    L.append("### Results")
    L.append("")
    L.append("| Threshold | MC hits | MC rate | Expected (165 triples) | Observed |")
    L.append("|-----------|---------|---------|------------------------|----------|")
    L.append(f"| |Q-2/3| < 0.001 | {mc_hits_001} | {mc_rate_001*100:.2f}% | {expected_001:.1f} | {count_001} |")
    L.append(f"| |Q-2/3| < 0.0007 | {mc_hits_0007} | {mc_rate_0007*100:.2f}% | {expected_0007:.1f} | {count_0007} |")
    L.append("")

    if expected_001 > 0:
        ratio_001 = count_001 / expected_001
        L.append(f"Observed/Expected ratio at 0.001 threshold: {ratio_001:.2f}")
    if expected_0007 > 0:
        ratio_0007 = count_0007 / expected_0007
        L.append(f"Observed/Expected ratio at 0.0007 threshold: {ratio_0007:.2f}")
    elif expected_0007 == 0 and count_0007 > 0:
        L.append(f"Observed/Expected ratio at 0.0007 threshold: {count_0007}/0 (no MC hits at this threshold)")
    L.append("")

    # All hits within 0.005 for reference
    hits_005 = [r for r in results if r['dev'] < 0.005]
    L.append("## Complete List: All Configurations with |Q - 2/3| < 0.005")
    L.append("")
    L.append("| Rank | Triple | Q | |Q - 2/3| | Dev % |")
    L.append("|------|--------|---|----------|-------|")
    for r in hits_005:
        rank = results.index(r) + 1
        L.append(f"| {rank} | {r['label']} | {r['Q']:.8f} | {r['dev']:.8f} | {r['dev_pct']:.4f}% |")
    L.append("")

    # Write
    outpath = "/home/codexssh/phys3/results/meson_koide_scan.md"
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, 'w') as f:
        f.write('\n'.join(L) + '\n')

    print(f"Report written to {outpath}")

if __name__ == "__main__":
    main()
