#!/usr/bin/env python3
"""
Regression checks for calculations/rg_running_qcd.py (issue #3).

Language: Python 3 (numpy, scipy). Run from the repository root:
    python3 calculations/check_rg_running_qcd.py
Exit status 0 when every check passes.

Reference values come from CRunDec 3.1 (Herren-Steinhauser,
arXiv:1703.03751) through the `rundec` Python wrapper 0.7. They are stored
below, so the script runs without CRunDec. When `import rundec` works, the
references are also recomputed live and compared with the stored ones.

CRunDec's AlphasExact integrates the RGE with a fixed Runge-Kutta step
(1e-4 in ln mu) even at one loop, so its values carry a relative error of
about 1e-9; the analytic one-loop solution here is exact. Tolerances:
  * fixed-n_f rows printed in issue #3 (9 decimals): |diff| <= 1.5e-9
  * all CRunDec comparisons: relative 5e-9 per value; 1e-8 on a
    round-trip residual (two legs of CRunDec integration error)
  * identities, round trips, ratio invariance: relative 1e-12
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rg_running_qcd as rg  # noqa: E402

MZ, A_MZ = rg.MZ, rg.ALPHA_S_MZ
FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAIL.append(name)


def rel(a, b):
    return abs(a - b) / abs(b)


# ----------------------------------------------------------------------
# Stored CRunDec 3.1 references (alpha_s, m(mu1)/m(mu0)), alpha_s(mu0) = 0.1179
# Cases: (mu0, nf0, mu1, nf1); thresholds at 1.275, 4.18, 162.5 GeV,
# MS-bar decoupling (DecAsDownMS/UpMS, DecMqDownMS/UpMS) at mu_h = m_h.
# ----------------------------------------------------------------------
ISSUE_ROWS = {  # as printed in issue #3
    10.0: (0.172868568, 1.220995902),
    20.0: (0.150817929, 1.137088712),
}
OLD_ROWS = {  # 'Repository' columns of issue #3 (former code)
    10.0: (0.131177213, 1.057255098),
    20.0: (0.126702692, 1.038283244),
}
CASES = {
    'down_1cross': (MZ, 5, 3.0, 4),
    'up_1cross': (3.0, 4, 10.0, 5),
    'at_b_nf5': (MZ, 5, 4.18, 5),
    'at_b_nf4': (MZ, 5, 4.18, 4),
    'down_2cross': (MZ, 5, 1.0, 3),
    'up_top': (MZ, 5, 500.0, 6),
}
CRUNDEC = {
    ('fixed5', 10.0): (0.17286856783003698, 1.2209959021372585),
    ('fixed5', 20.0): (0.15081792883603365, 1.1370887116922113),
    ('down_1cross', 1): (0.23361897444372776, 1.4229189760405996),
    ('up_1cross', 1): (0.1001398793467416, 0.9202838826245617),
    ('at_b_nf5', 1): (0.211846290459217, 1.3576452092023406),
    ('at_b_nf4', 1): (0.211846290459217, 1.3576452092023406),
    ('down_2cross', 1): (0.3574470649538047, 1.7379110018920867),
    ('up_top', 1): (0.09579560688857491, 0.8916579735293795),
    ('down_1cross', 2): (0.25067632253058963, 1.5434155210119274),
    ('up_1cross', 2): (0.099445159526491, 0.9112533776037285),
    ('at_b_nf5', 2): (0.22320859938942364, 1.4478130728735983),
    ('at_b_nf4', 2): (0.22320859938942364, 1.4478130728735983),
    ('down_2cross', 2): (0.4428785478864057, 2.126990202514348),
    ('up_top', 2): (0.09516147991136176, 0.8804290488925136),
    ('down_1cross', 3): (0.2523517797547945, 1.5607447520894466),
    ('up_1cross', 3): (0.09939314178711924, 0.9103467267814406),
    ('at_b_nf5', 3): (0.22396183249889934, 1.45799495402153),
    ('at_b_nf4', 3): (0.224135725651176, 1.4595215023569612),
    ('down_2cross', 3): (0.4631604805020409, 2.2429006864839938),
    ('up_top', 3): (0.09513764062341147, 0.8795114889334336),
}


def crundec_live():
    """Recompute CRUNDEC with the rundec wrapper; None if unavailable."""
    try:
        import rundec
    except ImportError:
        return None
    cr = rundec.CRunDec()

    def segs(mu0, nf0, mu1, nf1):
        out, mu, nf = [], mu0, nf0
        while nf != nf1:
            mth = rg.THRESHOLDS[nf + 1] if nf1 > nf else rg.THRESHOLDS[nf]
            out.append((mu, mth, nf, 1 if nf1 > nf else -1))
            nf += 1 if nf1 > nf else -1
            mu = mth
        out.append((mu, mu1, nf, 0))
        return out

    def path(mu0, nf0, mu1, nf1, L):
        a, m = A_MZ, 1.0
        for s, e, nf, st in segs(mu0, nf0, mu1, nf1):
            ae = cr.AlphasExact(a, s, e, nf, L) if s != e else a
            m = cr.mMS2mMS(m, a, ae, nf, L)
            a = ae
            if st < 0:
                m = cr.DecMqDownMS(m, a, e, e, nf - 1, L)
                a = cr.DecAsDownMS(a, e, e, nf - 1, L)
            elif st > 0:
                m = cr.DecMqUpMS(m, a, e, e, nf, L)
                a = cr.DecAsUpMS(a, e, e, nf, L)
        return a, m

    live = {}
    for mu in (10.0, 20.0):
        a = cr.AlphasExact(A_MZ, MZ, mu, 5, 1)
        live[('fixed5', mu)] = (a, cr.mMS2mMS(1.0, A_MZ, a, 5, 1))
    for L in (1, 2, 3):
        for k, c in CASES.items():
            live[(k, L)] = path(*c, L)
    return live


# ----------------------------------------------------------------------
print("1. Fixed n_f = 5, one loop, alpha_s(M_Z) = 0.1179, no thresholds")
for mu, (a_ref, m_ref) in ISSUE_ROWS.items():
    a = rg.alpha_s_fixed_nf(A_MZ, MZ, mu, 5, 1)
    m = rg.mass_factor_fixed_nf(A_MZ, a, 5, 1)
    check(f"alpha_s({mu:g} GeV) vs issue row", abs(a - a_ref) <= 1.5e-9, f"{a:.10f} vs {a_ref}")
    check(f"m({mu:g})/m(M_Z) vs issue row", abs(m - m_ref) <= 1.5e-9, f"{m:.10f} vs {m_ref}")
    a_old = rg.alpha_s_1loop_old(mu, MZ, A_MZ, 5)
    m_old = (a_old / A_MZ) ** rg.gamma_m(5)
    check(f"former code reproduces the issue's repository column at {mu:g} GeV",
          abs(a_old - OLD_ROWS[mu][0]) <= 1.5e-9 and abs(m_old - OLD_ROWS[mu][1]) <= 1.5e-9,
          f"{a_old:.9f}, {m_old:.9f}")
    check(f"former code fails the absolute check at {mu:g} GeV", abs(a_old - a_ref) > 1e-3)
check("identity at the starting scale: alpha_s",
      rg.alpha_s_fixed_nf(A_MZ, MZ, MZ, 5, 1) == A_MZ)
check("identity at the starting scale: mass factor",
      rg.mass_factor_fixed_nf(A_MZ, A_MZ, 5, 1) == 1.0)
check("identity of run_mass at the starting scale (3 loops)",
      rel(rg.run_mass(93.4, 2.0, 2.0, loops=3), 93.4) < 1e-15)
# one-loop analytic vs numerical integration of the same RGE
b0 = rg.beta_coeffs(5)[0]
from scipy.integrate import solve_ivp  # noqa: E402
sol = solve_ivp(lambda _, y: [-b0 * y[0]**2], (0, np.log(10.0**2 / MZ**2)), [A_MZ / np.pi],
                method='DOP853', rtol=1e-13, atol=1e-16)
check("one-loop analytic solution equals numerical integration",
      rel(sol.y[0, -1] * np.pi, rg.alpha_s_fixed_nf(A_MZ, MZ, 10.0, 5, 1)) < 1e-12)

# ----------------------------------------------------------------------
print("\n2. Thresholds: segment-by-segment running and MS-bar matching vs CRunDec")
live = crundec_live()
if live is None:
    print("  (rundec not importable: comparing with stored CRunDec 3.1 values only)")
else:
    worst = max(max(rel(live[k][0], v[0]), rel(live[k][1], v[1])) for k, v in CRUNDEC.items())
    check("live CRunDec reproduces the stored references", worst < 1e-14, f"max rel {worst:.1e}")
for L in (1, 2, 3):
    for k, (mu0, nf0, mu1, nf1) in CASES.items():
        a = rg.alpha_s_path(A_MZ, mu0, nf0, mu1, nf1, L)
        m = rg.run_mass_path(1.0, A_MZ, mu0, nf0, mu1, nf1, L)
        a_ref, m_ref = CRUNDEC[(k, L)]
        check(f"{L}-loop {k:12s} alpha_s and mass",
              rel(a, a_ref) < 5e-9 and rel(m, m_ref) < 5e-9,
              f"rel {rel(a, a_ref):.1e}, {rel(m, m_ref):.1e}")
check("active_nf at a threshold endpoint includes that quark",
      rg.active_nf(4.18) == 5 and rg.active_nf(1.275) == 4 and rg.active_nf(162.5) == 6)
check("default endpoint at m_b equals the explicit n_f = 5 value",
      rg.alpha_s_at(4.18, loops=3) == rg.alpha_s_at(4.18, 5, loops=3))
# Round trips. With continuous matching (1, 2 loops) they close exactly.
# At 3 loops the truncated up- and down-decoupling relations are inverse only
# through O(a^2); the residual is O(a^4) and must equal CRunDec's own
# round-trip residual (stored: alpha_s back at M_Z, mass ratio back at M_Z,
# alpha_s back at 1 GeV after M_Z -> 500 GeV -> 1 GeV).
ROUNDTRIP_3LOOP = (0.11789940812393357, 0.9999825531488528, 0.4631602388026455)
for L in (1, 2, 3):
    a1 = rg.alpha_s_at(1.0, loops=L)
    back = rg.alpha_s_path(a1, 1.0, 3, MZ, 5, L)
    m1 = rg.run_mass(1.0, MZ, 1.0, loops=L)
    mback = rg.run_mass_path(m1, a1, 1.0, 3, MZ, 5, L)
    a5 = rg.alpha_s_at(500.0, loops=L)
    down = rg.alpha_s_path(a5, 500.0, 6, 1.0, 3, L)
    if L < 3:
        check(f"{L}-loop round trip M_Z -> 1 GeV -> M_Z (two thresholds) closes",
              rel(back, A_MZ) < 1e-12 and rel(mback, 1.0) < 1e-12,
              f"rel {rel(back, A_MZ):.1e}, {rel(mback, 1.0):.1e}")
        check(f"{L}-loop round trip M_Z -> 500 GeV -> 1 GeV (three thresholds) closes",
              rel(down, a1) < 1e-12, f"rel {rel(down, a1):.1e}")
    else:
        ra, rm, rd = ROUNDTRIP_3LOOP
        a1_crd = CRUNDEC[('down_2cross', 3)][0]
        res = (back / A_MZ - 1, mback - 1, down / a1 - 1)
        res_crd = (ra / A_MZ - 1, rm - 1, rd / a1_crd - 1)
        check("3-loop round trip M_Z -> 1 GeV -> M_Z: residual O(a^4), equal to CRunDec's",
              all(abs(x - y) < 1e-8 for x, y in zip(res[:2], res_crd[:2])) and abs(res[0]) < 1e-4,
              f"residual {res[0]:.3e} vs {res_crd[0]:.3e} (alpha), {res[1]:.3e} vs {res_crd[1]:.3e} (mass)")
        check("3-loop round trip M_Z -> 500 GeV -> 1 GeV: residual equal to CRunDec's",
              abs(res[2] - res_crd[2]) < 1e-8, f"residual {res[2]:.3e} vs {res_crd[2]:.3e}")

# ----------------------------------------------------------------------
print("\n3. Flavour ratios")
trip = (93.4, 1275.0, 4180.0)
Q0 = rg.koide_Q(trip)
for (mu0, mu1, L) in ((5.0, 10.0, 1), (5.0, 40.0, 2), (5.0, MZ, 3)):
    run = [rg.run_mass(m, mu0, mu1, 5, 5, L) for m in trip]
    check(f"common-scale ratio invariance, n_f = 5, {L}-loop, {mu0:g} -> {mu1:g} GeV",
          abs(rg.koide_Q(run) - Q0) < 1e-12, f"{abs(rg.koide_Q(run) - Q0):.1e}")
run = [rg.run_mass(m, 2.0, MZ, loops=3) for m in trip]
check("common-scale ratio invariance across m_b with 3-loop matching",
      abs(rg.koide_Q(run) - Q0) < 1e-12)
a5_old = rg.alpha_s_1loop_old(5.0, MZ, A_MZ, 5)
fac_old = (rg.alpha_s_1loop_old(10.0, 5.0, a5_old, 5) / a5_old) ** rg.gamma_m(5)
run_old = [m * fac_old for m in trip]
check("the ratio control alone cannot see the former extra pi",
      abs(rg.koide_Q(run_old) - Q0) < 1e-12)
# regression of the recalibrated table (1-loop, eft prescription)
t1 = rg.table_values(rg.masses_at(2.0, 1, 'eft'))
tz = rg.table_values(rg.masses_at(MZ, 1, 'eft'))
check("delta(-s,c,b) at a common scale, 1-loop", abs(t1['(-s,c,b)'] - 3.1316) < 1e-3,
      f"{t1['(-s,c,b)']:.4f}%")
check("delta(c,b,t) at a common scale, 1-loop", abs(t1['(c,b,t)'] - 6.1183) < 1e-3,
      f"{t1['(c,b,t)']:.4f}%")
check("eft common-scale table is independent of the common scale",
      all(abs(t1[k] - tz[k]) < 1e-6 for k in t1), "2 GeV vs M_Z")
t3 = rg.table_values(rg.masses_at(2.0, 3, 'eft'))
check("3-loop common-scale values", abs(t3['(-s,c,b)'] - 4.5715) < 1e-3
      and abs(t3['(c,b,t)'] - 7.9316) < 1e-3, f"{t3['(-s,c,b)']:.4f}%, {t3['(c,b,t)']:.4f}%")

print()
if FAIL:
    print(f"{len(FAIL)} check(s) FAILED")
    sys.exit(1)
print("all checks passed")
