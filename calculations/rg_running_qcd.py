#!/usr/bin/env python3
"""
QCD running of MS-bar quark masses and the common-scale flavour table.

Language: Python 3 (numpy, scipy).

Library part (importable; used by check_rg_running_qcd.py):
  * alpha_s_1loop          fixed-n_f one-loop coupling (corrected: one factor of pi)
  * alpha_s_fixed_nf       fixed-n_f coupling at 1, 2 or 3 loops
  * mass_factor_fixed_nf   m(mu1)/m(mu0) at fixed n_f, RunDec c-function
  * alpha_s_at             coupling at any scale with segment-by-segment n_f
  * run_mass               MS-bar mass between any two scales, same path
  * active_nf              threshold convention

Conventions (all stated once here):
  * a = alpha_s/pi; d a / d ln mu^2 = -sum_i beta_i a^(i+2);
    d ln m / d ln mu^2 = -sum_i gamma_i a^(i+1)  (RunDec normalisation,
    Chetyrkin-Kuhn-Steinhauser, hep-ph/0004189; Herren-Steinhauser,
    arXiv:1703.03751).
  * L-loop running uses beta_0..beta_{L-1}, gamma_0..gamma_{L-1}; the mass
    factor is c(a1)/c(a0) with c expanded to order a^(L-1) (RunDec eq. 10),
    exactly as CRunDec's mMS2mMS.
  * Matching: at mu_h = m_h(m_h) (MS-bar heavy mass) with L-loop decoupling,
    i.e. terms through a^(L-1): alpha^(nl) = alpha^(nl+1) (1 + 11/72 a^2),
    m^(nl) = m^(nl+1) (1 + 89/432 a^2). At one loop both are continuous.
  * Active flavours: n_f(mu) counts the heavy flavours with mu_h <= mu, so a
    point exactly at a threshold belongs to the theory in which that quark is
    active (m_c(m_c) is a 4-flavour, m_b(m_b) a 5-flavour, m_t(m_t) a
    6-flavour quantity). Any call can override n_f at either endpoint.
  * Decoupled heavy flavour (a heavy quark's mass below its own threshold):
    'eft'  - the mass parameter is run segment by segment in the active
             theory, with the same matching as a light mass (default);
    'full' - the mass is run in the theory where the quark is active, with
             the coupling of that theory continued below threshold.
    Both are prescriptions; their difference is quoted as a convention
    uncertainty.

Report part (python3 calculations/rg_running_qcd.py): the flavour table at
the published reference scales and at a common scale, before/after the
one-loop fix, at 1/2/3 loops, both prescriptions, with uncertainties.
"""

import numpy as np
from scipy.integrate import solve_ivp

ZETA3 = 1.2020569031595942

# Reference coupling: PDG-style MS-bar alpha_s^(5)(M_Z)
MZ = 91.1876          # GeV
ALPHA_S_MZ = 0.1179
ALPHA_S_MZ_ERR = 0.0009

# Matching scales mu_h = m_h(m_h) (GeV), keyed by n_f above the threshold
THRESHOLDS = {4: 1.275, 5: 4.18, 6: 162.5}

# PDG MS-bar inputs at published reference scales: (mass MeV, scale GeV,
# minus error MeV, plus error MeV). Light quarks at 2 GeV, heavy at m_q(m_q).
MASSES_PDG = {
    'u': (2.16, 2.0, 0.26, 0.49),
    'd': (4.67, 2.0, 0.17, 0.48),
    's': (93.4, 2.0, 3.4, 8.6),
    'c': (1275.0, 1.275, 25.0, 25.0),
    'b': (4180.0, 4.18, 20.0, 30.0),
    't': (162500.0, 162.5, 1500.0, 2100.0),
}
OWN_NF = {'u': 3, 'd': 3, 's': 3, 'c': 4, 'b': 5, 't': 6}


# ----------------------------------------------------------------------
# Coefficients
# ----------------------------------------------------------------------
def beta_coeffs(nf):
    return (
        (11.0 - 2.0 * nf / 3.0) / 4.0,
        (102.0 - 38.0 * nf / 3.0) / 16.0,
        (2857.0 / 2.0 - 5033.0 * nf / 18.0 + 325.0 * nf**2 / 54.0) / 64.0,
    )


def gamma_coeffs(nf):
    return (
        1.0,
        (202.0 / 3.0 - 20.0 * nf / 9.0) / 16.0,
        (1249.0 + (-2216.0 / 27.0 - 160.0 * ZETA3 / 3.0) * nf
         - 140.0 * nf**2 / 81.0) / 64.0,
    )


def gamma_m(nf):
    """One-loop exponent gamma_0/beta_0 = 12/(33 - 2 n_f) of alpha_s."""
    return 12.0 / (33 - 2 * nf)


def active_nf(mu):
    """Active flavours at mu (GeV); a point exactly at mu_h has the quark active."""
    nf = 3
    for n in (4, 5, 6):
        if mu >= THRESHOLDS[n]:
            nf = n
    return nf


# ----------------------------------------------------------------------
# Fixed-n_f running
# ----------------------------------------------------------------------
def alpha_s_1loop(mu, mu0, alpha0, nf):
    """One-loop alpha_s at fixed n_f:
    1/alpha(mu) = 1/alpha0 + b0 ln(mu^2/mu0^2),  b0 = (33 - 2 n_f)/(12 pi).
    (The pre-2026-11 version divided b0 by pi a second time.)"""
    b0 = (33 - 2 * nf) / (12 * np.pi)
    return alpha0 / (1 + alpha0 * b0 * np.log(mu**2 / mu0**2))


def alpha_s_1loop_old(mu, mu0, alpha0, nf):
    """The former implementation, kept only to document the error."""
    b0 = (33 - 2 * nf) / (12 * np.pi)
    return alpha0 / (1 + alpha0 * b0 * np.log(mu**2 / mu0**2) / np.pi)


def alpha_s_fixed_nf(alpha0, mu0, mu1, nf, loops=1):
    """alpha_s^(nf)(mu1) from alpha_s^(nf)(mu0); exact solution of the
    truncated RGE (analytic at one loop, numerical otherwise)."""
    if mu0 == mu1:
        return alpha0
    if loops == 1:
        return alpha_s_1loop(mu1, mu0, alpha0, nf)
    b = beta_coeffs(nf)[:loops]

    def rhs(_, y):
        a = y[0]
        return [-sum(bi * a**(i + 2) for i, bi in enumerate(b))]

    sol = solve_ivp(rhs, (0.0, np.log(mu1**2 / mu0**2)), [alpha0 / np.pi],
                    method='DOP853', rtol=1e-13, atol=1e-16)
    return sol.y[0, -1] * np.pi


def c_function(alpha, nf, loops=1):
    """RunDec c-function truncated at order a^(loops-1)."""
    a = alpha / np.pi
    b, g = beta_coeffs(nf), gamma_coeffs(nf)
    B1, B2 = b[1] / b[0], b[2] / b[0]
    C0, C1, C2 = g[0] / b[0], g[1] / b[0], g[2] / b[0]
    terms = [1.0,
             (C1 - B1 * C0) * a,
             0.5 * ((C1 - B1 * C0)**2 + C2 - B1 * C1 + B1**2 * C0 - B2 * C0) * a**2]
    return a**C0 * sum(terms[:loops])


def mass_factor_fixed_nf(alpha0, alpha1, nf, loops=1):
    """m(mu1)/m(mu0) at fixed n_f given alpha_s at both scales."""
    return c_function(alpha1, nf, loops) / c_function(alpha0, nf, loops)


# ----------------------------------------------------------------------
# Decoupling at mu_h = m_h(m_h)
# ----------------------------------------------------------------------
def dec_alpha(alpha, down, loops):
    """alpha^(nl) <-> alpha^(nl+1) at mu = m_h(m_h); terms through a^(loops-1)."""
    if loops < 3:
        return alpha
    a = alpha / np.pi
    return alpha * (1 + 11.0 / 72.0 * a**2) if down else alpha * (1 - 11.0 / 72.0 * a**2)


def dec_mass(m, alpha, down, loops):
    """m^(nl) <-> m^(nl+1) at mu = m_h(m_h); alpha in the theory being left."""
    if loops < 3:
        return m
    a = alpha / np.pi
    return m * (1 + 89.0 / 432.0 * a**2) if down else m * (1 - 89.0 / 432.0 * a**2)


# ----------------------------------------------------------------------
# Segment-by-segment running
# ----------------------------------------------------------------------
def _path(mu0, nf0, mu1, nf1):
    """Segments [(mu_start, mu_end, nf)] and crossings between them."""
    segs, mu, nf = [], mu0, nf0
    while nf != nf1:
        if nf1 > nf:
            mth = THRESHOLDS[nf + 1]
            segs.append((mu, mth, nf, +1))
            nf += 1
        else:
            mth = THRESHOLDS[nf]
            segs.append((mu, mth, nf, -1))
            nf -= 1
        mu = mth
    segs.append((mu, mu1, nf, 0))
    return segs


def alpha_s_path(alpha0, mu0, nf0, mu1, nf1, loops=1):
    """alpha_s^(nf1)(mu1) from alpha_s^(nf0)(mu0) through the active theories."""
    alpha = alpha0
    for start, end, nf, step in _path(mu0, nf0, mu1, nf1):
        alpha = alpha_s_fixed_nf(alpha, start, end, nf, loops)
        if step:
            alpha = dec_alpha(alpha, down=(step < 0), loops=loops)
    return alpha


def alpha_s_at(mu, nf=None, loops=1, alpha_mz=ALPHA_S_MZ):
    """alpha_s^(nf)(mu) from alpha_s^(5)(M_Z); nf defaults to active_nf(mu)."""
    nf = active_nf(mu) if nf is None else nf
    return alpha_s_path(alpha_mz, MZ, 5, mu, nf, loops)


def run_mass(m0, mu0, mu1, nf0=None, nf1=None, loops=1, alpha_mz=ALPHA_S_MZ):
    """MS-bar mass m^(nf1)(mu1) from m^(nf0)(mu0), segment by segment,
    with alpha_s obtained from alpha_s^(5)(M_Z) along the same rules."""
    nf0 = active_nf(mu0) if nf0 is None else nf0
    nf1 = active_nf(mu1) if nf1 is None else nf1
    alpha = alpha_s_at(mu0, nf0, loops, alpha_mz)
    return run_mass_path(m0, alpha, mu0, nf0, mu1, nf1, loops)


def run_mass_path(m0, alpha0, mu0, nf0, mu1, nf1, loops=1):
    """m^(nf1)(mu1) from m^(nf0)(mu0) given alpha_s^(nf0)(mu0)."""
    alpha, m = alpha0, m0
    for start, end, nf, step in _path(mu0, nf0, mu1, nf1):
        a_end = alpha_s_fixed_nf(alpha, start, end, nf, loops)
        m *= mass_factor_fixed_nf(alpha, a_end, nf, loops)
        alpha = a_end
        if step:
            m = dec_mass(m, alpha, down=(step < 0), loops=loops)
            alpha = dec_alpha(alpha, down=(step < 0), loops=loops)
    return m


def run_quark(name, mu1, loops=1, prescription='eft', alpha_mz=ALPHA_S_MZ, m0=None):
    """Mass of quark `name` at mu1 starting from its PDG reference point."""
    mref, muref = MASSES_PDG[name][:2]
    m0 = mref if m0 is None else m0
    nf0 = active_nf(muref)
    if prescription == 'full' and mu1 < THRESHOLDS.get(OWN_NF[name], 0.0):
        # stay in the theory where the quark is active
        return run_mass(m0, muref, mu1, nf0, OWN_NF[name], loops, alpha_mz)
    return run_mass(m0, muref, mu1, nf0, None, loops, alpha_mz)


# ----------------------------------------------------------------------
# Former (pre-fix) common-target helper, kept to document 'before' numbers
# ----------------------------------------------------------------------
def _alpha_old(mu):
    alpha, mu_ref = ALPHA_S_MZ, MZ
    if mu >= MZ:
        if mu > 162.5:
            return alpha_s_1loop_old(mu, 162.5, alpha_s_1loop_old(162.5, MZ, alpha, 5), 6)
        return alpha_s_1loop_old(mu, MZ, alpha, 5)
    for th, nfa in ((4.18, 5), (1.275, 4)):
        if mu >= th:
            return alpha_s_1loop_old(mu, mu_ref, alpha, nfa)
        alpha, mu_ref = alpha_s_1loop_old(th, mu_ref, alpha, nfa), th
    return alpha_s_1loop_old(mu, mu_ref, alpha, 3)


def run_mass_old(m0, mu0, mu1):
    """The former run_mass_with_thresholds (extra pi, n_f+1 below crossings)."""
    m, mu = m0, mu0
    nf_of = lambda x: 6 if x >= 162.5 else 5 if x >= 4.18 else 4 if x >= 1.275 else 3
    if mu1 > mu:
        for th, nfb, _ in [(1.275, 3, 4), (4.18, 4, 5), (162.5, 5, 6)]:
            if mu < th <= mu1:
                m *= (_alpha_old(th - 1e-6) / _alpha_old(mu)) ** gamma_m(nfb + 1)
                mu = th
            elif mu < th:
                break
    else:
        for th, nfa, _ in [(162.5, 6, 5), (4.18, 5, 4), (1.275, 4, 3)]:
            if mu > th >= mu1:
                m *= (_alpha_old(th + 1e-6) / _alpha_old(mu)) ** gamma_m(nfa)
                mu = th
    return m * (_alpha_old(mu1) / _alpha_old(mu)) ** gamma_m(nf_of(mu1))


# ----------------------------------------------------------------------
# Flavour observables
# ----------------------------------------------------------------------
def koide_Q(masses, signs=(1, 1, 1)):
    """Q = sum m / (sum s_k sqrt m_k)^2; Q = 2/3 on the Koide surface."""
    m = np.asarray(masses, dtype=float)
    return m.sum() / np.dot(signs, np.sqrt(m))**2


def koide_dev(masses, signs=(1, 1, 1)):
    """delta = <z^2>/(2 z0^2) - 1 = (3/2) Q - 1, in percent."""
    return 100.0 * (1.5 * koide_Q(masses, signs) - 1.0)


TRIPLES = {
    '(-s,c,b)': (('s', 'c', 'b'), (-1, 1, 1)),
    '(c,b,t)': (('c', 'b', 't'), (1, 1, 1)),
    '(0,u+d,s)': None,
}


def table_values(mset):
    out = {}
    out['(-s,c,b)'] = koide_dev([mset['s'], mset['c'], mset['b']], (-1, 1, 1))
    out['(c,b,t)'] = koide_dev([mset['c'], mset['b'], mset['t']])
    out['(0,u+d,s)'] = 100.0 * (1.5 * (mset['u'] + mset['d'] + mset['s'])
                                / (np.sqrt(mset['u'] + mset['d']) + np.sqrt(mset['s']))**2 - 1.0)
    out['sqrt(ms/mc)'] = np.sqrt(mset['s'] / mset['c'])
    return out


def masses_at(mu, loops=1, prescription='eft', alpha_mz=ALPHA_S_MZ, inputs=None):
    inputs = inputs or {}
    return {q: run_quark(q, mu, loops, prescription, alpha_mz, inputs.get(q))
            for q in MASSES_PDG}


def affected_predictions(loops=1):
    """Mass relations re-evaluated with all masses at one common scale
    (eft prescription; the common scale drops out). Returned masses are
    converted back to the published reference points with the same running."""
    from scipy.optimize import brentq
    M = masses_at(2.0, loops, 'eft')
    Rc = M['c'] / MASSES_PDG['c'][0]
    Rb = M['b'] / MASSES_PDG['b'][0]
    Rt = M['t'] / MASSES_PDG['t'][0]
    seed_mc = (2 + np.sqrt(3))**2 * M['s'] / Rc
    v0_mb = (3 * np.sqrt(M['s']) + np.sqrt(M['c']))**2 / Rb
    cbt_mt = brentq(lambda mt: koide_Q([M['c'], M['b'], mt]) - 2 / 3, 1e4, 1e7) / Rt
    dual = koide_Q([1 / M['d'], 1 / M['s'], 1 / M['b']])
    return {'seed m_c(m_c)': seed_mc, 'v0-doubling m_b(m_b)': v0_mb,
            'Koide(c,b,t) m_t(m_t)': cbt_mt, 'dual Q(1/m_d,1/m_s,1/m_b)': dual}


def main():
    from scipy.optimize import brentq
    print("=" * 78)
    print("QCD running of MS-bar quark masses and the common-scale flavour table")
    print("=" * 78)
    print(f"alpha_s^(5)(M_Z) = {ALPHA_S_MZ} +- {ALPHA_S_MZ_ERR},  M_Z = {MZ} GeV")
    print("matching at mu_h = m_h(m_h):", {k: v for k, v in THRESHOLDS.items()})
    print("inputs (MS-bar, MeV @ GeV):",
          {q: f"{v[0]}@{v[1]}" for q, v in MASSES_PDG.items()})

    print("\n--- Fixed n_f = 5, one loop, from M_Z (compare CRunDec rows in issue #3) ---")
    for mu in (10.0, 20.0):
        a_new = alpha_s_fixed_nf(ALPHA_S_MZ, MZ, mu, 5, 1)
        a_old = alpha_s_1loop_old(mu, MZ, ALPHA_S_MZ, 5)
        print(f"  mu = {mu:4.0f} GeV: alpha_s = {a_new:.9f} (old {a_old:.9f}); "
              f"m(mu)/m(MZ) = {mass_factor_fixed_nf(ALPHA_S_MZ, a_new, 5, 1):.9f} "
              f"(old {(a_old/ALPHA_S_MZ)**gamma_m(5):.9f})")

    pdg = {q: v[0] for q, v in MASSES_PDG.items()}
    t_pdg = table_values(pdg)
    print("\n--- Flavour table, delta = (3/2)Q - 1 in % ---")
    keys = ['(-s,c,b)', '(c,b,t)', '(0,u+d,s)', 'sqrt(ms/mc)']
    hdr = f"  {'masses':44s}" + "".join(f"{k:>13s}" for k in keys)
    print(hdr)
    print(f"  {'published reference scales (mixed)':44s}"
          + "".join(f"{t_pdg[k]:13.4f}" for k in keys))
    old2 = {q: run_mass_old(v[0], v[1], 2.0) for q, v in MASSES_PDG.items()}
    t_old = table_values(old2)
    print(f"  {'common 2 GeV, former code (before)':44s}"
          + "".join(f"{t_old[k]:13.4f}" for k in keys))
    rows = {}
    for loops in (1, 2, 3):
        for presc in ('eft', 'full'):
            for mu in (2.0, MZ):
                t = table_values(masses_at(mu, loops, presc))
                rows[(loops, presc, mu)] = t
                lab = f"common {mu:g} GeV, {loops}-loop, {presc}"
                print(f"  {lab:44s}" + "".join(f"{t[k]:13.4f}" for k in keys))

    # uncertainties on the 1-loop eft common-scale values at 2 GeV
    base = rows[(1, 'eft', 2.0)]
    print("\n--- Uncertainties on the common-scale values (2 GeV, 1-loop, eft) ---")
    for k in keys:
        loop_spread = max(abs(rows[(l, 'eft', 2.0)][k] - base[k]) for l in (2, 3))
        presc = abs(rows[(1, 'full', 2.0)][k] - base[k])
        da = max(abs(table_values(masses_at(2.0, 1, 'eft', ALPHA_S_MZ + s * ALPHA_S_MZ_ERR))[k]
                     - base[k]) for s in (-1, 1))
        inp = 0.0
        for q, (m, _, lo, hi) in MASSES_PDG.items():
            d = max(abs(table_values(masses_at(2.0, 1, 'eft', inputs={q: m + e}))[k] - base[k])
                    for e in (-lo, hi))
            inp += d**2
        print(f"  {k:12s} value {base[k]:9.4f}   loop order {loop_spread:.4f}   "
              f"alpha_s {da:.4f}   heavy-quark prescription {presc:.4f}   "
              f"input masses {np.sqrt(inp):.4f}")

    print("\n--- Mass relations: published reference scales vs one common scale ---")
    mixed = {'seed m_c(m_c)': (2 + np.sqrt(3))**2 * pdg['s'],
             'v0-doubling m_b(m_b)': (3 * np.sqrt(pdg['s']) + np.sqrt(pdg['c']))**2,
             'Koide(c,b,t) m_t(m_t)': brentq(lambda mt: koide_Q([pdg['c'], pdg['b'], mt]) - 2 / 3,
                                             1e4, 1e7),
             'dual Q(1/m_d,1/m_s,1/m_b)': koide_Q([1 / pdg['d'], 1 / pdg['s'], 1 / pdg['b']])}
    aff = {L: affected_predictions(L) for L in (1, 2, 3)}
    print(f"  {'relation':28s}{'mixed scales':>14s}{'common 1-loop':>15s}"
          f"{'2-loop':>11s}{'3-loop':>11s}   (masses in MeV)")
    for k in mixed:
        print(f"  {k:28s}{mixed[k]:14.4f}" + "".join(f"{aff[L][k]:13.4f}" for L in (1, 2, 3)))

    print("\n--- Control: three active masses run from a common scale ---")
    trip = {'s': 93.4, 'c': 1275.0, 'b': 4180.0}
    Q0 = koide_Q(list(trip.values()))
    for mu1, loops in ((10.0, 1), (40.0, 2), (91.1876, 3)):
        run = [run_mass(m, 5.0, mu1, 5, 5, loops) for m in trip.values()]
        print(f"  5 GeV -> {mu1:7.3f} GeV, n_f = 5, {loops}-loop: "
              f"|Q - Q0| = {abs(koide_Q(run) - Q0):.1e}")
    run = [run_mass(m, 2.0, 91.1876, None, None, 3) for m in trip.values()]
    print(f"  2 GeV -> M_Z across m_b, 3-loop with matching: |Q - Q0| = "
          f"{abs(koide_Q(run) - Q0):.1e}")


if __name__ == "__main__":
    main()
