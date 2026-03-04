"""
staged_koide.py

Koide parametrization analysis for three fermion mass triples.

Parametrization:
    m_i = M0 * (1 + sqrt(2)*cos(2*pi*(i-1)/3 + delta))^2   for i = 1, 2, 3

Any (M0, delta) automatically gives Q = (m1+m2+m3)/(sum_i s_i*sqrt(m_i))^2 = 2/3
when s_i = +1 for all i. With sign flips, the Koide formula uses signed square roots
sigma_i = s_i * sqrt(m_i).

In the parametrization, sigma_i = sqrt(M0) * (1 + sqrt(2)*cos(2*pi*(i-1)/3 + delta))
and the Koide Q = 2/3 automatically by the algebra of the cosine vectors.

Definitions:
    sigma_i   = s_i * sqrt(m_i)           (signed square roots)
    v0        = (1/3) * sum_i sigma_i      (vacuum charge, per Rivero convention)
    M0        = (v0 / (1/3 * sum f_i))^2  but sum f_i = 3 always, so v0 = sqrt(M0)
              => M0 = v0^2 / 1 ... wait: sum_i (1+sqrt2*cos(theta_i)) = 3,
                 so v0 = (sum sigma_i)/3 = sqrt(M0).
    Thus: M0 = v0^2, where v0 = (1/3) sum_i sigma_i.

Seed condition: mass slot k vanishes when f_k(delta) = 0, i.e.
    1 + sqrt(2)*cos(2*pi*k/3 + delta) = 0
    => cos(2*pi*k/3 + delta) = -1/sqrt(2)
    => 2*pi*k/3 + delta = 3*pi/4 + 2*pi*n  or  5*pi/4 + 2*pi*n
    => delta_seed = 3*pi/4 - 2*pi*k/3  (mod 2*pi)  [principal]
       delta_seed = 5*pi/4 - 2*pi*k/3  (mod 2*pi)  [second branch]

v0-doubling (per v0_doubling.md):
    v0(seed)  = (1/3) * (sigma_j + sigma_k)        where j,k are nonzero slots at seed
    v0(full)  = (1/3) * (sigma_1 + sigma_2 + sigma_3)  physical value
    Ratio v0(full)/v0(seed) — reported as ~2 for the quark bloom.
"""

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from itertools import permutations as _perms

# ──────────────────────────────────────────────────────────────────────────────
# Physical masses (MeV), PDG values
# ──────────────────────────────────────────────────────────────────────────────

M = {
    "e":   0.51099895,
    "mu":  105.6583755,
    "tau": 1776.86,
    "u":   2.16,
    "d":   4.67,
    "s":   93.4,
    "c":   1270.0,
    "b":   4180.0,
    "t":   172760.0,
}

# ──────────────────────────────────────────────────────────────────────────────
# Core functions
# ──────────────────────────────────────────────────────────────────────────────

def f_vec(delta):
    """Return f = (f0, f1, f2) where f_k = 1 + sqrt(2)*cos(2*pi*k/3 + delta)."""
    sq2 = np.sqrt(2)
    return np.array([
        1 + sq2 * np.cos(delta),
        1 + sq2 * np.cos(delta + 2*np.pi/3),
        1 + sq2 * np.cos(delta + 4*np.pi/3),
    ])


def koide_model_masses(M0, delta):
    """Return (m0, m1, m2) = M0 * f_k^2."""
    return M0 * f_vec(delta)**2


def koide_Q(m_vals, signs):
    """Q = sum(m_i) / (sum(s_i * sqrt(m_i)))^2."""
    m = np.array(m_vals, dtype=float)
    s = np.array(signs, dtype=float)
    num = np.sum(m)
    denom = np.sum(s * np.sqrt(m))
    return num / denom**2


def sigma_vec(m_vals, signs):
    """Signed square roots: sigma_i = s_i * sqrt(m_i)."""
    return np.array(signs, dtype=float) * np.sqrt(np.array(m_vals, dtype=float))


def v0_from_sigma(sig):
    """v0 = (1/3) * sum(sigma_i). Since sum f_k = 3, we have sqrt(M0) = v0."""
    return np.sum(sig) / 3.0


# ──────────────────────────────────────────────────────────────────────────────
# Fitting: given sigma_phys (signed sqrt masses in physical order),
# find (M0, delta, permutation) that minimizes reconstruction error.
#
# For each permutation (physical -> model slot assignment) and a given delta:
#   sp[k] = sigma_phys[perm[k]]  (sigma in model slot order)
#   Model:  sp[k] = sqrt(M0) * f_k(delta)
#   Optimal sqrt(M0) = (sp . f) / (f . f)   [least squares]
#   Residual = ||sp - sqM0 * f||^2 / ||sp||^2   [fractional]
# ──────────────────────────────────────────────────────────────────────────────

def fit_koide(m_vals, signs):
    """
    Fit Koide parametrization to physical masses with sign assignments.
    Returns (M0, delta, perm_list, rms_frac_error) where perm[k] = physical index
    assigned to model slot k (0-indexed).
    """
    sig_phys = sigma_vec(m_vals, signs)

    best = None
    best_res = np.inf

    for perm in _perms(range(3)):
        sp = sig_phys[list(perm)]  # sigma in model slot order

        def residual(delta):
            f = f_vec(delta)
            ff = np.dot(f, f)
            if ff < 1e-30:
                return 1e30
            sqM0 = np.dot(sp, f) / ff
            if sqM0 <= 0:
                return 1e30
            return np.sum((sp - sqM0 * f)**2) / np.dot(sp, sp)

        # Dense scan over [0, 2*pi)
        n = 7200
        ds = np.linspace(0, 2*np.pi, n, endpoint=False)
        rs = np.array([residual(d) for d in ds])
        idx = np.argmin(rs)

        # Refine around minimum
        d_lo = ds[max(0, idx-2)]
        d_hi = ds[min(n-1, idx+2)]
        opt = minimize_scalar(residual, bounds=(d_lo, d_hi), method='bounded',
                              options={'xatol': 1e-14})
        d_opt = opt.x % (2*np.pi)
        r_opt = opt.fun

        f = f_vec(d_opt)
        ff = np.dot(f, f)
        sqM0 = np.dot(sp, f) / ff
        M0 = sqM0**2

        if r_opt < best_res and sqM0 > 0:
            best_res = r_opt
            best = (M0, d_opt, list(perm), np.sqrt(r_opt))

    return best


# ──────────────────────────────────────────────────────────────────────────────
# Seed angles
# ──────────────────────────────────────────────────────────────────────────────

def all_seed_info():
    """
    All 6 seed (delta, zero_slot) pairs in [0, 2*pi).
    Zero slot k: f_k(delta_seed) = 0.
    """
    seeds = []
    for k in range(3):
        for branch_angle in [3*np.pi/4, 5*np.pi/4]:
            ds = (branch_angle - 2*np.pi*k/3) % (2*np.pi)
            fk = 1 + np.sqrt(2)*np.cos(2*np.pi*k/3 + ds)
            seeds.append((ds, k, fk))
    return seeds


def find_nearest_seed(delta_full):
    """
    Return (delta_seed, zero_slot) for the seed angle closest to delta_full.
    """
    seeds = all_seed_info()
    best_ds, best_slot = None, None
    best_dist = np.inf
    for (ds, slot, check) in seeds:
        dist = abs((ds - delta_full + np.pi) % (2*np.pi) - np.pi)
        if dist < best_dist:
            best_dist = dist
            best_ds, best_slot = ds, slot
    return best_ds, best_slot


# ──────────────────────────────────────────────────────────────────────────────
# v0 computation: seed vs full
#
# v0(full) = (1/3) * sum_i sigma_i  (physical signed sqrt masses)
#
# v0(seed): at the seed, one sigma vanishes. The seed v0 is defined as
#   v0(seed) = (1/3) * sum_{i: nonzero} sigma_i
# where the "nonzero" sigma values are the physical values of those two masses.
# This is the definition used in v0_doubling.md.
#
# The v0-doubling ratio r = v0(full) / v0(seed).
# ──────────────────────────────────────────────────────────────────────────────

def compute_v0_pair(sig_phys, perm, zero_slot):
    """
    sig_phys: signed sqrt masses in physical order.
    perm: model->physical index map.
    zero_slot: model slot that vanishes at seed.
    Returns (v0_full, v0_seed, ratio).
    """
    v0_full = np.sum(sig_phys) / 3.0

    # v0_seed: sum of nonzero physical sigmas / 3
    nonzero_phys_indices = [perm[k] for k in range(3) if k != zero_slot]
    v0_seed = np.sum(sig_phys[nonzero_phys_indices]) / 3.0

    ratio = v0_full / v0_seed if abs(v0_seed) > 1e-30 else None
    return v0_full, v0_seed, ratio


# ──────────────────────────────────────────────────────────────────────────────
# Overlap prediction
# ──────────────────────────────────────────────────────────────────────────────

def solve_for_mb_given_mc(ms, mc, mt, constraint):
    """
    constraint = 'scb': Q(-s,c,b)=2/3, return mb.
    constraint = 'cbt': Q(c,b,t)=2/3, return mb.
    Returns list of positive solutions for mb.
    """
    if constraint == 'scb':
        # (-sqrt(ms) + sqrt(mc) + sqrt(mb))^2 = (3/2)*(ms + mc + mb)
        # Let A = -sqrt(ms)+sqrt(mc), y=sqrt(mb):
        # (A+y)^2 = (3/2)(ms+mc+y^2)
        # A^2+2Ay+y^2 = (3/2)(ms+mc) + (3/2)y^2
        # -(1/2)y^2 + 2Ay + A^2 - (3/2)(ms+mc) = 0
        # y^2 - 4Ay - 2A^2 + 3(ms+mc) = 0
        A = -np.sqrt(ms) + np.sqrt(mc)
        disc = 4*A**2 - (3*(ms+mc) - 2*A**2)
        # More carefully: discriminant of y^2 - 4Ay + (3(ms+mc) - 2A^2) = 0... let me redo
        # From Q=2/3: 3*sum = 2*(sum_sigma)^2
        # 3*(ms+mc+mb) = 2*(-sqrt(ms)+sqrt(mc)+sqrt(mb))^2
        # 3*(ms+mc+mb) = 2*(A+y)^2 where A=-sqrt(ms)+sqrt(mc), y=sqrt(mb)
        # 3*ms + 3*mc + 3*y^2 = 2*A^2 + 4*A*y + 2*y^2
        # y^2 - 4*A*y + 3*(ms+mc) - 2*A^2 = 0
        disc = 16*A**2 - 4*(3*(ms+mc) - 2*A**2)
        # = 16A^2 - 12(ms+mc) + 8A^2 = 24A^2 - 12(ms+mc)
        disc = 24*A**2 - 12*(ms+mc)
        if disc < 0:
            return []
        sq = np.sqrt(disc)
        return [y**2 for y in [(4*A + sq)/2, (4*A - sq)/2] if y > 0]

    elif constraint == 'cbt':
        # (sqrt(mc)+sqrt(mb)+sqrt(mt))^2 = (3/2)*(mc+mb+mt)
        # A = sqrt(mc)+sqrt(mt), y=sqrt(mb):
        # 3*(mc+mt+y^2) = 2*(A+y)^2
        # y^2 - 4Ay + 3*(mc+mt) - 2*A^2 = 0
        A = np.sqrt(mc) + np.sqrt(mt)
        disc = 24*A**2 - 12*(mc+mt)
        if disc < 0:
            return []
        sq = np.sqrt(disc)
        return [y**2 for y in [(4*A + sq)/2, (4*A - sq)/2] if y > 0]

    return []


def overlap_prediction(ms, mt):
    """
    Find (mc, mb) satisfying both Q(-s,c,b)=2/3 and Q(c,b,t)=2/3.
    Returns list of (mc, mb) solutions.
    """
    def diff(mc):
        sols_scb = solve_for_mb_given_mc(ms, mc, None, 'scb')
        sols_cbt = solve_for_mb_given_mc(None, mc, mt,  'cbt')
        # Physical: mb > mc (from scb), mc < mb < mt (from cbt)
        mb_scb = max([x for x in sols_scb if x > mc], default=None)
        mb_cbt_cands = [x for x in sols_cbt if mc < x < mt]
        mb_cbt = min(mb_cbt_cands) if mb_cbt_cands else None
        if mb_scb is None or mb_cbt is None:
            return np.nan
        return mb_scb - mb_cbt

    mc_scan = np.linspace(200, mt/2, 100000)
    d_arr = np.array([diff(mc) for mc in mc_scan])
    solutions = []
    for i in range(len(d_arr)-1):
        if np.isfinite(d_arr[i]) and np.isfinite(d_arr[i+1]) and d_arr[i]*d_arr[i+1] < 0:
            mc_sol = brentq(diff, mc_scan[i], mc_scan[i+1], xtol=1e-8)
            mb_list = solve_for_mb_given_mc(ms, mc_sol, None, 'scb')
            mb_sol = max([x for x in mb_list if x > mc_sol], default=None)
            if mb_sol:
                solutions.append((mc_sol, mb_sol))
    return solutions


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    out = []

    def p(s=""):
        out.append(s + "\n")
        print(s)

    p("# Staged Koide Analysis")
    p()
    p("Parametrization:  m_i = M0 * (1 + sqrt(2)*cos(2*pi*(i-1)/3 + delta))^2")
    p("This automatically gives Q = 2/3 for any (M0, delta).")
    p()
    p("Definitions:")
    p("  sigma_i  = s_i * sqrt(m_i)  (signed square root)")
    p("  v0       = (1/3) * sum_i sigma_i   (vacuum charge; equals sqrt(M0) in the parametrization)")
    p("  M0       = v0^2  (since sum_k f_k = 3 always)")
    p()
    p("Seed: mass slot k vanishes when f_k(delta) = 0, i.e.")
    p("  delta_seed = 3*pi/4 - 2*pi*k/3  (mod 2*pi)  [principal branch]")
    p("  delta_seed = 5*pi/4 - 2*pi*k/3  (mod 2*pi)  [second branch]")
    p()
    p("v0-doubling definition (from v0_doubling.md):")
    p("  v0(seed) = (1/3) * sum_{nonzero slots} sigma_i")
    p("  v0(full) = (1/3) * sum_i sigma_i  (all three)")
    p("  Ratio v0(full)/v0(seed) -- the v0-doubling ratio.")

    # ── Define triples ────────────────────────────────────────────────────────

    triples = [
        dict(
            name="(e, mu, tau)",
            label="emu_tau",
            masses=[M["e"], M["mu"], M["tau"]],
            signs=(+1, +1, +1),
            names=["e", "mu", "tau"],
        ),
        dict(
            name="(-s, c, b)",
            label="scb",
            masses=[M["s"], M["c"], M["b"]],
            signs=(-1, +1, +1),
            names=["s", "c", "b"],
        ),
        dict(
            name="(c, b, t)",
            label="cbt",
            masses=[M["c"], M["b"], M["t"]],
            signs=(+1, +1, +1),
            names=["c", "b", "t"],
        ),
    ]

    # ── Part 1 ────────────────────────────────────────────────────────────────

    p()
    p("## Part 1: Fitted Parameters for Each Triple")

    results = {}

    for tr in triples:
        name  = tr["name"]
        label = tr["label"]
        mv    = tr["masses"]
        signs = tr["signs"]
        pn    = tr["names"]

        p()
        p(f"### Triple {name}")
        p()
        p(f"Physical masses (MeV):  {pn[0]} = {mv[0]},  {pn[1]} = {mv[1]},  {pn[2]} = {mv[2]}")

        sig_phys = sigma_vec(mv, signs)
        p(f"Signed sqrt masses (sigma):  ({sig_phys[0]:.6f}, {sig_phys[1]:.6f}, {sig_phys[2]:.6f})  sqrt(MeV)")

        Q_dir = koide_Q(mv, signs)
        p(f"Koide Q = {Q_dir:.8f}   (exact 2/3 = {2/3:.8f})")

        v0_full = v0_from_sigma(sig_phys)
        M0_v0 = v0_full**2
        p(f"v0_full = (1/3)*sum(sigma) = {v0_full:.6f}  sqrt(MeV)")
        p(f"M0 from v0^2 = {M0_v0:.4f}  MeV  (valid for exact Q=2/3)")

        # Fit parametrization
        M0f, df, perm, rms_err = fit_koide(mv, signs)
        v0f = np.sqrt(M0f)   # In parametrization, v0 = sqrt(M0)

        p()
        p("Fitted (M0, delta) by least-squares on sigma values:")
        p(f"  M0_full    = {M0f:.6f} MeV")
        p(f"  v0_full (from sqrt(M0_full)) = {v0f:.6f}  sqrt(MeV)")
        p(f"  delta_full = {df:.8f} rad  =  {np.degrees(df):.6f} deg  =  {df/np.pi:.8f} pi")
        p(f"  Index mapping: model slot (0,1,2) -> physical index {tuple(perm)}")
        p(f"  RMS fractional reconstruction error = {rms_err*100:.4f}%")

        # Reconstruction check
        f_full = f_vec(df)
        m_recon = M0f * f_full**2
        p()
        p("Reconstruction check (model slot order):")
        for k in range(3):
            phys_name = pn[perm[k]]
            phys_mass = mv[perm[k]]
            err = (m_recon[k] - phys_mass)/phys_mass * 100
            p(f"  slot {k}: {phys_name}  model = {m_recon[k]:.4f} MeV  physical = {phys_mass:.4f} MeV  err = {err:+.4f}%")

        # Seed
        ds, zero_slot = find_nearest_seed(df)
        zero_phys = pn[perm[zero_slot]]
        dbloom = df - ds
        # Wrap to (-pi/3, pi/3)
        dbloom_w = (dbloom + np.pi/3) % (2*np.pi/3) - np.pi/3

        p()
        p("Seed parameters:")
        p(f"  delta_seed = {ds:.8f} rad  =  {np.degrees(ds):.6f} deg  =  {ds/np.pi:.8f} pi")
        p(f"  Zero slot: model slot {zero_slot}  ->  physical mass {zero_phys}")
        p(f"  Bloom rotation: Delta_delta = delta_full - delta_seed = {dbloom:.6f} rad  =  {np.degrees(dbloom):.6f} deg")
        p(f"  Delta_delta / (2*pi/3) = {dbloom/(2*np.pi/3):.6f}")

        # v0 computation (the physical definition from v0_doubling.md)
        v0_full_phys, v0_seed_phys, v0_ratio = compute_v0_pair(sig_phys, perm, zero_slot)
        p()
        p("v0 values (physical definition: (1/3)*sum sigma):")
        p(f"  v0_full = (1/3)*({'+'.join(f'{s:.4f}' for s in sig_phys)}) = {v0_full_phys:.6f}  sqrt(MeV)")
        nonzero_slots = [k for k in range(3) if k != zero_slot]
        nonzero_sig = [sig_phys[perm[k]] for k in nonzero_slots]
        nonzero_names = [pn[perm[k]] for k in nonzero_slots]
        p(f"  v0_seed = (1/3)*(sigma_{nonzero_names[0]} + sigma_{nonzero_names[1]}) = (1/3)*({nonzero_sig[0]:.4f} + {nonzero_sig[1]:.4f}) = {v0_seed_phys:.6f}  sqrt(MeV)")
        if v0_ratio is not None:
            p(f"  v0 ratio: v0_full / v0_seed = {v0_ratio:.6f}   (ideal: 2.0000)")

        # Seed masses: what the parametrization gives at delta_seed with M0 from v0_seed
        M0_seed = v0_seed_phys**2
        f_seed = f_vec(ds)
        m_seed = M0_seed * f_seed**2
        p()
        p(f"Seed triple masses (M0_seed = v0_seed^2 = {M0_seed:.4f} MeV, delta_seed):")
        for k in range(3):
            phys_name = pn[perm[k]]
            if k == zero_slot:
                p(f"  {phys_name}_seed = {m_seed[k]:.3e} MeV  (zero slot)")
            else:
                m_phys_k = mv[perm[k]]
                p(f"  {phys_name}_seed = {m_seed[k]:.4f} MeV   physical = {m_phys_k:.4f} MeV   ratio = {m_seed[k]/m_phys_k:.6f}")

        results[label] = dict(
            M0_full=M0f, delta_full=df, v0_full_param=v0f,
            v0_full=v0_full_phys, v0_seed=v0_seed_phys, v0_ratio=v0_ratio,
            M0_seed=M0_seed, v0_seed_M0=np.sqrt(M0_seed),
            delta_seed=ds, delta_bloom=dbloom, delta_bloom_w=dbloom_w,
            zero_slot=zero_slot, zero_mass_name=zero_phys,
            perm=perm, Q_direct=Q_dir,
            names=pn, masses=mv, signs=signs, sig_phys=sig_phys,
            m_seed=m_seed, f_seed=f_seed,
        )

        p()

    # ── Part 2: Staged picture ─────────────────────────────────────────────

    p()
    p("## Part 2: The Staged Picture")
    p()
    p("Hypothesis: quark mass spectrum arises in stages:")
    p("  Stage 1 (EW scale): (c,b,t) triple; top Yukawa dominates. Seed has c_seed = 0.")
    p("  Stage 2 (chiral):  (-s,c,b) triple; s acquires mass. Seed has s_seed = 0.")
    p("  Stage 3 (bloom):   s acquires its physical mass, rotating delta from delta_seed.")
    p()

    for stage, label, desc in [
        (1, "cbt",     "(c, b, t) — EW-determined"),
        (2, "scb",     "(-s, c, b) — chiral"),
        (3, "emu_tau", "(e, mu, tau) — leptonic"),
    ]:
        r = results[label]
        p(f"### Stage {stage}: {desc}")
        p(f"Zero mass at seed: {r['zero_mass_name']}")
        p(f"delta_seed = {r['delta_seed']:.6f} rad = {r['delta_seed']/np.pi:.6f} pi = {np.degrees(r['delta_seed']):.4f} deg")
        p(f"delta_full = {r['delta_full']:.6f} rad = {r['delta_full']/np.pi:.6f} pi = {np.degrees(r['delta_full']):.4f} deg")
        p(f"Bloom Delta_delta = {r['delta_bloom']:.6f} rad = {np.degrees(r['delta_bloom']):.4f} deg")
        p()
        p("Seed triple masses:")
        mv, perm, pn = r["masses"], r["perm"], r["names"]
        for k in range(3):
            phys_name = pn[perm[k]]
            m_phys = mv[perm[k]]
            if k == r["zero_slot"]:
                p(f"  {phys_name}_seed = 0  (zero slot; physical mass = {m_phys:.4f} MeV)")
            else:
                m_seed_k = r["m_seed"][k]
                p(f"  {phys_name}_seed = {m_seed_k:.4f} MeV   physical = {m_phys:.4f} MeV   ratio = {m_seed_k/m_phys:.6f}")
        p()

    # ── Part 3: Overlap prediction ──────────────────────────────────────────

    p()
    p("## Part 3: Overlap Prediction from Dual Q=2/3 Constraints")
    p()
    p("Given m_s = 93.4 MeV and m_t = 172760 MeV as inputs,")
    p("find (m_c, m_b) satisfying Q(-s,c,b) = 2/3 AND Q(c,b,t) = 2/3 simultaneously.")
    p()

    ms_in = M["s"]
    mt_in = M["t"]
    sols = overlap_prediction(ms_in, mt_in)

    if sols:
        for i, (mc_p, mb_p) in enumerate(sols):
            p(f"Solution {i+1}:")
            p(f"  m_c (predicted) = {mc_p:.4f} MeV   physical = {M['c']:.4f} MeV   deviation = {(mc_p-M['c'])/M['c']*100:+.3f}%")
            p(f"  m_b (predicted) = {mb_p:.4f} MeV   physical = {M['b']:.4f} MeV   deviation = {(mb_p-M['b'])/M['b']*100:+.3f}%")
            Qscb = koide_Q([ms_in, mc_p, mb_p], (-1, +1, +1))
            Qcbt = koide_Q([mc_p, mb_p, mt_in], (+1, +1, +1))
            p(f"  Verification: Q(-s,c,b) = {Qscb:.10f}   Q(c,b,t) = {Qcbt:.10f}")
            p(f"                (both should equal 2/3 = {2/3:.10f})")
            p()
    else:
        p("No solution found in scanned range.")

    # ── Part 4: v0-doubling ─────────────────────────────────────────────────

    p()
    p("## Part 4: v0-Doubling")
    p()
    p("v0 = (1/3)*sum(sigma_i) where sigma_i are the signed square roots.")
    p("v0(seed) = (1/3)*sum of the two nonzero sigma values at the seed.")
    p("v0-doubling hypothesis: v0(full)/v0(seed) = 2.")
    p()

    for tr in triples:
        r = results[tr["label"]]
        p(f"### Triple {tr['name']}")
        sig = r["sig_phys"]
        perm, pn = r["perm"], r["names"]
        nonzero = [k for k in range(3) if k != r["zero_slot"]]
        nsigs = [sig[perm[k]] for k in nonzero]
        nnames = [pn[perm[k]] for k in nonzero]
        p(f"  sigma values: {pn[0]}={sig[0]:.6f}, {pn[1]}={sig[1]:.6f}, {pn[2]}={sig[2]:.6f}  sqrt(MeV)")
        p(f"  v0_full = (1/3)*({sig[0]:.4f}+{sig[1]:.4f}+{sig[2]:.4f}) = {r['v0_full']:.6f}  sqrt(MeV)")
        p(f"  v0_seed = (1/3)*(sigma_{nnames[0]}+sigma_{nnames[1]}) = (1/3)*({nsigs[0]:.4f}+{nsigs[1]:.4f}) = {r['v0_seed']:.6f}  sqrt(MeV)")
        if r["v0_ratio"] is not None:
            dev = abs(r["v0_ratio"] - 2.0)*100/2.0
            p(f"  v0_full/v0_seed = {r['v0_ratio']:.6f}   ideal: 2.000   deviation: {dev:.3f}%")
        p()

    p("### Cross-check: explicit v0-doubling formula for quarks")
    p()
    p("From v0_doubling.md: the quark bloom has")
    p("  v0(seed) = (sqrt(m_s) + sqrt(m_c))/3")
    p("  v0(full) = (-sqrt(m_s) + sqrt(m_c) + sqrt(m_b))/3")
    p()
    sq_s = np.sqrt(M["s"])
    sq_c = np.sqrt(M["c"])
    sq_b = np.sqrt(M["b"])
    v0_seed_explicit = (sq_s + sq_c)/3
    v0_full_explicit = (-sq_s + sq_c + sq_b)/3
    ratio_explicit = v0_full_explicit / v0_seed_explicit
    p(f"  sqrt(m_s) = {sq_s:.6f}, sqrt(m_c) = {sq_c:.6f}, sqrt(m_b) = {sq_b:.6f}  sqrt(MeV)")
    p(f"  v0(seed) = ({sq_s:.6f} + {sq_c:.6f})/3 = {v0_seed_explicit:.6f}")
    p(f"  v0(full) = (-{sq_s:.6f} + {sq_c:.6f} + {sq_b:.6f})/3 = {v0_full_explicit:.6f}")
    p(f"  Ratio v0(full)/v0(seed) = {ratio_explicit:.6f}   (reported value: 2.0005)")
    p()
    p("v0-doubling prediction for m_b: set ratio = 2 exactly.")
    p("  2*v0(seed) = v0(full):")
    p("  2*(sqrt(m_s)+sqrt(m_c))/3 = (-sqrt(m_s)+sqrt(m_c)+sqrt(m_b))/3")
    p("  2*sqrt(m_s)+2*sqrt(m_c) = -sqrt(m_s)+sqrt(m_c)+sqrt(m_b)")
    p("  sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c)")
    sq_b_pred = 3*sq_s + sq_c
    mb_pred_v0 = sq_b_pred**2
    p(f"  sqrt(m_b) = 3*{sq_s:.6f} + {sq_c:.6f} = {sq_b_pred:.6f}")
    p(f"  m_b (predicted) = {mb_pred_v0:.4f} MeV   physical = {M['b']:.4f} MeV   deviation = {(mb_pred_v0-M['b'])/M['b']*100:+.4f}%")

    p()
    p("### Cross-consistency: (c,b,t) seed b-mass vs physical b")
    p()
    cbt_r = results["cbt"]
    perm_cbt = cbt_r["perm"]
    pn_cbt = cbt_r["names"]
    p("The (c,b,t) seed sets c=0. Seed masses are computed from M0_seed = v0_seed^2.")
    p(f"v0_seed(c,b,t) = {cbt_r['v0_seed']:.6f}  sqrt(MeV)")
    p(f"M0_seed(c,b,t) = {cbt_r['M0_seed']:.4f} MeV")
    p()
    for k in range(3):
        phys_name = pn_cbt[perm_cbt[k]]
        m_seed_k = cbt_r["m_seed"][k]
        if k == cbt_r["zero_slot"]:
            p(f"  {phys_name}_seed = {m_seed_k:.3e} MeV  (zero)")
        else:
            m_phys_k = cbt_r["masses"][perm_cbt[k]]
            p(f"  {phys_name}_seed = {m_seed_k:.4f} MeV   physical = {m_phys_k:.4f} MeV   ratio = {m_seed_k/m_phys_k:.6f}")

    # ── Part 5: Summary tables ──────────────────────────────────────────────

    p()
    p()
    p("## Part 5: Summary Tables")
    p()

    p("### Table 1: Koide Q Values")
    p()
    p("| Triple | Signs | Q (direct) | Q - 2/3 |")
    p("|--------|-------|-----------|---------|")
    sign_strs = {"emu_tau": "(+,+,+)", "scb": "(-,+,+)", "cbt": "(+,+,+)"}
    for tr in triples:
        r = results[tr["label"]]
        p(f"| {tr['name']} | {sign_strs[tr['label']]} | {r['Q_direct']:.8f} | {r['Q_direct']-2/3:+.2e} |")

    p()
    p("### Table 2: Full Triple Parameters")
    p()
    p("| Triple | M0_full (MeV) | delta_full / pi | v0_full (sqrt(MeV)) |")
    p("|--------|--------------|----------------|---------------------|")
    for tr in triples:
        r = results[tr["label"]]
        p(f"| {tr['name']} | {r['M0_full']:.4f} | {r['delta_full']/np.pi:.6f} | {r['v0_full']:.6f} |")

    p()
    p("### Table 3: Seed Parameters")
    p()
    p("| Triple | Zero mass | delta_seed / pi | v0_seed (sqrt(MeV)) | M0_seed (MeV) |")
    p("|--------|-----------|----------------|--------------------|--------------| ")
    for tr in triples:
        r = results[tr["label"]]
        p(f"| {tr['name']} | {r['zero_mass_name']} | {r['delta_seed']/np.pi:.6f} | {r['v0_seed']:.6f} | {r['M0_seed']:.4f} |")

    p()
    p("### Table 4: Bloom Rotations")
    p()
    p("| Triple | Delta_delta (rad) | Delta_delta (deg) | Delta_delta / (2*pi/3) |")
    p("|--------|------------------|------------------|------------------------|")
    for tr in triples:
        r = results[tr["label"]]
        p(f"| {tr['name']} | {r['delta_bloom']:.6f} | {np.degrees(r['delta_bloom']):.4f} | {r['delta_bloom']/(2*np.pi/3):.6f} |")

    p()
    p("### Table 5: v0-Doubling Ratios")
    p()
    p("v0 = (1/3)*sum(sigma_i).  v0_seed = (1/3)*sum of nonzero sigma values at seed.")
    p()
    p("| Triple | Zero mass | v0_full | v0_seed | v0_full/v0_seed | Deviation from 2 |")
    p("|--------|-----------|--------|--------|----------------|-----------------|")
    for tr in triples:
        r = results[tr["label"]]
        if r["v0_ratio"] is not None:
            dev = abs(r["v0_ratio"] - 2)*100/2
            p(f"| {tr['name']} | {r['zero_mass_name']} | {r['v0_full']:.5f} | {r['v0_seed']:.5f} | {r['v0_ratio']:.6f} | {dev:.3f}% |")

    p()
    p("### Table 6: Overlap Prediction")
    p()
    p(f"Inputs: m_s = {M['s']} MeV,  m_t = {M['t']} MeV")
    p()
    p("| Mass | Predicted (MeV) | Physical (MeV) | Deviation |")
    p("|------|----------------|---------------|-----------|")
    if sols:
        mc_p, mb_p = sols[0]
        p(f"| m_c | {mc_p:.3f} | {M['c']:.3f} | {(mc_p-M['c'])/M['c']*100:+.3f}% |")
        p(f"| m_b | {mb_p:.3f} | {M['b']:.3f} | {(mb_p-M['b'])/M['b']*100:+.3f}% |")

    p()
    p("### Table 7: Seed vs Physical Mass Comparison")
    p()
    p("Seed masses are computed from M0_seed = v0_seed^2 and delta_seed.")
    p("The nonzero seed masses are NOT equal to the physical masses because")
    p("f_k(delta_seed) differs from f_k(delta_full).")
    p()
    p("| Triple | Mass | Seed (MeV) | Physical (MeV) | seed/physical |")
    p("|--------|------|-----------|---------------|---------------|")
    for tr in triples:
        r = results[tr["label"]]
        for k in range(3):
            phys_name = r["names"][r["perm"][k]]
            m_phys = r["masses"][r["perm"][k]]
            m_seed_k = r["m_seed"][k]
            if k == r["zero_slot"]:
                p(f"| {tr['name']} | {phys_name} | 0 | {m_phys:.4f} | — |")
            else:
                p(f"| {tr['name']} | {phys_name} | {m_seed_k:.4f} | {m_phys:.4f} | {m_seed_k/m_phys:.6f} |")

    p()
    p("---")
    p("All masses in MeV.  Exact Koide condition: Q = 2/3 = 0.666...")
    p()

    # ── Appendix ──────────────────────────────────────────────────────────────

    p()
    p("## Appendix: Seed Angle Reference")
    p()
    p("Principal seed angles delta_seed = 3*pi/4 - 2*pi*k/3 (mod 2*pi):")
    for k in range(3):
        ds = (3*np.pi/4 - 2*np.pi*k/3) % (2*np.pi)
        fk = 1 + np.sqrt(2)*np.cos(2*np.pi*k/3 + ds)
        p(f"  k={k}: delta_seed = {ds:.8f} rad = {ds/np.pi:.8f} pi = {np.degrees(ds):.4f} deg,  f_k = {fk:.2e}")
    p()
    p("Second branch delta_seed = 5*pi/4 - 2*pi*k/3 (mod 2*pi):")
    for k in range(3):
        ds = (5*np.pi/4 - 2*np.pi*k/3) % (2*np.pi)
        fk = 1 + np.sqrt(2)*np.cos(2*np.pi*k/3 + ds)
        p(f"  k={k}: delta_seed = {ds:.8f} rad = {ds/np.pi:.8f} pi = {np.degrees(ds):.4f} deg,  f_k = {fk:.2e}")
    p()
    p("Note: delta_seed/pi values are exact fractions:")
    p("  Principal branch: 3/4, 3/4-2/3=1/12, 3/4-4/3=-7/12 = 5/12 (mod 1) ... checking:")
    for k in range(3):
        frac = (3/4 - 2*k/3) % 1
        p(f"  k={k}: {frac} pi  (3/4 - {k}*2/3 mod 1 = {frac:.10f})")
    p()
    p("  Second branch: 5/4, 5/4-2/3=7/12, 5/4-4/3=-1/12 = 11/12 ... checking:")
    for k in range(3):
        frac = (5/4 - 2*k/3) % 1
        p(f"  k={k}: {frac} pi  (5/4 - {k}*2/3 mod 1 = {frac:.10f})")

    text = "".join(out)
    md_path = "/home/codexssh/phys3/results/staged_koide.md"
    with open(md_path, "w") as fh:
        fh.write(text)
    print(f"\n[Written to {md_path}]")


if __name__ == "__main__":
    main()
