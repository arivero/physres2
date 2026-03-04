import numpy as np
from scipy.optimize import minimize

# ============================================================
# Helper functions
# ============================================================

def Q_unsigned(m1, m2, m3):
    """Standard Koide Q with unsigned (positive) square roots."""
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return s**2 / (m1 + m2 + m3)

def Q_signed(s1, s2, s3):
    """Koide Q with signed square roots.  m_i = s_i^2."""
    return (s1 + s2 + s3)**2 / (s1**2 + s2**2 + s3**2)

def solve_m3_branches(m1, m2):
    """
    Given Q(m1, m2, m3) = 3/2 with all positive sqrt, solve for sqrt(m3).

    The formula:
      sqrt(m3) = (sqrt(m1)+sqrt(m2)) * [2 - sqrt(3 + 6*sqrt(m1*m2)/(sqrt(m1)+sqrt(m2))^2)]

    Four branches from sign choices:
      outer_sign * (sqrt(m1)+sqrt(m2)) * (2 - inner_sign * sqrt(...))

    Returns list of dicts with outer, inner signs, sqrt_m3, m3.
    """
    s1 = np.sqrt(m1)
    s2 = np.sqrt(m2)
    S = s1 + s2
    ratio = 6.0 * s1 * s2 / S**2   # = 6*sqrt(m1*m2)/(sqrt(m1)+sqrt(m2))^2
    inner_arg = 3.0 + ratio
    inner_sqrt = np.sqrt(inner_arg)
    branches = []
    for outer_sign in [+1, -1]:
        for inner_sign in [+1, -1]:
            sqrt_m3 = outer_sign * S * (2.0 - inner_sign * inner_sqrt)
            m3 = sqrt_m3**2
            branches.append({
                'outer': outer_sign,
                'inner': inner_sign,
                'sqrt_m3': sqrt_m3,
                'm3': m3
            })
    return branches

def print_branches(branches, label=""):
    print("\n" + "="*80)
    print("  " + label)
    print("="*80)
    print("  {:>8s}  {:>6s}  {:>6s}  {:>18s}  {:>18s}".format(
        "Branch", "outer", "inner", "sqrt(m3)", "m3 (GeV)"))
    print("  " + "-"*70)
    for i, b in enumerate(branches):
        print("  {:>8d}  {:>+6d}  {:>+6d}  {:>18.10f}  {:>18.10f}".format(
            i+1, b['outer'], b['inner'], b['sqrt_m3'], b['m3']))


# ============================================================
# CHAIN A: Descending from top, bottom
# ============================================================
print("\n" + "#"*80)
print("#  CHAIN A: Descending chain starting from m_t, m_b")
print("#"*80)

m_t = 172.69   # GeV
m_b = 4.18     # GeV

# PDG reference masses (GeV)
m_c_pdg = 1.27
m_s_pdg = 0.0934   # 93.4 MeV
m_d_pdg = 0.00467  # 4.67 MeV
m_u_pdg = 0.00216  # 2.16 MeV

# ----------------------------------------------------------
# Step 1: Q(m_b, m_?, m_t) = 3/2   =>  solve for m_?
# The formula takes sqrt(m1)=sqrt(m_b), sqrt(m2)=sqrt(m_t)
# and solves for sqrt(m3) = sqrt(m_?)
# ----------------------------------------------------------
branches_1 = solve_m3_branches(m_b, m_t)
print_branches(branches_1, "Step 1: Q(m_b, m_?, m_t) = 3/2, solving for m_?")

# Back-substitution: use signed Q with s1=sqrt(m_b), s2=sqrt(m_t), s3=sqrt_m3
s_b = np.sqrt(m_b)
s_t = np.sqrt(m_t)
print("\n  Back-substitution (signed Q, with s_b={:.8f}, s_t={:.8f}):".format(s_b, s_t))
for i, b in enumerate(branches_1):
    q_signed = Q_signed(s_b, s_t, b['sqrt_m3'])
    q_unsigned = Q_unsigned(m_b, m_t, b['m3'])
    resid_s = abs(q_signed - 1.5)
    resid_u = abs(q_unsigned - 1.5)
    dev = abs(b['m3'] - m_c_pdg) / m_c_pdg * 100
    print("    Br {}: m_? = {:.8f} GeV, Q_signed = {:.15f}, |res| = {:.2e}, "
          "Q_unsigned = {:.15f}, |res| = {:.2e}, dev(m_c) = {:.2f}%".format(
        i+1, b['m3'], q_signed, resid_s, q_unsigned, resid_u, dev))

diffs = [abs(b['m3'] - m_c_pdg) for b in branches_1]
best_1 = np.argmin(diffs)
m_c_pred = branches_1[best_1]['m3']
sqrt_m_c = branches_1[best_1]['sqrt_m3']
print("\n  >>> Best match to m_c = {:.4f} GeV: Branch {}, m_? = {:.8f} GeV "
      "(sqrt = {:+.8f})".format(m_c_pdg, best_1+1, m_c_pred, sqrt_m_c))
print("      Deviation: {:.4f}%".format(abs(m_c_pred - m_c_pdg)/m_c_pdg*100))

# ----------------------------------------------------------
# Step 2: Q(m_c_pred, m_?, m_b) = 3/2
# ----------------------------------------------------------
branches_2 = solve_m3_branches(m_c_pred, m_b)
print_branches(branches_2, "Step 2: Q(m_c={:.8f}, m_?, m_b={:.2f}) = 3/2".format(m_c_pred, m_b))

s_c = np.sqrt(m_c_pred)
print("\n  Back-substitution (signed Q, with s_c={:.8f}, s_b={:.8f}):".format(s_c, s_b))
for i, b in enumerate(branches_2):
    q_signed = Q_signed(s_c, s_b, b['sqrt_m3'])
    q_unsigned = Q_unsigned(m_c_pred, m_b, b['m3'])
    resid_s = abs(q_signed - 1.5)
    resid_u = abs(q_unsigned - 1.5)
    dev_s = abs(b['m3'] - m_s_pdg) / m_s_pdg * 100
    print("    Br {}: sqrt = {:+.10f}, m_? = {:.8f} GeV, Q_signed = {:.15f}, |res| = {:.2e}, "
          "Q_unsigned = {:.15f}, |res| = {:.2e}, dev(m_s) = {:.2f}%".format(
        i+1, b['sqrt_m3'], b['m3'], q_signed, resid_s, q_unsigned, resid_u, dev_s))

# The branch with negative sqrt nearest to m_s
diffs_s = [abs(b['m3'] - m_s_pdg) for b in branches_2]
best_2 = np.argmin(diffs_s)
m_s_pred = branches_2[best_2]['m3']
sqrt_m_s = branches_2[best_2]['sqrt_m3']
print("\n  >>> Best match to m_s = {:.4f} GeV: Branch {}, m_? = {:.8f} GeV "
      "(sqrt = {:+.10f})".format(m_s_pdg, best_2+1, m_s_pred, sqrt_m_s))
print("      sqrt is {}".format('NEGATIVE' if sqrt_m_s < 0 else 'positive'))
print("      Deviation: {:.4f}%".format(abs(m_s_pred - m_s_pdg)/m_s_pdg*100))

# Verify: Q_signed(s_c, sqrt_m_s, s_b) should be 3/2
q_check2 = Q_signed(s_c, sqrt_m_s, s_b)
print("      Q_signed(s_c, s_s, s_b) = {:.15f}, |res| = {:.2e}".format(q_check2, abs(q_check2-1.5)))

# ----------------------------------------------------------
# Step 3: Q(m_s_pred, m_?, m_c_pred) = 3/2
# ----------------------------------------------------------
branches_3 = solve_m3_branches(m_s_pred, m_c_pred)
print_branches(branches_3, "Step 3: Q(m_s={:.8f}, m_?, m_c={:.8f}) = 3/2".format(m_s_pred, m_c_pred))

s_s = np.sqrt(m_s_pred)
print("\n  Back-substitution (signed Q, with s_s={:.8f}, s_c={:.8f}):".format(s_s, s_c))
for i, b in enumerate(branches_3):
    q_signed = Q_signed(s_s, s_c, b['sqrt_m3'])
    q_unsigned = Q_unsigned(m_s_pred, m_c_pred, b['m3'])
    resid_s = abs(q_signed - 1.5)
    resid_u = abs(q_unsigned - 1.5)
    dev_d = abs(b['m3'] - m_d_pdg) / m_d_pdg * 100
    dev_u = abs(b['m3'] - m_u_pdg) / m_u_pdg * 100
    print("    Br {}: sqrt = {:+.12f}, m = {:.10e} GeV = {:.6f} MeV, "
          "Q_s = {:.15f} |r| = {:.2e}, Q_u = {:.15f} |r| = {:.2e}, "
          "dev_d={:.2f}%, dev_u={:.2f}%".format(
        i+1, b['sqrt_m3'], b['m3'], b['m3']*1000,
        q_signed, resid_s, q_unsigned, resid_u, dev_d, dev_u))

diffs_d = [abs(b['m3'] - m_d_pdg) for b in branches_3]
diffs_u = [abs(b['m3'] - m_u_pdg) for b in branches_3]
best_3d = np.argmin(diffs_d)
best_3u = np.argmin(diffs_u)
m_d_pred = branches_3[best_3d]['m3']
sqrt_m_d = branches_3[best_3d]['sqrt_m3']
m_u_pred_val = branches_3[best_3u]['m3']
sqrt_m_u_val = branches_3[best_3u]['sqrt_m3']

print("\n  >>> Nearest to m_d = {:.5f} GeV ({:.2f} MeV):".format(m_d_pdg, m_d_pdg*1000))
print("      Branch {}, m = {:.10e} GeV = {:.6f} MeV (sqrt = {:+.10f})".format(
    best_3d+1, m_d_pred, m_d_pred*1000, sqrt_m_d))
print("      Deviation: {:.4f}%".format(abs(m_d_pred - m_d_pdg)/m_d_pdg*100))

print("  >>> Nearest to m_u = {:.5f} GeV ({:.2f} MeV):".format(m_u_pdg, m_u_pdg*1000))
print("      Branch {}, m = {:.10e} GeV = {:.6f} MeV (sqrt = {:+.10f})".format(
    best_3u+1, m_u_pred_val, m_u_pred_val*1000, sqrt_m_u_val))
print("      Deviation: {:.4f}%".format(abs(m_u_pred_val - m_u_pdg)/m_u_pdg*100))

# ----------------------------------------------------------
# Step 4: Continue one more step
# ----------------------------------------------------------
m_step3 = m_d_pred
sqrt_step3 = sqrt_m_d
branches_4 = solve_m3_branches(m_step3, m_s_pred)
print_branches(branches_4, "Step 4: Q(m_?3={:.10f}, m_?, m_s={:.8f}) = 3/2".format(m_step3, m_s_pred))

s_3 = np.sqrt(m_step3)
print("\n  Back-substitution (signed Q, with s_3={:.10f}, s_s={:.8f}):".format(s_3, s_s))
for i, b in enumerate(branches_4):
    q_signed = Q_signed(s_3, s_s, b['sqrt_m3'])
    q_unsigned = Q_unsigned(m_step3, m_s_pred, b['m3'])
    resid_s = abs(q_signed - 1.5)
    resid_u = abs(q_unsigned - 1.5)
    print("    Br {}: sqrt = {:+.12f}, m = {:.10e} GeV = {:.8f} MeV, "
          "Q_s = {:.15f} |r| = {:.2e}, Q_u = {:.15f} |r| = {:.2e}".format(
        i+1, b['sqrt_m3'], b['m3'], b['m3']*1000,
        q_signed, resid_s, q_unsigned, resid_u))


# ============================================================
# CHAIN B: Ascending, lepton-seeded
# ============================================================
print("\n\n" + "#"*80)
print("#  CHAIN B: Lepton-seeded ascending chain")
print("#"*80)

# PDG charged lepton masses (GeV) - more precise values
m_e = 0.00051099895000   # 0.51099895 MeV
m_mu = 0.1056583755      # 105.6583755 MeV
m_tau = 1.77686           # 1776.86 MeV

print("\n  Charged lepton masses (GeV):")
print("    m_e   = {:.11f}".format(m_e))
print("    m_mu  = {:.10f}".format(m_mu))
print("    m_tau = {:.5f}".format(m_tau))

Q_lep = Q_unsigned(m_e, m_mu, m_tau)
print("\n  Q_unsigned(m_e, m_mu, m_tau) = {:.12f}".format(Q_lep))
print("  |Q - 3/2| = {:.6e}".format(abs(Q_lep - 1.5)))

# Parametrization: sqrt(m_k) = sqrt(M0) * (1 + sqrt(2)*cos(2*pi*k/3 + delta0))
# Sum of sqrt(m_k) over k=0,1,2 = 3*sqrt(M0) since sum of cos(2*pi*k/3 + d) = 0
s_e = np.sqrt(m_e)
s_mu = np.sqrt(m_mu)
s_tau = np.sqrt(m_tau)
S_lep = s_e + s_mu + s_tau
M0 = (S_lep / 3.0)**2
sM0 = np.sqrt(M0)

print("\n  sqrt(m_e)   = {:.12f}".format(s_e))
print("  sqrt(m_mu)  = {:.12f}".format(s_mu))
print("  sqrt(m_tau) = {:.12f}".format(s_tau))
print("  Sum         = {:.12f}".format(S_lep))
print("  M0 = (Sum/3)^2 = {:.12f} GeV".format(M0))
print("  sqrt(M0)    = {:.12f}".format(sM0))

# Extract delta0 from the electron equation
# sqrt(m_e)/sqrt(M0) = 1 + sqrt(2)*cos(delta0)
# cos(delta0) = (sqrt(m_e)/sqrt(M0) - 1)/sqrt(2)
r_e = s_e / sM0
cos_d0 = (r_e - 1.0) / np.sqrt(2)
print("\n  r_e = sqrt(m_e)/sqrt(M0) = {:.12f}".format(r_e))
print("  cos(delta0) = {:.12f}".format(cos_d0))

# delta0 from arccos (principal value)
delta0 = np.arccos(cos_d0)
print("  delta0 = arccos(cos_d0) = {:.12f} rad = {:.8f} deg".format(delta0, np.degrees(delta0)))

# Verify all three
print("\n  Verification with delta0 = {:.10f}:".format(delta0))
fit_ok = True
for k, (name, m_pdg) in enumerate(zip(['e', 'mu', 'tau'], [m_e, m_mu, m_tau])):
    sqrt_pred = sM0 * (1.0 + np.sqrt(2) * np.cos(2*np.pi*k/3.0 + delta0))
    m_pred = sqrt_pred**2
    actual = np.sqrt(m_pdg)
    diff = abs(sqrt_pred - actual)
    dev = abs(m_pred - m_pdg)/m_pdg * 100
    print("    k={} ({:>3s}): sqrt_pred = {:.12f}, sqrt_pdg = {:.12f}, diff = {:.6e}, "
          "m_pred = {:.12f}, dev = {:.8f}%".format(k, name, sqrt_pred, actual, diff, m_pred, dev))

# The fit is not perfect because Q(leptons) != exactly 3/2.
# The system is overconstrained (3 eqns, 2 unknowns), and Q != 3/2 exactly
# means there is a small residual.
#
# We can do a least-squares fit.
print("\n  --- Least-squares fit for M0, delta0 ---")

def cost_lep(params):
    M, d = params
    sM = np.sqrt(abs(M))
    err = 0.0
    for k, m in enumerate([m_e, m_mu, m_tau]):
        pred = sM * (1 + np.sqrt(2)*np.cos(2*np.pi*k/3 + d))
        err += (pred - np.sqrt(m))**2
    return err

best_res = None
best_val = np.inf
for d_init in np.linspace(0, 2*np.pi, 30):
    res = minimize(cost_lep, [M0, d_init], method='Nelder-Mead',
                   options={'xatol':1e-14, 'fatol':1e-20, 'maxiter':100000})
    if res.fun < best_val:
        best_val = res.fun
        best_res = res

M0_fit = best_res.x[0]
delta0_fit = best_res.x[1]
sM0_fit = np.sqrt(M0_fit)
print("  M0_fit     = {:.12f} GeV  ({:.8f} MeV)".format(M0_fit, M0_fit*1000))
print("  delta0_fit = {:.12f} rad = {:.8f} deg".format(delta0_fit, np.degrees(delta0_fit)))
print("  Fit SSE    = {:.4e}".format(best_val))

print("\n  Fit verification:")
for k, (name, m_pdg) in enumerate(zip(['e', 'mu', 'tau'], [m_e, m_mu, m_tau])):
    sqrt_pred = sM0_fit * (1.0 + np.sqrt(2) * np.cos(2*np.pi*k/3.0 + delta0_fit))
    m_pred = sqrt_pred**2
    dev = abs(m_pred - m_pdg)/m_pdg * 100
    print("    k={} ({:>3s}): m_pred = {:.12f} GeV, m_pdg = {:.12f} GeV, dev = {:.8f}%".format(
        k, name, m_pred, m_pdg, dev))

# Also report the analytical delta0 (exact for the M0 derived from the sum)
print("\n  Note: analytical M0 (from sum constraint) = {:.12f}".format(M0))
print("  Note: analytical delta0 (from electron)   = {:.12f}".format(delta0))
print("  These differ slightly from the fit because Q(leptons) = {:.12f} != 3/2 exactly.".format(Q_lep))

# For the second triplet, use the fitted values
print("\n\n  ====== SECOND TRIPLET: M1 = 3*M0, delta1 = 3*delta0 ======")
M1 = 3.0 * M0_fit
delta1 = 3.0 * delta0_fit
sM1 = np.sqrt(M1)

print("  M1     = 3 * M0_fit = {:.12f} GeV".format(M1))
print("  delta1 = 3 * delta0_fit = {:.12f} rad = {:.8f} deg".format(delta1, np.degrees(delta1)))
print("  sqrt(M1) = {:.12f}".format(sM1))

quark_pred_masses = []
quark_pred_sqrts = []
print("\n  Predicted masses from second triplet:")
for k in range(3):
    sqrt_pred = sM1 * (1.0 + np.sqrt(2) * np.cos(2*np.pi*k/3.0 + delta1))
    m_pred = sqrt_pred**2
    quark_pred_masses.append(m_pred)
    quark_pred_sqrts.append(sqrt_pred)
    print("    k={}: sqrt(m) = {:+.12f}, m = {:.10f} GeV = {:.6f} MeV".format(
        k, sqrt_pred, m_pred, m_pred*1000))

# Check Q for this triplet
Q_q_unsigned = Q_unsigned(quark_pred_masses[0], quark_pred_masses[1], quark_pred_masses[2])
Q_q_signed = Q_signed(quark_pred_sqrts[0], quark_pred_sqrts[1], quark_pred_sqrts[2])
print("\n  Q_unsigned(triplet) = {:.12f}, |res| = {:.2e}".format(Q_q_unsigned, abs(Q_q_unsigned-1.5)))
print("  Q_signed(triplet)   = {:.12f}, |res| = {:.2e}".format(Q_q_signed, abs(Q_q_signed-1.5)))

# Sort and compare to PDG s, c, b
pred_sorted = sorted(quark_pred_masses)
pdg_quarks = [m_s_pdg, m_c_pdg, m_b]
pdg_qnames = ['s', 'c', 'b']

print("\n  Comparison to PDG quarks (sorted by mass):")
print("  {:>10s}  {:>16s}  {:>12s}  {:>10s}".format("Quark", "Predicted (GeV)", "PDG (GeV)", "Dev %"))
print("  " + "-"*55)
for p, pdg, name in zip(pred_sorted, pdg_quarks, pdg_qnames):
    dev = abs(p - pdg)/pdg * 100
    print("  {:>10s}  {:>16.8f}  {:>12.8f}  {:>10.4f}".format(name, p, pdg, dev))


# ============================================================
# COMPREHENSIVE FINAL TABLE
# ============================================================
print("\n\n" + "="*100)
print("  COMPREHENSIVE FINAL TABLE")
print("="*100)

print("\n  ---- CHAIN A: Descending from (m_t = 172.69 GeV, m_b = 4.18 GeV) ----")
print()

# Recompute everything for the summary
print("  Step 1: Q(m_b, m_?, m_t) = 3/2  -->  m_? (candidate for m_c)")
print("  {:<8s} {:<8s} {:<8s} {:>16s} {:>16s} {:>10s} {:>14s}".format(
    "Branch", "outer", "inner", "sqrt(m_?)", "m_? (GeV)", "Dev(m_c)%", "|Q_s-3/2|"))
print("  " + "-"*85)
for i, b in enumerate(branches_1):
    qs = Q_signed(s_b, s_t, b['sqrt_m3'])
    dev = abs(b['m3'] - m_c_pdg)/m_c_pdg*100
    tag = " <-- BEST" if i == best_1 else ""
    print("  {:<8d} {:<+8d} {:<+8d} {:>16.8f} {:>16.8f} {:>10.2f} {:>14.2e}{}".format(
        i+1, b['outer'], b['inner'], b['sqrt_m3'], b['m3'], dev, abs(qs-1.5), tag))

print("\n  Step 2: Q(m_c_pred, m_?, m_b) = 3/2  -->  m_? (candidate for m_s)")
print("  {:<8s} {:<8s} {:<8s} {:>16s} {:>16s} {:>10s} {:>14s}".format(
    "Branch", "outer", "inner", "sqrt(m_?)", "m_? (GeV)", "Dev(m_s)%", "|Q_s-3/2|"))
print("  " + "-"*85)
for i, b in enumerate(branches_2):
    qs = Q_signed(s_c, s_b, b['sqrt_m3'])
    dev = abs(b['m3'] - m_s_pdg)/m_s_pdg*100
    tag = " <-- BEST (sqrt<0)" if i == best_2 else ""
    print("  {:<8d} {:<+8d} {:<+8d} {:>+16.8f} {:>16.8f} {:>10.2f} {:>14.2e}{}".format(
        i+1, b['outer'], b['inner'], b['sqrt_m3'], b['m3'], dev, abs(qs-1.5), tag))

print("\n  Step 3: Q(m_s_pred, m_?, m_c_pred) = 3/2  -->  m_? (candidate for m_d or m_u)")
print("  {:<8s} {:<8s} {:<8s} {:>16s} {:>16s} {:>10s} {:>10s} {:>14s}".format(
    "Branch", "outer", "inner", "sqrt(m_?)", "m_? (MeV)", "Dev(d)%", "Dev(u)%", "|Q_s-3/2|"))
print("  " + "-"*95)
for i, b in enumerate(branches_3):
    qs = Q_signed(s_s, s_c, b['sqrt_m3'])
    dev_d = abs(b['m3'] - m_d_pdg)/m_d_pdg*100
    dev_u = abs(b['m3'] - m_u_pdg)/m_u_pdg*100
    tag = ""
    if i == best_3d:
        tag += " <-- near d"
    if i == best_3u:
        tag += " <-- near u"
    print("  {:<8d} {:<+8d} {:<+8d} {:>+16.10f} {:>16.8f} {:>10.2f} {:>10.2f} {:>14.2e}{}".format(
        i+1, b['outer'], b['inner'], b['sqrt_m3'], b['m3']*1000, dev_d, dev_u, abs(qs-1.5), tag))

print("\n  Step 4: Q(m_d_pred, m_?, m_s_pred) = 3/2  -->  m_? (continuation)")
print("  {:<8s} {:<8s} {:<8s} {:>16s} {:>16s} {:>14s}".format(
    "Branch", "outer", "inner", "sqrt(m_?)", "m_? (MeV)", "|Q_s-3/2|"))
print("  " + "-"*75)
for i, b in enumerate(branches_4):
    qs = Q_signed(s_3, s_s, b['sqrt_m3'])
    print("  {:<8d} {:<+8d} {:<+8d} {:>+16.10f} {:>16.8f} {:>14.2e}".format(
        i+1, b['outer'], b['inner'], b['sqrt_m3'], b['m3']*1000, abs(qs-1.5)))


print("\n\n  ---- CHAIN B: Lepton-seeded ----")
print("  Fit parameters:")
print("    M0     = {:.12f} GeV".format(M0_fit))
print("    delta0 = {:.12f} rad = {:.8f} deg".format(delta0_fit, np.degrees(delta0_fit)))
print("    M1     = {:.12f} GeV (= 3*M0)".format(M1))
print("    delta1 = {:.12f} rad = {:.8f} deg (= 3*delta0)".format(delta1, np.degrees(delta1)))
print()
print("  Lepton fit:")
print("  {:>10s}  {:>16s}  {:>16s}  {:>10s}".format("Particle", "Pred (GeV)", "PDG (GeV)", "Dev %"))
print("  " + "-"*55)
for k, (name, m_pdg) in enumerate(zip(['e', 'mu', 'tau'], [m_e, m_mu, m_tau])):
    sqrt_pred = sM0_fit * (1.0 + np.sqrt(2) * np.cos(2*np.pi*k/3.0 + delta0_fit))
    m_pred = sqrt_pred**2
    dev = abs(m_pred - m_pdg)/m_pdg * 100
    print("  {:>10s}  {:>16.10f}  {:>16.10f}  {:>10.6f}".format(name, m_pred, m_pdg, dev))

print()
print("  Quark predictions (second triplet):")
print("  {:>10s}  {:>16s}  {:>16s}  {:>10s}".format("Particle", "Pred (GeV)", "PDG (GeV)", "Dev %"))
print("  " + "-"*55)
for p, pdg, name in zip(pred_sorted, pdg_quarks, pdg_qnames):
    dev = abs(p - pdg)/pdg * 100
    print("  {:>10s}  {:>16.8f}  {:>16.8f}  {:>10.4f}".format(name, p, pdg, dev))
print("  Q_signed(s,c,b predicted) = {:.12f}, |res| = {:.2e}".format(Q_q_signed, abs(Q_q_signed-1.5)))
print("  Q_unsigned(s,c,b predicted) = {:.12f}, |res| = {:.2e}".format(Q_q_unsigned, abs(Q_q_unsigned-1.5)))

print("\n  DONE.")
