import numpy as np

def Q_signed(s1, s2, s3):
    return (s1 + s2 + s3)**2 / (s1**2 + s2**2 + s3**2)

def Q_unsigned(m1, m2, m3):
    return (np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3))**2 / (m1 + m2 + m3)

# =================================================================
# PDG REFERENCE VALUES
# =================================================================
m_t = 172.69
m_b = 4.18
m_c_pdg = 1.27
m_s_pdg = 0.0934
m_d_pdg = 0.00467
m_u_pdg = 0.00216
m_e = 0.00051099895
m_mu = 0.1056583755
m_tau = 1.77686

print("=" * 90)
print("  COMPLETE RESULTS: KOIDE Q=3/2 CHAINS AND LEPTON-SEEDED PREDICTIONS")
print("=" * 90)

# =================================================================
# CHAIN A
# =================================================================
print("\n" + "=" * 90)
print("  CHAIN A: DESCENDING FROM (m_t = 172.69 GeV, m_b = 4.18 GeV)")
print("=" * 90)

def solve_and_report(m1, m2, step_label, m1_name, m2_name, pdg_targets):
    """Solve Q(m1, m2, m3) = 3/2 for m3, report all branches."""
    s1 = np.sqrt(m1)
    s2 = np.sqrt(m2)
    S = s1 + s2
    ratio = 6 * s1 * s2 / S**2
    inner = np.sqrt(3 + ratio)

    # Two physically meaningful branches (outer_sign = +1)
    solutions = []
    for inner_sign, label in [(+1, "(+,+)"), (-1, "(+,-)")]:
        s3 = S * (2 - inner_sign * inner)
        m3 = s3**2
        q = Q_signed(s1, s2, s3)
        solutions.append((s3, m3, q, label))

    # Also include the outer_sign=-1 branches
    for inner_sign, label in [(+1, "(-,+)"), (-1, "(-,-)")]:
        s3 = -S * (2 - inner_sign * inner)
        m3 = s3**2
        q = Q_signed(s1, s2, s3)
        solutions.append((s3, m3, q, label))

    print("\n  {}".format(step_label))
    print("  Input: {}={:.8f} GeV, {}={:.8f} GeV".format(m1_name, m1, m2_name, m2))
    print("  {:>10s}  {:>16s}  {:>16s}  {:>14s}".format(
        "Branch", "sqrt(m3)", "m3 (GeV)", "|Q_s - 3/2|"))
    print("  " + "-" * 62)

    for s3, m3, q, label in solutions:
        resid = abs(q - 1.5)
        targets_str = ""
        for tname, tval in pdg_targets:
            dev = abs(m3 - tval) / tval * 100
            targets_str += "  dev({})={:.2f}%".format(tname, dev)
        print("  {:>10s}  {:>+16.10f}  {:>16.10f}  {:>14.2e}{}".format(
            label, s3, m3, resid, targets_str))

    return solutions

# Step 1
sols_1 = solve_and_report(m_b, m_t,
    "STEP 1: Q(m_b, m_?, m_t) = 3/2  -->  m_? ~ m_c",
    "m_b", "m_t", [("m_c", m_c_pdg)])

# Best match to m_c: branch (+,+), m = 1.3570 GeV
m_c_pred = sols_1[0][1]  # first solution
s_c_pred = sols_1[0][0]
print("\n  ==> Selected: m_c_pred = {:.8f} GeV (sqrt = {:+.8f})".format(m_c_pred, s_c_pred))
print("      PDG m_c = {:.4f} GeV, deviation = {:.4f}%".format(
    m_c_pdg, abs(m_c_pred - m_c_pdg) / m_c_pdg * 100))

# Step 2
sols_2 = solve_and_report(m_c_pred, m_b,
    "STEP 2: Q(m_c_pred, m_?, m_b) = 3/2  -->  m_? ~ m_s",
    "m_c_pred", "m_b", [("m_s", m_s_pdg)])

# Branch (+,+) naturally gives negative sqrt!
m_s_pred = sols_2[0][1]
s_s_pred = sols_2[0][0]  # This is negative
print("\n  ==> Selected: m_s_pred = {:.8f} GeV = {:.4f} MeV (sqrt = {:+.10f})".format(
    m_s_pred, m_s_pred * 1000, s_s_pred))
print("      sqrt(m_s) is NEGATIVE (inner sqrt > 2, so (2 - inner) < 0)")
print("      PDG m_s = {:.4f} MeV, deviation = {:.4f}%".format(
    m_s_pdg * 1000, abs(m_s_pred - m_s_pdg) / m_s_pdg * 100))

# Step 3
sols_3 = solve_and_report(m_s_pred, m_c_pred,
    "STEP 3: Q(m_s_pred, m_?, m_c_pred) = 3/2  -->  m_? ~ m_d or m_u?",
    "m_s_pred", "m_c_pred", [("m_d", m_d_pdg), ("m_u", m_u_pdg)])

m_light_pred = sols_3[0][1]
s_light_pred = sols_3[0][0]
print("\n  ==> Branch (+,+): m = {:.6f} MeV".format(m_light_pred * 1000))
print("      This is {:.4f} MeV, much smaller than m_d = {:.2f} MeV or m_u = {:.2f} MeV".format(
    m_light_pred * 1000, m_d_pdg * 1000, m_u_pdg * 1000))
print("      Deviation from m_d: {:.2f}%, from m_u: {:.2f}%".format(
    abs(m_light_pred - m_d_pdg) / m_d_pdg * 100,
    abs(m_light_pred - m_u_pdg) / m_u_pdg * 100))
print("      NOTE: The chain does not reach m_d or m_u in a single step.")

# Step 4
sols_4 = solve_and_report(m_light_pred, m_s_pred,
    "STEP 4: Q(m_step3, m_?, m_s_pred) = 3/2  -->  continuation",
    "m_step3", "m_s_pred", [("m_d", m_d_pdg), ("m_c", m_c_pdg)])

print("\n  ==> Branch (+,+): m = {:.6f} MeV (~ {:.2f} MeV)".format(
    sols_4[0][1] * 1000, sols_4[0][1] * 1000))
print("      Branch (+,-): m = {:.6f} GeV -- echoes back to m_c_pred!".format(sols_4[1][1]))

# =================================================================
# CHAIN B
# =================================================================
print("\n\n" + "=" * 90)
print("  CHAIN B: LEPTON-SEEDED ASCENDING CHAIN")
print("=" * 90)

s_e = np.sqrt(m_e)
s_mu = np.sqrt(m_mu)
s_tau = np.sqrt(m_tau)
S_lep = s_e + s_mu + s_tau
M0_analytic = (S_lep / 3.0)**2
sM0 = np.sqrt(M0_analytic)

Q_lep = Q_unsigned(m_e, m_mu, m_tau)

print("\n  Lepton masses and Koide check:")
print("    m_e   = {:.11f} GeV  ({:.8f} MeV)".format(m_e, m_e * 1000))
print("    m_mu  = {:.10f} GeV  ({:.6f} MeV)".format(m_mu, m_mu * 1000))
print("    m_tau = {:.5f} GeV  ({:.2f} MeV)".format(m_tau, m_tau * 1000))
print("    Q(e, mu, tau) = {:.12f}  |Q - 3/2| = {:.6e}".format(Q_lep, abs(Q_lep - 1.5)))

# Parametrization fit
# Analytical: M0 from sum, delta0 from electron
cos_d0 = (s_e / sM0 - 1) / np.sqrt(2)
delta0_analytic = np.arccos(cos_d0)

print("\n  Parametrization: sqrt(m_k) = sqrt(M0) * (1 + sqrt(2)*cos(2*pi*k/3 + delta0))")
print("    M0 (analytic, from sum)     = {:.12f} GeV".format(M0_analytic))
print("    delta0 (analytic, from m_e) = {:.12f} rad = {:.8f} deg".format(
    delta0_analytic, np.degrees(delta0_analytic)))

# Least-squares fit
from scipy.optimize import minimize
def cost(params):
    M, d = params
    sM = np.sqrt(abs(M))
    return sum((sM * (1 + np.sqrt(2)*np.cos(2*np.pi*k/3 + d)) - np.sqrt(m))**2
               for k, m in enumerate([m_e, m_mu, m_tau]))

best = None
best_f = np.inf
for d0 in np.linspace(0, 2*np.pi, 30):
    r = minimize(cost, [M0_analytic, d0], method='Nelder-Mead',
                 options={'xatol': 1e-14, 'fatol': 1e-20, 'maxiter': 100000})
    if r.fun < best_f:
        best_f = r.fun
        best = r

M0 = best.x[0]
delta0 = best.x[1]
sM0 = np.sqrt(M0)

print("\n    M0 (least-squares fit)      = {:.12f} GeV".format(M0))
print("    delta0 (least-squares fit)  = {:.12f} rad = {:.8f} deg".format(
    delta0, np.degrees(delta0)))
print("    Fit SSE = {:.4e}".format(best_f))

# Angle analysis
d0_mod = delta0 % (2 * np.pi / 3)
print("\n    delta0 mod (2*pi/3) = {:.12f} rad".format(d0_mod))
print("    2/9                 = {:.12f}".format(2.0 / 9.0))
print("    Difference          = {:.6e} rad".format(abs(d0_mod - 2.0 / 9.0)))
print("    ==> delta0 = 2*pi/3 + 2/9 to high precision!")

print("\n  Lepton fit quality:")
print("  {:>10s}  {:>14s}  {:>14s}  {:>10s}".format("Particle", "Pred (MeV)", "PDG (MeV)", "Dev %"))
print("  " + "-" * 52)
for k, (name, m_pdg) in enumerate(zip(['e', 'mu', 'tau'], [m_e, m_mu, m_tau])):
    s_pred = sM0 * (1 + np.sqrt(2) * np.cos(2 * np.pi * k / 3 + delta0))
    m_pred = s_pred**2
    dev = abs(m_pred - m_pdg) / m_pdg * 100
    print("  {:>10s}  {:>14.8f}  {:>14.8f}  {:>10.6f}".format(
        name, m_pred * 1000, m_pdg * 1000, dev))

# Second triplet
M1 = 3 * M0
delta1 = 3 * delta0
sM1 = np.sqrt(M1)

print("\n  Second triplet: M1 = 3*M0 = {:.10f} GeV, delta1 = 3*delta0 = {:.10f} rad".format(M1, delta1))
print("  delta1 mod 2*pi = {:.10f} rad = {:.6f} deg".format(delta1 % (2*np.pi), np.degrees(delta1 % (2*np.pi))))

quark_sqrts = []
quark_masses = []
for k in range(3):
    s = sM1 * (1 + np.sqrt(2) * np.cos(2 * np.pi * k / 3 + delta1))
    quark_sqrts.append(s)
    quark_masses.append(s**2)

# Sort by mass
idx_sorted = np.argsort(quark_masses)
pdg_q = [m_s_pdg, m_c_pdg, m_b]
pdg_qn = ['s', 'c', 'b']

print("\n  Quark predictions from second triplet:")
print("  {:>5s}  {:>3s}  {:>14s}  {:>14s}  {:>14s}  {:>10s}".format(
    "k", "ID", "sqrt(m)", "Pred (GeV)", "PDG (GeV)", "Dev %"))
print("  " + "-" * 65)

for rank, (iq, (pname, pval)) in enumerate(zip(idx_sorted, zip(pdg_qn, pdg_q))):
    s = quark_sqrts[iq]
    m = quark_masses[iq]
    dev = abs(m - pval) / pval * 100
    print("  {:>5d}  {:>3s}  {:>+14.10f}  {:>14.8f}  {:>14.8f}  {:>10.4f}".format(
        iq, pname, s, m, pval, dev))

# Q check
Q_s = Q_signed(quark_sqrts[0], quark_sqrts[1], quark_sqrts[2])
Q_u = Q_unsigned(quark_masses[0], quark_masses[1], quark_masses[2])
print("\n  Q_signed(s, c, b) = {:.15f}  |res| = {:.2e}".format(Q_s, abs(Q_s - 1.5)))
print("  Q_unsigned(s, c, b) = {:.15f}  |res| = {:.2e}".format(Q_u, abs(Q_u - 1.5)))
print("  Note: Q_signed = 3/2 by construction. Q_unsigned != 3/2 because sqrt(m_s) < 0.")


# =================================================================
# UNIFIED COMPARISON TABLE
# =================================================================
print("\n\n" + "=" * 90)
print("  UNIFIED COMPARISON: ALL PREDICTIONS vs PDG")
print("=" * 90)

print("\n  {:>12s}  {:>18s}  {:>12s}  {:>10s}  {:>8s}  {:>14s}".format(
    "Particle", "Predicted", "PDG", "Dev %", "Source", "|Q_s - 3/2|"))
print("  " + "-" * 80)

# Chain A predictions
q1_check = Q_signed(np.sqrt(m_b), np.sqrt(m_t), s_c_pred)
print("  {:>12s}  {:>14.8f} GeV  {:>8.4f} GeV  {:>10.4f}  {:>8s}  {:>14.2e}".format(
    "charm (c)", m_c_pred, m_c_pdg, abs(m_c_pred-m_c_pdg)/m_c_pdg*100, "Chain A", abs(q1_check-1.5)))

q2_check = Q_signed(np.sqrt(m_c_pred), np.sqrt(m_b), s_s_pred)
print("  {:>12s}  {:>11.4f} MeV  {:>8.1f} MeV  {:>10.4f}  {:>8s}  {:>14.2e}".format(
    "strange (s)", m_s_pred*1000, m_s_pdg*1000, abs(m_s_pred-m_s_pdg)/m_s_pdg*100, "Chain A", abs(q2_check-1.5)))

q3_check = Q_signed(np.sqrt(m_s_pred), np.sqrt(m_c_pred), s_light_pred)
print("  {:>12s}  {:>11.6f} MeV  {:>8.2f} MeV  {:>10.2f}  {:>8s}  {:>14.2e}".format(
    "light (?)", m_light_pred*1000, m_d_pdg*1000, abs(m_light_pred-m_d_pdg)/m_d_pdg*100, "Chain A", abs(q3_check-1.5)))

print()
# Chain B predictions
for iq, (pname, pval) in zip(idx_sorted, zip(pdg_qn, pdg_q)):
    m = quark_masses[iq]
    dev = abs(m - pval)/pval*100
    if pname == 's':
        print("  {:>12s}  {:>11.4f} MeV  {:>8.1f} MeV  {:>10.4f}  {:>8s}  {:>14.2e}".format(
            "strange (s)", m*1000, pval*1000, dev, "Chain B", abs(Q_s-1.5)))
    else:
        print("  {:>12s}  {:>14.8f} GeV  {:>8.4f} GeV  {:>10.4f}  {:>8s}  {:>14.2e}".format(
            pname + " quark", m, pval, dev, "Chain B", abs(Q_s-1.5)))


# =================================================================
# BACK-SUBSTITUTION RESIDUALS SUMMARY
# =================================================================
print("\n\n" + "=" * 90)
print("  BACK-SUBSTITUTION RESIDUALS (Self-Check)")
print("=" * 90)
print("  All residuals |Q_signed - 3/2| should be < 10^-8.")
print()

residuals = [
    ("Chain A Step 1", "Q(b, c_pred, t)", abs(q1_check - 1.5)),
    ("Chain A Step 2", "Q(c_pred, s_pred, b)", abs(q2_check - 1.5)),
    ("Chain A Step 3", "Q(s_pred, light, c_pred)", abs(q3_check - 1.5)),
    ("Chain A Step 4", "Q(light, 5.33MeV, s_pred)", abs(Q_signed(np.sqrt(m_light_pred), np.sqrt(m_s_pred), sols_4[0][0]) - 1.5)),
    ("Chain B leptons", "Q(e, mu, tau)", abs(Q_lep - 1.5)),
    ("Chain B quarks", "Q_signed(s, c, b)", abs(Q_s - 1.5)),
]

for name, desc, resid in residuals:
    status = "PASS" if resid < 1e-8 else "NOTE"
    print("  {:>20s}  {:>30s}  |res| = {:.2e}  {}".format(name, desc, resid, status))

print("\n  Note: Chain B lepton Q residual reflects the empirical deviation of")
print("  charged lepton masses from exact Koide relation (|Q-3/2| ~ 1.4e-5).")
print("  This is not a computational error but rather the well-known near-miss.")

print("\n" + "=" * 90)
print("  COMPUTATION COMPLETE")
print("=" * 90)
