import numpy as np
from scipy.optimize import minimize

###############################################################################
# TASK 1
###############################################################################
print("=" * 70)
print("TASK 1: The 7 ppm vs 33 ppm question")
print("=" * 70)

m_e = 0.51099895
m_mu = 105.6583755
m_tau = 1776.86
masses = [m_e, m_mu, m_tau]
sq = [np.sqrt(m) for m in masses]
S = sum(sq)
M0 = (S / 3) ** 2
sqM0 = np.sqrt(M0)

print(f"sqrt(m_e)  = {sq[0]:.10f}")
print(f"sqrt(m_mu) = {sq[1]:.10f}")
print(f"sqrt(m_tau)= {sq[2]:.10f}")
print(f"M0 = {M0:.10f} MeV, sqrt(M0) = {sqM0:.10f}")

best_d = None
best_c = 1e30
for d0i in np.linspace(0, 2 * np.pi, 300):
    r = minimize(
        lambda x: sum(
            (sqM0 * (1 + np.sqrt(2) * np.cos(2 * np.pi * k / 3 + x[0])) - sq[k]) ** 2
            for k in range(3)
        ),
        [d0i],
        method="Nelder-Mead",
        options={"xatol": 1e-15, "fatol": 1e-30, "maxiter": 100000},
    )
    if r.fun < best_c:
        best_c = r.fun
        best_d = r.x[0]

d0 = best_d % (2 * np.pi)
per = 2 * np.pi / 3
dm = d0 % per
res = dm - 2.0 / 9.0
n = int(d0 // per)

print(f"\ndelta0 (mod 2pi)   = {d0:.15f} rad")
print(f"Fit SSR            = {best_c:.6e}")

for k in range(3):
    p = sqM0 * (1 + np.sqrt(2) * np.cos(2 * np.pi * k / 3 + d0))
    print(f"  k={k}: pred={p:.10f}, act={sq[k]:.10f}, diff={p - sq[k]:.2e}")

print(f"\ndelta0 mod (2pi/3) = {dm:.15f} rad")
print(f"2/9                = {2.0 / 9.0:.15f} rad")
print(f"Residual           = {res:.15e} rad")
print(f"\nresidual / (2pi/3) = {res / per:.15e}  =>  {res / per * 1e6:.4f} ppm")
print(f"residual / (2/9)   = {res / (2.0 / 9.0):.15e}  =>  {res / (2.0 / 9.0) * 1e6:.4f} ppm")
print(f"\ndelta0 = {n} * (2pi/3) + 2/9 + {res:.15e}")

###############################################################################
# TASK 2
###############################################################################
print("\n" + "=" * 70)
print("TASK 2: Error propagation for Koide Q")
print("=" * 70)

dm_e = 0.00000015
dm_mu = 0.0000023
dm_tau = 0.12


def Q_val(m1, m2, m3):
    Sv = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return Sv ** 2 / (m1 + m2 + m3)


Q0 = Q_val(m_e, m_mu, m_tau)
print(f"Q                  = {Q0:.15f}")
print(f"|Q - 3/2|          = {abs(Q0 - 1.5):.15e}")

Sv = np.sqrt(m_e) + np.sqrt(m_mu) + np.sqrt(m_tau)
Sigma = m_e + m_mu + m_tau
dms = [dm_e, dm_mu, dm_tau]
partials = []
for i, m in enumerate(masses):
    dQdm = Sv * (Sigma / np.sqrt(m) - Sv) / Sigma ** 2
    partials.append(dQdm)

dQ = np.sqrt(sum((p * d) ** 2 for p, d in zip(partials, dms)))
print(f"dQ                 = {dQ:.15e}")
for i, (p, d) in enumerate(zip(partials, dms)):
    print(f"  |dQ/dm_{i}|*dm_{i} = {abs(p * d):.6e}")
nsigma = abs(Q0 - 1.5) / dQ
print(f"|Q - 3/2| / dQ     = {nsigma:.6f} sigma")

###############################################################################
# TASK 3
###############################################################################
print("\n" + "=" * 70)
print("TASK 3: Chain stall value")
print("=" * 70)

m_s = 92.274
m_c = 1356.8
sqrt_s = np.sqrt(m_s)
sqrt_c = np.sqrt(m_c)
print(f"m_s={m_s}, m_c={m_c}")
print(f"sqrt(m_s)={sqrt_s:.10f}, sqrt(m_c)={sqrt_c:.10f}")

task3_results = []
for sign_q, sign_c in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
    A = sqrt_s + sign_c * sqrt_c
    B = m_s + m_c
    b_coeff = -4 * A * sign_q
    c_coeff = -2 * A ** 2 + 3 * B
    disc = b_coeff ** 2 - 4 * c_coeff
    sq_l = "+" if sign_q > 0 else "-"
    sc_l = "+" if sign_c > 0 else "-"
    label = f"({sq_l},{sc_l})"
    print(f"\nBranch {label}: A={A:.6f}, disc={disc:.6f}")
    if disc < 0:
        print("  No real solutions")
        continue
    for pm in [1, -1]:
        x = (-b_coeff + pm * np.sqrt(disc)) / 2
        pm_l = "+" if pm == 1 else "-"
        if x > 0:
            m_q = x ** 2
            Sval = sqrt_s + sign_c * sqrt_c + sign_q * x
            Q_check = Sval ** 2 / (m_s + m_q + m_c)
            print(f"  root({pm_l}): x={x:.10f}, m_?={m_q:.6f} MeV, Q={Q_check:.15f}")
            task3_results.append((label, pm_l, m_q, Q_check))
        elif abs(x) < 1e-10:
            print(f"  root({pm_l}): x~0")
        else:
            print(f"  root({pm_l}): x={x:.6f} negative, skip")

print("\n--- Task 3 Summary ---")
stall_mass = None
for label, pm, m_q, Q in task3_results:
    flag = " <-- STALL" if m_q < 0.1 else ""
    print(f"  {label} root({pm}): m_? = {m_q:.6f} MeV{flag}")
    if m_q < 0.1:
        stall_mass = m_q

###############################################################################
# TASK 4
###############################################################################
print("\n" + "=" * 70)
print("TASK 4: Cyclic echo")
print("=" * 70)

m_stall = stall_mass
m_s2 = 92.274
sqrt_st = np.sqrt(m_stall)
sqrt_s2 = np.sqrt(m_s2)
print(f"m_stall={m_stall:.6f}, m_s={m_s2}")

task4_results = []
for sign_q, sign_s in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
    A = sqrt_st + sign_s * sqrt_s2
    B = m_stall + m_s2
    b_coeff = -4 * A * sign_q
    c_coeff = -2 * A ** 2 + 3 * B
    disc = b_coeff ** 2 - 4 * c_coeff
    if disc < 0:
        continue
    for pm in [1, -1]:
        x = (-b_coeff + pm * np.sqrt(disc)) / 2
        if x > 0:
            m_q = x ** 2
            Sval = sqrt_st + sign_s * sqrt_s2 + sign_q * x
            Q_check = Sval ** 2 / (m_stall + m_q + m_s2)
            sq_l = "+" if sign_q > 0 else "-"
            ss_l = "+" if sign_s > 0 else "-"
            pm_l = "+" if pm == 1 else "-"
            label = f"({sq_l},{ss_l})"
            print(f"  {label} root({pm_l}): m_?={m_q:.6f} MeV, Q={Q_check:.15f}")
            task4_results.append((label, pm_l, m_q))

print("\n--- Check for m_c = 1356.8 ---")
for label, pm, m_q in task4_results:
    if abs(m_q - 1356.8) < 200:
        r4 = m_q - 1356.8
        print(f"  {label} root({pm}): m_?={m_q:.6f}, residual={r4:.6f} MeV, frac={r4/1356.8:.6e}")

###############################################################################
# TASK 5
###############################################################################
print("\n" + "=" * 70)
print("TASK 5: M_W from Casimir formula")
print("=" * 70)

x1 = np.sqrt(3) - 1
x12 = (-3 + np.sqrt(57)) / 8
M_Z = 91.1876
dM_Z = 0.0021
print(f"x(s=1)={x1:.15f}")
print(f"x(s=1/2)={x12:.15f}")

m_sc = M_Z / np.sqrt(x1)
M_W = m_sc * np.sqrt(x12)
ratio = np.sqrt(x12 / x1)
dM_W = dM_Z * ratio
print(f"m={m_sc:.10f} GeV")
print(f"M_W={M_W:.10f} GeV")
print(f"M_W/M_Z={ratio:.15f}")
print(f"dM_W={dM_W:.10f} GeV")

M_W_pdg = 80.3692
dM_W_pdg = 0.0133
diff = M_W - M_W_pdg
sig_comb = np.sqrt(dM_W ** 2 + dM_W_pdg ** 2)
print(f"M_W(PDG)={M_W_pdg} +/- {dM_W_pdg} GeV")
print(f"diff={diff:.6f} GeV")
print(f"combined_sigma={sig_comb:.6f} GeV")
print(f"tension={abs(diff)/sig_comb:.2f} sigma (combined)")
print(f"tension_pred_only={abs(diff)/dM_W:.2f} sigma")
print(f"tension_pdg_only={abs(diff)/dM_W_pdg:.2f} sigma")
