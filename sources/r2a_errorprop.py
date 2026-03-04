#!/usr/bin/env python3
"""R2-A substitute: Error propagation computations"""
import numpy as np

print("=" * 60)
print("TASK: Error propagation on Q for (e, μ, τ)")
print("=" * 60)

# PDG 2023 values with uncertainties (MeV)
m_e = 0.51099895; dm_e = 0.00000015
m_mu = 105.6583755; dm_mu = 0.0000023
m_tau = 1776.86; dm_tau = 0.12

def Q(m1, m2, m3):
    return (np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3))**2 / (m1 + m2 + m3)

Q_central = Q(m_e, m_mu, m_tau)
print(f"Q(e,μ,τ) = {Q_central:.15f}")
print(f"|Q - 3/2| = {abs(Q_central - 1.5):.6e}")

# Numerical partial derivatives
eps = 1e-8
dQ_dme = (Q(m_e*(1+eps), m_mu, m_tau) - Q(m_e*(1-eps), m_mu, m_tau)) / (2*eps*m_e)
dQ_dmu = (Q(m_e, m_mu*(1+eps), m_tau) - Q(m_e, m_mu*(1-eps), m_tau)) / (2*eps*m_mu)
dQ_dtau = (Q(m_e, m_mu, m_tau*(1+eps)) - Q(m_e, m_mu, m_tau*(1-eps))) / (2*eps*m_tau)

print(f"\nPartial derivatives:")
print(f"  dQ/dm_e   = {dQ_dme:.6e} / MeV")
print(f"  dQ/dm_μ   = {dQ_dmu:.6e} / MeV")
print(f"  dQ/dm_τ   = {dQ_dtau:.6e} / MeV")

# Error propagation
dQ = np.sqrt((dQ_dme * dm_e)**2 + (dQ_dmu * dm_mu)**2 + (dQ_dtau * dm_tau)**2)
print(f"\nδQ = {dQ:.6e}")
print(f"|Q - 3/2| / δQ = {abs(Q_central - 1.5) / dQ:.2f} sigma")

# Breakdown
print(f"\nContributions to δQ:")
print(f"  From m_e:  {abs(dQ_dme * dm_e):.6e} ({(dQ_dme*dm_e)**2/dQ**2*100:.1f}%)")
print(f"  From m_μ:  {abs(dQ_dmu * dm_mu):.6e} ({(dQ_dmu*dm_mu)**2/dQ**2*100:.1f}%)")
print(f"  From m_τ:  {abs(dQ_dtau * dm_tau):.6e} ({(dQ_dtau*dm_tau)**2/dQ**2*100:.1f}%)")

# ============================================================
print("\n" + "=" * 60)
print("TASK: M_W prediction from Casimir formula")
print("=" * 60)

# Casimir eigenvalues
x_half = (-3 + np.sqrt(57)) / 8  # s = 1/2
x_one = np.sqrt(3) - 1           # s = 1

print(f"x(1/2) = {x_half:.15f}")
print(f"x(1)   = {x_one:.15f}")

# R ratio
R = 1 - x_half / x_one
print(f"R = 1 - x(1/2)/x(1) = {R:.15f}")

# Verify closed form
R_closed = (np.sqrt(19) - 3) * (np.sqrt(19) - np.sqrt(3)) / 16
print(f"R = (√19-3)(√19-√3)/16 = {R_closed:.15f}")
print(f"|difference| = {abs(R - R_closed):.2e}")

# Mass predictions
M_Z = 91.1876  # GeV
dM_Z = 0.0021  # GeV

m_scale = M_Z / np.sqrt(x_one)
M_W_pred = m_scale * np.sqrt(x_half)

print(f"\nm = M_Z / √x(1) = {m_scale:.6f} GeV")
print(f"M_W = m × √x(1/2) = {M_W_pred:.6f} GeV")

# Error propagation: M_W = M_Z * sqrt(x_half / x_one)
# dM_W/dM_Z = sqrt(x_half/x_one) = M_W/M_Z
ratio = np.sqrt(x_half / x_one)
dM_W = ratio * dM_Z
print(f"\nM_W = M_Z × √(x(1/2)/x(1)) = M_Z × {ratio:.10f}")
print(f"δM_W = {ratio:.6f} × δM_Z = {dM_W:.6f} GeV")

# Compare to PDG
M_W_pdg = 80.3692  # GeV (PDG 2023 world average)
dM_W_pdg = 0.0133  # GeV

tension_pred = abs(M_W_pred - M_W_pdg) / dM_W_pdg
tension_combined = abs(M_W_pred - M_W_pdg) / np.sqrt(dM_W_pdg**2 + dM_W**2)

print(f"\nM_W (predicted)  = {M_W_pred:.4f} ± {dM_W:.4f} GeV")
print(f"M_W (PDG)        = {M_W_pdg:.4f} ± {dM_W_pdg:.4f} GeV")
print(f"Difference       = {M_W_pred - M_W_pdg:.4f} GeV = {(M_W_pred - M_W_pdg)*1000:.1f} MeV")
print(f"Tension (PDG err) = {tension_pred:.2f}σ")
print(f"Tension (combined) = {tension_combined:.2f}σ")

# Also compare to CDF-II controversial value
M_W_cdf = 80.4335
dM_W_cdf = 0.0094
tension_cdf = abs(M_W_pred - M_W_cdf) / np.sqrt(dM_W_cdf**2 + dM_W**2)
print(f"\nM_W (CDF-II)     = {M_W_cdf:.4f} ± {dM_W_cdf:.4f} GeV")
print(f"Tension (CDF-II) = {tension_cdf:.2f}σ")

# The fourth eigenvalue (|M(1/2,-)|)
x_half_neg = (-3 - np.sqrt(57)) / 8
M_half_neg = m_scale * np.sqrt(abs(x_half_neg))
print(f"\n|M(1/2,-)| = {M_half_neg:.4f} GeV")
print(f"Compare to M_H = 125.25 ± 0.17 GeV: deviation = {abs(M_half_neg - 125.25)/125.25*100:.2f}%")

# Cyclic echo verification
print("\n" + "=" * 60)
print("TASK: Cyclic echo verification")
print("=" * 60)

# From chain: stall at 0.03465 MeV
# Q(m_stall, m_?, m_s) = 3/2 where m_s = 92.166 MeV
# Does branch (+,-) give m_c = 1357 MeV?

m_stall = 0.034650  # MeV (from R2-B Task 3)
m_s = 92.166141  # MeV (from chain, not PDG)

# The quadratic in sqrt(m2):
# (s1*sqrt(m1) + s2*sqrt(m2) + s3*sqrt(m3))^2 = 3/2 * (m1 + m2 + m3)
# Try all sign/branch combos
print(f"Input: m_stall = {m_stall:.6f} MeV, m_s = {m_s:.6f} MeV")
for s1 in [1, -1]:
    for s3 in [1, -1]:
        a = s1*np.sqrt(m_stall) + s3*np.sqrt(m_s)
        disc = 6*a**2 - 3*(m_stall + m_s)
        if disc < 0:
            continue
        sqrt_disc = np.sqrt(disc)
        for s2 in [1, -1]:
            for pm in [1, -1]:
                u = 2*a*s2 + pm*sqrt_disc
                if u < 0:
                    continue
                m2 = u**2
                Q_val = (s1*np.sqrt(m_stall) + s2*np.sqrt(m2) + s3*np.sqrt(m_s))**2 / (m_stall + m2 + m_s)
                if abs(Q_val - 1.5) < 1e-8:
                    signs = f"({'+' if s1>0 else '-'},{'+' if s2>0 else '-'},{'+' if s3>0 else '-'})"
                    echo_flag = ""
                    if abs(m2 - 1356.96) < 2:
                        echo_flag = " ** ECHOES m_c **"
                    elif abs(m2 - 4180) < 50:
                        echo_flag = " ** ECHOES m_b **"
                    print(f"  signs={signs} pm={'+' if pm>0 else '-'}: m = {m2:.4f} MeV = {m2/1000:.4f} GeV, |Q-3/2|={abs(Q_val-1.5):.2e}{echo_flag}")

print("\n" + "=" * 60)
print("ALL TASKS COMPLETE")
print("=" * 60)
