"""
Refinements: higher-precision Taylor coefficients and detailed analysis at t=sqrt(3).
"""
from mpmath import mp, mpf, sqrt, log, diff, pi, taylor

mp.dps = 80

def mass_spectrum(t, y):
    t = mpf(t); y = mpf(y)
    disc_f = sqrt(t**2 + 1)
    fm_plus = (disc_f + t)**2
    fm_minus = (disc_f - t)**2
    disc_sigma = sqrt((2*t**2 + y)**2 + 4*t**2)
    sig_plus = 2*t**2 + y + 1 + disc_sigma
    sig_minus = 2*t**2 + y + 1 - disc_sigma
    disc_pi = sqrt((2*t**2 - y)**2 + 4*t**2)
    pi_plus = 2*t**2 - y + 1 + disc_pi
    pi_minus = 2*t**2 - y + 1 - disc_pi
    return [sig_plus, sig_minus, pi_plus, pi_minus], [fm_plus, fm_minus]

def vcw(t, y):
    t = mpf(t); y = mpf(y)
    scalars, fermions = mass_spectrum(t, y)
    result = mpf(0)
    for msq in scalars:
        if msq > 0:
            result += msq**2 * log(msq)
    for msq in fermions:
        if msq > 0:
            result -= 2 * msq**2 * log(msq)
    return result

y = mpf('0.25')
t3 = sqrt(mpf(3))

# ---- Task 3 refined: derivatives at sqrt(3) ----
print("TASK 3 REFINED: derivatives at t = sqrt(3)")
print("=" * 60)

def f(t): return vcw(t, mpf('0.25'))

# Use diff with method='quad' for higher accuracy
v0 = f(t3)
dv1 = diff(f, t3, 1)
dv2 = diff(f, t3, 2)
dv3 = diff(f, t3, 3)

print(f"V(sqrt(3))   = {mp.nstr(v0, 25)}")
print(f"V'(sqrt(3))  = {mp.nstr(dv1, 20)}")
print(f"V''(sqrt(3)) = {mp.nstr(dv2, 20)}")
print(f"V'''(sqrt(3))= {mp.nstr(dv3, 20)}")

# ---- Task 4 refined: Taylor with more points ----
print("\nTASK 4 REFINED: Taylor coefficients")
print("=" * 60)

from mpmath import matrix, lu_solve

# Use 5 points for 5 coefficients: c0, c2, c4, c6, c8
# V(h) ~ c0 + c2*h^2 + c4*h^4 + c6*h^6 + c8*h^8
# (V(h) - V(0))/h^2 = c2 + c4*h^2 + c6*h^4 + c8*h^6

eps = mpf('1e-25')
v_at_0 = vcw(eps, y)
print(f"V(0) = {mp.nstr(v_at_0, 20)}")

# Analytic V(0) for y=1/4
v0_exact = mpf(9)/4 * log(mpf(3)/2) + mpf(1)/4 * log(mpf(1)/2)
print(f"V(0) analytic = {mp.nstr(v0_exact, 20)}")

h_pts = [mpf(s) for s in ['0.0005', '0.001', '0.002', '0.004', '0.008']]
r_pts = [(vcw(h, y) - v0_exact) / h**2 for h in h_pts]

A = matrix(5, 5)
for i in range(5):
    for j in range(5):
        A[i,j] = h_pts[i]**(2*j)
b = matrix(r_pts)
sol = lu_solve(A, b)

names = ['c2', 'c4', 'c6', 'c8', 'c10']
for i, name in enumerate(names):
    print(f"  {name} = {mp.nstr(sol[i], 15)}")

# Verification: reconstruct V at h=0.003
h_test = mpf('0.003')
v_recon = v0_exact + h_test**2 * sum(sol[j] * h_test**(2*j) for j in range(5))
v_direct = vcw(h_test, y)
print(f"\nVerification at h=0.003:")
print(f"  Reconstructed: {mp.nstr(v_recon, 20)}")
print(f"  Direct:        {mp.nstr(v_direct, 20)}")
print(f"  Difference:    {mp.nstr(abs(v_recon - v_direct), 5)}")

# ---- Supertrace check ----
print("\nSUPERTRACE CHECK at t=sqrt(3), y=1/4:")
print("=" * 60)
scalars, fermions = mass_spectrum(t3, y)
str_m2 = sum(scalars) - 2*sum(fermions)
str_m4 = sum(s**2 for s in scalars) - 2*sum(f**2 for f in fermions)
print(f"  STr[M^2] = {mp.nstr(str_m2, 15)}")
print(f"  STr[M^4] = {mp.nstr(str_m4, 15)}")

# Include the zero modes: 2 scalar zeros + 1 fermion zero
# Total DOF: 6 real scalars, 3 Dirac fermions
# STr[M^0] = 6 - 2*3*2 = 6 - 12 ... wait, Weyl fermions
# For Weyl: STr = sum_scalars - 2*sum_fermions (the factor 2 for Dirac = 2 Weyl)
# With zeros: STr[M^2] += 2*0 - 2*0 = same
print(f"  (STr[M^2] should be 2y*m^2 = {mp.nstr(2*y, 8)} in our units)")

# ---- y=0 special case analysis ----
print("\ny=0 SPECIAL CASE:")
print("=" * 60)
scalars0, fermions0 = mass_spectrum(t3, mpf(0))
print(f"  Scalars: {[mp.nstr(s, 15) for s in scalars0]}")
print(f"  Fermions: {[mp.nstr(f, 15) for f in fermions0]}")
print(f"  Note: at y=0, sigma and pi sectors are degenerate.")
print(f"  Each scalar pair = fermion pair -> perfect cancellation -> V_CW = 0")
