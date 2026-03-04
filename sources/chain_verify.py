import numpy as np

print("="*70)
print("  ADDITIONAL VERIFICATIONS AND NOTES")
print("="*70)

# Verify the formula derivation
# Q = (s1+s2+s3)^2 / (s1^2+s2^2+s3^2) = 3/2
# => 2*(s1+s2+s3)^2 = 3*(s1^2+s2^2+s3^2)
# => 2*s1^2+2*s2^2+2*s3^2+4*s1*s2+4*s1*s3+4*s2*s3 = 3*s1^2+3*s2^2+3*s3^2
# => s3^2 - 4*(s1+s2)*s3 + (s1^2+s2^2 - 4*s1*s2) = 0
# quadratic in s3:
# s3 = [4*(s1+s2) +/- sqrt(16*(s1+s2)^2 - 4*(s1^2+s2^2-4*s1*s2))] / 2
# s3 = 2*(s1+s2) +/- sqrt(4*(s1+s2)^2 - (s1^2+s2^2-4*s1*s2))
# disc = 4*s1^2+8*s1*s2+4*s2^2 - s1^2-s2^2+4*s1*s2
#      = 3*s1^2+12*s1*s2+3*s2^2
#      = 3*(s1+s2)^2 + 6*s1*s2
#
# So s3 = 2*(s1+s2) +/- sqrt(3*(s1+s2)^2 + 6*s1*s2)
#        = (s1+s2) * [2 +/- sqrt(3 + 6*s1*s2/(s1+s2)^2)]
#
# This matches the formula given. Note s3 here is sqrt(m3) and can be + or -.

# Check: for m_b=4.18, m_t=172.69
m_b = 4.18
m_t = 172.69
s_b = np.sqrt(m_b)
s_t = np.sqrt(m_t)
S = s_b + s_t

ratio = 6*s_b*s_t/S**2
inner = np.sqrt(3 + ratio)

print("\n  Verification of quadratic formula:")
print("  s_b = {:.10f}, s_t = {:.10f}".format(s_b, s_t))
print("  S = s_b+s_t = {:.10f}".format(S))
print("  ratio = 6*s_b*s_t/S^2 = {:.10f}".format(ratio))
print("  inner = sqrt(3+ratio) = {:.10f}".format(inner))
print("  S*(2-inner) = {:.10f}  (branch +,+)".format(S*(2-inner)))
print("  S*(2+inner) = {:.10f}  (branch +,-)".format(S*(2+inner)))

# Branch 1: s3 = S*(2-inner) = 1.16488716...
s3_1 = S*(2-inner)
q1 = (s_b+s_t+s3_1)**2 / (s_b**2+s_t**2+s3_1**2)
print("  Q_signed for branch (+,+): {:.15f}".format(q1))

# Branch 2: s3 = S*(2+inner) = 59.577...
s3_2 = S*(2+inner)
q2 = (s_b+s_t+s3_2)**2 / (s_b**2+s_t**2+s3_2**2)
print("  Q_signed for branch (+,-): {:.15f}".format(q2))

# Now the NEGATIVE branches mean s3 is negative
# Branch 3: s3 = -S*(2-inner) = -1.16488716...
s3_3 = -S*(2-inner)
q3 = (s_b+s_t+s3_3)**2 / (s_b**2+s_t**2+s3_3**2)
print("  Q_signed for branch (-,+): {:.15f}  -- does NOT give 3/2".format(q3))

# Branch 4: s3 = -S*(2+inner) = -59.577...
s3_4 = -S*(2+inner)
q4 = (s_b+s_t+s3_4)**2 / (s_b**2+s_t**2+s3_4**2)
print("  Q_signed for branch (-,-): {:.15f}  -- does NOT give 3/2".format(q4))

print("\n  CONCLUSION: Only outer_sign=+1 branches satisfy Q_signed = 3/2.")
print("  The outer_sign=-1 branches only satisfy Q_unsigned = 3/2 (with |sqrt|).")
print("  For the Koide formula with SIGNED sqrt, there are exactly TWO solutions")
print("  (inner_sign = +1 or -1), not four.")

# Now re-examine Step 2 carefully
# Q(m_c_pred, m_?, m_b) = 3/2 where m_c_pred = 1.35696209
m_c_pred = 1.35696209
s_c = np.sqrt(m_c_pred)
# We need s_? such that Q_signed(s_c, s_b, s_?) = 3/2
# Using the formula with m1=m_c_pred, m2=m_b:
# s_? = (s_c+s_b)*(2 +/- sqrt(3 + 6*s_c*s_b/(s_c+s_b)^2))

S2 = s_c + s_b
ratio2 = 6*s_c*s_b/S2**2
inner2 = np.sqrt(3+ratio2)

s_s_plus = S2*(2-inner2)
s_s_minus = S2*(2+inner2)

print("\n\n  Step 2 re-examination:")
print("  s_c = {:.10f}, s_b = {:.10f}".format(s_c, s_b))
print("  S = {:.10f}".format(S2))
print("  Branch (+): s_? = {:.10f}".format(s_s_plus))
print("  Branch (-): s_? = {:.10f}".format(s_s_minus))
print("  m_?(+) = {:.10f} GeV".format(s_s_plus**2))
print("  m_?(-) = {:.10f} GeV".format(s_s_minus**2))

# Check: is s_s_plus negative?
print("  s_s_plus is {}".format('NEGATIVE' if s_s_plus < 0 else 'positive'))

# Verify Q
q_plus = (s_c + s_b + s_s_plus)**2 / (s_c**2 + s_b**2 + s_s_plus**2)
q_minus = (s_c + s_b + s_s_minus)**2 / (s_c**2 + s_b**2 + s_s_minus**2)
print("  Q_signed(+): {:.15f}".format(q_plus))
print("  Q_signed(-): {:.15f}".format(q_minus))

print("\n  So in Step 2, the solution s_? = {:.10f} is NATURALLY negative".format(s_s_plus))
print("  from the formula with inner_sign=+1 (no extra sign flip needed).")
print("  This is a KEY insight: the formula (s1+s2)*(2-sqrt(3+...)) can naturally")
print("  give a negative sqrt(m) when the inner sqrt exceeds 2.")

# Let's check: is inner2 > 2?
print("  inner = sqrt(3+ratio) = {:.10f}".format(inner2))
print("  inner > 2? {}  (inner - 2 = {:.10f})".format(inner2 > 2, inner2 - 2))
# inner2 ~ 1.905... so 2-inner2 ~ 0.094 > 0 ... wait
# Let me recheck
print("  2-inner = {:.10f}".format(2-inner2))
print("  So s_? = S*(2-inner) = {:.10f}*{:.10f} = {:.10f}".format(S2, 2-inner2, S2*(2-inner2)))
# That gives -0.30358..., so 2-inner2 < 0, meaning inner2 > 2
print("  inner2 = {:.10f} > 2, so (2-inner) < 0 and s_? < 0".format(inner2))

# Now for Step 3
m_s_pred = s_s_plus**2
s_s = np.sqrt(m_s_pred)  # positive sqrt of the mass
# But the "signed" sqrt for s in the Koide relation is s_s_plus which is negative!
# When we continue the chain, we use m_s_pred and m_c_pred as inputs
# and solve for s_next from (sqrt(m_s_pred) + sqrt(m_c_pred))*(2 - sqrt(...))
# This uses the POSITIVE square roots of the masses as inputs to the formula.
# The "sign" information is lost when we square.

print("\n\n  Step 3: chain continuation")
print("  Using m_s = {:.10f} (from s_s = {:.10f})".format(m_s_pred, s_s_plus))
print("  The formula uses sqrt(m_s) = {:.10f} (positive)".format(s_s))
print("  and sqrt(m_c) = {:.10f}".format(s_c))

S3 = s_s + s_c
ratio3 = 6*s_s*s_c/S3**2
inner3 = np.sqrt(3+ratio3)
s_next_p = S3*(2-inner3)
s_next_m = S3*(2+inner3)
print("  inner = {:.10f}".format(inner3))
print("  Branch (+): s_next = {:.12f}, m = {:.10e} GeV = {:.6f} MeV".format(
    s_next_p, s_next_p**2, s_next_p**2*1000))
print("  Branch (-): s_next = {:.12f}, m = {:.10e} GeV = {:.6f} MeV".format(
    s_next_m, s_next_m**2, s_next_m**2*1000))
print("  Both are much smaller than m_d = 4.67 MeV or m_u = 2.16 MeV")
print("  Branch (+): {:.4f} MeV vs m_d = 4.67 MeV, m_u = 2.16 MeV".format(s_next_p**2*1000))

# Step 4
m_next = s_next_p**2
s_next = np.sqrt(m_next)
S4 = s_next + s_s
ratio4 = 6*s_next*s_s/S4**2
inner4 = np.sqrt(3+ratio4)
s4_p = S4*(2-inner4)
s4_m = S4*(2+inner4)
print("\n  Step 4:")
print("  Branch (+): s = {:.12f}, m = {:.10e} GeV = {:.8f} MeV".format(
    s4_p, s4_p**2, s4_p**2*1000))
print("  Branch (-): s = {:.12f}, m = {:.10e} GeV = {:.8f} MeV".format(
    s4_m, s4_m**2, s4_m**2*1000))

# Note that m_c appears again!
print("  Branch (-) gives m = {:.6f} GeV which equals m_c_pred = {:.6f} GeV!".format(
    s4_m**2, m_c_pred))
print("  This is the 'echo' of the original m_c.")

# Now check: in Chain B, the Q for the quark triplet must use signed sqrts
print("\n\n  === Chain B: Q check with signed square roots ===")
m_e = 0.00051099895
m_mu = 0.1056583755
m_tau = 1.77686

M0_fit = 0.313838229066
delta0_fit = 2.316624733889
sM0 = np.sqrt(M0_fit)

M1 = 3*M0_fit
delta1 = 3*delta0_fit
sM1 = np.sqrt(M1)

sqrts = []
masses = []
for k in range(3):
    s = sM1 * (1 + np.sqrt(2)*np.cos(2*np.pi*k/3 + delta1))
    sqrts.append(s)
    masses.append(s**2)

print("  Signed sqrts: {:.10f}, {:.10f}, {:.10f}".format(*sqrts))
print("  Masses: {:.10f}, {:.10f}, {:.10f}".format(*masses))

Q_s = (sqrts[0]+sqrts[1]+sqrts[2])**2 / (masses[0]+masses[1]+masses[2])
Q_u = (np.sqrt(masses[0])+np.sqrt(masses[1])+np.sqrt(masses[2]))**2 / sum(masses)
print("  Q_signed = {:.15f}".format(Q_s))
print("  Q_unsigned = {:.15f}".format(Q_u))
print("  The signed Q = 3/2 BY CONSTRUCTION (the parametrization guarantees it).")
print("  The unsigned Q != 3/2 because one sqrt is negative (k=1: {:.10f})".format(sqrts[1]))

# Note about the delta0 value
print("\n\n  === Interesting numerical coincidence ===")
print("  delta0 = {:.10f} rad".format(delta0_fit))
print("  2/9 = {:.10f}".format(2.0/9.0))
print("  delta0 - 2/9 = {:.10f}".format(delta0_fit - 2.0/9.0))
print("  pi*3/4 = {:.10f}".format(np.pi*3/4))

# Koide's original waterfall
print("\n\n  === Koide's original angle ===")
# In Koide's notation: delta = 2/27 is sometimes cited
# Let's check: the angle in the Koide formula for leptons
# sqrt(m_k) = sqrt(M)*(1 + sqrt(2)*cos(2*pi*k/3 + delta))
# delta for leptons is approximately 0.2222 radians = 2/9
# Or in some conventions delta = 2/27 with a different phase convention

# Actually, various sources use different conventions.
# Let me compute delta0 mod (2*pi/3) to see if there's a simpler form
delta0_mod = delta0_fit % (2*np.pi/3)
print("  delta0 mod (2*pi/3) = {:.10f}".format(delta0_mod))
print("  2*pi/3 = {:.10f}".format(2*np.pi/3))

# The classic Koide angle
# Some references use: sqrt(m_k) = sqrt(M)*(1 + sqrt(2)*cos(theta_k + delta))
# where theta_k = 0, 2*pi/3, 4*pi/3
# This is the same as k=0,1,2 with 2*pi*k/3

# Rivero's convention might differ. Let's just report what we have.
print("\n  Summary of angle:")
print("  delta0 = {:.10f} rad".format(delta0_fit))
print("  = {:.10f} * pi".format(delta0_fit/np.pi))
print("  = {:.6f} degrees".format(np.degrees(delta0_fit)))

print("\n  delta1 = 3*delta0 = {:.10f} rad".format(delta1))
print("  delta1 mod 2*pi = {:.10f} rad = {:.6f} deg".format(
    delta1 % (2*np.pi), np.degrees(delta1 % (2*np.pi))))
