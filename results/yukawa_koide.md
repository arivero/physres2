# Yukawa Coupling Matrices from the Koide Parametrization (tan beta = 1)

## Setup

The Koide parametrization writes three signed square roots as:

    sqrt(m_k) = sqrt(M0) * (1 + sqrt(2) * cos(2*pi*k/3 + delta)),  k = 0, 1, 2

where M0 = (sum signed sqrt(m_k))^2 / 9 and delta is the Koide phase.

In a two-Higgs-doublet model with tan beta = v_u/v_d = 1,
both vevs equal v/sqrt(2) and the Yukawa couplings are

    y_k = 2 m_k / v

for all fermion types (up, down, charged lepton). Here v = 246.22 GeV.

The Koide parametrization transfers to Yukawa couplings:

    sqrt(y_k) = sqrt(y0) * (1 + sqrt(2) * cos(2*pi*k/3 + delta))

with y0 = 2*M0/v. The Koide phase delta is identical for masses and Yukawa couplings.

PDG masses used (MeV): m_e = 0.51100, m_mu = 105.658, m_tau = 1776.86,
m_s = 93.4, m_c = 1270.0, m_b = 4180.0, m_t = 172760.0.

## Koide Parameters and Yukawa Couplings

### (e, mu, tau)

| Parameter | Value |
|-----------|-------|
| M0 (MeV) | 313.840920 |
| delta (rad) | 2.31662426 |
| delta (deg) | 132.7328 |
| Q (signed) | 0.6666608187 |
| y0 = 2*M0/v | 2.549272e-03 |

| Fermion | m (MeV) | y = 2m/v | sqrt(y) (signed) |
|---------|---------|----------|-------------------|
| e | 0.5110 | 4.150759e-06 | +2.037341e-03 |
| mu | 105.6580 | 8.582406e-04 | +2.929574e-02 |
| tau | 1776.8600 | 1.443311e-02 | +1.201379e-01 |

Verification: (sum signed sqrt(y))^2 / 9 = 2.549272e-03 = y0 = 2.549272e-03

### (-s, c, b)

| Parameter | Value |
|-----------|-------|
| M0 (MeV) | 912.555714 |
| delta (rad) | 2.74384326 |
| delta (deg) | 157.2106 |
| Q (signed) | 0.6749542234 |
| y0 = 2*M0/v | 7.412523e-03 |

| Fermion | m (MeV) | y = 2m/v | sqrt(y) (signed) |
|---------|---------|----------|-------------------|
| s | 93.4000 | 7.586711e-04 | -2.754398e-02 |
| c | 1270.0000 | 1.031598e-02 | +1.015676e-01 |
| b | 4180.0000 | 3.395338e-02 | +1.842644e-01 |

Verification: (sum signed sqrt(y))^2 / 9 = 7.412523e-03 = y0 = 7.412523e-03

### (c, b, t)

| Parameter | Value |
|-----------|-------|
| M0 (MeV) | 29576.439062 |
| delta (rad) | 2.16303821 |
| delta (deg) | 123.9330 |
| Q (signed) | 0.6694893550 |
| y0 = 2*M0/v | 2.402440e-01 |

| Fermion | m (MeV) | y = 2m/v | sqrt(y) (signed) |
|---------|---------|----------|-------------------|
| c | 1270.0000 | 1.031598e-02 | +1.015676e-01 |
| b | 4180.0000 | 3.395338e-02 | +1.842644e-01 |
| t | 172760.0000 | 1.403298e+00 | +1.184609e+00 |

Verification: (sum signed sqrt(y))^2 / 9 = 2.402440e-01 = y0 = 2.402440e-01

## Cartan Subalgebra Decomposition

The diagonal matrix sqrt(Y) = diag(sqrt(y_0), sqrt(y_1), sqrt(y_2)) decomposes as:

    sqrt(Y) = sqrt(y0) * (I + (sqrt(6)/2) * (n1 * lambda3 + n2 * lambda8))

where lambda3 = diag(1, -1, 0) and lambda8 = diag(1, 1, -2)/sqrt(3) are
the diagonal Gell-Mann matrices, and n_hat = (n1, n2) = (cos theta, sin theta)
is a unit vector in the Cartan plane.

The Cartan angle theta is related to the Koide phase by an exact identity:

    theta = pi/6 - delta

This follows from the trigonometric identities:

    n1 propto cos(delta) - cos(2*pi/3 + delta) = sqrt(3) * sin(pi/3 + delta)
    n2 propto cos(delta) + cos(2*pi/3 + delta) - 2*cos(4*pi/3 + delta) = 3*cos(pi/3 + delta)

so theta = arctan2(cos(pi/3 + delta), sin(pi/3 + delta)) = pi/2 - pi/3 - delta = pi/6 - delta.

### (e, mu, tau)

| Cartan parameter | Value |
|------------------|-------|
| n_hat | (-0.220405, -0.975409) |
| n_norm | 0.99999123 |
| theta (deg) | -102.7328 |
| pi/6 - delta (deg) | -102.7328 |

### (-s, c, b)

| Cartan parameter | Value |
|------------------|-------|
| n_hat | (-0.604747, -0.796418) |
| n_norm | 1.01235501 |
| theta (deg) | -127.2106 |
| pi/6 - delta (deg) | -127.2106 |

### (c, b, t)

| Cartan parameter | Value |
|------------------|-------|
| n_hat | (-0.068589, -0.997645) |
| n_norm | 1.00422511 |
| theta (deg) | -93.9330 |
| pi/6 - delta (deg) | -93.9330 |

## Seed Yukawa Matrices and Bloom Rotations

At delta_seed = 3*pi/4 = 135 deg, the k=0 coefficient vanishes:
1 + sqrt(2)*cos(3*pi/4) = 1 - 1 = 0, so m_0 = y_0 = 0.

The bloom rotation Delta_delta = delta_phys - 3*pi/4 parametrizes how far
the physical triple has rotated from the seed configuration.

### (e, mu, tau)

| Bloom parameter | Value |
|-----------------|-------|
| delta_phys (deg) | 132.7328 |
| delta_seed (deg) | 135.0000 |
| Delta_delta (deg) | -2.2672 |
| Delta_delta (rad) | -0.03957023 |

Seed Yukawa couplings at delta = 3*pi/4 (with same M0 = 313.8409 MeV):

| k | Seed m (MeV) | Seed y | Physical m (MeV) | Physical y |
|---|--------------|--------|------------------|------------|
| 0 | 0.0000e+00 | 0.0000e+00 | 0.5110 | 4.1508e-06 |
| 1 | 1.2614e+02 | 1.0246e-03 | 105.6580 | 8.5824e-04 |
| 2 | 1.7569e+03 | 1.4271e-02 | 1776.8600 | 1.4433e-02 |

### (-s, c, b)

| Bloom parameter | Value |
|-----------------|-------|
| delta_phys (deg) | 157.2106 |
| delta_seed (deg) | 135.0000 |
| Delta_delta (deg) | 22.2106 |
| Delta_delta (rad) | 0.38764877 |

Seed Yukawa couplings at delta = 3*pi/4 (with same M0 = 912.5557 MeV):

| k | Seed m (MeV) | Seed y | Physical m (MeV) | Physical y |
|---|--------------|--------|------------------|------------|
| 0 | 0.0000e+00 | 0.0000e+00 | 93.4000 | 7.5867e-04 |
| 1 | 3.6678e+02 | 2.9793e-03 | 1270.0000 | 1.0316e-02 |
| 2 | 5.1086e+03 | 4.1496e-02 | 4180.0000 | 3.3953e-02 |

### (c, b, t)

| Bloom parameter | Value |
|-----------------|-------|
| delta_phys (deg) | 123.9330 |
| delta_seed (deg) | 135.0000 |
| Delta_delta (deg) | -11.0670 |
| Delta_delta (rad) | -0.19315628 |

Seed Yukawa couplings at delta = 3*pi/4 (with same M0 = 29576.4391 MeV):

| k | Seed m (MeV) | Seed y | Physical m (MeV) | Physical y |
|---|--------------|--------|------------------|------------|
| 0 | 0.0000e+00 | 0.0000e+00 | 1270.0000 | 1.0316e-02 |
| 1 | 1.1887e+04 | 9.6560e-02 | 4180.0000 | 3.3953e-02 |
| 2 | 1.6557e+05 | 1.3449e+00 | 172760.0000 | 1.4033e+00 |

## Bloom as Cartan Plane Rotation

Since theta = pi/6 - delta, shifting delta by Delta_delta shifts theta by -Delta_delta.
In the Cartan plane, the bloom acts as a rigid rotation:

    (n1, n2) -> R(-Delta_delta) * (n1, n2)

where R(alpha) = [[cos alpha, -sin alpha], [sin alpha, cos alpha]].

Equivalently, the traceless part of sqrt(Y) rotates in the (lambda3, lambda8) plane
without changing its norm (which is fixed at sqrt(6)/2 for Q = 2/3) or the trace
(which is 3*sqrt(y0) and set by M0).

The three masses after a bloom of Delta_delta from the seed are:

    m_k(M0, Delta) = M0 * (1 + sqrt(2) * cos(2*pi*k/3 + 3*pi/4 + Delta))^2

At Delta = 0 (seed): m_0 = 0, m_1 = M0*(1+sqrt(2)*cos(2*pi/3+3*pi/4))^2,
m_2 = M0*(1+sqrt(2)*cos(4*pi/3+3*pi/4))^2.

Summary of bloom rotations:

| Triple | Delta_delta (deg) | Delta_delta (rad) | Cartan rotation (deg) |
|--------|-------------------|-------------------|-----------------------|
| (e, mu, tau) | -2.2672 | -0.039570 | 2.2672 |
| (-s, c, b) | 22.2106 | 0.387649 | -22.2106 |
| (c, b, t) | -11.0670 | -0.193156 | 11.0670 |

## Numerical Verification

The Koide condition Q = sum(y_k) / (sum signed sqrt(y_k))^2 = 2/3
is algebraically identical to Q = sum(m_k) / (sum signed sqrt(m_k))^2 = 2/3
because y_k = (2/v)*m_k and the factor (2/v) cancels between numerator and denominator:

    Q = sum((2/v)*m) / (sum(sqrt(2/v)*sqrt(m)))^2 = (2/v)*sum(m) / ((2/v)*(sum sqrt(m))^2) = sum(m)/(sum sqrt(m))^2

| Triple | Q (signed) | |Q - 2/3| |
|--------|------------|-----------|
| (e, mu, tau) | 0.6666608187 | 5.85e-06 |
| (-s, c, b) | 0.6749542234 | 8.29e-03 |
| (c, b, t) | 0.6694893550 | 2.82e-03 |

For (e, mu, tau), Q is 2/3 to 6e-6 (0.001%), confirming the near-exact Koide relation.
For (-s, c, b) and (c, b, t), the deviations (1.2% and 0.4%) reflect how far the
PDG masses lie from the nearest Koide curve. The algebraic identity Q_signed = 2/3
holds exactly for the *fitted* parametrization values, not for the raw PDG masses.

## Summary

| Triple | M0 (MeV) | y0 | delta (deg) | Delta_delta (deg) | theta (deg) |
|--------|----------|-----|-------------|-------------------|-------------|
| (e, mu, tau) | 313.8409 | 2.5493e-03 | 132.7328 | -2.2672 | -102.7328 |
| (-s, c, b) | 912.5557 | 7.4125e-03 | 157.2106 | 22.2106 | -127.2106 |
| (c, b, t) | 29576.4391 | 2.4024e-01 | 123.9330 | -11.0670 | -93.9330 |

Key results:

1. The Koide parametrization transfers directly from masses to Yukawa couplings
   via y0 = 2*M0/v. The phase delta is identical.

2. The Cartan subalgebra decomposition of sqrt(Y) has an exact relationship
   theta = pi/6 - delta between the Cartan plane angle and the Koide phase.

3. The bloom from the seed (delta = 3*pi/4, one mass zero) to the physical
   triple is a rigid rotation in the Cartan plane by angle -Delta_delta.

4. The Q = 2/3 condition is invariant under the rescaling m -> y = 2m/v,
   as it must be since Q is a ratio of homogeneous-degree-1 expressions.

5. At tan beta = 1 the top Yukawa coupling is y_t = 2*m_t/v = 1.403,
   which is O(1) as expected. The full range spans y_e = 4.15e-06
   to y_t = 1.4033, a hierarchy of 3e+05.