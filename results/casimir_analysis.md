# Casimir Eigenvalue Equation Analysis

## The equation

$$x^2 + C_2 \, x - C_2 = 0, \qquad C_2 = s(s+1)$$

### Structural constraints

The most general quadratic eigenvalue equation x^2 + f(s)*x + g(s) = 0 with f and g polynomials in the Casimir C_2 = s(s+1) is constrained by three requirements:

1. **f(0) = 0**: Zero spin (zero coupling) must give zero mass. Since f is a polynomial in C_2 and C_2(0) = 0, this is automatic for f = f_1 C_2 + f_2 C_2^2 + ...
2. **g(0) = 0**: The vacuum equation (s=0) must be x^2 = 0 with a double root at zero. Same vanishing condition.
3. **g_1 = -f_1**: Coefficient antisymmetry at linear order in C_2.

At leading order in C_2, these fix f = C_2, g = -C_2 (up to overall normalization), giving

$$x^2 + C_2 \, x - C_2 = 0$$

## Roots

### s = 1/2 (fermion)

- C_2 = 3/4
- Equation: x^2 + (3/4)x - 3/4 = 0, or equivalently 4x^2 + 3x - 3 = 0
- Discriminant: 9 + 48 = 57 = 3 * 19
- x_+ = (-3 + sqrt(57))/8 = 0.568729304408844
- x_- = (-3 - sqrt(57))/8 = -1.318729304408844

### s = 1 (vector boson)

- C_2 = 2
- Equation: x^2 + 2x - 2 = 0
- Discriminant: 4 + 8 = 12 = 4 * 3
- x_+ = -1 + sqrt(3) = 0.732050807568877
- x_- = -1 - sqrt(3) = -2.732050807568877

## Electroweak prediction

$$R = 1 - \frac{x_+(1/2)}{x_+(1)} = 0.2231013223008663...$$

| Quantity | Predicted | Experimental | Deviation |
|----------|-----------|-------------|-----------|
| sin^2(theta_W) | 0.2231013223 | 0.22306 +/- 0.00031 (PDG on-shell) | 0.13 sigma |
| M_W (GeV) | 80.3744 | 80.3692 +/- 0.0133 (PDG) | 0.39 sigma |
| M_W (GeV) | 80.3744 | 80.4335 +/- 0.0094 (CDF-II) | 6.3 sigma |

The predicted M_W disfavors the CDF-II measurement at 6.3 sigma.

## Algebraic identity with de Vries angle

The de Vries angle is defined as

$$R_{\rm dV} = \frac{(\sqrt{19}-3)(\sqrt{19}-\sqrt{3})}{16}$$

**This is algebraically identical to the Casimir ratio R.** The identity has been verified symbolically (sympy) and to 50 decimal digits (mpmath).

### Proof

$$R = 1 - \frac{-3+\sqrt{57}}{8(\sqrt{3}-1)}$$

Rationalize the denominator by multiplying by (sqrt(3)+1)/(sqrt(3)+1):

$$R = 1 - \frac{(-3+\sqrt{57})(\sqrt{3}+1)}{8 \cdot 2} = 1 - \frac{-3\sqrt{3}-3+3\sqrt{19}+\sqrt{57}}{16}$$

where we used sqrt(57)*sqrt(3) = sqrt(171) = 3*sqrt(19).

$$R = \frac{16 + 3\sqrt{3} + 3 - 3\sqrt{19} - \sqrt{57}}{16} = \frac{19 + 3\sqrt{3} - 3\sqrt{19} - \sqrt{57}}{16}$$

Expanding the de Vries formula:

$$R_{\rm dV} = \frac{19 - \sqrt{57} - 3\sqrt{19} + 3\sqrt{3}}{16}$$

These are term-by-term identical. QED.

### Significance

The number 19 arises from the s=1/2 discriminant: 57/3 = 19. The number 3 appears in both discriminants (57 = 3*19 for s=1/2; 12 = 4*3 for s=1). The de Vries formula, originally derived from mass relation considerations, decomposes the Casimir discriminant 57 = 3 * 19 into its prime factors and recombines them through sqrt(19) and sqrt(3). The Casimir eigenvalue equation provides a representation-theoretic derivation of the same algebraic number.

## Exact algebraic form

$$R = \frac{19 + 3\sqrt{3} - 3\sqrt{19} - \sqrt{57}}{16}$$

Equivalently:

$$\cos^2\theta_W = 1 - R = \frac{-3 + 3\sqrt{19} + \sqrt{57} - 3\sqrt{3}}{16}$$

## Search for alternative forms

Searched 582 polynomial equations including:
- Quadratics x^2 + f(s)*x + g(s) = 0 with 10 different Casimir-like functions for f, g
- Cubics x^3 + f(s)*x + g(s) = 0 (depressed cubics)
- Quartics x^4 + f(s)*x + g(s) = 0
- Rational coefficient variations x^2 + (p/q)C_2*x + (r/t)C_2 = 0 with |p|,|r| <= 3
- Alternative spin assignments (s1, s2) from {0.5, 1, 1.5, 2, 2.5, 3}

### Results

Of the 13 forms found within 1% of sin^2(theta_W):

- **11 are equivalent** to the original equation x^2 + C_2*x - C_2 = 0 (they are the same equation with rational coefficients that simplify to 1 and -1, e.g. (2/2)C_2 = C_2). All give R = 0.2231013223 (0.019% from PDG central value).

- **2 genuinely different quartic forms** appear, but with worse agreement:
  - x^4 + s*x - s(s+1) = 0: R = 0.2243 (0.54% deviation)
  - x^4 + s(s+1)/2*x - s(2s+1) = 0: R = 0.2252 (0.94% deviation)

No alternative spin pair (s1, s2) with the original equation gives a value within 1%.

**The original equation x^2 + C_2*x - C_2 = 0 with (s=1/2, s=1) is unique** among the forms tested in matching the observed weak mixing angle to sub-percent accuracy.

## Representation-theoretic interpretation

### 2x2 matrix realization

The equation x^2 + C_2 x - C_2 = 0 is the characteristic equation of the companion matrix

$$M = \begin{pmatrix} 0 & C_2 \\ 1 & -C_2 \end{pmatrix}$$

with tr(M) = -C_2 and det(M) = -C_2. The negative determinant guarantees one positive and one negative eigenvalue.

### Seesaw structure

The matrix M can be decomposed as

$$M = \begin{pmatrix} 0 & C_2 \\ 1 & 0 \end{pmatrix} + \begin{pmatrix} 0 & 0 \\ 0 & -C_2 \end{pmatrix}$$

This has the structure of a seesaw: off-diagonal mixing of order C_2 with a diagonal mass term -C_2. The positive eigenvalue x_+ < C_2 is the "light" state, analogous to the light neutrino in the seesaw mechanism.

### Self-consistency condition

The equation can be rewritten as

$$\frac{x^2}{1-x} = C_2$$

This defines x as a self-consistent solution where the "squared coupling" x^2 equals C_2 times the "complement" (1-x). The function x^2/(1-x) is monotonically increasing for 0 < x < 1, so there is exactly one positive root in this range for each C_2 > 0.

### Golden ratio connection

For C_2 = 1, the equation becomes x^2 + x - 1 = 0, the golden ratio equation with positive root x = (sqrt(5)-1)/2. The spin value at which C_2 = s(s+1) = 1 is itself s = (sqrt(5)-1)/2, the golden ratio.

### Lorentz group

The Casimir C_2 = s(s+1) is the quadratic Casimir of the SU(2) rotation subgroup of the Lorentz group. The two values s = 1/2 (Weyl fermion representations (1/2,0) and (0,1/2)) and s = 1 (vector representation (1/2,1/2)) are precisely the spin content of electroweak matter and gauge bosons. The ratio of positive eigenvalues at these two spins gives the weak mixing angle.

### Open question

Whether the equation x^2 + C_2 x - C_2 = 0 arises as the eigenvalue equation of a specific operator in the representation theory of the Poincare group or some extension (e.g., as a mass-shell condition modified by spin-dependent terms, or from a quadratic Casimir relation in a larger algebra containing both SU(2)_L and U(1)_Y) remains to be established. The seesaw-like structure of the companion matrix and the involvement of the SU(2) Casimir are suggestive, but a first-principles derivation from known representation theory has not been identified.
