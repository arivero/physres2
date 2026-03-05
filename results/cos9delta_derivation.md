# How (det M)³/Λ¹⁸ Generates a δ-Dependent Potential

## Key Formula
f(δ) = ∏_k [1 + √2 cos(δ + 2πk/3)] = −1/2 + cos(3δ)/√2

This is EXACT (verified by SymPy).

## Fourier Decomposition of [f(δ)]⁶
det M = z₀⁶ [f(δ)]², so (det M)³ = z₀¹⁸ [f(δ)]⁶

| Harmonic | Coefficient | Relative amplitude |
|----------|------------|-------------------|
| const | 41/64 = 0.641 | — |
| cos 3δ | −51√2/64 = −1.127 | 1.00 |
| cos 6δ | 195/256 = 0.762 | 0.68 |
| cos 9δ | −35√2/128 = −0.387 | 0.34 |
| cos 12δ | 9/64 = 0.141 | 0.12 |
| cos 15δ | −3√2/128 = −0.033 | 0.03 |
| cos 18δ | 1/256 = 0.004 | 0.003 |

## Result: cos(9δ) is NOT dominant
- cos(3δ) is 3× larger than cos(9δ)
- Actual potential [f(δ)]⁶ has 6 local minima, not 9
- For pure Z₉ structure, lower harmonics must cancel

## Implication
The statement "(det M)³ generates a cos(9δ) potential" is misleading.
It generates a potential with ALL harmonics cos(3nδ), n=0,...,6.
The cos(9δ) term IS the unique Z₉-periodic component, but is subdominant.
