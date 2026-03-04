# Bootstrapping generations

**Author:** Alejandro Rivero (arivero@unizar.es)
**Affiliation:** BIFI - Universidad de Zaragoza

## Abstract

A supersymmetric version of Chew's old democratic bootstrap argument predicts the existence of three generations of quarks and leptons with one quark, of type "up", more massive than the other five.

---

## Historical Background

The origin of the "dual models" traces to Geoffrey Chew's advocacy in the early sixties for the study of the S-matrix, looking for a unique consistent solution that "bootstraps" itself. A guiding principle was "nuclear democracy"; no particle should be claimed to be elementary. Each particle was providing the force needed to bind the others.

At the start of the seventies, it was noticed that dual models with fermions seem to require a special symmetry, supersymmetry, relating bosonic and fermionic degrees of freedom. Schwarz proposed to interpret this need as a "quark-gluon dual model".

Given the stringy character of these models, it seems reasonable to consider that the gluon part, bosonic, is constituted by a string terminated in a pair of quarks. With this imagery, we can try a version of the democratic bootstrap. We require that scalar leptons and scalar quarks are composed of pairs of quarks, bound by some gluonic string. Still, we keep agnostic about if quarks and leptons are composed of themselves or other particles, or if they are truly elementary.

Each possible pair of "string terminations" will produce one kind of scalar, as a triplet of colour in the case of pairs quark-quark, and as a colour singlet in the case of quark-antiquark; i.e., we are taking the triplet in $3 \times 3 = 3 + 6$ and the singlet in $3 \times \bar 3 = 1 + 8$ when combining the colour charges of the constituents. As for the electric charge, it is simply the sum of the charges of the two elements of the pair.

Following an emergent tradition, we will call squarks and sleptons to the scalar particles, and we will refer to the requirement of matching between the number of composite scalar particles and the number of degrees of freedom of the standard model fermions as the **sBootstrap**.

---

## The sBootstrap Condition

Let's first taste the sBootstrap condition at cask strength: $N$ generations of quarks of types $u, d$ should produce $N^2$ scalars of charge $+1/3$, but we have $2N$ fermionic degrees of freedom with such charge, so $N=2$ for the matching in the "down antiquark" sector. And they also produce, from pairs of two down-type quarks, a total of $N(N+1)/2$ scalars of charge $-2/3$, again to be matched to $2N$ degrees of freedom in the "up antiquark" sector. So $N=3$ for the matching in this sector. The mechanism fails, which is no surprise.

We dilute the condition by requesting that not every quark can terminate the string. We postulate that only some subset of "light quarks" are in the terminations of the string. With three parameters — number of generations $N$, number of light down quarks $r$, and number of light up quarks $s$ — the matching conditions are:

$$\begin{align}
r s &= 2 N \\
{r (r+1) \over 2} &= 2N
\end{align}$$

respectively for $+1/3$ and $-2/3$ particles. The integer solution is that $N$ must be half of an hexagonal number:

$$2 N = \cancel{1}, 6, \cancel{15}, 28, \cancel{45}, 66, \cancel{91}, 120, \ldots$$

So the **smallest admissible solution is $N=3$**.

---

## Asymptotic Freedom

The requirement of asymptotic freedom in the beta function of QCD puts an upper limit, asking the number of flavours $n_f$ to be such that $2 n_f < 33$. If this limitation applies to all flavours, light or heavy, then the solution is **unique**.

**Result: N=3 generations of quarks with r=3 light down quarks but only s=2 light up quarks.**

If the beta function is built only with the light quarks, we have a short list:

| N | s | r | Total light flavours |
|---|---|---|----|
| 3 | 2 | 3 | 5 |
| 14 | 4 | 7 | 11 |
| 33 | 6 | 11 | 17 |

Note that we always have $r = 2s-1$.

---

## Uniqueness from Slepton Sector

For neutral leptons, using $SU(r+s)$ to classify colour neutral composites (not $U(r+s)$, which gives an odd number of neutral states with no solution):

**1) With both left and right handed neutrinos:**

$$r^2 + s^2 -1 = 4N$$

**2) With only left handed neutrinos:**

$$r^2 + s^2 -1 = 2N$$

The second option can be discarded: combined with the previous equations, it does not produce an integer solution. **The sBootstrap can be extended to the lepton sector only if we use the extra degrees of freedom that come with the addition of right-handed neutrinos.**

Furthermore, the extension produces a system of equations with **an unique integer solution, the one with $N=3$**.

---

## Group-Theoretic Perspective

The scalars are produced from 3 "down quarks" and 2 "up quarks" according the decomposition of a flavour group $SU(5)$ to $SU(3) \times SU(2)$.

Sleptons from **24** of flavour group:

$$24=(1,1)+(3,1)+(2,3)+(2, \bar 3)+(1,8)$$

Squarks from the sextets of **15**:

$$15=(3,1)+(2,3)+(1,6)$$

And similarly the anti-squarks.

---

## Conclusion

- The sBootstrap requisites predict that the number of fermion families must be **three**.
- The solution uses $N=3$ generations with $r=3$ light down-type quarks and $s=2$ light up-type quarks.
- The construction requires right-handed neutrinos.
- The solution group-theoretically corresponds to the decomposition of $SU(5) \to SU(3) \times SU(2)$.

Having a third generation discovered with one massive "up type" quark, this should be considered a striking postdiction from the theory of dual models, or string theory how they call it nowadays.
