# Bootstrapping charges via Supersymmetry

**Author:** Alejandro Rivero (arivero@unizar.es)

## Abstract

A self-referencing postulate in the Chan-Paton factors of a string allows to fix some of the freedom in the choosing of standard-model-like vacua.

---

## The View from Above: SO(32) Decomposition

Consider $SO(32)$. We can factor it as flavour times colour, $G^F \times SU(3)$, following the procedure of Gell-Mann, Ramond and Slansky (case 4, equation 2.18):

$$\begin{align}
(32 \times 32)_A = &[1,1,1^c] + \mathbf{(1,24,1^c)} + (1,1,1^c) + (1,1,8^c) \\
&+ (2,5,3^c) + (2,\bar 5, \bar 3^c) \\
&+ \mathbf{[1,15,\bar3^c]} + \mathbf{[1, \bar {15}, 3^c]} \\
&+ [1,10,6^c] + [1, \bar{10}, \bar 6^c]
\end{align}$$

where the derived flavour group is $G^F=SO(2) \times U(5) \times U(1)$. The same structure plus an extra $(1,28,8^c)$ can be derived branching down via:

$$SO(32) \supset SO(30) \times U(1) \supset SU(15) \times U(1) \times U(1) \supset SU(5) \times SU(3) \times U(1) \times U(1)$$

from where, following the branching rules, we get:

$$\begin{align}
\mathbf{496} = &(24, 8^c)_{0,0} \oplus (1,1^c)_{0,0} \oplus \mathbf{(24, 1^c)_{0,0}} \oplus (1,1^c)_{0,0} \\
&\oplus (1, 8^c)_{0,0} \oplus (5, 3^c)_{2,2} \oplus (5, 3^c)_{2,−2} \oplus (5, 3^c)_{-2,2} \oplus (5, 3^c)_{-2,−2} \\
&\oplus \mathbf{(15, 3^c)_{0,4} \oplus (15, 3^c)_{0,−4}} \oplus (10, 6^c)_{0,4} \oplus (10, 6^c)_{0,−4}
\end{align}$$

The $\mathbf{24}$ descends as:

$$24 = (1, 1)(0) + (1, 3)(0) + (3, 2)(5) + (\bar 3, 2)(−5) + (8, 1)(0)$$

and with $Q=\frac 15 Y$ it looks as three generations of scalar leptons. Furthermore, the $\mathbf{15}$ descends as:

$$15 = (1, 3)(−6) + (3, 2)(−1) + (6, 1)(4)$$

and extending the $U(1)$ charge assignment the total $\mathbf{15+\bar{15}+24}$ can be interpreted as three generations of scalar quarks and leptons:

| irrep | N | $Y_1$ | $Y_2$ | $Q = \frac{1}{30}Y_1 - \frac{1}{5}Y_2$ |
|-------|---|-------|-------|-----|
| (6,1) | 6 | 4 | 4 | -2/3 |
| (3,2) | 6 | 4 | -1 | +1/3 |
| (1,3) | 3 | 4 | -6 | +4/3 |
| $(\bar 6,1)$ | 6 | -4 | -4 | +2/3 |
| $(\bar 3,2)$ | 6 | -4 | 1 | -1/3 |
| $(1,\bar 3)$ | 3 | -4 | 6 | -4/3 |
| $(\bar 3,2)$ | 6 | 0 | -5 | +1 |
| (3,2) | 6 | 0 | 5 | -1 |
| (1,3) | 12 | 0 | 0 | 0 |

This is possibly the most straightforward way to obtain three generations of standard model colour and electric charge out of a string motivated group. This is probably left unnoticed because one wants to get also the chiral charge, and then one looks for spinors in complex representations by exploring groups $SO(4n+2)$, which excludes $SO(32)$.

---

## The View from Below: Self-Referential Bootstrap

The goal is to point out a possible principle, in the frame of open superstrings, to bootstrap this particular group and charge from a self-referential argument.

Consider the $54=\mathbf{15+\bar{15}+24}$ as produced by an unoriented bosonic string with $SO(10)$ or $Sp(10)$ Chan-Paton labels. We need unorientability of the string to produce all the 54 states forming a set of mesons, quarks and diquarks. We could then discard the unoriented states to reduce to a theory of mesons with SU(5) group—which in this case would be the three generations of scalar leptons.

The postulate: there exists a supersymmetric theory containing this same set of scalar states, and there is a nuclear democracy in the labeling—the Chan-Paton charges are a subset of the quarks. We call them the light quarks or, indistinctly, the preons of the theory.

So, consider $r$ quarks of $+2/3$ charge and $s$ quarks of charge $-1/3$. Our postulate is that they combine pairwise to form consistent generations of squarks and sleptons.

### Bootstrap constraints

**Condition 1** (diquark balance): The combination must build the same number of "up type" and "down type" diquarks:

$$rs = s(s+1)/2 = 2n \Rightarrow s= 2r-1$$

**Condition 2** (general counting): The total number of combinations must be an integer multiple of $rs$:

$$rs + {s^2+s \over 2} + {r^2 + r \over 2} + (r+s)^2 -1 = K rs$$

yielding $K=9$. This implies one generation is composed of $(K-1)/2$ Dirac-like tuples plus one extra state (either neutral or charged with +4/3).

Either requirement fixes $s=2r-1$ and then the allowed groups are $SO(2s+2r)=SO(6r-2)$, with $r$ even. The smallest possible group is SO(10).

### Uniqueness condition

From the table of solutions:

| $r$ | $s$ | $n_g$ | $\frac{r(r+1)}{2}$ | excess neutral |
|-----|-----|-------|-----|-----|
| 2 | 3 | 3 | 3 | 0 |
| 4 | 7 | 14 | 10 | 4 |
| 6 | 11 | 33 | 21 | 12 |

The simplest extra postulate is to ask for **absence of excess neutral bosons** (number of neutral bosons = number of charged ones = $2rs$). With this extra requisite, the solution is unique:

$$r=2, \quad s=3, \quad n_g=3, \quad G=SO(10)$$

---

## Colour Recovery and SO(32)

Colour can be recovered by assigning SU(3) singlet to the 24 of SU(5) and SU(3) triplets/anti-triplets to the $15$ and $\bar{15}$:

$$\mathbf{?} \longrightarrow (15,3) + (\bar{15}, \bar 3) + (24, 1)$$

Self-referentially, $N=3$ colour is derived from:

$$\mathbf{3 \times 3 = 6 + 3}$$

so a composite 3 is composed of itself. The mesonic strings go to:

$$\mathbf{3 \times \bar 3 = 8 + 1}$$

The process is: $\text{oriented} \to \text{unoriented} \to \text{reoriented}$

Colouring expands the 10-plet $dsb \; uc \; \bar d\bar s\bar b \; \bar u\bar c$ to a 30-plet $d_{rgb}s_{rgb}\ldots$, going from $SO(10)$ to $SO(30)$, and we are only one "uncolored preon" out of reaching $SO(32)$.

### Green-Schwarz counting

Adding $U(3)$ instead of $SU(3)$:

$$\begin{align}
&(15,3)+(15,6)+(\bar{15},\bar 3)+(\bar{15},6)+(24,8)+(24,1)+(1,8)+(1,1) \\
&= (54+1,9) = 495
\end{align}$$

A super-bootstrapped three-family structure supplemented with three colours has about the right size to potentially reach the Green-Schwarz cancellation.

---

## Complete State Listing

**From the $24+1$ (meson/slepton sector):**
- Charge +1, irrep $(3,2)$: $u\bar d, u\bar s, u\bar b, c\bar d, c\bar s, c\bar b$
- Charge -1, irrep $(\bar 3,2)$: $\bar u d, \bar u s, \bar u b, \bar c d, \bar c s, \bar c b$
- Charge 0: $u\bar u, c\bar c, u\bar c, c\bar u, d\bar s, \ldots$

**From the $\bar{15} \oplus 15$ (diquark/squark sector):**
- Charge $\pm 4/3$: $uu, cc, cu$ and conjugates
- Charge $\pm 2/3$: $\bar d\bar d, \bar s\bar s, \bar b\bar b, \ldots$ and conjugates
- Charge $\pm 1/3$: $ud, us, ub, \ldots$ and conjugates

---

## Conclusion

The postulates:
1. Preons are of two kinds by electric charge: "up" (+2/3) and "down" (-1/3)
2. In the diquark sector: number of "up" diquarks = number of "down" diquarks
3. In the meson sector: number of charged mesons = number of neutral mesons

**Main result**: There is a consistent assignment of charges, and the only possible number of generations is three. Adding colour compatible with the construction isolates $SU(3)$ and strongly signals the full group must be $SO(32)$.

It is interesting that the number of quarks involved in the construction is one less than the number of quarks obtained, and that the one odd out is of "up" type — suggesting the process forces a mass for the top quark separated from the other five.
