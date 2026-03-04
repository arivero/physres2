# Bootstrapping some structure of the Standard Model via Supersymmetry

**Author:** Alejandro Rivero (arivero@unizar.es)

## Abstract

Aiming to select special sections of the string landscape, we suggest a simple self-referencing postulate in the Chan-Paton charges that fixes some of the freedom of model building.

**Note:** This is a draft to clarify this idea and eventually incorporate it to the literature in some letter journal. The initial proposal is ten years old, so probably I am assuming a lot of details not written, and writing a lot of details too trivial for readers, so any feedback is specially welcome.

---

## Introduction

Bootstrapping, as it is known, refers to the possibility of determining some parameters—usually of a scattering process—via some set of consistency conditions and a self-referential mechanism—usually a duality. This was expected to generate some non-linearity pinpointing special places of the space of parameters. Back in 1970, such program was still considered unable to pinpoint the quantum numbers of a model, as Chew himself remarked in an article in Physics Today:

> "The unfriendly question raised most often ... is how self-consistency can possibly be expected to generate "internal quantum numbers" such as hypercharge and baryon number..., but how can you bootstrap a symmetry? A conceivable response is that symmetries (or the associated quantum numbers) are related to particle multiplicities, and the non-linear unitarity condition responds to the number of different particles."

Nowadays the progress on the analysis of anomalies makes us more knowledgeable: we know that we need complete generations of the model, that there is a restriction in the total charge, and that if we want to use more sophisticated models, as in string theory, we are limited in our options to choose a gauge group. Today it is more likely to fix a gauge group via consistency conditions than to bootstrap mass relations or coupling constants, but even here the path has been barred: from the bottom we are puzzled by the number of generations, and from the top by a landscape of possible options.

In this note we want to point out a possible principle, in the frame of open superstrings, with a self-referential argument not far from the old idea of nuclear democracy, where all the particles were considered equally elementary and composite. Here it is going to be enough to demand all the bosonic particles to be composite, in the sense of having a Chan-Paton group, and in a way such that the labeling references to the same model.

## SO(10) Chan-Paton Structure

Let's start by noticing that an unoriented bosonic string with SO(10) Chan-Paton labels contains a precursor of the susy scalars of the standard model, already with three generations. The 54 of SO(10) goes down to *flavour* $SU(5) \times U_1(1)$:

$$54 = {15}(4) + \bar{15}(−4) + {24} (0)$$

And each representation goes down to (flavour) $SU(3) \times SU(2) \times U_2(1)$:

$$\begin{align}
15 &= (1, 3)(−6) + (3, 2)(−1) + (6, 1)(4) \\
24 &= (1, 1)(0) + (1, 3)(0) + (3, 2)(5) + (\bar 3, 2)(−5) + (8, 1)(0)
\end{align}$$

In this model we look only to electric charge. Setting the two $U(1)$ charges to $Q_1+Q_2=-\frac 16$, $Q_2=-1/5$ we can produce a total combination $Q$ as in the following table. Note that $N$ here is the number of states.

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
| (1,1) | — | — | — | — |
| (1,3) | 12 | 0 | 0 | 0 |
| (8,1) | — | — | — | — |

The table contains, as promised, three generations of scalars with the right electric charge, plus three "half-generations" of an extra quark. To propose an identification for the components of the $10$ of $SO(10)$ we can also branch the 10-plet of the group first to $SU(5)$, with $10 = 5 (2) + \bar 5 (-2)$ and then to $SU(3) \times SU(2)$, via $5 = (3, 1)(2) + (1,2) (-3)$. We see that the 10-plet has three elements of $Q=\frac {2}{30} - \frac {2}{5} = -\frac 13$, two elements of $Q=\frac {2}{30} + \frac {3}{5} = \frac {23}$, and then the corresponding opposite charges. Abusing a bit on traditional notation, we can call $d,s,b$ to the elements of the triplet and $u,c$ to the doublet, and draw the whole construction as showing that the bosons are actually being generated from pairing "quarks" or "antiquarks".

## Self-Referential Postulate

We need to use both these "quarks" and the opposite "antiquarks" as labels of the Chan-Paton group, and to request unorientability to produce the 54 states (in principle, we discard the singlet) that reproduce a theory of mesons, quarks and diquarks. Of course we could then again discard the unoriented states to reduce to a theory of mesons with SU(5) group, which in this case would be the three generations of scalar leptons.

Looking for a postulate for this structure, the idea that we propose is to assume that there exists a supersymmetric theory containing this same set of scalar states, and that there is a nuclear democracy in the labeling: the Chan-Paton charges of this theory are a subset of the quarks. We can call them the light quarks or, indistinctly, the preons of the theory.

So, consider $r$ quarks of $+2/3$ charge and $s$ quarks of charge $-1/3$. Our postulate is that they combine pairwise—say, at the extremes of an open string—to form consistent generations of squarks and sleptons.

## Bootstrap Constraints

We can formulate this requisite in two ways:

- We can ask that the combination must build the same number of "up type" and "down type" diquarks. This is, that $rs = s(s+1)/2 = 2n$ and so $s= 2 r -1$. Plus, we can ask $r$ to be even, to make sure we build a pair number of scalars of each charge (we want to be able to promote them to a supersymmetric theory in the future).

- Or, we can go for a more general condition, that we will see implies the former one: we can ask that the total number of combinations must be an integer multiple of a single set, this is a multiple of $rs$. In this case we look for integer positive solutions of:

$$rs + {s^2+s \over 2} + {r^2 + r \over 2} + (r+s)^2 -1 = K rs$$

and then $K=9$. This in turns implies that one generation is composed of $(K-1)/2$ Dirac-like tuples—including the neutrino—plus one extra state, with its antiparticle. This extra state can be either neutral or charged with +4/3.

In any case, either requirement fixes $s=2 r -1$ and then the allowed groups are $SO(2s+2r)=SO(6r -2)$, with $r$ even. The smallest possible group is SO(10), but bigger groups are possible.

We need an extra postulate to fix $n_g=3$. To look for it, lets examine a table of solutions, with the corresponding total of "standard model generations" of scalars and its extra content:

| $r$ | $s$ | $n_g$ | $\frac{r(r+1)}{2}$ | $\frac{r^2+s^2 -1}{2} - 2 n_g$ |
|-----|-----|-------|-----|-----|
| 2 | 3 | 3 | 3 | 0 |
| 4 | 7 | 14 | 10 | 4 |
| 6 | 11 | 33 | 21 | 12 |
| 8 | 15 | 60 | 36 | 24 |

From the table, it seems that the simplest extra postulate is to ask for absence of excess neutral bosons. An equivalent requisite is to ask that the number of neutral bosons be equal to the number of charged ones, given that the latter is always $2 r s$ as in the "quark" sector.

## Color Recovery

As for colour, it can be recovered if one considers to assign a SU(3) singlet to the 24 of SU(5) and SU(3) triplets and anti-triplets to the $15$ and $\bar 15$. The decomposition should then be promoted to something under $SO(10) \otimes SU(3)$ such as:

$$\bf ? \longrightarrow (15,3) + (\bar 15, \bar 3) + (24, 1)$$

In the sense of being strings, and thus related to a large $SU(N)$ limit, all the states are colored; the triplet or singlet state being associated to the planarity of the string. The pairing particle/antiparticle is used to recover an oriented sector. What we need is some argument to pinpoint exactly $SU(3)$, either as a preferred labeling or extracting it out of a $SU(3 + N)$ in large $N$.

We argue that a bootstrap of $N=3$ is derived from the fact that:

$$\bf 3 \times 3 = 6 + 3$$

so in some sense a composite $3$ is composed of itself. The mesonic strings, on the other hand, go to the singlet from:

$$\bf 3 \times \bar 3 = 8 + 1$$

and the colouring paints the matrix as desired.

After adding colour, it can be seen more clearly the process:

$$\text{oriented} \to \text{unoriented} \to \text{reoriented}$$

that allow to extract out a set of mesons with Chan-Paton flavour charges under $SU(r+s)$. This technique is used by Armori et al under the market name "planar orientifold" for the calculation of meson masses in a large $N$ theory generating an "Hadronic Susy". Note that the diquark, unoriented part, can not survive in the large $N$, making already a sort of confinement.

Amusingly, colouring the scheme has expanded the original 10-plet to a 30-plet, and so from $SO(10)$ to $SO(30)$, and we are really only one "uncolored preon" out of reaching $SO(32)$.

## SO(32) Decomposition

Even if we had no results about consistency of gauge theories in strings, we could have considered, by symmetry reasons, to leave the colour boxes to overflow the diagonal by adding singlets, this is, to consider $U(3)$ instead of $SU(3)$ and the $\bf 54+1$ of $SO(10)$ so that the decomposition now should be something as:

$$\begin{align}
\bf
&(15,3)+(15,6)+(\bar {15},\bar 3)+(\bar {15},6)+(24,8)+(24,1)+(1,8)+(1,1) \\
&= 45+90+45+90+192+24+8+1 \\
&= 495
\end{align}$$

So we see that a super-bootstrapped three family structure supplemented with three colours has about the right size to potentially reach the Green-Schwarz cancellation. Note that honestly we are still in the world of the bosonic string; we can take this number as sheer coincidence or as signal that the "sBoostrap" postulate actually has some susy in it. In any case, it isolates a pathway for model building with exactly three generations.

Actually the sum is a bit of luck; a more accurate decomposition, for instance, is to branch out from $SO(32)$ to $SU(5)\times SU(3)$ via say $SO(30)$ and $SU(15)$. Such pathway decomposes:

$$\begin{align}
496 &= (435)_0 ⊕ (30)_2 ⊕ (30)_{−2} ⊕ (1)_0 \\
&= (224)_{0,0} ⊕ (105)_{0,4} ⊕ (105)_{0,−4} ⊕ (1)_{0,0} \\
&\quad ⊕ (15)_{2,2} ⊕ (15)_{2,−2} ⊕ (15)_{-2,2} ⊕ (15)_{-2,−2} \\
&= (24, 8)_{0,0} ⊕ (24, 1)_{0,0} ⊕ (1, 8)_{0,0} ⊕ \\
&\quad (15, 3)_{0,4} ⊕ (10, 6)_{0,4} ⊕ (15, 3)_{0,−4} ⊕ (10, 6)_{0,−4} ⊕ (1)_{0,0} \\
&\quad ⊕ (5, 3)_{2,2} ⊕ (5, 3)_{2,−2} ⊕ (5, 3)_{-2,2} ⊕ (5, 3)_{-2,−2}
\end{align}$$

In the current context, as we are not including the weak force in the model, there is no value on locating a concrete descent. Also, we ignore if there is some more direct way to relate $SO(10)$ to $SO(32) = SO(2^{D/2})$; it is a bit intriguing that adding colour seems to be enough.

## Conclusion

In conclusion, let's review what we have got. We have started from the postulate that a bosonic unoriented string could be constructed such that:

- The "preons" labeling the string are of two kinds according electric charge: "up" and "down".
- In the diquark sector, the number of "up" diquarks is equal to the number of "down" diquarks.
- In the meson sector, the number of charged mesons is equal to the number of neutral mesons.

We call this arrangement "a supersymmetrical bootstrap" because in a field theory of quarks and mesons the preons should be a subset of the quarks of the theory, assuming that such quarks are recovered via susy from the bosonic sector.

The main conclusion is that there is a consistent assignment of charges, and the only possible number of generations is three.

Then we add colouring compatible with the construction: a meson must be a singlet of colour, a diquark must be in the same representation of colour than a (anti-)quark. This isolates $SU(3)$ and strongly signals that the full group must be $SO(32)$.

The complete program to build such model should imply to start from this bosonic "hadron-diquark" string model, upgrade to superstrings then going to D=10, and then analyze the branes corresponding to Chan-Paton matrices. Note that it would still be true that removing coloured sectors the hadronic mesons of the standard model are recovered, so as a minimum the approach meets the goals of the original program of 1970!

It is interesting that the number of quarks actually involved in the construction is one less than the number of quarks obtained in the process, and that the one odd out is of "up" type. So it could be expected that the process of building the model allows, or even forces, a mass for the top quark separated from the other five.
