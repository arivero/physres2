# The number of generations in supersymmetric models with preons

**Author:** Alejandro Rivero
**Affiliation:** BIFI - Universidad de Zaragoza, Mariano Esquillor, Edificio I + D - 50018 Zaragoza (Spain)
**Email:** arivero@unizar.es

## Abstract

If the scalar sector of the supersymmetric standard model is built with different pairs of entities (preons) having the same colour and electric charge assignment that common quarks, then the model has exactly three generations.

---

## Introduction

The goal of this brief note is to show the result claimed in the abstract, that has not been mentioned—up to my knowledge—in the published literature, and then proceed further to speculate about the entities that could do the role of such preons.

Our rules, very general, are:

- Each scalar, either squark or slepton, gets its electric and colour charge from a different pair of preons, following a model of $SU(n_f)$ flavour. A useful visualization is to think of such pairs as possible terminations of an open string, but other models of composites are equally acceptable.

- Colour comes from the singlet of $\bf 3 \times \bar 3 = 1 + 8$ in the case of pairs preon-antipreon and from the triplet of $\bf 3 \times 3 = 3 + 6$ in the case of preon-preon pairs (and similarly for antipreon-antipreon). As for the electric charge, it is simply the sum of the charges of the two elements of the pair.

- The preon flavours must be the same that the standard model quarks, i.e. either +2/3 "up-type" or -1/3 "down type".

The question to solve is how many preons do we need to produce exactly an integer number $N$ of generations of squarks and sleptons, with "generation" understood in the usual way, the scalar matter that pairs—via susy—with one generation of standard model fermions.

## The matching problem

Let's first check that there is no solution if we use the same number and kind of flavours in the preon sector that in the reproduced supersymetric model. $N$ 'generations' of such preons of types $u, d$ should produce $N^2$ scalars of charge $+1/3$. And by supersymmetry with the $N$ down anti-quarks, we have $2 N$ anti-squarks with such charge. So $N=2$ for the matching in the "down antiquark" sector.

And they also produce, from pairs of two down-type preons, a total of $N (N+1)/2$ scalars of charge $-2/3$, again to be matched to $2 N$ degrees of freedom in the "up antiquark" sector. So $N=3$ for the matching in this sector and the solutions do not coincide. The mechanism fails, which is no a big surprise.

## General solution with r and s preons

We need thus consider a more general solution, where we allow for a number $r$ of preons of down-type, say $r$, and a number $s$ of preons of up-type. The matching of combinations and squarks, again for $N$ generations in the supersymetrized standard model, is now:

$$\begin{align}
r s &= 2 N \\
{r (r+1) \over 2} &= 2N
\end{align}$$

respectively for $+1/3$ and $-2/3$ particles. The integer solution of the system is that $N$ must be half of an hexagonal number:

$$2 N = \cancel { \, 1}, 6, \cancel {15}, 28, \cancel {45}, 66, \cancel {91}, 120, ...$$

So $N=3$ is the smallest admissible solution, but it is not unique yet.

## Asymptotic freedom constraint

In the particular case that preons were assumed to be fermions bound by a SU(3) interaction, we could consider the requirement of asymptotic freedom in the beta function to put an upper limit to the number of solutions, asking the number of preon flavours $n_f=r+s$ to be such that $2 n_f < 33$.

In this case we have still a very short list of solutions:

| N | s | r | Total $n_f$ |
|---|---|---|----|
| 3 | 2 | 3 | 5 |
| 14 | 4 | 7 | 11 |
| 33 | 6 | 11 | 17 |

Where we have included the third solution only because 17 flavours are very near of the theoretical bound of 16.5. Note that we always have $r = 2 s -1$.

To get uniqueness along this way, we should have more flavours in the bag. It could happen if a particular model of compositeness were demanding, as we did in the initial example, to have the same number of generations in the preons that in the initial theory. In such case it should be needed to justify why only $r$ and $s$ preons enter in the composition; for instance they could be massless and the rest very massive, so they could not be at the ends of a string, or they should disintegrate before arriving to a bound state—a phenomena that we already know from the top quark in the SM.

In any case, we can still search for another source of uniqueness looking at the slepton sector, which we have not considered until now.

## The slepton sector

For the charged leptons, we have nothing new. Every charged scalar lepton is a composite of a preon and antipreon of different type, so $ r s = 2 N$.

For neutral leptons, we have some options. On the Standard Model side, we could consider the known degrees of freedom of neutrinos, only as left-handed particles, or add a right-handed neutrino in each generation.

On the composition side, note that in the initial requeriments we have fixed the question of what group must be use to classify the colour neutral composites. We could try $U(r+s)$, having an extra neutral state, but then we had an odd number of neutral states and no solution. Then we use $SU(r+s)$ and the total number of neutrals will be one less than in the $U(r+s)$ case.

So we have finally two possible equations:

**1) With both left and right handed neutrinos:**

$$r^2 + s^2 -1 = 4N$$

**2) With only left handed neutrinos:**

$$r^2 + s^2 -1 = 2N$$

And we can discard the second option: when combined with the previous equations, it does not produce an integer solution. Interestingly, we can build the slepton sector only if we use the extra degrees of freedom that come with the addition of right-handed neutrinos. With it, the match to sleptons produces a system of equations with an unique integer solution, the one with $N=3$, which was the result we were aiming to.

So, even if we do not want to use QCD beta-function as an argument, we also have an unique solution when including, as we should, all the sfermions in our scheme. And it is a very interesting one: we have $N=3$ generations of quarks and leptons, and the solution uses $r=3$ preons of type "down" and $s=2$ preons of type "up".

## Group-theoretic perspective

Looking the solution group-theoretically, we can say that the scalars are produced according the decomposition of a flavour group $SU(5)$ to $SU(3) \times SU(2)$.

We get the sleptons extracted from a **24** of this flavour group:

$$24=(1,1)+(3,1)+(2,3)+(2, \bar 3)+(1,8)$$

And the squarks from the two sextets that appear in the decomposition of **15**:

$$15=(3,1)+(2,3)+(1,6)$$

And similarly the anti-squarks.

## The sBootstrap condition

Finally, let's speculate on some principle that could produce the initial rules. Let's assume that the preons are fermions. Three fermions of charge $-2/3$ and two fermions of charge $+1/3$ fail the condition of anomaly matching, and they need to be extended with more particles, that will not participate in the construction of the scalar sector.

The most obvious extension is the traditional, non susy, standard model, where each generation is anomaly-free by itself. So we can use three generations of preons similar to the SM, with "leptonic preons" not having the "strong" interaction needed to participate in the composites. With this extension, the only issue is that we use $s=2$ "up preons" and not three; so we must produce this effect by asking that the "top" preon is special. For instance we can follow the argument we have invoked before and tell that all the preons could be massless except the "top" preon, so that it would not be allowed at the end of a "colour" string, while all the other five could be terminations of such string.

Having this extension, the second speculation follows, and makes the first one a lot more elegant: if they are so similar, could it be that the preons of the susy scalar sector and the quarks of the Standard Model are actually the same particles?

We can call this idea the **sBootstrap**, because if the scalar sector can produce itself from the fermion sector, and then the fermion sector can be produced from the scalar sector with the susy generator, the SSM matter bootstraps itself. It is not as democratic as the original idea from Chew, but it keeps some of its spirit.

And it is also in the spirit of the first models with supersymmetry, that were called "quark-gluon" or "fermion-meson" dual models, because at that time the scalar sector was for sure a composite (the strongly interacting mesons). Of course, having a third generation discovered with one massive "up type" quark, this should be considered a striking postdiction from the theory of dual models, but it was a very bold step at that time. At most, this early antecendent justifies our imaginery of a colour-binding string as the excuse to use only pairs of preons. And we can keep ourselves agnostic about if quarks and leptons are composed of themselves or other particles, or if they are a truly elementary sector.

Also, the sBootstrap condition resurrects our first source of uniqueness, asymptotic freedom, as now the number of generations is always $N$, and then the two solutions of $N=14$ and $N=33$ can be discarded even without using the matching of the number of sneutrinos. Still, it is reasonable to keep asking the reconstruction of all the scalar matter sector, and not only the coloured part.

## Conclusions

To conclude:

- A very general class of preon models implies the prediction of only 3 generations of the standard model, with the extra bonus of asking for right-handed neutrinos.

- And besides, the more speculative 'sBootstrap condition' also predicts that the top quark must be different of the others in some aspect related to its ability to form bound states. For a relativistic string model of binding, this condition can be that in some limit the top quark is massive and the rest are massless.
