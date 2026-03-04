# Koide equations for quark mass triplets

**Author:** Alejandro Rivero
**Email:** arivero@unizar.es
**Affiliation:** Institute for Biocomputation and Physics of Complex Systems (BIFI), University of Zaragoza, Mariano Esquillor, Edificio I + D, 50018 Zaragoza, Spain

## Abstract

We study Koide equation sequentially in a chain of mass triplets. We notice at least one not previously published triplet whose existence allows to build a quark mass hierarchy from fixed yukawas for top and up quarks. Also, the new triplets are used to build mass predictions either descending from the experimental values of top and bottom, or ascending from the original triplet of charged leptons.

**PACS:** 12.15Ff

---

## Section 1: Origin of Koide formula

An early intuition of the mass of the $d$ quark as the result of Cabibbo mixing with $s$ quark, $m_d\sim \theta_c^2 m_s$, guided the efforts of modellers in the seventies to derive either this relationship or a close one, namely:

$$\tan \theta_c = \sqrt {m_d\over m_s}$$

Particularly, Harari, Haut and Weyers used a discrete symmetry to produce not only the above equation but also two exact equations:

$$m_u=0; \quad {m_d \over m_s} = { {2 - \sqrt 3} \over {2 + \sqrt 3}}$$

His model was a sequential symmetry breaking, using the trivial and bidimensional irreducible representations of the permutation group $S_3$.

The work was promptly critiqued because at the end it was unclear how do they get the mass eigenvalues; they rotate the mass matrix to select a particular representation and then they disregard non diagonal terms. Besides they did not provide an explicit representation for the Higgs sector. As for $m_u=0$, it was not really a problem: other QCD mechanism could happen later, to give it a mass of order $m_d m_s / M_\Lambda$. Still, the result contributed to sustain an interest on the use of discrete symmetry to predict structures in the mass matrix.

Following explicitly this topic, in 1981 Koide attempted another approach, incorporating the symmetry in the context of a model of preons. Given that such preons were components both of quarks and leptons, the result was a formula for Cabibbo angle using the mass of charged leptons and, as a by-product, a formula linking their masses:

$${ (\sqrt m_1 + \sqrt m_2 + \sqrt m_3 )^2\over m_1 + m_2 + m_3 } = \frac 32$$

The formula predicts a tau mass of 1776.96894(7) MeV, perfectly matching the experimental value of $1776.82 \pm 0.16$ MeV. In 1981, the measurement was still $1783 \pm 4$ MeV, so in some sense Koide's work was a real prediction, even if its composite character was ruled out by experiments.

This formula, considered independently of an underlying model, is referred as "Koide formula" or "Koide equation", and a tuple of three masses fulfilling it is called a "Koide triplet" or "Koide tuple".

The masses derived from the Harari-Haut-Weyers work are a Koide triplet.

---

## Section 2: Modern Research on Koide Triplets

The use of a composite model can be substituted by an ad-hoc Higgs sector with discrete symmetries. Such models of "flavons" have been revisited by Koide in later work. In some proposals, more symmetry is added in order to protect Koide equation from renormalisation running.

Only with the SM, the equation is not protected, and its running under the renormalisation group has been studied for leptons as well as for some quark triplets. A detailed study is for instance Xing et al. Generically, it is noticed that the running up to GUT scale corrects the equation by less than a ten per cent in the most extreme case, and that in the case of charged leptons the equation works better for pole masses. The reason for the stability of the formula is that we are considering mass quotients and the correction becomes proportional, for the case of leptons, to $\sqrt{m_\mu/m_\tau}$, and similarly for triplets of quarks with same charge.

### Reformulations of Koide equation

Koide formula has had two relevant rewritings. First, Foot suggested an exact angle between the vector of square roots of the masses and the permutation invariant tuple:

$$(\sqrt m_1, \sqrt m_2, \sqrt m_3) \angle (1,1,1) = 45° $$

More recently, when generalised to neutrinos, it was seen that some fits with negative signs of the square root could be needed, and so nowadays a presentation more agnostic about such signs is preferred:

$$\sqrt m_k = M_0^{\frac 1 2} \left(1 + \sqrt 2 \, \cos \left({2 \pi \over 3} k + \delta_0\right)\right)$$

For the charged lepton, we have then three masses fitted with two parameters:

$$M_{e\mu\tau}=313.8 \text{ MeV}, \quad \delta_{e\mu\tau}=0.222$$

and as we move $\delta_0$, we can produce zero and negative values for some square root.

D. Look and M. D. Sheppeard proposed this formulation in a generalised way, allowing for any real value $\lambda_0$ instead of $\sqrt 2$. Sheppeard noticed that the angles $\delta_0$ have some extra regularity, compared to its value for leptons, when the match is applied to the up type $(uct)$ and the down type $(dbs)$ quarks: $\delta_{e\mu\tau}=3 \delta_{uct}$ and $\delta_{dsb} = 2 \delta_{uct}$.

This insight was revised by Zenczykowski. He notes that the value of $\delta_0$ can be extracted directly from:

$$\begin{align}
{1\over \sqrt 2}(\sqrt m_2 - \sqrt m_1) &= \sqrt{3M_0} \sin \delta_0  \\
{1\over \sqrt 6}(2 \sqrt m_3 -\sqrt m_2 - \sqrt m_1) &= \sqrt{3M_0} \cos \delta_0
\end{align}$$

and then he studies its experimental value independently of the whole formula, confirming:

$$\delta_{e\mu\tau}=3 \delta_{uct}$$

The above equations are related to the initial research on mixing angles, as Koide used the combination of the first two equations to produce Cabibbo angle, while a third condition:

$${1\over \sqrt 3}( \sqrt m_3 +\sqrt m_2 + \sqrt m_1) = \sqrt{3M_0}$$

joined to the previous two, implies Koide formula from a pythagorean condition on the LHS of the three equations:

$$(\text{eq1})^2+(\text{eq2})^2=(\text{eq3})^2$$

The above equations relate to the doublet and singlet irreps of $S_3$, where $\delta_0$ can be understood as a rotation of the basis for the bidimensional irrep.

Generalisations of Koide equation to neutrinos have been addressed in the literature but without a conclusive result. We do not address them in this paper. Last, extensions to six-quark equations were done by Goffinet and more recently by Kartavsev, but a definite method to build tuples with an arbitrary number of particles has not been built.

---

## Section 3: Recent Empirical Findings

Two years ago, Rodejohann and Zhang noticed that if one forgets the prejudice of using equal-charge quarks, then the triplet $(c,b,t)$ matches Koide equation accurately.

We can inquire if this success is repeatable and we can keep connecting triplets of quarks with alternating isospin. Note that already the first solution of Harari et al involved quarks of different charge.

To explore this possibility, we solve the Koide equation to get one of the masses as a function of the other two:

$$m_3= \left((\sqrt m_1 + \sqrt m_2 ) \left(2- \sqrt{3+6 {\sqrt{m_1 m_2} \over (\sqrt m_1+\sqrt m_2)^2}} \right)\right)^2$$

And then iterate a waterfall from experimental values of top and bottom:

$$\begin{align}
m_t &= 173.5 \text{ GeV} \\
m_b &= 4.18 \text{ GeV} \\
m_c(173.5,4.18) &= 1.374 \text{ GeV} \\
m_s(4.18,1.374) &= 94 \text{ MeV} \\
m_q(1.374,0.094) &= 0.030 \text{ MeV} \\
m_{q'}(0.094,0.000030) &= 5.5 \text{ MeV}
\end{align}$$

Experimental $m_c$ and $m_s$ are, respectively, $1.275 \pm 0.025$ GeV and $95 \pm 5$ MeV.

So a new Koide triplet, not detected previously, is apparent: $(s,c,b)$. The reason of the miss, besides the need of alternate isospin, is that the solution in this case has a negative sign for $\sqrt m_s$.

This has also the collateral effect that $-\sqrt m_s, \sqrt m_c, \sqrt m_b$ is almost orthogonal to the leptonic tuple $\sqrt m_\tau, \sqrt m_\mu, \sqrt m_e$. From the view of Foot's cone, they are in opposite generatrices.

That we can choose such negative signs in the search for Koide triplets can be inferred for the rewriting above, as they are produced forcefully for some values of the angle. Furthermore, the formalism allows to justify, in a model independent way, that the sign of the square root can be arbitrary with the only condition of keeping the same choosing in the three equations. Formally, their LHS can be seen as a row sum of the product matrix Q:

$$Q=
\begin{pmatrix}
0 & {1\over \sqrt 2} & - {1\over \sqrt 2 } \\
{2\over \sqrt 6} & - {1\over \sqrt 6} & - {1\over \sqrt 6} \\
{1\over \sqrt 3}& {1\over \sqrt 3}  & {1\over \sqrt 3}
\end{pmatrix}
\begin{pmatrix}
\sqrt{m_3} & 0 & 0 \\
0& \sqrt m_2 & 0\\
0&0 & \sqrt m_1
\end{pmatrix}$$

The matrix Q being constructed so that $QQ^+$ (and $Q^+Q$ too) is a mass matrix having $m_1,m_2,m_3$ as eigenvalues, independently of the choosing of signs (or phase) in the square root of the diagonal matrix. Of course, the argument to impose the pythagorean equation in the sum of each row will be model dependent.

In the obtained chain of equations, we infer that $q$ and $q'$ correspond respectively to the masses of up and down quarks, and so the other two triplets are $(usc)$ and $(dus)$.

A zero mass for the up quark allows to suspect that we are seeing some symmetry which is only slightly displaced. Plus, a situation where the up quark is massless at some stage of the breaking constitutes an escape from the QCD $\theta$ problem.

### Exhaustive search for Koide triplets

We can question if there is some other lost triplet between the 20 possible combinations. From the empirical masses we can check that only two triplets meet Koide equation (for either sign) with a 1% of tolerance: $(cbt)$ and $(scb)$. Next, $(usc)$, $(dsb)$ and $(dsc)$ can fit within a 10%, and the rest of the triplets go worse. If we use instead the GUT-level masses, the five better matches, within a 10% of tolerance, are still the same but now neither of them are in the 1% accuracy, again suggesting that the formula works better in the electroweak scale than at GUT energies.

Comparing the "waterfall" chain with this exploration, the disagreement comes only with the last step, $(dus)$, obviously because of the non zero mass of the up quark in the measured data. But the alternatives really only match Koide equation within a ten percent and break the pattern of chaining, that here was proceeding regularly (generation-wise and isospin-wise).

---

## Section 4: A prediction upstream

In the extra symmetry revised by Zenczykowski, the angle $\delta_0$ of charged leptons was three times the angle for up-type quarks. A similar proportionality is also happening here. For the triplet bottom-charm-strange, the fitted parameters are:

$$M_{scb}=941.2 \text{ MeV}, \quad \delta_{scb}=0.671$$

Here not only the angle, but also the mass seems related by a factor three.

Were we to accept the relationships $M_{scb}=3 M_{e\mu\tau}$ and $\delta_{scb}=3 \delta_{e\mu\tau}$, the set of Koide equations would be the most predictive semi-empirical formula for standard model masses: taking electron and muon as input, we can build the parameters of the lepton triplet and then seed them to produce the $(scb)$ triplet. The chain predicts good values for down, strange, charm, bottom, and top quarks, as well as for the tau lepton. For the top quark, with usual values of $m_e$ and $m_\mu$ as input, the predicted mass is 173.264 GeV, right in the target of current averages from Tevatron and LHC of the pole mass.

For reference, the complete comparison between experiment and prediction is:

| | exp. | pred. |
|---|---|---|
| $t$ | $173.5 \pm 1.0$ | 173.26 |
| $b$ | $4.18 \pm 0.03$ | 4.197 |
| $c$ | $1.275 \pm 0.025$ | 1.359 |
| $\tau$ | $1.77682 (16)$ | 1.776968 |
| $s$ | $95 \pm 5$ | 92.275 |
| $d$ | $\sim 4.8$ | 5.32 |
| $u$ | $\sim 2.3$ | 0.0356 |

Where we do not quote error widths in the prediction side because the input data, electron and muon mass, is known with a precision higher than the rest of the measurements, so all shown digits are exact. The table can be taken more as an illustration of the power of Koide triplets than as a working prediction; the experimental quark masses are not the pole masses, except for the top quark, and the relationships $M_{scb}=3 M_{e\mu\tau}$ and $\delta_{scb}=3 \delta_{e\mu\tau}$ are ad-hoc as explained.

---

## Section 5: Results with Simple Inputs

To understand the factors of three in the last section, it is instructive to examine a case without experimental inputs, considering instead that the Yukawa coupling of the top quark is exactly one and the Yukawa coupling of the up quark is exactly zero. The point of asking for $y_u=0$ is that there is the possibility of having some extra symmetry that becomes restored in this limit.

Recall from the standard model that when $y_t=1$, the mass of the top is 174.1 GeV, via $m_q=y_q \langle v \rangle /\sqrt 2$; we are going to use this proportionality instead of quoting the adimensional value of the yukawa coupling.

We should consider the four possible solutions in each step of the Koide equation. A fast exploration with a computer shows that it is possible to try different combinations of triplets to get the goal $y_u=0$, but that the best matching with data is the case of our sequential chain, and that it is also the case with the greatest gap between top and bottom. Thus we are going to use this fact as a postulate to select the solution, in two steps:

First, we ask the difference between the top and the next quark to be the greatest possible within Koide equation. This amounts to set the angle $\delta_{cbt}=0$. The masses of bottom and charm are then the same, 2.56 GeV. The next triplet has an angle, $\delta_{scb}=60°$ and predicts a mass for the strange quark of 150 MeV. And then the next triplet has an angle of $12.62°$ and predicts a small but non zero mass for the next quark: 0.810 MeV. All these steps are unambiguous.

So this case with a single input, $y_t=1$, already provides a good match with the standard model patterns. Now we proceed to smoothly move the parameters in order to approach this $y_u$ to zero, fixing it at the cost of losing exact angles in the upper mass triplets. We get:

$$\begin{align}
&m_t = 174.1 \text{ GeV} \\
&m_b = 3.64 \text{ GeV} \\
&m_c = 1.698 \text{ GeV} \\
&m_s = 121.95 \text{ MeV} \\
&m_u = 0 \text{ MeV} \\
&m_d = 8.75 \text{ MeV}
\end{align}$$

In the above table, we have chosen by hand a position for lepton masses. This has some value: contemplate simply the lower part of the table, written now as:

$$\begin{align}
m_b &= 4 M_{scb} \\
m_c &= {2 + \sqrt 3 \over 2} M_{scb} \\
m_s &= {2 - \sqrt 3 \over 2} M_{scb}\\
m_u &= 0 \\
\end{align}$$

It can be seen that in this case $M_{e\mu\tau}= M_{usc}= M_{scb}/3$, giving some ground to the ad-hoc use of this factor in the previous section; and also that the angles are one of them three times the other: $\delta_{e\mu\tau}=15°$ and $\delta_{scb}=45°$, then agreeing both with our observation and Zenczykowski's.

---

## Section 6: Discussion

A question about the validity of Koide triplets is if they are equations for the pole masses or for running masses, and if the later, for which scale should they be applied. Looking at the running, it can be argued that they can be used by modellers in any situation, as they are good at electroweak scale and still approximate enough at GUT scale, but that their better agreement seems to happen with pole masses (or, for quarks, for $m_q(m_q)$ values).

This was already known for the lepton triplet, and here we have provided an "upstream prediction" going from the e, $\mu$ pole mass as input to the experimental top mass, which is also a pole mass.

Consider $(s,c,b)$. For $\overline {MS}$ masses in current data, the quotient LHS/RHS of the Koide equation is 0.986. The value contemporary to earlier work was 0.974. Using running masses from earlier work, the quotient at $M_Z$ is 0.949, and at GUT scale it is 0.947.

Pole masses can be more natural on models where mixing is the origin of Koide equation, while running masses are appropriate for models where a symmetry breaking at some scale produces the mass pattern.

### The model question

There is the question of the model itself. We do not aim in this paper to propose such; there are some attempts in the literature to include more symmetry. Now, we can say some words about the fact that we seem to have not one, but four equations in the quark sector.

We believe that $S_4=V_4 \ltimes S_3$ is an interesting candidate for a bigger group. It is the group of rotations of a cube, and in fact we can classify our triplets by drawing such cube with quarks in the faces, such that $u,c,t$ meet in a vertex, with respective opposite faces $b,d,s$. The vertexes, arranged by oppositeness across the diagonal, are:

$$\begin{array}{cccc}
bds & usc & scb & cbt \\
uct & btd & tdb & dus
\end{array}$$

These four diagonal axis have associated $S_3$ subgroups, and each axis has at least one of the triplets we are interested on. It can be noticed that an axis is "overloaded" with the triplet $dus$, and this can become an argument to substitute it by $bds$. Some other substitutions are possible, and while we can not concrete the option without resorting to an explicit model, which is beyond this work, the point of being able to select here, from the 20 possible triplets, the ones having a best matching with data, is an argument to use $S_4$ at least for an initial organisation.

---

## Conclusions

Let's conclude by reviewing what we have found in this paper:

We have done a complete exploration of Koide quark triplets allowing for triplets with quarks of different charge. We have found at least one new exact tuple for Koide equation, for quarks $s,c,b$, and two more, approximate, for $u,s,c$ and $d,u,s$; even if the latter is not really new in the literature and we could be in doubt about substituting it for $d,s,b$.

We noticed that these tuples allow to build a chain of equations, starting from $c,b,t$ from earlier work; we have verified its empirical validity and we have suggested yet two other methods to reproduce the mass scales of the Standard Model fermions: either to fix the yukawas of top and up quark to the natural values 1 and 0 or, independently, to propose a proportionality between the parameters of the charged leptons tuple and those of the $s,c,b$ tuple. Both methods have some virtue:

If the last method is applied using as input the experimental values of two lepton masses, all the quark masses, except the up quark, coincide with the experimental measurements within one or two sigmas at most, and a very good prediction of the top quark pole mass is reached. And in the former method, with fixed $y_t=1$, and $y_u=0$, we could tell that the only input from experiment is Fermi constant (to transform yukawas to mass values) and the number of generations, and that Koide triplets provide all the necessary levels for mass.

---

## Acknowledgments

The author wants to acknowledge C. Brannen, M. Porter and J. Yablon by their support and comments.
