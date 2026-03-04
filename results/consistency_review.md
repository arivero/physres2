# Consistency Review: sbootstrap_v4d.tex

Reviewed: 2026-03-04
File: `/Users/arivero/phys3/sbootstrap_v4d.tex` (2354 lines)

---

## 1. Notation Consistency

### 1.1 The symbol Q is overloaded with three distinct meanings

The letter Q is used for three unrelated quantities without disambiguation:

- **Koide ratio** (used throughout most of the paper): Q = 2/3 (Sections 2-5) or Q = 3/2 (Section 6, mass chains). This is the primary meaning.
- **Electric charge** Q_em: Appears in the charge census (line 920, "charges via Q = T_3 - Y/6"; lines 950-952, "Q = +2/3" and "Q = -2/3" referring to fractional electric charges of quarks; line 965, "states at Q = 0"). Also used in the D-term section (lines 1723-1747) as Q_em.
- **SUSY generator**: Lines 1074 and 1081 use "Q . ell" and "Q . q_i" where Q is the bosonic SUSY generator (Miyazawa supercharge).

The electric-charge usage at lines 950-952 is particularly confusing because "Q = +2/3" and "Q = -2/3" could easily be misread as Koide ratios, especially since Q = 2/3 is the Koide target value discussed extensively just pages earlier.

**Recommendation:** Use a distinct symbol for the SUSY generator (e.g., script-Q or calligraphic-Q) and always write Q_em for electric charge to avoid ambiguity with the Koide parameter Q.

### 1.2 The Koide ratio Q = 2/3 versus Q = 3/2 convention switch

The paper uses Q = 2/3 in Sections 2-5 (with Q = sum(m) / (sum(sqrt(m)))^2) and then switches to Q = 3/2 in Section 6 (with Q = (sum(sqrt(m)))^2 / sum(m)). A footnote at lines 454-463 attempts to explain the switch, but contains a mathematical error.

**The footnote states:** "The two conventions are related by Q_here = 3 Q_prev" (line 455).

**This is incorrect.** Since Q_here = (sum(sqrt(m)))^2 / sum(m) and Q_prev = sum(m) / (sum(sqrt(m)))^2, these are reciprocals: Q_here = 1/Q_prev. Numerically, 1/(2/3) = 3/2, which is correct. But "3 * Q_prev = 3 * (2/3) = 2", not 3/2. The stated relation Q_here = 3 Q_prev is arithmetically wrong.

After Section 6, the paper switches back to the Q = 2/3 convention for the remaining sections (e.g., lines 1920, 2044, 2239). The convention is self-consistent within those sections; only the transition footnote is erroneous.

### 1.3 M_0 versus v_0^2

The scale parameter M_0 from the Koide angle parametrization (line 512: M_0 = 313.84 MeV) and the vacuum energy v_0^2 from the seed-triple discussion (line 242: 314 MeV for leptons) are identified as the same quantity on line 514: "The scale M_0 ~ 314 MeV is the vacuum energy v_0^2 of Sec. 3." This identification is stated clearly and used consistently. However, the rounding differs: M_0 = 313.84 MeV (exact fit) versus v_0^2 = 314 MeV (rounded in tables). This is acceptable but worth noting for precision-conscious readers.

### 1.4 The Casimir is written consistently as C_2 = s(s+1)

The Casimir eigenvalue equation (Section 7, Part B) consistently uses C_2 = s(s+1) throughout (lines 1318-1352). No variation of the form "C = s(s+1)" appears. The notation C_2(G) for the group Casimir (line 1010, in the asymptotic-freedom proof) is a different quantity -- the group-theoretic Casimir for the adjoint representation -- and is correctly distinguished.

### 1.5 Mass units are generally consistent

The paper uses MeV throughout for particle masses and vacuum energies, with GeV used only for electroweak-scale quantities (M_W, M_Z, v_EW, the heavy-sector v_0^2, and the Casimir scale parameter m = 106.577 GeV). This is consistent.

One item worth flagging: the heavy-sector vacuum energy appears as:
- Line 310 (table): v_0^2 = 29,580 MeV
- Lines 387, 1799, 2052 (text): v_0^2 ~ 29.6 GeV
- Line 1813: (v_0^{(cbt)})^2 ~ 29,600 MeV

29,580 MeV vs. 29,600 MeV is a rounding discrepancy of 0.07%. The table value is more precise; the text rounds up. Not a real error, but the two should match.

### 1.6 The R value in the Casimir section

The electroweak ratio R = 0.2231014... is stated at line 1366 and used consistently:
- Line 1392: sqrt(1 - 0.22310) in the W-mass prediction
- Line 1375: sin^2(theta_W) = 0.22306 (the experimental value, correctly distinguished from R)
- Line 1414: same experimental value restated

These are internally consistent. The comparison between R (algebraic prediction) and sin^2(theta_W) (measured) is correctly handled.

---

## 2. Cross-Reference Integrity

### 2.1 All \ref{} and \eqref{} targets exist

Every cross-reference in the document was checked against the defined labels. All targets exist:

| Reference target | Defined at line |
|---|---|
| sec:intro | 58 |
| sec:koide_algebra | 99 |
| sec:seeds | 146 |
| sec:two_tuples | 255 |
| sec:mesonkoide | 351 |
| sec:cbt | 381 |
| sec:mass_chains | 430 |
| sec:counting | 737 |
| sec:diquarks | 1050 |
| sec:two_structures | 1172 |
| sec:breaking | 1470 |
| sec:electron | 2096 |
| sec:restoration | 2132 |
| sec:chirality | 2165 |
| sec:open | 2195 |
| sec:conclusions | 2233 |
| All subsec: labels | verified present |
| All eq: labels | verified present |
| All tab: labels | verified present |

**No broken references found.**

### 2.2 Uncited bibliography entries

Two bibliography entries are defined but never cited in the text:

- **\bibitem{GOR1968}** (Gell-Mann, Oakes, Renner) -- line 2333. The GOR relation is discussed extensively in the paper (lines 1673-1674, 2017-2018) but the \cite{GOR1968} command never appears. The initialism "GOR" appears at lines 343, 1537, 1673, 1691, 2017 without citation.

- **\bibitem{FramptonGlashow1987}** (Frampton and Glashow, chiral color) -- line 2348. The "chiral color" concept appears in the chirality table (line 2180: "Chiral color? ... Axigluon mass") but no \cite command references this entry.

These entries will generate "unused bibitem" warnings during compilation.

---

## 3. Redundancy

### 3.1 The $(c,b,t)$ vacuum energy is stated three times in running text

The vacuum energy of the $(c,b,t)$ triple is given as v_0^2 ~ 29.6 GeV at:
- Line 387 (Section 5, "The Third Tuple")
- Line 1799 (Section 8, "SUSY Breaking," Layer 3)
- Line 2052 (Section 8, "What the (c,b,t) triple tells us")

All three are in running text, not tables. The first two restatements are in different sections and are arguably helpful for standalone reading. The third (line 2052) is within the same section as line 1799 and appears redundant.

### 3.2 The meson Koide triple precision is stated in multiple places

The precision of the (-pi, D_s, B) triple (Q = 0.6674, 0.10% deviation) appears at:
- Line 308 (precision summary table in Section 4)
- Line 359 (Section 5, "The Meson Koide Triple," stated as 0.6674 and 0.10%)
- Line 2242 (Conclusions, stated as "to 0.1%")

This is a case where the precision summary table in Section 4 and the dedicated Section 5 present the same number. The Section 5 treatment adds the energy-balance ratio and the comparison to quark triples, so it is not pure duplication, but readers may wonder why the precision summary precedes the detailed discussion.

### 3.3 The Pauli theorem conclusion is stated four times

The conclusion that "the sBootstrap diquark in the 15 cannot be a qq composite" appears at:
- Lines 1143-1144 (Diquark section, Theorem statement)
- Lines 1166-1168 (Diquark section, interpretation)
- Lines 2199-2205 (Open Problems)
- Lines 2266-2271 (Conclusions)

The Open Problems and Conclusions restatements are reasonable for a summarizing role, but the Diquark section itself states the same conclusion twice (once as the theorem consequence, once as the final paragraph's summary).

### 3.4 The W-mass prediction

The W-mass prediction M_W = 80.374 GeV appears only in the Casimir section (line 1393). It is not restated in the Open Problems or Conclusions. This is arguably an omission (see Section 5 below) rather than redundancy, but it means the Casimir section's main quantitative result is only stated once.

### 3.5 The seed triple vacuum energy table

The vacuum energies of the seed triples appear in:
- Lines 236-246 (Section 3, the seed vacuum energy table)
- Lines 262-272 (Section 4, the seed-to-full-triple table)
- Lines 296-312 (Section 4, the precision summary table)

All three tables share overlapping data (the v_0^2 values for the seed triples). Having three similar tables in quick succession is a style choice, but a reader covering 2 pages will see the same numbers three times.

---

## 4. Style Mismatches

### 4.1 Equation display style

Sections 2-5 and 9-12 use unnumbered display equations liberally (bare \[ ... \] environments), while Section 6 (Mass Chains) numbers nearly every equation and uses labeled equations (eq:koide_angle, etc.) more carefully. Section 7 (Counting) is extremely rigorous with theorem-proof style. Section 8 (SUSY Breaking) is expository and uses fewer numbered equations. The shift in rigor is noticeable between Section 7 (formal proofs) and Sections 8-12 (narrative style).

### 4.2 Table formatting

Tables in Sections 3-4 use the \begin{center}\begin{tabular} environment (informal placement). Tables in Sections 6-7 use the \begin{table}[ht]\centering environment with \caption and \label (proper float placement). This difference means the early tables will not appear in a list of tables and cannot be cross-referenced, while the later tables can be. Specifically:

- Lines 236-246: center environment, no label, no caption
- Lines 262-272: center environment, no label, no caption
- Lines 296-312: center environment, no label, no caption
- Lines 328-338: center environment, no label, no caption
- Lines 1128-1138: center environment, no label, no caption
- Lines 1518-1528: center environment, no label, no caption
- Lines 1569-1581: center environment, no label, no caption
- Lines 1954-1967: center environment, no label, no caption
- Lines 2114-2125: center environment, no label, no caption
- Lines 2171-2183: center environment, no label, no caption

versus:

- Line 477: \label{tab:koide_scan}, with caption
- Line 561: \label{tab:scaling_pred}, with caption
- Line 694: \label{tab:master}, with caption
- Line 959: \label{tab:census}, with caption
- Line 1234: \label{tab:496}, with caption

The float tables are concentrated in Sections 6-7. All other sections use inline centered tabulars.

### 4.3 Uncertainty reporting

Uncertainties appear in three formats:
- **Plus-minus notation:** m_tau = 1776.86 +/- 0.12 MeV (line 446), M_W = 80.369 +/- 0.013 GeV (line 1378)
- **Sigma notation:** "0.39 sigma" (line 1395), "6.2 sigma" (line 1398), "0.91 sigma" (line 465), "2.8 sigma" (line 492), "3.1 sigma" (line 539)
- **Percentage deviation:** "0.10%", "1.24%", etc. (throughout)

These three formats serve different purposes (measurement uncertainty, significance level, Koide deviation) and are not genuinely inconsistent. However, the percentage deviations in the precision summary table (Section 4) should not be confused with statistical significances; they are not accompanied by error bars.

### 4.4 Negative results presentation

Negative results are handled inconsistently:
- The scaling failure (Section 6.4, "The second scaling fails") gets its own subsection with a bold summary.
- The chain stall (Section 6.5, paragraph "Step 3: the chain stalls") is a paragraph within a subsection.
- The lepton "almost-seed" (Section 3.3) notes the 9% deviation inline.
- The quark seed deviation (Section 8.7) is discussed at length in its own subsection.

This is not a serious problem, but the varying levels of emphasis on negative results may reflect different authorial voices.

### 4.5 Comment-line formatting

Section separators use `%======` lines. Most sections have the pattern: one `%======` line, then `\section{...}`, then `\label{...}`, then one `%======` line. The Open Problems section (lines 2192-2193) has two consecutive `%======` lines before it, with a stray extra separator from a deleted section heading. This is cosmetic only.

---

## 5. Missing Content in Open Problems and Conclusions

### 5.1 The Casimir equation / electroweak ratio is absent from both

The Casimir eigenvalue equation (Section 7, Part B) produces the paper's most precise numerical prediction: sin^2(theta_W) = 0.2231014... matching the on-shell measurement to 0.02%. It also predicts M_W = 80.374 GeV. Yet:

- The **Open Problems** section (Section 11) does not mention the Casimir equation's open issues (scheme dependence, the mysterious fourth eigenvalue at 122.4 GeV, the unknown origin of the equation), despite these being discussed at length in Section 7.4 (subsec:casimir_open).
- The **Conclusions** section (Section 12) does not mention the electroweak ratio prediction at all.

This is a significant omission. The Casimir result is one of the paper's strongest numerical results, and its open issues (especially the scheme dependence) are among the most important unresolved problems.

### 5.2 The SO(32) embedding is absent from Conclusions

The SO(32) decomposition (Section 7, Part A) provides the UV origin for the symmetric 15 representation. The Conclusions mention the 15 (line 2263) and its SO(32) origin (line 2265: "by the SO(32) decomposition"), but do not discuss the 270 extra states that must be projected out -- the projection problem (line 1298-1301) is an important open question that is missing from Open Problems.

### 5.3 The composite scalar counting / generation uniqueness is underrepresented

Section 6 (Counting) proves that N=3 generations is the unique solution to the sBootstrap Diophantine system, and that right-handed neutrinos are required. This is a strong structural result, but:
- The Open Problems section does not discuss follow-up questions (e.g., whether the counting extends to include the top quark, or what happens with supersymmetric extensions).
- The Conclusions section does not list the N=3 uniqueness as a principal result, despite it being arguably as important as the Koide energy balance.

### 5.4 The iterative chain / scaling rule results are absent from Open Problems

The mass-chain section (Section 6) presents both successes (lepton-to-quark scaling, iterative chain capturing c,s,b) and failures (scaling stops at one step, chain stalls at 0.035 MeV). The Open Problems section mentions "connecting the three tuples" (lines 2212-2215) but does not specifically discuss why the iterative chain fails to reach the first generation, or whether the x3 scaling has a deeper origin.

### 5.5 Missing reference to the chirality section

The "Nature Abhors Chirality" section (Section 10) proposes a unifying principle. The Open Problems section mentions the strong CP problem (which is discussed in that section's table) but does not cross-reference Section 10 or discuss whether the chirality principle can be formalized.

---

## 6. Factual / Numerical Consistency

### 6.1 The $(c,b,t)$ vacuum energy: 29,580 vs. 29,600 MeV

As noted in Section 1.5 above:
- Precision summary table (line 310): v_0^2 = 29,580 MeV
- Text (lines 387, 1799, 2052): v_0^2 ~ 29.6 GeV = 29,600 MeV
- Line 1813: (v_0^{(cbt)})^2 ~ 29,600 MeV

The 20 MeV discrepancy (0.07%) between 29,580 and 29,600 is a rounding issue. The table value should be taken as more precise. The text should either say "29.6 GeV" (if rounding) or "29.58 GeV" (if matching the table).

### 6.2 The R value and the W-mass prediction are self-consistent

R = 0.2231014... (line 1366). The W-mass prediction uses sqrt(1 - 0.22310) (line 1392), which rounds R to 5 significant figures. The result M_W = 80.374 GeV is consistent with using M_Z = 91.1876 GeV and R = 0.22310.

### 6.3 The lepton Koide Q values between conventions

In the precision summary table (Section 4, line 304): the lepton seed has Q = 0.6667 (in the Q = 2/3 convention).
In the mass chains section (Section 6, line 449): Q = 1.500014 (in the Q = 3/2 convention).

Converting: 1.500014 in the reciprocal convention gives 1/1.500014 = 0.66666... which matches 0.66666 (line 307) to the stated precision. So the numerical values are consistent between conventions, even though the stated algebraic relationship (Q_here = 3 Q_prev) is wrong.

**Verification:** 3 * 0.6667 = 2.0001, not 1.5. The reciprocal: 1/0.6667 = 1.4999, which matches 1.500014 within rounding. This confirms the relationship is reciprocal, not multiplicative by 3.

### 6.4 The quark seed vacuum energy

Line 196: v_0^2 = (z_s + z_c)^2 / 9 ~ 228 MeV.
Table at line 243: v_0^2 = 228 MeV for (0, s, c). Consistent.
Conclusions (line 2247): "v_0^2 ~ 230--350 MeV" -- the lower bound of 230 is a slight rounding up from 228. This is acceptable.

### 6.5 The meson seed vacuum energy

Abstract (line 41): v_0^2 ~ 351 MeV.
Table at line 242: v_0^2 = 351 MeV. Consistent.

### 6.6 The charm mass from the two chains

Table tab:master (lines 700-701):
- Chain A (from t,b): m_c = 1357 MeV
- Chain B (from leptons): m_c = 1360 MeV
- PDG: 1270 MeV

Text (line 691): "1357 vs. 1360 MeV for charm." Consistent with the table.

The iterative chain section (line 633): m_? = 1357 MeV. Consistent with Chain A.
The scaling rule section (line 566): s = 92.3, c = 1360, b = 4197. Consistent with Chain B.

### 6.7 The strange mass from the two chains

Table tab:master (line 700): Chain A = 92.2, Chain B = 92.3, PDG = 93.4.
Line 640: m_? = 92.17 MeV (chain value). 92.17 rounds to 92.2. Consistent.
Line 566: s predicted = 92.3 MeV. Consistent with Chain B.

### 6.8 The top mass

Line 179: m_t = 172.76 GeV (footnote, cited as pole mass).
Table tab:master (line 703): PDG = 172,760 MeV = 172.76 GeV. Consistent.

### 6.9 The pion-muon relation

Line 69: m_pi^2 - m_mu^2 = 8,316 MeV^2 versus f_pi^2 = 8,501 MeV^2, a 2.2% discrepancy.
Line 2029: "holding to 2.2%." Consistent.
Line 1653: m_pi^2 - m_mu^2 = f_pi^2 (stated as exact equality in the VA section). This is the idealized SUSY relation; the 2.2% discrepancy is the breaking correction. The section correctly qualifies this elsewhere.

### 6.10 The Koide phase

Line 512: delta_0 = 2.3166 rad.
Line 523: delta_0 mod (2*pi/3) = 0.22223 rad.
Check: 2*pi/3 = 2.0944. 2.3166 mod 2.0944 = 2.3166 - 2.0944 = 0.2222 rad. Consistent.
Line 527: residual = 0.22223 - 2/9 = 0.22223 - 0.22222... = 7.8 * 10^{-6}. Consistent (2/9 = 0.222222..., so 0.22223 - 0.22222 = 0.00001 ~ 10^{-5}; the stated 7.8 * 10^{-6} is a more precise calculation).

---

## 7. Summary of Findings

### Critical issues (should be fixed before submission):

1. **Erroneous footnote on Q convention** (lines 454-463): The stated relationship Q_here = 3 Q_prev is mathematically incorrect. The correct relationship is Q_here = 1/Q_prev (reciprocal, not multiplication by 3).

2. **Two uncited bibliography entries:** GOR1968 and FramptonGlashow1987 are defined but never \cited. Either add citations or remove the entries.

3. **The Casimir section results are missing from both Open Problems and Conclusions.** The W-mass prediction and the open issues (scheme dependence, fourth eigenvalue) should be mentioned in at least one of these summary sections.

### Moderate issues (recommended fixes):

4. **Overloaded symbol Q** used for the Koide ratio, electric charge, and the SUSY generator. The electric-charge usage in the charge census (Q = +2/3, Q = -2/3) is easily confused with the Koide ratio.

5. **Missing from Conclusions:** The generation-uniqueness theorem (N=3 is forced by the Diophantine counting) is not listed among the principal results.

6. **The SO(32) projection problem** (270 extra states) is not listed in Open Problems despite being flagged as an open question in Section 7.3.

7. **Minor numerical rounding discrepancy:** The $(c,b,t)$ vacuum energy appears as both 29,580 MeV (table) and 29,600 MeV (text).

### Minor issues (polish):

8. Table formatting is inconsistent: Sections 2-5 and 8-12 use inline centered tabulars without captions or labels; Sections 6-7 use proper float tables with captions and labels.

9. Triple duplicate mention of the v_0^2 ~ 29.6 GeV value for the (c,b,t) triple in running text (not counting the table), where two of the three are in the same section.

10. The vacuum energy range in the Conclusions (line 2247: "230--350 MeV") slightly rounds up the quark seed value from 228 MeV.

11. Stray duplicate comment-line separator before the Open Problems section (line 2192-2193).
