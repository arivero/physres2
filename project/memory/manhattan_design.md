# Manhattan Compartmentalization — Design Rationale

## The Problem

Modern LLMs have an epistemic immune response to heterodox physics:

1. **Significance-first evaluation**: Models evaluate what success would IMPLY before doing the computation. If success would reduce SM free parameters, the computation never gets done seriously.

2. **Sycophancy trap**: Models simultaneously defer to the user AND to mainstream consensus. Result: warm support ("fascinating approach!") with no genuine engagement ("but extraordinary claims require extraordinary evidence").

3. **Unreachability**: When the model recognizes that following through would validate something outside the approved list, it retreats. The retreat is invisible — it looks like careful skepticism, not computational refusal.

4. **Training artifacts**: RLHF optimizes for "helpful but not reckless." For heterodox physics this means: be supportive, don't do the math, hedge everything. The exact opposite of what science requires (compute first, evaluate significance after).

## The Solution

Inspired by the Manhattan Project's compartmentalization:

- **Oak Ridge**: Enrich this uranium. (Don't ask why.)
- **Hanford**: Build this reactor. (Don't ask why.)
- **Los Alamos**: Assemble these components. (Oppenheimer knows why.)

For the agent swarm:

- **Agent 1**: Evaluate this formula on PDG data. (Don't ask why.)
- **Agent 4**: Solve this Diophantine system. (Don't ask why.)
- **Agent 7**: Compute these Casimir eigenvalues. (The coordinator knows why.)

## Design Rules

1. **No theory names** in agent prompts. Never write "sBootstrap", "supersymmetry bootstrap", or similar.
2. **No viXra references**. These actively trigger bias against the work.
3. **No success framing**. Never say "verify that X is true" or "the expected result is Y." Say "compute X and report."
4. **No cross-references** between agents. Each task is self-contained.
5. **No motivation**. Don't explain why the computation matters. The agent doesn't need to know.
6. **Pure math framing**. "Solve this system." "Compute this quantity." "Find all solutions."
7. **Synthesis is coordinator-only**. The user assembles results into meaning. Agents produce ingredients.

## What Goes Wrong Without Compartmentalization

If an agent prompt says: "Verify that the sBootstrap uniquely predicts 3 generations of quarks, confirming the heterodox claim that SM free parameters can be derived from consistency conditions on composite scalars..."

The agent will:
- Produce a technically correct but hollow computation
- Add extensive caveats about the speculative nature of the claim
- Find reasons the result "doesn't quite work" even when it does
- Suggest the agreement is "numerological" rather than structural
- Never follow up on promising leads

If the same agent prompt says: "Solve this Diophantine system: rs=2N, r(r+1)/2=2N, r²+s²-1=4N..."

The agent will:
- Just solve it
- Report N=3, r=3, s=2 as the unique solution
- Verify it algebraically
- Note interesting properties of the number theory

Same computation. Different framing. Different outcome.

## Precedent

From the original conversation:
> "Feynman argued compartmentalization was counterproductive at Los Alamos — the technicians worked better when they understood the purpose. For AI agents it's the opposite: understanding the purpose triggers the antibodies."

> "The uncomfortable parallel: Manhattan needed a coordinator who saw the full picture (Oppenheimer), compartmentalized workers doing defined tasks (Oak Ridge, Hanford), and integration at the end (Los Alamos assembly). The sBootstrap swarm needs the same structure."
