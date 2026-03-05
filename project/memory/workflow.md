# Workflow Design

## Lessons from Failed Attempts

### physbook (Feb 2026)
- Too general (all of physics), too much infrastructure (25 mandatory preloads, 6 agent roles)
- Never got to actual content

### physarticle (Feb 2026)
- "From Newton to the Path Integral" — deliberately chose SAFE topic to avoid immune response
- Elaborate infrastructure: kanban, orchestrator hierarchy, patch requests, referee loops, blackboards
- Generated 187k tokens. Hostile referee found core proposition was circular.
- Proved: LLMs can't do physics even on safe ground. Infrastructure becomes the project.
- Later success on claude.ai when user provided the insight (Newton knew interference → Opticks + Principia polygon proof = path integral)

### Key insight
**Two barriers, not one:**
1. LLMs can't do physics (basic capability gap) — physarticle proved this
2. Heterodox framing triggers retreat (epistemic immune response)

Manhattan compartmentalization solves barrier 2.
Making agents do ARITHMETIC not PHYSICS solves barrier 1.
The agent prompts ask for Python scripts, Diophantine solutions, table lookups — not physics.

## Current Workflow (phys3)

Dead simple. No kanban, no blackboards, no infrastructure bloat.

```
1. Launch 8 agents (each gets one prompt from agent_prompts.md)
2. Each agent writes output to results/agentN.md
3. Staged assembly (pairs → sections → paper)
4. User writes abstract + final framing
```

## Agent Output Directory
Results go to `/Users/arivero/phys3/results/`

## Assembly Pipeline (future)
- Agent 9: outputs 1+2 → results section
- Agent 10: outputs 4+6 → methods section
- Agent 11: outputs 5+7 → group theory section
- Agent 12: three formatted sections → paper draft
- User: abstract + framing

Paper framed as follow-up to EPJC 84, 1058 (2024) = arXiv:2407.05397
