#!/bin/bash
# Sanitizer: catch internal project labels leaking into paper .tex files
# Run after any agent writes to a .tex file.

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'
FAIL=0

for f in *.tex; do
  [ -f "$f" ] || continue

  # 1. Internal gap/agent labels (G12, M13, B4, R2-A, R3-C etc.)
  if grep -nE '\bG[0-9]{1,2}\b|\bM[0-9]{1,2}\b|\bB[0-9]{1,2}\b|\bR[0-9]-[A-F]\b' "$f" | grep -v '^[0-9]*:%'; then
    echo -e "${RED}FAIL${NC}: $f contains internal project labels (G/M/B/R-agent)"
    FAIL=1
  fi

  # 2. Agent infrastructure language
  if grep -niE 'Round [0-9]|coordinator|subagent|Manhattan|compartmental' "$f" | grep -v '^[0-9]*:%'; then
    echo -e "${RED}FAIL${NC}: $f contains agent infrastructure language"
    FAIL=1
  fi

  # 3. Orphan labels (defined but never referenced)
  LABELS=$(grep -oE '\\label\{[^}]+\}' "$f" | sed 's/\\label{//;s/}//' | sort)
  REFS=$(grep -oE '\\(eq)?ref\{[^}]+\}' "$f" | sed 's/\\ref{//;s/\\eqref{//;s/}//' | sort -u)
  ORPHANS=$(comm -23 <(echo "$LABELS") <(echo "$REFS"))
  N_ORPHAN=$(echo "$ORPHANS" | grep -c '[a-z]')
  N_TOTAL=$(echo "$LABELS" | grep -c '[a-z]')
  if [ "$N_ORPHAN" -gt 0 ]; then
    echo -e "${RED}INFO${NC}: $f has $N_ORPHAN/$N_TOTAL unreferenced labels"
  fi
done

if [ "$FAIL" -eq 0 ]; then
  echo -e "${GREEN}PASS${NC}: No internal labels found in .tex files"
fi
