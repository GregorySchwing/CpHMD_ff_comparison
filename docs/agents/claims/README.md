# Claim files

A claim file is the first commit of a claim: `docs/agents/claims/<issue number>.md`, added by the first commit on the branch `claude/<issue number>`, which is pushed from the default branch as docs/agents/PROTOCOL.md says. It names the account, the plan and the files the work expects to touch, so the board can show who holds what and what might collide, and so a released claim leaves its plan behind for the next account. The account's name is its label without the `account:` prefix (`main` for the label account:main); the board reads it from the `- Account:` line, and takes the list under "## Files to touch" for its Files column until the pull request has files of its own. The pull request that finishes the issue deletes the file. Files here are named by issue number only; this README is the one exception, and nothing reads it as a claim.

The template, with N the issue's number:

```
# Claim: #N

- Issue: https://github.com/GregorySchwing/CpHMD_ff_comparison/issues/N
- Account: main
- Branch: claude/N

## Plan

Two to five sentences: what will change and why, in the order it will land.

## Files to touch

- the paths this work expects to change, one per line
```

Keep it short; the pull request carries the detail. A claim is not a design document.
