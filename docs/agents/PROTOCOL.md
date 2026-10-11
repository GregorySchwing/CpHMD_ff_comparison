# How work is claimed

The queue is GitHub Issues on the repository. The board is the file STATE.md on the branch `board`, rendered from the issues, the branches, the claim files and the pull requests by the board workflow (`python3 -m sdlc_protocol board`, whose rules are sdlc_protocol.board), which rewrites that branch's one commit on every issue and pull request event, every branch made or deleted, every push to the default branch, and every six hours as a fallback; nobody edits it by hand. Read it with `git fetch origin board && git show origin/board:STATE.md`; a repository that follows this protocol also prints it from its session-start hook, so every session begins with it. A claim is the first push of a branch: GitHub refuses the second push of the same name, so two accounts can never hold the same issue. Everything below is done by the account doing the work, in the open, and the claim skill that ships with the template and with the sdlc plugin is the same steps as commands.

## Labels

- todo: the queue. An issue labelled todo and shown unclaimed on the board is free to take.
- claimed: an account has pushed its branch. The account's own label goes on with it.
- in-review: the pull request is ready; the owner merges it and the merge closes the issue.
- blocker: the issue stops everyone until it is done.
- account:main: the main account, the first one, which runs the worker routine and, once a second account joins, the coordinator. One label per account, account:<name>, where the name is the one on the claim file's Account line.
- epic, blocked and agent:<name>: team mode's (below).

sdlc_protocol.labels spells these names and the two prefixes once, with the colour and description of every label a repository makes up front; an agent:<name> label is made when team mode first needs it.

## Pick

Fetch, then read the board (`git fetch origin board && git show origin/board:STATE.md`). Take an open issue labelled todo that the Queue section shows as unclaimed and that no open blocker stands in the way of. An issue with a branch `claude/<issue>` on origin, or a file `docs/agents/claims/<issue>.md` on the default branch, is another account's, whatever its labels say.

## Claim

Push a new branch `claude/<issue number>` from the default branch whose first commit adds `docs/agents/claims/<issue number>.md`: the account, the plan, the files to touch, the issue's link. Nothing else is in that first commit. If the push is refused because the branch exists, another account has it: delete the local branch and pick another issue.

The claim file takes this shape, with N the issue's number. The board reads two parts of it (sdlc_protocol.claims): the `- Account:` line, whose name is the account's label without the `account:` prefix, and the list under `## Files to touch`, which fills the board's Files column until the pull request has files of its own.

```
# Claim: #N

- Issue: https://github.com/OWNER/NAME/issues/N
- Account: main
- Branch: claude/N

## Plan

Two to five sentences: what will change and why, in the order it will land.

## Files to touch

- the paths this work expects to change, one per line
```

## Announce

Add the label claimed and the account's label to the issue, and post the plan as a comment, so the board and the people reading the issue see the same thing the branch says.

## Work in the open

Open a draft pull request at once, from the branch to the default branch, with "Closes #N" in its body. Push often: a branch with no commit for 24 hours is released (below). When the board comments "main moved" on the pull request, the default branch has moved past the branch: rebase on it and push again with `--force-with-lease`. That is the one forced push the protocol allows, and only on your own claim branch.

## Finish

The pull request deletes its claim file, says "Closes #N", passes the repository's checks and is marked ready; the issue gets the label in-review. The owner merges; the merge closes the issue, and the next board no longer lists the claim. Delete the branch after the merge if the repository does not do it by itself. An issue closed without a merge loses its branch too: once no open pull request uses `claude/<issue>`, the repository's tidy workflow deletes it.

## Blockers

An issue labelled blocker stops everyone: work it stands in the way of waits until it is done. It is claimed like any other issue, by one account, and nobody else fixes it in parallel; a second account that sees it waits, or helps on the one pull request by review.

## Release

A claim with no commit for 24 hours is released by the coordinator: it comments on the issue, removes the label claimed, and leaves a checkpoint comment saying what the branch holds and what seems left, so the next account can pick the issue up or start over. The branch stays as it is; the next account reads it and starts its own from the default branch. sdlc_protocol.claims.release_due says which claims are due.

## Team mode

When a repository turns team mode on (the template does it with the repository variables CLAUDE_TEAM and CLAUDE_WORKER), an issue labelled epic is planned by a lead into sub-issues pre-assigned to members (label agent:<name>), whom this protocol binds like any account, with two differences: members bring the default branch in by merging it, not by rebasing, so they never force-push, and they land by auto-merge on a green check rather than by the owner's hand. The board's queue leaves out an epic and its members' tasks, so no other account takes them. The clerical part, writing the plan into the epic, opening the tasks and telling each member what to do next, is code, not Claude: sdlc_protocol.team, run as `python3 -m sdlc_protocol team`.

## Policy and worlds

The rules above are the same wherever they run. Two things change from one place to another, and sdlc_protocol.gate names them: a Policy, which says how the protocol behaves, and a World, which is where the branches live and which an adapter implements. The gate's decisions (claim, open a pull request, mark it ready, report its status, draw the board) are functions of the policy, the world and a queue of issues and pull requests, and they run no git and no HTTP themselves.

The policy's fields, with the values this document describes as the defaults:

- claim_prefix, `claude/`: what a claim branch's name starts with, before the issue's number.
- required_check, `check`: the check that must be green on the branch's head before the gate merges; empty where merging cleanly is the whole check.
- bring_main_in, `rebase`: how a branch takes in a default branch that has moved. `rebase` is an account's way, above, with the one `--force-with-lease`; `merge` is a team member's, merging the default branch in so nobody force-pushes. The "main moved" answer says the steps in the policy's own way (sdlc_protocol.gate.main_moved_instruction).
- who_merges, `owner`: who moves the default branch. `owner` is the repository's owner by hand, as above; `gate` is the gate itself, when a pull request is marked ready; `auto-merge` is GitHub, once the required check is green, as team mode does.
- idle_after, 24 hours: how long a claim may go without a commit before it is released.

Two presets name the ways of working above: ACCOUNT_POLICY (rebase, the owner merges) and TEAM_POLICY (merge, auto-merge). A policy also carries remote and default_branch, `origin` and `main` unless the place calls them otherwise, which are only the names its answers use, so an instruction reads `origin/main` or `team/main` as the place does.

A world answers six questions: the default branch's commit (main_sha), a branch's commit (branch_sha), how a branch stands to the default branch (standing: clean, behind, conflicting or unknown), the state of the required check on its head (check: green, red, pending or none), the merge itself (merge, which returns the default branch's new commit or raises why not: MainMoved, the default branch left as it was; MergeRefused, for a reason the branch's author fixes; GateFailed, the world could not act), and the time (now). What each adapter supplies:

- sdlc_protocol.adapters.github, a repository on GitHub through the REST API, with a workflow's GITHUB_TOKEN or an App installation token: the standing from the pull request's mergeable_state, told apart further by mergeable and by how far the head is behind its base where GitHub says only that something blocks the merge, and unknown for a branch with no open pull request; the check from the check runs on the branch's head that carry the required check's name; the merge through GitHub's merge API, pinned to the head the gate looked at, and only when who_merges is `gate`, any other policy being refused with NotAllowed; the clock as given, or the wall clock.
- sdlc_protocol.adapters.bare_git, a bare repository that agents share, as coordbench's sdlc arm does: the standing from git itself, in a scratch clone, with the conflicting files named and two more words a git world can tell apart (no branch, nothing new); no checks, so the check is always none; the merge as `git merge --no-ff` under the gate's own committer identity, pushed without force, a conflict or a refused push being MainMoved; the wall clock.
- sdlc_protocol.adapters.memory, a world held in memory for tests and simulations: standings and checks set by hand, a clock that moves only when told, so a simulation can hold a check pending for as long as it means CI to take, and a record of every merge.

A refusal carries the code its command exits with, so a script can tell one kind from another: 1 refused (claim first, push, commit something), 2 someone else's (another holds the claim, the issue or the pull request), 3 main moved (bring the default branch in and try again), 4 gate failed (the world could not act).

## Rules that hold everywhere

- One issue, one branch `claude/<issue>`, one pull request.
- Never push to another account's branch. Never force-push the default branch or any branch that is not your own claim. The push guard (sdlc_protocol.guard, run as a Claude Code hook through `python3 -m sdlc_protocol guard`) refuses those pushes before they run, and branch protection on GitHub stands behind it.
- Memory rides in the same pull request as the change it came from.
- The repository's checks pass before a pull request is marked ready.
- CI holds a `claude/<issue>` pull request to this protocol (sdlc_protocol.protocol_check, run as `python3 -m sdlc_protocol check`, part of the required check): the claim file is there while the pull request is a draft, gone once it is ready, the body closes its issue, and the branch lives in the repository itself, never in a fork.
