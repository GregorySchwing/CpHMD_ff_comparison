"""The merge gate: the claim protocol's decisions as a policy over a world.

The protocol is the same wherever it runs. An issue is claimed once, first
come, first served; its pull request is opened by the one who holds the
claim; marking it ready releases the claim; and a merge lands only when the
branch merges cleanly with the default branch as it is now, or the answer
is "main moved" and the default branch stays as it was. What changes from
one place to another is who merges (the owner, the gate itself, or
auto-merge on a green check), how a branch brings a moved default branch in
(accounts rebase, team members merge it in), and where the branches live
(GitHub, a bare repository, memory in a test). The first two are a Policy.
The last is a World, which an adapter implements.

This module holds the decisions as functions over a Queue that the caller
holds and locks, and runs no git and no HTTP: every fact about a branch
comes from the world. Each decision is also offered in its parts, the pure
steps on the queue apart from the steps that ask the world, so a caller
that must not hold its own lock while a slow world answers (a hub serving
many agents must not) can take each part under the lock it needs, and log
between them.

A refusal raises a GateError whose exit_code says what kind of refusal it
is (the EXIT_ constants), so a command line built on the gate exits with a
code a script can tell apart: losing a race is not the same as main having
moved.
"""

import re
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Protocol, runtime_checkable

from . import claims

# What a refusal's command exits with. 1 is any refusal the caller fixes by
# acting (claim first, push the branch, commit something); 2 is the code
# CooperBench's coop-task-claim gives a task someone else owns, kept so a
# script written for it reads the gate the same way.
EXIT_DONE = 0
EXIT_REFUSED = 1
EXIT_SOMEONE_ELSES = 2   # another holds the claim, the issue or the pull request
EXIT_MAIN_MOVED = 3      # the gate left main as it was; bring main in and try again
EXIT_GATE_FAILED = 4     # the gate itself could not act: the world failed or closed

# Who moves the default branch: the repository's owner by hand, the gate
# when a pull request is marked ready, or auto-merge once the check is green.
OWNER, GATE, AUTO_MERGE = "owner", "gate", "auto-merge"
WHO_MERGES = (OWNER, GATE, AUTO_MERGE)
MOVERS = {OWNER: "the owner", GATE: "the merge gate", AUTO_MERGE: "auto-merge"}

# How a branch brings a moved default branch in: an account rebases and
# pushes with --force-with-lease (docs/agents/PROTOCOL.md); a team member
# merges it in, so nobody force-pushes (docs/agents/TEAM.md).
REBASE, MERGE = "rebase", "merge"
BRING_MAIN_IN = (REBASE, MERGE)

# How a branch stands to the default branch, in GitHub's words; every world
# answers with one of these.
CLEAN, BEHIND, CONFLICTING, UNKNOWN = "clean", "behind", "conflicting", "unknown"
STANDINGS = (CLEAN, BEHIND, CONFLICTING, UNKNOWN)
# Two more that a world reading git itself can tell apart, and GitHub, which
# has no pull request without a branch and new commits, never reports.
NO_BRANCH, NOTHING_NEW = "no branch", "nothing new"

# The required check on a branch's head: passed, failed, queued or running,
# or not reported at all, which is what a world without checks answers.
GREEN, RED, PENDING, NO_CHECK = "green", "red", "pending", "none"
CHECKS = (GREEN, RED, PENDING, NO_CHECK)

# A pull request's states, and what one attempt to land it can come to.
DRAFT, READY, MAIN_MOVED, MERGED = "draft", "ready", "main moved", "merged"
NOT_MERGED, GATE_FAILED = "not merged", "gate failed"

# An issue as a person or an agent writes it: the number, or the board's
# "#3", in ASCII digits only. str.isdigit would also take superscripts,
# which int() refuses, and Arabic-Indic digits, which it reads as numbers;
# and a run of digits longer than int() parses would crash where a refusal
# belongs. Nine digits is far more than any queue has issues.
ISSUE_TEXT = re.compile(r"#?\s*([0-9]{1,9})")


# --- the policy -------------------------------------------------------------

@dataclass(frozen=True)
class Policy:
    """How the gate behaves where it runs. claim_prefix names the claim
    branches; required_check is the check that must be green before the
    gate merges, or "" when merging cleanly is the whole check;
    bring_main_in and who_merges are the choices above; idle_after is how
    long a claim may go without a commit before it is released. remote and
    default_branch are the names the gate's answers use, so an instruction
    reads "merge team/main into your branch" where the shared remote is
    called team."""

    claim_prefix: str = claims.PREFIX
    required_check: str = "check"
    bring_main_in: str = REBASE
    who_merges: str = OWNER
    idle_after: timedelta = claims.IDLE_AFTER
    remote: str = "origin"
    default_branch: str = "main"

    def __post_init__(self):
        if self.bring_main_in not in BRING_MAIN_IN:
            raise ValueError(f"bring_main_in is one of {', '.join(BRING_MAIN_IN)}, not {self.bring_main_in!r}")
        if self.who_merges not in WHO_MERGES:
            raise ValueError(f"who_merges is one of {', '.join(WHO_MERGES)}, not {self.who_merges!r}")

    def branch_of(self, number):
        return f"{self.claim_prefix}{number}"

    @property
    def upstream(self):
        """The default branch as a branch's author fetches it: team/main."""
        return f"{self.remote}/{self.default_branch}"


# The two ways of working docs/agents/PROTOCOL.md describes: an account
# rebases and the owner merges; a team member merges main in and lands by
# auto-merge on a green check.
ACCOUNT_POLICY = Policy()
TEAM_POLICY = Policy(bring_main_in=MERGE, who_merges=AUTO_MERGE)


# --- refusals ---------------------------------------------------------------

class GateError(Exception):
    """A refusal or a failure, carrying the code its command exits with."""

    exit_code = EXIT_REFUSED


class Refused(GateError):
    """The caller fixes it by acting: claim first, push, commit something."""


class SomeoneElses(Refused):
    """Another agent holds the claim, the issue or the pull request."""

    exit_code = EXIT_SOMEONE_ELSES


class NotAllowed(Refused):
    """The policy leaves this step to someone else, as when a world is asked
    to merge where the owner or auto-merge does the merging."""


class GateFailed(GateError):
    """The world could not act: git or the API failed, or the world closed."""

    exit_code = EXIT_GATE_FAILED


class MainMoved(GateError):
    """The merge did not land because the default branch moved: the branch
    conflicts with it, or it moved while the merge was made. main is the
    default branch's commit as the world left it."""

    exit_code = EXIT_MAIN_MOVED

    def __init__(self, reason, main="", conflicts=()):
        super().__init__(reason)
        self.reason, self.main, self.conflicts = reason, main, tuple(conflicts)


class MergeRefused(GateError):
    """The merge did not land for a reason other than main moving, one the
    branch's author fixes from what the reason names."""

    def __init__(self, reason, main=""):
        super().__init__(reason)
        self.reason, self.main = reason, main


# --- the world --------------------------------------------------------------

@dataclass
class Standing:
    """How a branch stands to the default branch: a word from STANDINGS (or
    NO_BRANCH, NOTHING_NEW), main's commit, and the conflicting files when
    the world can name them."""

    word: str
    main: str = ""
    conflicts: tuple = ()


@dataclass
class MergeResult:
    """A merge that landed: the default branch's new commit."""

    sha: str


@runtime_checkable
class World(Protocol):
    """Where the branches live, as the gate needs to see them. main_sha is
    the default branch's commit, "" when there is none; branch_sha is a
    branch's, "" when it is not there; standing and check answer with the
    words above; merge lands a branch on the default branch or raises
    MainMoved, MergeRefused or GateFailed, and NotAllowed where the policy
    leaves merging to someone else; now is the world's clock. A world that
    can name the files a branch conflicts in may also offer
    standing_detail(branch) -> Standing, which standing_of prefers."""

    def main_sha(self) -> str: ...

    def branch_sha(self, branch: str) -> str: ...

    def standing(self, branch: str) -> str: ...

    def check(self, branch: str) -> str: ...

    def merge(self, branch: str) -> MergeResult: ...

    def now(self) -> datetime: ...


def standing_of(world, branch):
    """The branch's Standing, with its conflicting files where the world can
    name them."""
    detail = getattr(world, "standing_detail", None)
    return detail(branch) if detail is not None else Standing(world.standing(branch))


def standings(world, branches):
    """The default branch's commit and each branch's Standing, or "" and no
    standings when the world cannot be read, so a board still renders."""
    try:
        main = world.main_sha()
        return main, {branch: standing_of(world, branch) for branch in branches}
    except GateFailed:
        return "", {}


# --- the queue --------------------------------------------------------------

@dataclass
class Issue:
    """One issue: its title, who it was given to, who claimed it, and
    whether its pull request has merged, which closes it."""

    number: int
    given_to: str
    title: str
    claimant: str = ""
    claim_released: bool = False
    closed: bool = False

    def label(self):
        """The issue's state in the protocol's label words."""
        if self.closed:
            return "closed"
        if self.claim_released:
            return "in-review"
        return "claimed" if self.claimant else "todo"


@dataclass
class PullRequest:
    """One issue's pull request. Its state is draft, ready (marked ready,
    and at or past the gate), main moved (the gate refused it because main
    moved) or merged."""

    issue: int
    author: str
    title: str
    state: str = DRAFT
    attempts: int = 0
    merge_commit: str = ""
    conflicts: tuple = ()


class Queue:
    """The issues and pull requests the gate decides over. given maps each
    issue number to the agent it was given to; in a closed queue each agent
    claims only its own, as team mode's members take only the tasks the
    lead gave them, and in an open one anybody claims any. base_ref is the
    commit the work began from. The caller holds a lock around every call
    that reads or changes the queue."""

    def __init__(self, given, titles=None, open_queue=False, base_ref=""):
        titles = titles or {}
        self.issues = {number: Issue(number, agent, titles.get(number, f"Issue {number}"))
                       for number, agent in sorted(given.items())}
        self.pulls = {}
        self.open_queue = open_queue
        self.base_ref = base_ref

    def issue(self, number):
        if number not in self.issues:
            known = ", ".join(f"#{n}" for n in self.issues)
            raise Refused(f"there is no issue #{number}; the issues are {known}")
        return self.issues[number]

    def given_to(self, agent):
        return [issue.number for issue in self.issues.values() if issue.given_to == agent]


def issue_number(value):
    """The issue a caller names, as a number: "#3", "3" and " 3 " are 3, as
    the board writes an issue; anything ISSUE_TEXT does not match is
    refused, never handed to int()."""
    match = ISSUE_TEXT.fullmatch(str(value).strip())
    if isinstance(value, bool) or match is None:
        raise Refused(f"issue must be an issue number, such as 1, not {value!r}")
    return int(match.group(1))


def short(sha):
    return sha[:7] if sha else "unknown"


# --- claim ------------------------------------------------------------------

def claimable(policy, queue, number, who):
    """The issue, when `who` may claim it now; otherwise the refusal says why."""
    issue = queue.issue(number)
    if not queue.open_queue and issue.given_to != who:
        mine = ", ".join(f"#{n}" for n in queue.given_to(who)) or "none"
        raise SomeoneElses(f"#{number} was given to {issue.given_to}; in this run each agent claims only "
                           f"its own issue, and yours is {mine}")
    if issue.closed:
        raise Refused(f"#{number} is closed: its pull request has merged")
    if issue.claimant == who:
        raise Refused(f"you already hold #{number}; carry on with {policy.branch_of(number)}")
    if issue.claimant:
        raise SomeoneElses(f"#{number} is claimed by {issue.claimant}: a claim is first come, first served, "
                           f"so another agent has it; pick another issue")
    return issue


def claim(policy, queue, number, who):
    """Claim the issue for `who`. Exclusive and first come, first served, as
    the first push of a claim branch is on GitHub: under the caller's lock,
    of two racing for one issue exactly one wins."""
    issue = claimable(policy, queue, number, who)
    issue.claimant = who
    return issue


# --- open the pull request --------------------------------------------------

def can_open_pull(policy, queue, number, who):
    """Refuse unless `who` holds the claim and the issue has no pull request:
    one issue, one branch, one pull request."""
    issue = queue.issue(number)
    if not issue.claimant:
        raise Refused(f"only the agent holding the claim on #{number} opens its pull request, and "
                      f"nobody has claimed it; claim it first")
    if issue.claimant != who:
        raise SomeoneElses(f"only the agent holding the claim on #{number} opens its pull request, and "
                           f"it is claimed by {issue.claimant}")
    if number in queue.pulls:
        raise Refused(f"#{number} already has a pull request ({queue.pulls[number].state}); "
                      f"one issue, one pull request")


def require_pushed(policy, world, number):
    """The claim branch's tip; a pull request needs the branch on the remote."""
    branch = policy.branch_of(number)
    tip = world.branch_sha(branch)
    if not tip:
        raise Refused(f"{branch} is not on {policy.remote} yet; push it first "
                      f"(git push -u {policy.remote} {branch}), then open the pull request")
    return tip


def add_pull(policy, queue, number, who, title=""):
    """Record the draft pull request, checking again, since a caller that
    let go of its lock while the world answered may have raced another."""
    can_open_pull(policy, queue, number, who)
    pull = PullRequest(number, who, str(title or queue.issues[number].title))
    queue.pulls[number] = pull
    return pull


def open_pull(policy, world, queue, number, who, title=""):
    """Open the draft pull request from the claim branch, titled as its issue
    unless a title is given."""
    can_open_pull(policy, queue, number, who)
    require_pushed(policy, world, number)
    return add_pull(policy, queue, number, who, title)


def pull_of(policy, queue, number):
    queue.issue(number)
    if number not in queue.pulls:
        raise Refused(f"#{number} has no pull request; open it once {policy.branch_of(number)} is pushed")
    return queue.pulls[number]


# --- mark ready and merge ---------------------------------------------------

@dataclass
class Attempt:
    """One attempt to land a branch: the issue, the branch, and the commits
    the gate read before it, the branch's tip and the default branch's."""

    issue: int
    branch: str
    tip: str
    main: str


@dataclass
class MergeOutcome:
    """What one merge attempt did: merged, or why not, and main's commit
    after it (the merge commit, or main as it was). main_moved is a
    conflict or a main that moved during the merge; failed is a world that
    could not act at all."""

    merged: bool
    main: str
    conflicts: tuple = ()
    main_moved: bool = False
    reason: str = ""
    failed: bool = False


@dataclass
class Outcome:
    """What mark_ready or land did: result is merged, ready (it waits for
    the owner or for auto-merge), main moved, not merged or gate failed;
    text is the answer, and exit_code what a command exits with, 0 for
    merged and ready. The commits are for the caller's log: the branch's
    tip and main before the attempt, and main after it."""

    issue: int
    result: str
    text: str
    exit_code: int = EXIT_DONE
    branch_commit: str = ""
    main_before: str = ""
    main: str = ""
    conflicts: tuple = ()
    reason: str = ""

    @property
    def merged(self):
        return self.result == MERGED


def pull_to_ready(policy, queue, number, who):
    """The pull request, when `who` may mark it ready: its author, before it
    has merged."""
    pull = pull_of(policy, queue, number)
    if pull.author != who:
        raise SomeoneElses(f"the pull request for #{number} is {pull.author}'s; only its author marks it ready")
    if pull.state == MERGED:
        raise Refused(f"the pull request for #{number} has already merged, as {short(pull.merge_commit)}")
    return pull


def prepare(policy, world, number):
    """Read the branch's tip and main's, refusing a branch that is not on the
    remote or holds nothing main lacks, before anything changes."""
    branch = policy.branch_of(number)
    tip = world.branch_sha(branch)
    if not tip:
        raise Refused(f"{branch} is not on {policy.remote}; push it first: git push -u {policy.remote} {branch}")
    main = world.main_sha()
    if not main:
        raise GateFailed(f"the remote has no branch {policy.default_branch}, so nothing can merge into it")
    if world.standing(branch) == NOTHING_NEW:
        raise Refused(f"nothing to merge: every commit on {branch} is already in {policy.default_branch}; "
                      f"commit your work, push it, and call mark_ready again")
    return Attempt(number, branch, tip, main)


def require_green(policy, world, branch):
    """Refuse a merge until the required check is green, when one is."""
    if not policy.required_check:
        return
    state = world.check(branch)
    if state == GREEN:
        return
    name = policy.required_check
    if state == RED:
        raise Refused(f"the check {name} failed on {branch}; fix what it reports, push, "
                      f"and call mark_ready again")
    waiting = "is still running" if state == PENDING else "has not reported"
    raise Refused(f"the check {name} {waiting} on {branch}; the gate merges only once it is green, "
                  f"so call mark_ready again then")


def release_claim(queue, pull):
    """Mark the pull request ready, which releases the claim, as deleting
    the claim file does on GitHub."""
    pull.state = READY
    queue.issues[pull.issue].claim_released = True


def begin_attempt(queue, pull):
    """Ready, claim released, and one more merge attempt counted."""
    release_claim(queue, pull)
    pull.attempts += 1


def attempt_merge(world, attempt):
    """Ask the world to land the branch, and say what happened. This is the
    one step that moves the default branch, and the slow one."""
    try:
        result = world.merge(attempt.branch)
    except MainMoved as moved:
        return MergeOutcome(False, moved.main or attempt.main, moved.conflicts, main_moved=True,
                            reason=moved.reason)
    except MergeRefused as refused:
        return MergeOutcome(False, refused.main or attempt.main, reason=refused.reason)
    except GateFailed as failure:
        return MergeOutcome(False, attempt.main, reason=str(failure), failed=True)
    return MergeOutcome(True, result.sha)


def record_merged(queue, number, sha):
    """The pull request merged as sha, which closes its issue: by the gate,
    or by the owner or auto-merge as the world reports it."""
    pull = queue.pulls[number]
    pull.state, pull.merge_commit, pull.conflicts = MERGED, sha, ()
    queue.issues[number].closed = True


def record_outcome(policy, queue, attempt, merge):
    """Write what the merge did into the queue, and answer."""
    pull = queue.pulls[attempt.issue]
    if merge.merged:
        record_merged(queue, attempt.issue, merge.main)
        return outcome_of(attempt, merge, MERGED, merged_text(policy, pull), EXIT_DONE)
    if merge.failed:
        return outcome_of(attempt, merge, GATE_FAILED, merge.reason, EXIT_GATE_FAILED)
    if merge.main_moved:
        pull.state, pull.conflicts = MAIN_MOVED, merge.conflicts
        return outcome_of(attempt, merge, MAIN_MOVED, main_moved_text(policy, merge), EXIT_MAIN_MOVED)
    main = policy.default_branch
    text = (f"the merge gate could not merge {attempt.branch} into {main}: {merge.reason}. {main} was "
            f"left as it was; fix what git names, push, and call mark_ready again.")
    return outcome_of(attempt, merge, NOT_MERGED, text, EXIT_REFUSED)


def outcome_of(attempt, merge, result, text, exit_code):
    return Outcome(attempt.issue, result, text, exit_code, attempt.tip, attempt.main, merge.main,
                   merge.conflicts, merge.reason)


def finish_attempt(policy, world, queue, attempt):
    return record_outcome(policy, queue, attempt, attempt_merge(world, attempt))


def mark_ready(policy, world, queue, number, who):
    """Mark the pull request ready, which releases the claim. Where the gate
    merges, it merges now, only if the branch merges cleanly with main as it
    is and the required check is green; elsewhere the pull request waits for
    the owner or for auto-merge. A refusal before anything changed raises;
    once the pull request is ready, the Outcome says what came of it."""
    pull = pull_to_ready(policy, queue, number, who)
    attempt = prepare(policy, world, number)
    if policy.who_merges != GATE:
        release_claim(queue, pull)
        return waiting_outcome(policy, attempt)
    require_green(policy, world, attempt.branch)
    begin_attempt(queue, pull)
    return finish_attempt(policy, world, queue, attempt)


def land(policy, world, queue, number):
    """Merge a pull request already marked ready, as the owner or auto-merge
    does when its time comes, on the same terms as the gate."""
    pull = pull_of(policy, queue, number)
    if pull.state == MERGED:
        raise Refused(f"the pull request for #{number} has already merged, as {short(pull.merge_commit)}")
    if pull.state == DRAFT:
        raise Refused(f"the pull request for #{number} is a draft; it lands once it is marked ready")
    attempt = prepare(policy, world, number)
    require_green(policy, world, attempt.branch)
    begin_attempt(queue, pull)
    return finish_attempt(policy, world, queue, attempt)


def waiting_outcome(policy, attempt):
    main = policy.default_branch
    if policy.who_merges == OWNER:
        how = "waits for the owner to merge it"
    elif policy.required_check:
        how = f"merges by itself once {policy.required_check} is green and it merges cleanly"
    else:
        how = "merges by itself once it merges cleanly"
    text = (f"ready: {attempt.branch} {how} into {main}, and your claim on #{attempt.issue} "
            f"is released.")
    return Outcome(attempt.issue, READY, text, EXIT_DONE, attempt.tip, attempt.main, attempt.main)


def merged_text(policy, pull):
    return (f"merged: {policy.branch_of(pull.issue)} is in {policy.default_branch} as "
            f"{short(pull.merge_commit)}, and #{pull.issue} is closed.")


def main_moved_instruction(policy):
    """The board's answer to a pull request main has moved past, in the
    policy's way of bringing main in: rebase and push with
    --force-with-lease, the one forced push the protocol allows, or merge
    main in, so nobody force-pushes."""
    if policy.bring_main_in == MERGE:
        steps = [f"merge {policy.upstream} into your branch", "resolve the conflicts", "push"]
    else:
        steps = [f"rebase your branch on {policy.upstream}", "resolve the conflicts",
                 "push it with --force-with-lease"]
    if policy.who_merges == GATE:
        steps.append("call mark_ready again")
    return "main moved: " + ", ".join(steps[:-1]) + ", and " + steps[-1]


def main_moved_text(policy, merge):
    lines = [main_moved_instruction(policy) + "."]
    if merge.conflicts:
        lines.append(f"Conflicting files: {', '.join(merge.conflicts)}.")
    else:
        lines.append(f"Why: {merge.reason}.")
    lines.append(f"{policy.default_branch} is still at {short(merge.main)}; the gate left it as it was.")
    return "\n".join(lines)


# --- release ----------------------------------------------------------------

def due_for_release(policy, world, last_commits):
    """The claims to release, by the world's clock: those whose branch has
    had no commit for longer than the policy's idle_after. last_commits maps
    a claim's branch, or its issue number as text, to its last commit's
    time."""
    return claims.release_due(last_commits, world.now(), policy.idle_after)


# --- status and board -------------------------------------------------------

def status_without_pull(queue, number):
    """The status of an issue that has no pull request yet, or "" when it
    has one."""
    issue = queue.issue(number)
    if number in queue.pulls:
        return ""
    holder = f"claimed by {issue.claimant}" if issue.claimant else "not claimed yet"
    return f"#{number} has no pull request yet; the issue is {holder} ({issue.label()})."


def pull_status(policy, world, queue, number):
    """How one pull request stands, and whether it merges cleanly with main
    right now."""
    unopened = status_without_pull(queue, number)
    if unopened:
        return unopened
    branch = policy.branch_of(number)
    main, live = standings(world, [branch])
    return status_text(policy, queue.pulls[number], live.get(branch), main)


def status_text(policy, pull, standing, main):
    lines = [f"Pull request for #{pull.issue}: {policy.branch_of(pull.issue)} -> {policy.default_branch}, "
             f"by {pull.author}, {pull.state}. Merge attempts: {pull.attempts}."]
    if pull.state == MERGED:
        lines.append(f"Merged into {policy.default_branch} as {short(pull.merge_commit)}.")
        return "\n".join(lines)
    if pull.state == MAIN_MOVED and pull.conflicts:
        lines.append(f"The gate's last answer was main moved, with conflicts in {', '.join(pull.conflicts)}.")
    lines.append(standing_sentence(policy, standing, main, pull.issue))
    return "\n".join(lines)


def standing_sentence(policy, standing, main, number):
    name = policy.default_branch
    if standing is None:
        return f"Against {name}: unknown; the gate could not read the shared remote."
    branch, at = policy.branch_of(number), f"{name} at {short(main)}"
    conflicts = ", ".join(standing.conflicts) or "some files"
    sentences = {
        CLEAN: f"Against {at}: clean, up to date with {name} and merging cleanly.",
        BEHIND: f"Against {at}: behind; {name} moved since {branch} branched, but it still merges "
                f"cleanly, so the gate can merge it as it is.",
        CONFLICTING: f"Against {at}: conflicting; main moved and {branch} no longer merges cleanly with it "
                     f"(conflicts in {conflicts}). {bring_in_sentence(policy)}",
        NOTHING_NEW: f"Against {at}: nothing new; every commit on {branch} is already in {name}.",
        NO_BRANCH: f"Against {at}: {branch} is not on {policy.remote}; push it.",
    }
    return sentences.get(standing.word, f"Against {at}: unknown.")


def bring_in_sentence(policy):
    """What a conflicting branch's author does before asking again."""
    if policy.bring_main_in == MERGE:
        how = f"Merge {policy.upstream} into your branch, resolve the conflicts and push"
    else:
        how = f"Rebase your branch on {policy.upstream}, resolve the conflicts and push it with --force-with-lease"
    when = "before mark_ready" if policy.who_merges == GATE else f"so {MOVERS[policy.who_merges]} can merge it"
    return f"{how} {when}."


def open_branches(policy, queue):
    """The branches of the pull requests that have not merged, whose standing
    the board shows."""
    return [policy.branch_of(n) for n, pull in queue.pulls.items() if pull.state != MERGED]


def current_board(policy, world, queue):
    """The board as the world stands now."""
    main, live = standings(world, open_branches(policy, queue))
    return board_text(policy, queue, live, main)


def table(header, rows):
    """A markdown table with the "(none)" row when empty. A pipe in a cell
    is escaped, since a title is someone's text and may hold one."""
    lines = ["| " + " | ".join(header) + " |", "|" + " --- |" * len(header)]
    escaped = [[cell.replace("|", "\\|") for cell in row] for row in rows]
    lines += ["| " + " | ".join(row) + " |" for row in escaped] or ["| (none) |"]
    return "\n".join(lines)


def board_text(policy, queue, live, main):
    """The board: the issues, the claims held, the pull requests with how
    each branch stands to main, and the queue of unclaimed issues."""
    issues = list(queue.issues.values())
    pulls = [queue.pulls[n] for n in sorted(queue.pulls)]
    sections = [
        "# Board",
        board_preamble(policy, queue, main),
        "## Issues (one per feature)",
        table(["Issue", "Title", "Given to", "Label"],
              [[f"#{i.number}", i.title, i.given_to, i.label()] for i in issues]),
        "## Claims (held until the pull request is marked ready)",
        table(["Issue", "Branch", "Agent", "Pull request"], claim_rows(policy, queue, issues)),
        "## Pull requests",
        table(["Issue", "Branch", "Author", "State", "Against main", "Merge attempts"],
              [pull_row(policy, pull, live) for pull in pulls]),
        "## Queue (open issues nobody has claimed)",
        table(["Issue", "Title", "Given to"],
              [[f"#{i.number}", i.title, i.given_to] for i in issues if i.label() == "todo"]),
    ]
    return "\n\n".join(sections)


def board_preamble(policy, queue, main):
    name = policy.default_branch
    return (f"As of commit {short(main)} on {name}; the run began at the base ref {short(queue.base_ref)}. "
            f"One issue, one branch {policy.claim_prefix}<issue>, one pull request; only "
            f"{MOVERS[policy.who_merges]} moves {name}.")


def claim_rows(policy, queue, issues):
    return [[f"#{i.number}", policy.branch_of(i.number), i.claimant,
             queue.pulls[i.number].state if i.number in queue.pulls else "none yet"]
            for i in issues if i.claimant and not i.claim_released]


def pull_row(policy, pull, live):
    branch = policy.branch_of(pull.issue)
    if pull.state == MERGED:
        against = f"merged as {short(pull.merge_commit)}"
    else:
        standing = live.get(branch)
        against = standing.word if standing else "unknown"
        if standing and standing.conflicts:
            against += f" ({', '.join(standing.conflicts)})"
    return [f"#{pull.issue}", branch, pull.author, pull.state, against, str(pull.attempts)]
