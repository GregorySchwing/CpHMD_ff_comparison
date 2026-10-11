"""Team mode's clerical work, done by code so no Claude usage goes to it.

Team mode runs CooperBench's best-scoring configuration on GitHub: a lead
turns an epic into a plan and pre-assigned sub-issues, members work them
in parallel on claim branches, and the lead integrates at the end. Claude
does the thinking; this module does the rest, three ways, each run by a
workflow step through the command line (cli.py's "team" subcommand):

  seed   read the lead's plan (the JSON claude-code-action returns under
         --json-schema), write it into the epic, open one sub-issue per
         task labelled todo and agent:<owner>, and link each under the epic.
  route  after an issue event, say what the lead does next:
         "plan <epic>", "integrate <epic>", "unblock <task>" or "none".
  next   for one member, say what it does next: "carry <issue>" (its own
         claimed task needs work), "take <issue>" (its pre-assigned task,
         else its share of the unclaimed blockers and todos) or "none".

Every call to GitHub goes through a Client, the small interface below that
adapters/github.py implements over REST and a test implements in memory,
so the rules here can be pinned without a network. The answers are the
one-line strings the workflows read; what the clerk cannot do is a
TeamError, whose message never holds a token.
"""

import json
import re
from typing import Protocol, runtime_checkable

from .claims import PREFIX as CLAIM_PREFIX, number_of
# The label names and the two prefixes are spelled once, in labels.py: a
# team member names itself with agent:<name>; any other account, and the
# plain worker that runs as the main account, with account:<name>.
from .labels import (ACCOUNT_PREFIX, AGENT_PREFIX, BLOCKED as BLOCKED_LABEL,
                     BLOCKER as BLOCKER_LABEL, EPIC as EPIC_LABEL, TODO as TODO_LABEL,
                     agent_label as agent_of, holders, label_names)

# The one check the ruleset requires and a member can fix; the board's own
# job and Claude's runs also leave check runs, which no member can mend.
REQUIRED_CHECK = "check"
PLAN_MARKER = "<!-- team:plan -->"
TASKS_LINE = re.compile(r"<!-- team:tasks ([0-9 ]*) -->")
EPIC_LINE = re.compile(r"^Epic: #(\d+)\s*$", re.MULTILINE)
# A ready pull request with one of these states has fallen behind its
# base or conflicts with it; any other state waits on checks or a merge.
NEEDS_REBASE = {"behind", "dirty"}
FAILED = {"failure", "timed_out"}
AGENT_COLOUR = "c5def5"


class TeamError(Exception):
    """What the clerk could not do; the message never holds the token."""


class ClientError(Exception):
    """A call GitHub refused or did not answer. An adapter raises this, or a
    subclass of it, so the clerk can tell the refusals it tolerates (a
    sub-issue link that failed) from a bug in its own code; the message
    never holds the token."""


@runtime_checkable
class Client(Protocol):
    """The few calls the clerk makes, as one object handed in by the caller.

    Issues, pull requests and check runs are the dictionaries GitHub's REST
    API returns, so the clerk reads the same fields whoever made them: an
    issue has "number", "id", "title", "body", "state" and "labels" (each a
    dictionary with a "name"); a pull request has "draft", "mergeable_state"
    and "head" with its "sha"; a check run has "name" and "conclusion".
    A call GitHub refuses raises ClientError.
    """

    def issue(self, number: int) -> dict:
        """One issue, open or closed, by number."""

    def open_issues(self) -> list[dict]:
        """Every open issue of the repository."""

    def branch_names(self) -> list[str]:
        """The names of the repository's branches."""

    def open_pull(self, branch: str) -> dict | None:
        """The open pull request from the branch, in the detail that holds
        "mergeable_state" and the head's "sha", or None when there is none."""

    def check_runs(self, sha: str) -> list[dict]:
        """The check runs on one commit."""

    def create_label(self, name: str, colour: str) -> None:
        """Make the label; one that exists already is left as it is, since
        GitHub answers a second creation with 422 and that is not a failure."""

    def create_issue(self, title: str, body: str, labels: list[str]) -> dict:
        """Open an issue and return it, with its "number" and "id"."""

    def add_sub_issue(self, parent: int, child_id: int) -> None:
        """Hang the issue with that id under the parent issue's number."""

    def set_body(self, number: int, body: str) -> None:
        """Replace an issue's body."""

    def remove_label(self, number: int, name: str) -> None:
        """Take a label off an issue; one the issue does not carry is fine."""


# --- label rules --------------------------------------------------------------
# label_names, agent_of (labels.agent_label) and holders come from labels.py
# and are re-exported here, so a reader of the clerk finds them by the names
# the clerk uses them under.

def is_mine(issue, agent):
    """An issue this member was given or claimed, and nobody else holds:
    a second account can claim a member's task off the queue before the
    member gets to it, and its branch is then that account's."""
    return holders(issue) == {agent}


def epic_of(issue):
    """The epic a task belongs to, from the "Epic: #N" line seed writes."""
    found = EPIC_LINE.search(issue.get("body") or "")
    return int(found.group(1)) if found else None


def task_numbers(epic):
    """The task issues seed opened for an epic, from its machine line."""
    found = TASKS_LINE.search(epic.get("body") or "")
    return [int(n) for n in found.group(1).split()] if found else []


def open_issues(client):
    """Every open issue (pull requests, which GitHub's issue list mixes in
    and an adapter may pass on, left out)."""
    return [issue for issue in client.open_issues() if "pull_request" not in issue]


def claim_branches(client):
    """The issue numbers that have a claude/<number> branch on the remote."""
    numbers = (number_of(name) for name in client.branch_names())
    return {number for number in numbers if number is not None}


# --- route ------------------------------------------------------------------

def route(client, event, number, label=None):
    """What the lead does after an issue event: plan a newly labelled epic,
    help a member that labelled its task blocked, integrate an epic whose
    last task just closed, or nothing. A dispatch with the epic's own
    number and the event closed integrates it too, which is the way to
    recover an integration whose event was lost."""
    issue = client.issue(number)
    is_open = issue.get("state") == "open"
    if event == "labeled" and label == EPIC_LABEL:
        planned = PLAN_MARKER in (issue.get("body") or "")
        return f"plan {number}" if is_open and not planned else "none"
    if event == "labeled" and label == BLOCKED_LABEL:
        waiting = is_open and BLOCKED_LABEL in label_names(issue) and epic_of(issue)
        return f"unblock {number}" if waiting else "none"
    if event == "closed":
        return integration_due(client, epic_of(issue) or (number if task_numbers(issue) else None))
    return "none"


def integration_due(client, epic_number):
    """"integrate <epic>" when every task of an open epic is closed and the
    lead has not yet claimed the integration, or its claim needs work: two
    tasks closing together must not start two integrations."""
    if epic_number is None:
        return "none"
    epic = client.issue(epic_number)
    tasks = task_numbers(epic)
    if epic.get("state") != "open" or not tasks:
        return "none"
    if any(client.issue(n).get("state") != "closed" for n in tasks):
        return "none"
    if epic_number in claim_branches(client) and not claim_needs_work(client, epic_number):
        return "none"
    return f"integrate {epic_number}"


# --- next -------------------------------------------------------------------

def claim_needs_work(client, number):
    """Whether an open claim has something for its holder to do: no pull
    request yet, a draft, a base it has fallen behind or conflicts with, or
    a failed required check. One that is ready and only waiting for its
    checks or its merge needs nothing, and a run spent on it would be usage
    wasted."""
    pull = client.open_pull(f"{CLAIM_PREFIX}{number}")
    if not pull:
        return True
    if pull.get("draft") or pull.get("mergeable_state") in NEEDS_REBASE:
        return True
    runs = client.check_runs(pull["head"]["sha"]) or []
    return any(run.get("name") == REQUIRED_CHECK and run.get("conclusion") in FAILED
               for run in runs)


class Snapshot:
    """One look at the open issues and the claim branches, shared by every
    member's answer so the members of one run divide the work alike."""

    def __init__(self, client):
        self.client = client
        self.issues = sorted(open_issues(client), key=lambda issue: issue["number"])
        self.claimed = claim_branches(client)
        self.needs_work = {}

    def claim_needs_work(self, number):
        if number not in self.needs_work:
            self.needs_work[number] = claim_needs_work(self.client, number)
        return self.needs_work[number]

    def free(self, issue):
        return issue["number"] not in self.claimed and BLOCKED_LABEL not in label_names(issue)

    def own_move(self, agent):
        """A member's own work: a claim of its own that needs work, else
        the task the lead gave it; None when it has neither."""
        mine = [i for i in self.issues if is_mine(i, agent) and BLOCKED_LABEL not in label_names(i)]
        for issue in mine:
            if issue["number"] in self.claimed and self.claim_needs_work(issue["number"]):
                return f"carry {issue['number']}"
        for issue in mine:
            if self.free(issue) and TODO_LABEL in label_names(issue) and EPIC_LABEL not in label_names(issue):
                return f"take {issue['number']}"
        return None

    def pool(self):
        """The shared work: unclaimed issues labelled blocker or todo that
        no member was given and that are not epics, blockers first, since a
        blocker stops everyone (spec/PROTOCOL.md)."""
        def shared(issue):
            labels = label_names(issue)
            return (self.free(issue) and not holders(issue) and EPIC_LABEL not in labels
                    and (TODO_LABEL in labels or BLOCKER_LABEL in labels))
        return sorted((i for i in self.issues if shared(i)),
                      key=lambda i: (BLOCKER_LABEL not in label_names(i), i["number"]))


def next_for(client, agent, team=None):
    """The member's next move: its own claims first, then its pre-assigned
    task, then its share of the shared work. Every member of a run asks at
    the same moment, so the shared work is dealt out: the members with no
    work of their own, in team order, take the first, second, third issue,
    and one with no issue left for it waits, rather than start a session
    that would race a peer for the same claim."""
    team = list(team or [agent])
    if agent not in team:
        team.append(agent)
    snap = Snapshot(client)
    own = snap.own_move(agent)
    if own:
        return own
    idle = [member for member in team if member == agent or not snap.own_move(member)]
    pool = snap.pool()
    slot = idle.index(agent)
    return f"take {pool[slot]['number']}" if slot < len(pool) else "none"


def next_move(client, agent, team=None, label=None):
    """The member's next move as the worker workflow asks for it: next_for,
    and when the move is a take and the member named its label, the label
    is made first, so the session that follows can put it on the issue (a
    label gh cannot find makes `gh issue edit --add-label` add none)."""
    answer = next_for(client, agent, team)
    if label and answer.startswith("take "):
        ensure_label(client, label)
    return answer


# --- seed -------------------------------------------------------------------

def owners_for(tasks, agents):
    """Each task's owner: the one the lead named when it is on the team;
    a task whose owner is not goes first to a member no task names yet, so
    no member sits idle while another holds two, then to the members in
    turn."""
    named = [task.get("owner") if task.get("owner") in agents else None for task in tasks]
    spare = [agent for agent in agents if agent not in named]
    owners, turn = [], 0
    for owner in named:
        if owner is None and spare:
            owner = spare.pop(0)
        elif owner is None:
            owner, turn = agents[turn % len(agents)], turn + 1
        owners.append(owner)
    return owners


def task_body(epic_number, task):
    files = "\n".join(f"- {path}" for path in task.get("files", [])) or "- (the lead named none)"
    return (f"Epic: #{epic_number}\n\n{task.get('body', '').strip()}\n\n"
            f"## Files to touch\n\n{files}\n\n"
            "Stay inside these files; another member owns the rest of the epic "
            "(docs/agents/TEAM.md).\n")


def ensure_label(client, name):
    """Make a member's label; one that exists already is fine, and the
    client treats it so. A label gh cannot find makes `gh issue edit
    --add-label` add none of the labels it was given, so the label is made
    before it is needed."""
    client.create_label(name, AGENT_COLOUR)


def open_task(client, epic_number, task, owner):
    """Open one task issue and hang it under the epic as a sub-issue; the
    link is a convenience for people, so its failure is not one."""
    label = AGENT_PREFIX + owner
    ensure_label(client, label)
    made = client.create_issue(task["title"], task_body(epic_number, task), [TODO_LABEL, label])
    try:
        client.add_sub_issue(epic_number, made["id"])
    except ClientError:
        pass
    return made["number"]


def plan_section(plan, made):
    """What seed appends to the epic: the lead's plan, the task list, and
    the machine line route reads."""
    lines = [f"- [ ] #{number} {title} (agent:{owner})" for number, title, owner in made]
    numbers = " ".join(str(number) for number, _, _ in made)
    return (f"\n\n{PLAN_MARKER}\n## Plan, by the lead\n\n{plan.strip()}\n\n"
            "## Tasks\n\n" + "\n".join(lines) + f"\n\n<!-- team:tasks {numbers} -->\n")


def seed(client, epic_number, plan, agents):
    """Open the epic's tasks and write the plan into it, once."""
    epic = client.issue(epic_number)
    if PLAN_MARKER in (epic.get("body") or ""):
        return "already planned"
    tasks = plan.get("tasks") or []
    if not tasks or not agents:
        raise TeamError("the plan has no tasks, or the team has no members")
    made = []
    for task, owner in zip(tasks, owners_for(tasks, agents)):
        made.append((open_task(client, epic_number, task, owner), task["title"], owner))
    body = (epic.get("body") or "") + plan_section(plan.get("plan", ""), made)
    client.set_body(epic_number, body)
    # An epic opened with the Task form carries todo, which would put it on
    # the board's queue for any account to take whole; its tasks are the work.
    client.remove_label(epic_number, TODO_LABEL)
    return "seeded " + " ".join(str(number) for number, _, _ in made)


# --- the team variable ----------------------------------------------------------

def team_members(text):
    """The members from the CLAUDE_TEAM variable's JSON list."""
    try:
        members = json.loads(text or "[]")
    except ValueError:
        raise TeamError("CLAUDE_TEAM is not a JSON list") from None
    if not isinstance(members, list) or not all(isinstance(m, str) and m for m in members):
        raise TeamError("CLAUDE_TEAM must be a JSON list of names, as in [\"w1\",\"w2\"]")
    return members
