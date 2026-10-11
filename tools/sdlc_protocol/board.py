"""The board: what a repository's open work looks like, as one page.

The board is the file STATE.md on the branch "board": the open claims, the
open pull requests and how each stands to the default branch, the queue of
unclaimed issues, the blockers and the pull requests merged last, so every
session and every account reads the same page instead of asking the API.
This module is the pure part. A Snapshot holds what GitHub knows, filled by
an adapter or read back from JSON; render turns it into the page, with no
wall-clock value and every row sorted, so the same state gives the same
bytes and a workflow writes the branch only when something moved; and
comments_due says which pull requests the board would speak to, and what
it would say, for the adapter to post once each. Nothing here calls GitHub.
"""

import json
import unicodedata
from dataclasses import dataclass, field
from pathlib import Path
from typing import NamedTuple

from . import labels
from .claims import (
    IDLE_AFTER, by_number, claim_account, claim_files, is_claim_branch, is_idle, key_of,
)

RECENT_MERGED = 5
# Where a workflow writes the board when nothing else is named: a file in
# the current directory, since the workflow checks out the branch "board".
STATE_PATH = Path("STATE.md")
NO_CLAIM_FILE = "no claim file yet"
# A claim branch that outlived its issue: the row stays, and says what to do.
ISSUE_CLOSED = "issue closed"
DELETE_BRANCH = "delete the branch"
NO_PULL_YET = "(no pull request yet)"

MAIN_MOVED_MARKER = "<!-- sdlc-board: main-moved {sha} -->"
IDLE_MARKER = "<!-- sdlc-board: idle -->"
IDLE_TEXT = (
    "idle: no change for more than a day; if this work is parked, say so "
    "here, or the coordinator may release the claim."
)

CLAIM_HEADING = "## Claims (open claude/<issue> branches and docs/agents/claims/<issue>.md)"
CLAIM_COLUMNS = ("Issue", "Branch", "Account", "PR", "Files it touches")
PULL_HEADING = "## Open pull requests"
PULL_COLUMNS = ("PR", "Title", "Branch", "State (clean, behind, conflicting, blocked, unknown)", "Last update")
QUEUE_HEADING = "## Queue (open issues labelled todo, unclaimed)"
ISSUE_COLUMNS = ("Issue", "Title", "Labels")
BLOCKER_HEADING = "## Blockers (open issues labelled blocker)"
MERGED_HEADING = f"## Recently merged (the last {RECENT_MERGED})"
MERGED_COLUMNS = ("PR", "Title", "Branch", "Merged")

# The five words the board uses for how a pull request stands to the
# default branch, and GitHub's mergeable_state words that mean each.
STATE_WORDS = ("clean", "behind", "conflicting", "blocked", "unknown")
GITHUB_STATES = {"has_hooks": "clean", "unstable": "clean", "dirty": "conflicting"}


# --- the snapshot -----------------------------------------------------------

@dataclass
class Issue:
    """An open issue as the queue and the blockers show it."""

    number: int
    title: str = ""
    labels: list = field(default_factory=list)
    updated_at: str = ""

    @classmethod
    def from_dict(cls, data):
        return cls(
            number=int(data["number"]),
            title=data.get("title") or "",
            labels=labels.label_names(data),
            updated_at=data.get("updated_at") or "",
        )

    def to_dict(self):
        return {"number": self.number, "title": self.title, "labels": list(self.labels), "updated_at": self.updated_at}


@dataclass
class Pull:
    """An open pull request with its head branch, GitHub's word for how it
    stands to the default branch (mergeable_state, as the API gives it) and
    the paths it changes. The dict shape is the API's, head and all, so a
    pull request detail from GitHub is a fixture as it is."""

    number: int
    title: str = ""
    branch: str = ""
    head_sha: str = ""
    mergeable_state: str | None = None
    updated_at: str = ""
    files: list = field(default_factory=list)

    @property
    def state(self):
        """The board's word for the pull request's standing."""
        return state_word(self.mergeable_state)

    @classmethod
    def from_dict(cls, data):
        head = data.get("head") or {}
        return cls(
            number=int(data["number"]),
            title=data.get("title") or "",
            branch=head.get("ref") or "",
            head_sha=head.get("sha") or "",
            mergeable_state=data.get("mergeable_state"),
            updated_at=data.get("updated_at") or "",
            files=[file_name(f) for f in data.get("files") or []],
        )

    def to_dict(self):
        return {
            "number": self.number,
            "title": self.title,
            "head": {"ref": self.branch, "sha": self.head_sha},
            "mergeable_state": self.mergeable_state,
            "updated_at": self.updated_at,
            "files": list(self.files),
        }


@dataclass
class Merged:
    """A pull request that landed: enough for the recently merged table."""

    number: int
    title: str = ""
    branch: str = ""
    merged_at: str = ""

    @classmethod
    def from_dict(cls, data):
        return cls(
            number=int(data["number"]),
            title=data.get("title") or "",
            branch=(data.get("head") or {}).get("ref") or "",
            merged_at=data.get("merged_at") or "",
        )

    def to_dict(self):
        return {"number": self.number, "title": self.title, "head": {"ref": self.branch}, "merged_at": self.merged_at}


@dataclass
class ClaimIssue:
    """The issue behind a claim branch: its state and its label names, so
    the board can say that a branch outlived its issue, and can name the
    account from the issue's label once the claim file is gone."""

    state: str = ""
    labels: list = field(default_factory=list)

    @classmethod
    def from_dict(cls, data):
        return cls(state=data.get("state") or "", labels=labels.label_names(data))

    def to_dict(self):
        return {"state": self.state, "labels": list(self.labels)}


@dataclass
class Snapshot:
    """Everything the board renders, in the shape a --dry-run fixture has:
    the default branch's commit and name, the open issues (pull requests
    left out), the open pull requests with their standing and files, the
    branch names, the claim files by issue key (None for a branch with no
    file yet), the issue behind each claim branch, and the merges last."""

    sha: str = ""
    default_branch: str = "main"
    issues: list = field(default_factory=list)
    pulls: list = field(default_factory=list)
    branches: list = field(default_factory=list)
    claims: dict = field(default_factory=dict)
    claim_issues: dict = field(default_factory=dict)
    merged: list = field(default_factory=list)

    @classmethod
    def from_dict(cls, data):
        """A snapshot from its JSON shape, which is also what the API gives:
        labels may be objects or names, branches objects or names."""
        return cls(
            sha=data.get("sha") or "",
            default_branch=data.get("default_branch") or "main",
            issues=[Issue.from_dict(i) for i in data.get("issues") or []],
            pulls=[Pull.from_dict(p) for p in data.get("pulls") or []],
            branches=[branch_name(b) for b in data.get("branches") or []],
            claims={str(k): v for k, v in (data.get("claims") or {}).items()},
            claim_issues={str(k): ClaimIssue.from_dict(v or {}) for k, v in (data.get("claim_issues") or {}).items()},
            merged=[Merged.from_dict(m) for m in data.get("merged") or []],
        )

    def to_dict(self):
        return {
            "sha": self.sha,
            "default_branch": self.default_branch,
            "issues": [i.to_dict() for i in self.issues],
            "pulls": [p.to_dict() for p in self.pulls],
            "branches": list(self.branches),
            "claims": dict(self.claims),
            "claim_issues": {k: v.to_dict() for k, v in self.claim_issues.items()},
            "merged": [m.to_dict() for m in self.merged],
        }

    @classmethod
    def from_json(cls, text):
        return cls.from_dict(json.loads(text))

    def to_json(self):
        return json.dumps(self.to_dict(), indent=2) + "\n"


def branch_name(branch):
    """A branch as a name, whether the API's object or a fixture's string."""
    return branch if isinstance(branch, str) else branch["name"]


def file_name(changed):
    """A changed file as a path, whether the API's object from the pull
    request's files or a fixture's string."""
    return changed if isinstance(changed, str) else changed["filename"]


def as_snapshot(value):
    """A Snapshot from one, or from the dict shape a fixture or an API
    adapter holds, so render takes either."""
    return value if isinstance(value, Snapshot) else Snapshot.from_dict(value)


def as_pull(value):
    return value if isinstance(value, Pull) else Pull.from_dict(value)


# --- the rules --------------------------------------------------------------

def state_word(mergeable_state):
    """One of the five words the board uses for a pull request, from
    GitHub's mergeable_state. GitHub's has_hooks and unstable are mergeable
    states (the checks are what they are), so they read as clean; dirty is
    conflicting; blocked is a computed state of its own; unknown is kept for
    a state GitHub has not computed yet. The board's own words pass through,
    so a snapshot that already speaks them renders the same."""
    if mergeable_state in STATE_WORDS:
        return mergeable_state
    return GITHUB_STATES.get(mergeable_state, "unknown")


def pr_state(pull):
    """The board's word for a pull request, given as a Pull or as the dict
    the API returns."""
    return state_word(pull.mergeable_state if isinstance(pull, Pull) else pull.get("mergeable_state"))


def short_sha(sha):
    return (sha or "")[:7] or "unknown"


def issue_ref(key):
    return f"#{key}" if key.isdigit() else key


def claim_branches(branches):
    return [branch_name(b) for b in branches if is_claim_branch(branch_name(b))]


def claimed_numbers(branches):
    """The issue numbers that have a claim branch."""
    return {int(key_of(b)) for b in claim_branches(branches) if key_of(b).isdigit()}


def pulls_behind_main(pulls):
    """The pull requests the main-moved rule speaks to: behind or
    conflicting, and only those."""
    return [p for p in pulls if pr_state(p) in ("behind", "conflicting")]


def idle_pulls(pulls, now, idle_after=IDLE_AFTER):
    """The pull requests with no change for more than a day."""
    return [p for p in pulls if pull_is_idle(as_pull(p), now, idle_after)]


def pull_is_idle(pull, now, idle_after=IDLE_AFTER):
    """Whether one pull request has waited too long. One with no update
    time is not idle rather than failing the run: there is nothing to
    measure, and GitHub always gives one."""
    return bool(pull.updated_at) and is_idle(pull.updated_at, now, idle_after)


def main_moved_text(state, default_branch, sha7):
    how = "behind" if state == "behind" else "conflicting with"
    return (
        f"main moved: this branch is {how} {default_branch} at {sha7}; "
        f"rebase or merge {default_branch} before the next push."
    )


class DueComment(NamedTuple):
    """A comment the board owes a pull request: posted only when no comment
    on it carries the marker, which is the adapter's job."""

    number: int
    marker: str
    text: str

    @property
    def kind(self):
        return "idle" if self.marker == IDLE_MARKER else "main-moved"


def comments_due(snapshot, event_name, ref_name, now, idle_after=IDLE_AFTER):
    """What the board says after this event, as (pull number, marker, text).
    After a push to the default branch, every pull request behind it or in
    conflict with it hears so, once per commit of the default branch (the
    marker carries the short sha); a pull request nobody touched for a day
    is asked once whether it is parked. Main-moved comments come first,
    then idle ones, each in the order the snapshot lists the pull requests,
    as the template's board posted them."""
    snap = as_snapshot(snapshot)
    due = []
    if event_name == "push" and ref_name == snap.default_branch:
        sha7 = short_sha(snap.sha)
        marker = MAIN_MOVED_MARKER.format(sha=sha7)
        for pull in pulls_behind_main(snap.pulls):
            due.append(DueComment(pull.number, marker, main_moved_text(pull.state, snap.default_branch, sha7)))
    for pull in idle_pulls(snap.pulls, now, idle_after):
        due.append(DueComment(pull.number, IDLE_MARKER, IDLE_TEXT))
    return due


# --- the page ---------------------------------------------------------------

def cell(value):
    """A table cell: one line, ASCII, pipes escaped, so the board stays a
    valid table and passes a repository's own ASCII check."""
    text = str(value).replace("\r", " ").replace("\n", " ").replace("|", "\\|")
    return unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode().strip()


def table(columns, rows):
    """A Markdown table; an empty one keeps its header and says (none), so
    the layout is the same whatever the repository holds."""
    lines = ["| " + " | ".join(columns) + " |", "|" + " --- |" * len(columns)]
    for row in rows or [["(none)"] + [""] * (len(columns) - 1)]:
        lines.append("| " + " | ".join(cell(c) for c in row) + " |")
    return lines


def claim_row(key, branch, text, pull, issue):
    """One claim as the board shows it. The account is the claim file's, or
    the issue's label once the finishing pull request has deleted the file.
    A branch whose issue is closed is a leftover: the row says so in place
    of the pull request and asks for the branch to go in place of the files."""
    closed = issue is not None and issue.state == "closed"
    if pull:
        files = sorted(pull.files)
    else:
        files = claim_files(text) if text else []
    if pull:
        pull_cell = f"#{pull.number}"
    elif closed:
        pull_cell = ISSUE_CLOSED
    else:
        pull_cell = NO_PULL_YET
    if files:
        files_cell = ", ".join(files)
    elif closed and not pull:
        files_cell = DELETE_BRANCH
    else:
        files_cell = "-"
    if text:
        account = claim_account(text)
    else:
        account = labels.account_label({"labels": issue.labels} if issue else None) or NO_CLAIM_FILE
    return [issue_ref(key), branch, account, pull_cell, files_cell]


def claim_rows(snap):
    pulls_by_branch = {p.branch: p for p in snap.pulls}
    rows = []
    for branch in sorted(claim_branches(snap.branches), key=lambda b: by_number(key_of(b))):
        key = key_of(branch)
        rows.append(claim_row(key, branch, snap.claims.get(key), pulls_by_branch.get(branch), snap.claim_issues.get(key)))
    return rows


def pull_rows(snap):
    return [
        [f"#{p.number}", p.title, p.branch, p.state, p.updated_at]
        for p in sorted(snap.pulls, key=lambda p: p.number)
    ]


def issue_rows(snap, label, unclaimed_only):
    """The issues carrying a label, in number order; the queue leaves out
    an issue with the label claimed, a claim branch, or team mode's labels."""
    claimed = claimed_numbers(snap.branches)
    rows = []
    for issue in sorted(snap.issues, key=lambda i: i.number):
        names = sorted(issue.labels)
        if label not in names:
            continue
        if unclaimed_only and (labels.CLAIMED in names or issue.number in claimed or labels.reserved_for_team(names)):
            continue
        rows.append([f"#{issue.number}", issue.title, ", ".join(names)])
    return rows


def merged_rows(snap):
    return [
        [f"#{m.number}", m.title, m.branch, m.merged_at]
        for m in sorted(snap.merged, key=lambda m: m.merged_at, reverse=True)
    ]


def render(snapshot, default_branch=None):
    """STATE.md for a snapshot: the layout the protocol's STATE.md describes.
    The default branch named here wins over the snapshot's own, for a
    workflow that knows it from the event; None takes the snapshot's."""
    snap = as_snapshot(snapshot)
    branch = cell(default_branch or snap.default_branch or "main")
    lines = [
        "# Board",
        "",
        f"As of commit {cell(short_sha(snap.sha))} on {branch}. "
        "Rewritten by .github/workflows/board.yml; do not edit by hand.",
    ]
    sections = (
        (CLAIM_HEADING, CLAIM_COLUMNS, claim_rows(snap)),
        (PULL_HEADING, PULL_COLUMNS, pull_rows(snap)),
        (QUEUE_HEADING, ISSUE_COLUMNS, issue_rows(snap, labels.TODO, True)),
        (BLOCKER_HEADING, ISSUE_COLUMNS, issue_rows(snap, labels.BLOCKER, False)),
        (MERGED_HEADING, MERGED_COLUMNS, merged_rows(snap)),
    )
    for heading, columns, rows in sections:
        lines += ["", heading, ""] + table(columns, rows)
    return "\n".join(lines) + "\n"


def write_if_changed(path, text):
    """Write the board only when its bytes differ, so a caller comparing
    files sees no change when nothing moved; True when it was written. The
    one thing in this module that touches a file, kept here because the
    rendered page is ASCII by construction and this is where that is known."""
    path = Path(path)
    data = text.encode("ascii")
    if path.is_file() and path.read_bytes() == data:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return True
