"""The claim: a branch named for an issue and the file that describes it.

A claim is the first push of the branch claude/<issue number>. GitHub
refuses the second push of the same name, so the branch is the lock, and
its first commit adds docs/agents/claims/<issue number>.md, which names
the account, the plan and the files the work expects to touch. The board
reads that file until the pull request has files of its own, and the
coordinator releases a claim whose branch has had no commit for a day.
This module holds those conventions and the parsing, and nothing that
talks to git or to GitHub.
"""

import re
from datetime import datetime, timedelta, timezone

PREFIX = "claude/"
CLAIM_DIR = "docs/agents/claims"
# A claim with no commit for this long is released by the coordinator,
# and the board asks a pull request idle this long whether it is parked.
IDLE_AFTER = timedelta(hours=24)

ACCOUNT_LINE = re.compile(r"^\W*account\W*:\s*(.+?)\s*$", re.IGNORECASE)
# A note in parentheses after the account, like the claim template's
# "(label account:main)", is for the reader, not the board.
TRAILING_NOTE = re.compile(r"\s*\(.*\)\s*$")
FILES_HEADING = re.compile(r"^#+\s*files to touch\s*$", re.IGNORECASE)
LIST_ITEM = re.compile(r"^\s*[-*]\s+(.+?)\s*$")
UNNAMED = "unnamed"


def branch_of(number, prefix=PREFIX):
    """The claim branch for an issue: claude/<number>. The prefix is an
    argument so a gate policy that names another one uses the same rule."""
    return f"{prefix}{number}"


def is_claim_branch(branch, prefix=PREFIX):
    return branch.startswith(prefix)


def key_of(branch, prefix=PREFIX):
    """What follows claude/ in a claim branch: the issue number by protocol,
    kept as text because a hand-made branch may put something else there
    and the board still shows it. None for a branch that is no claim."""
    return branch[len(prefix):] if is_claim_branch(branch, prefix) else None


def number_of(branch, prefix=PREFIX):
    """The issue number a claim branch names, or None when the branch is no
    claim or what follows claude/ is not a number. Decimal digits only, as
    the template's check matched them with \\d: a superscript passes
    isdigit() but not int(), and a branch name must never raise."""
    key = key_of(branch, prefix)
    return int(key) if key is not None and key.isdecimal() else None


def claim_path(number):
    """The claim file's path inside the repository, as text with slashes."""
    return f"{CLAIM_DIR}/{number}.md"


def by_number(value):
    """A sort key that puts numbers in order and anything else after them,
    so claude/9 comes before claude/100 and claude/x after both."""
    text = str(value)
    return (0, int(text), "") if text.isdigit() else (1, 0, text)


def claim_account(text):
    """The account a claim file names, from a line like "Account: main" in
    any of the Markdown spellings a hand might give it, without a note in
    parentheses after it; "unnamed" when no line does."""
    for line in text.splitlines():
        match = ACCOUNT_LINE.match(line)
        if match:
            return TRAILING_NOTE.sub("", match.group(1)).strip("`*_ ")
    return UNNAMED


def claim_files(text):
    """The paths a claim file lists under its "Files to touch" heading, one
    list line each, up to the next heading, sorted: what the board shows
    for a claim whose pull request does not exist yet."""
    files = []
    inside = False
    for line in text.splitlines():
        if FILES_HEADING.match(line):
            inside = True
            continue
        if inside and line.startswith("#"):
            break
        match = LIST_ITEM.match(line) if inside else None
        if match:
            files.append(match.group(1).strip("`"))
    return sorted(files)


def parse_time(value):
    """An ISO 8601 time as GitHub or git writes it, as an aware datetime; a
    datetime given is returned as it is. A time without a zone is UTC."""
    if isinstance(value, datetime):
        moment = value
    else:
        text = value.strip()
        if text.endswith("Z"):
            text = text[:-1] + "+00:00"
        moment = datetime.fromisoformat(text)
    if moment.tzinfo is None:
        moment = moment.replace(tzinfo=timezone.utc)
    return moment


def is_idle(last_change, now, idle_after=IDLE_AFTER):
    """Whether more than idle_after has passed since the last change: the
    one rule behind the board's idle comment and the coordinator's release."""
    return now - parse_time(last_change) > idle_after


def release_due(claims, now, idle_after=IDLE_AFTER):
    """The claims the coordinator releases: those whose branch has had no
    commit for more than idle_after. The claims are a mapping, or pairs,
    from the claim's key (its branch, or the issue number) to the time of
    the branch's last commit, as a datetime or as GitHub writes it; the
    answer is the keys due, in issue order, so the comments the coordinator
    posts come out the same for the same state."""
    pairs = claims.items() if hasattr(claims, "items") else claims
    due = [key for key, last_commit in pairs if is_idle(last_commit, now, idle_after)]
    return sorted(due, key=issue_order)


def issue_order(key):
    """A sort key for a claim given as its branch or as its issue number."""
    text = str(key)
    return by_number(key_of(text) if is_claim_branch(text) else text)
