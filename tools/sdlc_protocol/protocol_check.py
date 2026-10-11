"""Hold a claude/<issue> pull request to the claim protocol.

The protocol says what a claim looks like; this module makes the checkable
parts of it a list of problems, for CI to fail on, so a pull request that
forgot them cannot merge however the session that made it was feeling:

- A claim branch belongs to the repository, never a fork.
- The body closes its own issue, so the merge closes it.
- A draft still carries docs/agents/claims/<issue>.md: the claim stays
  visible exactly as long as it holds.
- A pull request marked ready has deleted it: finishing lets go.

Any other head branch, dependabot's or a person's, is not a claim and
passes untouched. The pull request is the object GitHub's event carries;
whether the claim file exists in the checkout is the caller's to say, so
nothing here reads a file.
"""

import re

from .claims import claim_path, number_of

# The forms GitHub links a pull request to an issue with; the protocol
# writes "Closes #N" and the rest are the same promise.
CLOSES = r"(?i)\b(?:close[sd]?|fix(?:e[sd])?|resolve[sd]?)\s+#%d\b"
PROTOCOL = "docs/agents/PROTOCOL.md"


def head_branch(pull):
    return ((pull or {}).get("head") or {}).get("ref") or ""


def claim_number(pull):
    """The issue a pull request claims, from its head branch; None when the
    branch is no claim."""
    return number_of(head_branch(pull))


def repository_of(side):
    """The full name of the repository a pull request's head or base is in."""
    return ((side or {}).get("repo") or {}).get("full_name")


def closes_its_issue(pull, number):
    return re.search(CLOSES % number, pull.get("body") or "") is not None


def problems_for(pull, claim_file_exists):
    """Every way this pull request breaks the protocol; [] when it is no
    claim at all or keeps the rules. claim_file_exists says whether
    docs/agents/claims/<issue>.md is in the checkout: a bool, or a callable
    given the issue number, for a caller that would rather not parse the
    branch itself."""
    branch = head_branch(pull)
    number = claim_number(pull)
    if number is None:
        return []
    if callable(claim_file_exists):
        claim_file_exists = claim_file_exists(number)
    problems = []
    if repository_of(pull.get("head")) != repository_of(pull.get("base")):
        problems.append(f"{branch}: a claim branch lives in this repository, not in a fork "
                        f"({PROTOCOL}, Claim)")
    if not closes_its_issue(pull, number):
        problems.append(f"the body does not close #{number}, so the merge would leave the issue open; "
                        f'write "Closes #{number}" ({PROTOCOL}, Finish)')
    path = claim_path(number)
    if pull.get("draft"):
        if not claim_file_exists:
            problems.append(f"{path} is missing: a draft carries its claim file ({PROTOCOL}, Claim)")
    elif claim_file_exists:
        problems.append(f"{path} is still here: a pull request marked ready has deleted its claim file "
                        f"({PROTOCOL}, Finish)")
    return problems
