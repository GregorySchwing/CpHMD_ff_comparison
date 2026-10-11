"""The labels the claim protocol reads and writes.

GitHub Issues is the queue, and labels are how an issue says where it
stands: todo in the queue, claimed by a branch, in-review with a pull
request waiting, blocker stopping everyone. Two prefixes name who holds
work: account:<name> for a Claude account, agent:<name> for a team member.
The board, the protocol check and the team clerk all read the same names,
so they are spelled once, here, and every consumer imports them.
"""

TODO = "todo"
CLAIMED = "claimed"
IN_REVIEW = "in-review"
BLOCKER = "blocker"
EPIC = "epic"
BLOCKED = "blocked"

# An account's label and a team member's label: the name after the colon
# is who holds the work.
ACCOUNT_PREFIX = "account:"
AGENT_PREFIX = "agent:"
MAIN_ACCOUNT = "main"

# The labels a repository needs for GitHub Issues to be the queue, and the
# two team mode adds: name, colour (six hex digits, as the API takes them),
# description. A repository is set up by making every one of these.
LABELS = (
    (TODO, "0e8a16", "In the queue, open for a claim"),
    (CLAIMED, "fbca04", "A claude/<issue> branch holds it"),
    (IN_REVIEW, "1d76db", "A pull request is open and waiting"),
    (BLOCKER, "b60205", "Stops everyone until it is settled"),
    (ACCOUNT_PREFIX + MAIN_ACCOUNT, "5319e7", "Work for the main Claude account"),
    (EPIC, "8b5cf6", "A spec the team lead plans into tasks (docs/agents/TEAM.md)"),
    (BLOCKED, "e99695", "A member waits on a peer or the lead (docs/agents/TEAM.md)"),
)


def label_names(issue):
    """The names of an issue's labels, sorted, whether the API gave label
    objects or a fixture gave plain strings. Sorted, so a board rendered
    from the same issue twice shows the labels in the same order."""
    return sorted(label if isinstance(label, str) else label["name"]
                  for label in (issue or {}).get("labels") or [])


def names_with_prefix(labels, prefix):
    """What follows the prefix in every label that carries it, in order."""
    return [name[len(prefix):] for name in labels if name.startswith(prefix)]


def account_label(issue):
    """The account an issue's label names (account:main gives main), or
    None when no label does."""
    accounts = names_with_prefix(label_names(issue), ACCOUNT_PREFIX)
    return accounts[0] if accounts else None


def agent_label(issue):
    """The team member an issue's label names (agent:w1 gives w1), or None."""
    agents = names_with_prefix(label_names(issue), AGENT_PREFIX)
    return agents[0] if agents else None


def holders(issue):
    """Everyone an issue's labels name, member or account, as a set."""
    names = label_names(issue)
    return set(names_with_prefix(names, ACCOUNT_PREFIX)) | set(names_with_prefix(names, AGENT_PREFIX))


def reserved_for_team(labels):
    """Whether these labels keep an issue off every account's queue: an
    epic, which team mode's lead splits into tasks, or a task the lead
    gave a member (label agent:<name>), which only that member takes."""
    return EPIC in labels or any(name.startswith(AGENT_PREFIX) for name in labels)
