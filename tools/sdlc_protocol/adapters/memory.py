"""A world held in memory, for tests and simulations of the merge gate.

The gate's decisions need only a few facts about the branches: their
commits, how each stands to main, the state of the required check, and
what a merge did. This world keeps those facts in dictionaries that a test
sets by hand, so a decision can be checked against every standing and
every check state without git or GitHub, and a simulation can hold a check
pending for as long as it means CI to take. Commits are made-up shas,
numbered so the same steps give the same shas.

A merge lands a branch that stands clean or behind, as a gate that merges
only cleanly would; after it, the merged branch has nothing new, and every
other branch that stood clean is behind, since main moved past it. A
conflicting branch is MainMoved, and fail_next_merge makes the next merge
raise whatever a test needs, such as a push refused because main moved.
"""

import hashlib
import threading
from datetime import datetime, timedelta, timezone

from ..gate import (
    BEHIND, CLEAN, CONFLICTING, NO_BRANCH, NO_CHECK, NOTHING_NEW, UNKNOWN, MainMoved, MergeRefused,
    MergeResult, Standing,
)

START = datetime(2026, 1, 1, tzinfo=timezone.utc)


class MemoryWorld:
    """Branches, standings and checks set by hand, and a clock that moves
    only when told to."""

    def __init__(self, main="main", now=START):
        self.main = main
        self.lock = threading.RLock()
        self.commits = 0
        self.main_commit = self.new_sha()
        self.branches = {}
        self.words = {}
        self.conflicts = {}
        self.checks = {}
        self.clock = now
        self.merges = []
        self.next_failure = None

    def __enter__(self):
        self.lock.acquire()
        return self

    def __exit__(self, *exc):
        self.lock.release()

    # Setting the world up.

    def new_sha(self):
        self.commits += 1
        return hashlib.sha1(f"commit {self.commits}".encode("ascii")).hexdigest()

    def push(self, branch, standing=CLEAN, conflicts=()):
        """Push a new commit to the branch, standing as given to main."""
        self.branches[branch] = self.new_sha()
        self.set_standing(branch, standing, conflicts)
        return self.branches[branch]

    def set_standing(self, branch, word, conflicts=()):
        self.words[branch] = word
        self.conflicts[branch] = tuple(conflicts)

    def set_check(self, branch, state):
        self.checks[branch] = state

    def move_main(self):
        """Someone else moved main: every branch that stood clean is behind."""
        self.main_commit = self.new_sha()
        self.fall_behind()
        return self.main_commit

    def fail_next_merge(self, error):
        """The next merge raises error instead of landing."""
        self.next_failure = error

    def advance(self, **delta):
        """Move the clock on, by timedelta's keywords: advance(hours=25)."""
        self.clock += timedelta(**delta)
        return self.clock

    # The World the gate decides over.

    def main_sha(self):
        return self.main_commit

    def branch_sha(self, branch):
        return self.branches.get(branch, "")

    def standing(self, branch):
        if branch not in self.branches:
            return NO_BRANCH
        return self.words.get(branch, UNKNOWN)

    def standing_detail(self, branch):
        return Standing(self.standing(branch), self.main_commit, self.conflicts.get(branch, ()))

    def check(self, branch):
        return self.checks.get(branch, NO_CHECK)

    def merge(self, branch):
        """Land a clean or behind branch on main, or raise why not."""
        if self.next_failure is not None:
            error, self.next_failure = self.next_failure, None
            raise error
        word = self.standing(branch)
        if word == CONFLICTING:
            raise MainMoved("the branch conflicts with main", self.main_commit, self.conflicts.get(branch, ()))
        if word not in (CLEAN, BEHIND):
            raise MergeRefused(f"{branch} stands {word} to {self.main}", self.main_commit)
        return self.land(branch)

    def now(self):
        return self.clock

    # What a merge does to the rest of the world.

    def land(self, branch):
        self.main_commit = self.new_sha()
        self.merges.append((branch, self.main_commit))
        self.set_standing(branch, NOTHING_NEW)
        self.fall_behind()
        return MergeResult(self.main_commit)

    def fall_behind(self):
        for branch, word in self.words.items():
            if word == CLEAN:
                self.words[branch] = BEHIND
