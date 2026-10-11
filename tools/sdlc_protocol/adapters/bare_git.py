"""The world over a bare repository: the gate's merges made with git itself.

This is the coordbench sdlc arm's merge gate, moved here so the arm and
anything else that shares a bare repository between agents decide through
sdlc_protocol.gate. The world keeps one scratch clone of the remote, made
at first use in a temporary directory, runs git through subprocess with a
timeout, and never force-pushes, so the default branch only moves forward.
A merge is `git merge --no-ff --no-edit` onto main as the world last
fetched it, made with the gate's own committer identity and pushed without
force: a conflict, or a push refused because main moved meanwhile, is
MainMoved and leaves main as it was.

Every call is made holding the world's lock: `with world:` or `with
world.lock:`, or a single call, which takes the lock itself. Within one
`with world:` the world fetches once, at the first read, so one decision
sees one state of the remote, as the arm's gate fetched once per call;
outside it each read fetches. close frees the scratch clone for good, and
a later call is refused rather than cloning again.
"""

import os
import re
import subprocess
import tempfile
import threading
from datetime import datetime, timezone

from ..gate import (
    BEHIND, CLEAN, CONFLICTING, NO_BRANCH, NO_CHECK, NOTHING_NEW, UNKNOWN, GateFailed, MainMoved,
    MergeOutcome, MergeRefused, MergeResult, Standing, standings,
)

GIT_TIMEOUT_SECONDS = 300

# The identity the coordbench sdlc arm's gate has always merged with, kept
# so a merge made here is the same commit the arm's gate would have made.
GATE_NAME = "coordbench merge gate"
GATE_EMAIL = "gate@coordbench.local"

ISSUE_TAIL = re.compile(r"[0-9]+")


def git_options(name=GATE_NAME, email=GATE_EMAIL):
    """The options every git command runs with. The scratch clone is the
    gate's own: a signing key, a hook or a pager in the host's git
    configuration must not decide whether a merge lands."""
    return ("-c", f"user.name={name}", "-c", f"user.email={email}",
            "-c", "commit.gpgsign=false", "-c", f"core.hooksPath={os.devnull}")


GIT_OPTIONS = git_options()


def git_environment():
    """The host's environment, with git told never to wait for a password
    and to speak English, since the gate quotes its messages."""
    return {**os.environ, "GIT_TERMINAL_PROMPT": "0", "LC_ALL": "C"}


def run_git(args, cwd=None, options=GIT_OPTIONS):
    """One git command with the gate's options and a timeout. A git that
    hangs or is missing is GateFailed, which the caller can read."""
    try:
        return subprocess.run(["git", *options, *args], cwd=cwd, capture_output=True, text=True,
                              errors="replace", timeout=GIT_TIMEOUT_SECONDS, env=git_environment(),
                              stdin=subprocess.DEVNULL)
    except subprocess.TimeoutExpired:
        raise GateFailed(f"git {args[0]} took longer than {GIT_TIMEOUT_SECONDS} seconds and was stopped")
    except FileNotFoundError:
        raise GateFailed("git is not installed where the merge gate runs")


def last_line(result):
    lines = [line.strip() for line in (result.stderr + "\n" + result.stdout).splitlines() if line.strip()]
    return lines[-1] if lines else f"git exited with {result.returncode}"


def remote_main(remote, main="main"):
    """main's commit on the remote, or "" when it cannot be read."""
    try:
        result = run_git(["ls-remote", remote, f"refs/heads/{main}"])
    except GateFailed:
        return ""
    words = result.stdout.split()
    return words[0] if result.returncode == 0 and words else ""


def merge_message(branch):
    """The merge commit's message: a claim branch closes its issue."""
    tail = branch.rpartition("/")[2]
    return f"Merge {branch}: closes #{tail}" if ISSUE_TAIL.fullmatch(tail) else f"Merge {branch}"


class BareGitWorld:
    """The world over the bare repository at `remote` (a path or a URL git
    can fetch from), whose default branch is `main`."""

    def __init__(self, remote, main="main", name=GATE_NAME, email=GATE_EMAIL):
        self.remote = remote
        self.main = main
        self.main_ref = f"refs/remotes/origin/{main}"
        self.options = git_options(name, email)
        self.lock = threading.RLock()
        self.scratch = None
        self.closed = False
        self.depth = 0
        self.fresh = False

    # Holding the world: one fetch per hold.

    def __enter__(self):
        self.lock.acquire()
        self.depth += 1
        return self

    def __exit__(self, *exc):
        self.depth -= 1
        if self.depth == 0:
            self.fresh = False
        self.lock.release()

    def refresh(self):
        """Fetch, unless this hold of the world has fetched already."""
        if not self.fresh:
            self.fetch()
            self.fresh = True

    # The World the gate decides over.

    def main_sha(self):
        """main's commit, "" when the remote has no main."""
        with self:
            self.refresh()
            return self.commit_of(self.main_ref)

    def branch_sha(self, branch):
        with self:
            self.refresh()
            return self.branch_tip(branch)

    def standing(self, branch):
        return self.standing_detail(branch).word

    def standing_detail(self, branch):
        """How the branch stands to main, naming the conflicting files."""
        with self:
            self.refresh()
            return self.read_standing(branch)

    def check(self, branch):
        """A bare repository runs no checks."""
        return NO_CHECK

    def merge(self, branch):
        """Land the branch on main, or raise why not."""
        with self:
            self.refresh()
            tip, main = self.merge_inputs(branch)
            try:
                outcome = self.merge_commits(tip, main, merge_message(branch))
            finally:
                # The push moved main, or may have: the next read fetches.
                self.fresh = False
            return result_of(outcome)

    def now(self):
        return datetime.now(timezone.utc)

    def standings(self, branches):
        """main's commit and each branch's Standing, from one fetch."""
        with self:
            return standings(self, branches)

    # The scratch clone.

    @property
    def repo(self):
        return os.path.join(self.scratch.name, "repo")

    def git(self, *args, check=True):
        result = run_git(list(args), cwd=self.repo, options=self.options)
        if check and result.returncode != 0:
            raise GateFailed(f"the merge gate's git {args[0]} failed: {last_line(result)}")
        return result

    def ensure_clone(self):
        """Clone the remote without checking anything out; every merge
        checks out the main it merges into."""
        if self.scratch is not None:
            return
        if self.closed:
            raise GateFailed("the run has closed, and the merge gate with it: nothing more merges into "
                             f"{self.main}; submit what you have")
        scratch = tempfile.TemporaryDirectory(prefix="sdlc-gate-")
        result = run_git(["clone", "--quiet", "--no-checkout", self.remote, os.path.join(scratch.name, "repo")],
                         options=self.options)
        if result.returncode != 0:
            scratch.cleanup()
            raise GateFailed(f"the merge gate cannot clone the shared remote: {last_line(result)}")
        self.scratch = scratch

    def close(self):
        """Free the scratch clone for good. A call that reaches the world
        later is refused rather than cloning again, since nothing would free
        a second clone; one in progress finishes first, holding the lock."""
        with self.lock:
            self.closed = True
            self.fresh = False
            if self.scratch is not None:
                self.scratch.cleanup()
                self.scratch = None

    def fetch(self):
        self.ensure_clone()
        self.git("fetch", "--quiet", "--prune", "origin")

    def commit_of(self, ref):
        result = self.git("rev-parse", "--verify", "--quiet", f"{ref}^{{commit}}", check=False)
        return result.stdout.strip() if result.returncode == 0 else ""

    def branch_tip(self, branch):
        return self.commit_of(f"refs/remotes/origin/{branch}")

    def main_tip(self):
        main = self.commit_of(self.main_ref)
        if not main:
            raise GateFailed(f"the shared remote has no branch {self.main}; the run was not set up as "
                             f"the gate needs")
        return main

    def is_ancestor(self, older, newer):
        return self.git("merge-base", "--is-ancestor", older, newer, check=False).returncode == 0

    def tidy(self):
        """Leave the clone with no merge in progress and nothing stray, so a
        refused or interrupted merge never leaks into the next one."""
        self.git("merge", "--abort", check=False)
        self.git("reset", "--hard", "--quiet", check=False)
        self.git("clean", "-ffdq", check=False)

    def unmerged_paths(self):
        result = self.git("diff", "--name-only", "--diff-filter=U", check=False)
        return tuple(sorted(set(result.stdout.splitlines())))

    # The merge and the standing, as the arm's gate made them.

    def merge_inputs(self, branch):
        """The branch's tip and main's from the last fetch, refusing a branch
        that is not on the remote or holds nothing main lacks: the gate's
        decisions check both first, so this only answers a race."""
        main = self.main_tip()
        tip = self.branch_tip(branch)
        if not tip:
            raise MergeRefused(f"{branch} is not on the shared remote", main)
        if self.is_ancestor(tip, main):
            raise MergeRefused(f"nothing to merge: every commit on {branch} is already in {self.main}", main)
        return tip, main

    def merge_commits(self, tip, main, message):
        """Merge tip into main as `git merge --no-ff --no-edit` does and push
        main without force; a conflict or a refused push leaves main as it
        was on the remote."""
        self.tidy()
        self.git("checkout", "--quiet", "--force", "--detach", main)
        attempt = self.git("merge", "--no-ff", "--no-edit", "-m", message, tip, check=False)
        if attempt.returncode != 0:
            conflicts = self.unmerged_paths()
            self.tidy()
            return MergeOutcome(False, main, conflicts, main_moved=bool(conflicts), reason=last_line(attempt))
        merged = self.commit_of("HEAD")
        push = self.git("push", "--quiet", "origin", f"{merged}:refs/heads/{self.main}", check=False)
        if push.returncode != 0:
            # Someone moved main on the remote while the gate merged; the
            # answer is the same as for a conflict, since bringing the new
            # main in is again what the branch needs.
            return MergeOutcome(False, main, main_moved=True,
                                reason=f"the push of {self.main} was refused: {last_line(push)}")
        return MergeOutcome(True, merged)

    def read_standing(self, branch):
        """How the branch stands to main as last fetched."""
        main = self.commit_of(self.main_ref)
        tip = self.branch_tip(branch)
        if not tip or not main:
            return Standing(NO_BRANCH if not tip else UNKNOWN, main)
        if self.is_ancestor(tip, main):
            return Standing(NOTHING_NEW, main)
        if self.is_ancestor(main, tip):
            return Standing(CLEAN, main)
        probe = self.git("merge-tree", "--write-tree", "--name-only", "--no-messages", main, tip, check=False)
        if probe.returncode == 0:
            return Standing(BEHIND, main)
        if probe.returncode == 1:
            return Standing(CONFLICTING, main, tuple(sorted(set(probe.stdout.splitlines()[1:]))))
        return Standing(UNKNOWN, main)


def result_of(outcome):
    """The arm's MergeOutcome as the World answers: a MergeResult, or the
    exception that says why the merge did not land."""
    if outcome.merged:
        return MergeResult(outcome.main)
    if outcome.main_moved:
        raise MainMoved(outcome.reason, outcome.main, outcome.conflicts)
    raise MergeRefused(outcome.reason, outcome.main)
