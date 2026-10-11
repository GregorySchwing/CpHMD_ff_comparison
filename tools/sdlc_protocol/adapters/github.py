"""The protocol over GitHub's REST API, through urllib.

Three parts of the package read GitHub, and they read it the same way, so
the calls live once, here. The board collects a Snapshot and posts the
comments board.comments_due names; the team clerk calls the few methods
team.Client lists; and the merge gate asks GitHubWorld, gate.World over one
repository, how a branch stands, whether its required check passed, and
for the merge itself, which goes through GitHub's merge API and only when
the policy says the gate merges.

The token is an argument or GITHUB_TOKEN: a workflow's own token and an
App installation token are both sent as a bearer token. It travels once
per call, in the Authorization header, and never in a message: an error
names the method, the path and the status, and keeps GitHub's own words
apart, on the exception, for a caller that wants to say why.
"""

import base64
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone

from .. import board, claims, gate, labels, team
from ..gate import GateFailed, MainMoved, MergeRefused, MergeResult, NotAllowed

API = "https://api.github.com"
API_VERSION = "2022-11-28"
USER_AGENT = "sdlc-board"
TIMEOUT = 60
PER_PAGE = 100
# GitHub computes mergeability after the first ask; one more ask after a
# pause is enough, and "unknown" is the honest word when it still has not.
MERGEABLE_RETRY_SECONDS = 3
RECENT_MERGED = board.RECENT_MERGED

# GitHub's mergeable_state words that say only that something besides the
# branch stops the merge, a rule, a review or the draft itself: how the
# branch stands to its base is then read from mergeable and a comparison.
OPAQUE_STATES = {"blocked", "draft"}
# A required check run that finished with one of these let the merge go;
# any other finish (failure, timed_out, cancelled, action_required, stale)
# did not.
PASSING = {"success", "neutral", "skipped"}


class GitHubError(team.ClientError):
    """A call GitHub refused or did not answer. The message names the
    method, the path and the status, never the token; status is the HTTP
    status (None when there was no answer) and detail GitHub's own message
    from the answer's body, "" when it gave none."""

    def __init__(self, message, status=None, detail=""):
        super().__init__(message)
        self.status = status
        self.detail = detail


def next_link(headers):
    """The URL of the next page named by a Link header, or None."""
    link = (headers.get("Link") if headers is not None else None) or ""
    for part in link.split(","):
        pieces = part.split(";")
        if len(pieces) >= 2 and pieces[1].strip() == 'rel="next"':
            return pieces[0].strip().strip("<>")
    return None


def error_detail(problem):
    """GitHub's message from a refused call's JSON body, or "" when the
    body is empty, unreadable or not GitHub's shape."""
    try:
        body = json.loads(problem.read() or b"null")
    except (AttributeError, OSError, TypeError, ValueError):
        return ""
    return str(body.get("message") or "") if isinstance(body, dict) else ""


def quoted(name):
    """A branch or label name made safe for a URL path, slashes kept, since
    GitHub reads a ref's slashes as part of its name."""
    return urllib.parse.quote(name, safe="/")


class GitHub:
    """The REST calls the package makes, through urllib, and the methods
    team.Client lists.

    The opener and the sleep are arguments so a test can hand in fakes; the
    token is the argument, else GITHUB_TOKEN, and with neither there is no
    client, since every call would be refused.
    """

    def __init__(self, token, repo, opener=None, sleep=time.sleep):
        self.token = token or os.environ.get("GITHUB_TOKEN")
        if not self.token:
            raise GitHubError("no token: pass one, or set GITHUB_TOKEN")
        self.repo = repo
        self.opener = opener or urllib.request.build_opener()
        self.sleep = sleep

    def headers(self):
        return {
            "Authorization": f"Bearer {self.token}",
            "Accept": "application/vnd.github+json",
            "X-GitHub-Api-Version": API_VERSION,
            "User-Agent": USER_AGENT,
        }

    def url(self, path, params=None):
        return API + path + ("?" + urllib.parse.urlencode(params) if params else "")

    def request(self, method, url, body=None, missing_ok=False):
        """One call; the decoded JSON and the answer's headers. A 404 is
        (None, headers) when missing_ok says absence is an answer."""
        data = json.dumps(body).encode("utf-8") if body is not None else None
        headers = self.headers()
        if data is not None:
            headers["Content-Type"] = "application/json"
        call = urllib.request.Request(url, data=data, method=method, headers=headers)
        where = f"{method} {urllib.parse.urlsplit(url).path}"
        try:
            with self.opener.open(call, timeout=TIMEOUT) as answer:
                raw = answer.read()
                return (json.loads(raw) if raw else None), answer.headers
        except urllib.error.HTTPError as problem:
            if problem.code == 404 and missing_ok:
                return None, problem.headers
            raise GitHubError(f"{where}: HTTP {problem.code}", problem.code, error_detail(problem)) from None
        except OSError as problem:
            raise GitHubError(f"{where}: {getattr(problem, 'reason', problem)}") from None
        except ValueError:
            raise GitHubError(f"{where}: the answer was not JSON") from None

    def get(self, path, params=None, missing_ok=False):
        data, _ = self.request("GET", self.url(path, params), missing_ok=missing_ok)
        return data

    def pages(self, path, params=None):
        """A list endpoint page by page, following the Link header, so a
        caller that has seen enough can stop asking."""
        params = dict(params or {})
        params["per_page"] = PER_PAGE
        url = self.url(path, params)
        while url:
            data, headers = self.request("GET", url)
            yield data or []
            url = next_link(headers)

    def get_all(self, path, params=None):
        """Every item of a list endpoint."""
        return [item for page in self.pages(path, params) for item in page]

    def post(self, path, body):
        data, _ = self.request("POST", self.url(path), body=body)
        return data

    def put(self, path, body):
        data, _ = self.request("PUT", self.url(path), body=body)
        return data

    def patch(self, path, body):
        data, _ = self.request("PATCH", self.url(path), body=body)
        return data

    def delete(self, path, missing_ok=False):
        data, _ = self.request("DELETE", self.url(path), missing_ok=missing_ok)
        return data

    # --- team.Client ---------------------------------------------------------

    def issue(self, number):
        return self.get(f"/repos/{self.repo}/issues/{number}")

    def open_issues(self):
        """Every open issue, with the pull requests GitHub mixes into the
        list left in, as the API gives them; the clerk and the board each
        leave them out."""
        return self.get_all(f"/repos/{self.repo}/issues", {"state": "open"})

    def branch_names(self):
        return [branch["name"] for branch in self.get_all(f"/repos/{self.repo}/branches")]

    def open_pull(self, branch):
        """The open pull request from a branch of this repository, in the
        detail that holds mergeable_state and the head's sha, or None."""
        owner = self.repo.split("/")[0]
        found = self.get(f"/repos/{self.repo}/pulls", {"state": "open", "head": f"{owner}:{branch}"}) or []
        if not found:
            return None
        return self.get(f"/repos/{self.repo}/pulls/{found[0]['number']}")

    def check_runs(self, sha):
        """Every check run on one commit, page by page: the endpoint wraps
        each page's runs in an object, and a busy commit has more than one
        page of them."""
        runs = []
        for page in self.pages(f"/repos/{self.repo}/commits/{sha}/check-runs"):
            if isinstance(page, dict):
                runs.extend(page.get("check_runs") or [])
        return runs

    def create_label(self, name, colour):
        """Make the label; GitHub answers 422 for one that exists, which is
        the label being there, not a failure."""
        try:
            self.post(f"/repos/{self.repo}/labels", {"name": name, "color": colour})
        except GitHubError as problem:
            if problem.status != 422:
                raise

    def create_issue(self, title, body, labels):
        return self.post(f"/repos/{self.repo}/issues", {"title": title, "body": body, "labels": labels})

    def add_sub_issue(self, parent, child_id):
        self.post(f"/repos/{self.repo}/issues/{parent}/sub_issues", {"sub_issue_id": child_id})

    def set_body(self, number, body):
        self.patch(f"/repos/{self.repo}/issues/{number}", {"body": body})

    def remove_label(self, number, name):
        self.delete(f"/repos/{self.repo}/issues/{number}/labels/{urllib.parse.quote(name, safe='')}",
                    missing_ok=True)


# --- the board ----------------------------------------------------------------

def settled(api, detail):
    """A pull request's detail with mergeable_state settled: asked once
    more after a pause when GitHub had not computed it, and "unknown" when
    it still has not."""
    if detail.get("mergeable") is None:
        api.sleep(MERGEABLE_RETRY_SECONDS)
        detail = api.get(f"/repos/{api.repo}/pulls/{detail['number']}")
    if detail.get("mergeable") is None:
        detail["mergeable_state"] = "unknown"
    return detail


def pull_with_detail(api, number):
    """One pull request with its mergeable_state settled and its file paths."""
    path = f"/repos/{api.repo}/pulls/{number}"
    detail = settled(api, api.get(path))
    detail["files"] = [f["filename"] for f in api.get_all(path + "/files")]
    return detail


def claim_text(api, branch):
    """The claim file on a claude/<issue> branch, or None when there is none yet."""
    name = board.branch_name(branch)
    found = api.get(
        f"/repos/{api.repo}/contents/{claims.CLAIM_DIR}/{claims.key_of(name)}.md",
        {"ref": name}, missing_ok=True,
    )
    if found is None:
        return None
    return base64.b64decode(found.get("content") or "").decode("utf-8", "replace")


def recent_merged(api):
    """The pull requests merged last, newest first, at most RECENT_MERGED,
    from the first page of closed pull requests by last update (a hundred).
    The list depends on the repository alone, so two runs in a row render
    the same one: a merge fires two runs (the pull request closing and the
    push to the default branch), and a list measured from the previous
    board's commit lost, in the second, the merge the first had reported."""
    for page in api.pages(f"/repos/{api.repo}/pulls", {"state": "closed", "sort": "updated", "direction": "desc"}):
        merged = [p for p in page if p.get("merged_at")]
        merged.sort(key=lambda p: p["merged_at"], reverse=True)
        return merged[:RECENT_MERGED]
    return []


def claim_issues(api, branches, open_issues):
    """The issue behind each claim branch, by its key: the state and the
    label names. An open issue comes from the list already read; any other
    takes one call, so the board can say that a branch outlived its issue.
    A branch whose key is not a number has none."""
    by_number = {str(i["number"]): i for i in open_issues}
    found = {}
    for branch in board.claim_branches(branches):
        key = claims.key_of(branch)
        if not key.isdigit():
            continue
        issue = by_number.get(key) or api.get(f"/repos/{api.repo}/issues/{key}", missing_ok=True)
        if issue is not None:
            found[key] = {"state": issue.get("state", ""), "labels": labels.label_names(issue)}
    return found


def collect(api, sha, default_branch):
    """Everything the board renders, in the shape a --dry-run fixture has,
    so a live run's answers can be saved and rendered again offline."""
    repo = api.repo
    issues = [i for i in api.get_all(f"/repos/{repo}/issues", {"state": "open"}) if "pull_request" not in i]
    pulls = [pull_with_detail(api, p["number"]) for p in api.get_all(f"/repos/{repo}/pulls", {"state": "open"})]
    branches = board.claim_branches(api.get_all(f"/repos/{repo}/branches"))
    return {
        "sha": sha,
        "default_branch": default_branch,
        "issues": issues,
        "pulls": pulls,
        "branches": branches,
        "claims": {claims.key_of(b): claim_text(api, b) for b in branches},
        "claim_issues": claim_issues(api, branches, issues),
        "merged": recent_merged(api),
    }


def snapshot(api, sha, default_branch):
    """The board's Snapshot of the repository as GitHub knows it now."""
    return board.Snapshot.from_dict(collect(api, sha, default_branch))


def comment_once(api, number, marker, text):
    """Post the comment unless one on the pull request carries the marker."""
    comments = api.get_all(f"/repos/{api.repo}/issues/{number}/comments")
    if any(marker in (c.get("body") or "") for c in comments):
        return False
    api.post(f"/repos/{api.repo}/issues/{number}/comments", {"body": marker + "\n" + text})
    return True


def post_comments(api, due):
    """Post each comment board.comments_due names, once each; the lines say
    which were posted, one each, for the workflow's log."""
    posted = []
    for number, marker, text in due:
        kind = "idle" if marker == board.IDLE_MARKER else "main-moved"
        if comment_once(api, number, marker, text):
            posted.append(f"PR #{number}: {kind} comment posted")
    return posted


# --- the gate's world ---------------------------------------------------------

def standing_word(pull, behind_by=None):
    """How a pull request's branch stands to its base, in the gate's four
    words, from the pull's mergeable_state as the board reads it. Where
    GitHub's word says only that something else blocks the merge (a rule,
    a review, a draft), mergeable tells a conflict apart, and behind_by,
    the commits of the base the head lacks, tells behind from clean; with
    neither known the answer is unknown."""
    if pull.get("mergeable_state") not in OPAQUE_STATES:
        return board.pr_state(pull)
    if pull.get("mergeable") is False:
        return gate.CONFLICTING
    if pull.get("mergeable") is None or behind_by is None:
        return gate.UNKNOWN
    return gate.BEHIND if behind_by > 0 else gate.CLEAN


def check_word(runs, name):
    """The required check's state on a commit, in the gate's words: none
    when no run has its name, red when a run of it finished without
    passing, pending while one is still running, green when every run of
    it passed. Red wins over pending, since a failed run already means
    the branch needs work."""
    mine = [run for run in runs if run.get("name") == name]
    if not mine:
        return gate.NO_CHECK
    if any(run.get("status") == "completed" and run.get("conclusion") not in PASSING for run in mine):
        return gate.RED
    if any(run.get("status") != "completed" for run in mine):
        return gate.PENDING
    return gate.GREEN


class GitHubWorld:
    """gate.World over one repository on GitHub.

    The policy says which check is required, who merges and which branch
    is the default one, so the world reads the same branch the gate's
    answers name; default_branch overrides it for a caller that knows
    better. The clock is an argument so a test can fix it.
    """

    def __init__(self, api, policy, default_branch=None, clock=None):
        self.api = api
        self.policy = policy
        self.default_branch = default_branch or policy.default_branch
        self.clock = clock or (lambda: datetime.now(timezone.utc))

    def main_sha(self):
        return self.branch_sha(self.default_branch)

    def branch_sha(self, branch):
        """The branch's commit, or "" when the repository has no such branch."""
        found = self.api.get(f"/repos/{self.api.repo}/git/ref/heads/{quoted(branch)}", missing_ok=True)
        return ((found or {}).get("object") or {}).get("sha") or ""

    def pull_of(self, branch):
        """The branch's open pull request, settled, or None."""
        found = self.api.open_pull(branch)
        return settled(self.api, found) if found else None

    def behind_by(self, pull):
        """How many commits of the pull request's base its head lacks, from
        GitHub's comparison of the two; None when GitHub cannot say."""
        base = (pull.get("base") or {}).get("ref") or self.default_branch
        head = (pull.get("head") or {}).get("sha")
        if not head:
            return None
        found = self.api.get(f"/repos/{self.api.repo}/compare/{quoted(base)}...{head}", missing_ok=True)
        return (found or {}).get("behind_by")

    def standing(self, branch):
        """Unknown for a branch with no open pull request: GitHub computes
        mergeability only for one, and the gate asks only about those."""
        pull = self.pull_of(branch)
        if pull is None:
            return gate.UNKNOWN
        if pull.get("mergeable_state") in OPAQUE_STATES and pull.get("mergeable"):
            return standing_word(pull, self.behind_by(pull))
        return standing_word(pull)

    def check(self, branch):
        sha = self.branch_sha(branch)
        return check_word(self.api.check_runs(sha), self.policy.required_check) if sha else gate.NO_CHECK

    def merge(self, branch):
        """Merge the branch's pull request through GitHub's merge API, with
        the head the gate looked at, so a push that lands meanwhile is
        refused rather than merged unseen. Only where the policy says the
        gate merges: elsewhere the owner or auto-merge does, and a merge
        here would go behind their back."""
        if self.policy.who_merges != gate.GATE:
            raise NotAllowed(f"{gate.MOVERS[self.policy.who_merges]} merges here, not the gate; "
                             f"{branch} waits for it")
        pull = self.api.open_pull(branch)
        if pull is None:
            raise MergeRefused(f"{branch} has no open pull request to merge; open one first")
        try:
            done = self.api.put(f"/repos/{self.api.repo}/pulls/{pull['number']}/merge",
                                {"sha": pull["head"]["sha"], "merge_method": "merge"})
        except GitHubError as problem:
            raise self.refusal(branch, pull, problem) from None
        return MergeResult(done["sha"])

    def refusal(self, branch, pull, problem):
        """What a refused merge means for the gate. GitHub answers 405 for
        a pull request it will not merge and 409 when the head moved: main
        moved when the branch, read again, no longer fits the default
        branch, and otherwise the author has something to fix (a review, a
        rule, a push to look at). Any other answer is the gate failing."""
        said = f": {problem.detail}" if problem.detail else ""
        if problem.status not in (405, 409):
            return GateFailed(f"GitHub did not merge #{pull['number']}: {problem}{said}")
        if self.standing(branch) in (gate.BEHIND, gate.CONFLICTING):
            return MainMoved(f"{branch} no longer merges cleanly into {self.default_branch}{said}")
        return MergeRefused(f"GitHub would not merge #{pull['number']}{said or ' as it stands'}")

    def now(self):
        return self.clock()
