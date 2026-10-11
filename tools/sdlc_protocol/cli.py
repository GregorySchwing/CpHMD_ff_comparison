"""The command line: python3 -m sdlc_protocol <board|check|team|guard|verify-vendored>.

Each subcommand is a thin shell over the modules and takes exactly the
arguments the template's tool for the same job took, so a workflow line
changes only from "python3 tools/board.py ..." to "python3 -m
sdlc_protocol board ...":

  board  render STATE.md from GitHub (--default-branch, --sha, --out,
         --render-only) and comment on the pull requests that wait, or
         render a JSON snapshot with --dry-run FIXTURE and no network
  check  hold the pull request in GITHUB_EVENT_PATH to the claim protocol
  team   [--repo OWNER/NAME] seed, route or next, team mode's clerk
  guard  the PreToolUse hook: the tool call's JSON on stdin, and exit 2
         with the reason on stderr to refuse a push
  verify-vendored VENDOR_DIR
         hold a vendored copy, VENDOR_DIR/sdlc_protocol, to the record
         VENDOR_DIR/VENDORED beside it (vendoring.py, which a copy also
         runs as a script of its own)

The exit codes are the old tools' too: 0 when done, 1 for a failure told
in one line on stderr (or, for check, the problems on stdout), 2 for a
push the guard refuses or arguments that cannot be read.

Only board and team import the GitHub adapter, and only when they run:
the guard answers before every Bash call of a session, and the start-up
time of urllib is no use to it.
"""

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

from . import board, claims, guard, protocol_check, team, vendoring

PROG = "python3 -m sdlc_protocol"
USAGE = f"usage: {PROG} {{board,check,team,guard,verify-vendored}} [arguments]; {PROG} <command> -h says more"


class CommandError(Exception):
    """What a subcommand cannot start without, in one line, never a token."""


def github_adapter():
    """The GitHub adapter, imported when a subcommand first calls GitHub."""
    from .adapters import github
    return github


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def environment(given):
    """The environment a subcommand reads: the one a test hands in, else
    the process's own."""
    return os.environ if given is None else given


# --- board --------------------------------------------------------------------

def board_arguments(argv):
    parser = argparse.ArgumentParser(
        prog=f"{PROG} board",
        description="Render the board (STATE.md) and comment on waiting pull requests.")
    parser.add_argument("--dry-run", metavar="FIXTURE", help="render this JSON snapshot instead of asking GitHub")
    parser.add_argument("--render-only", action="store_true", help="write the board and post no comment")
    parser.add_argument("--out", help=f"where to write the board (default: {board.STATE_PATH} in the current "
                                      "directory; stdout with --dry-run)")
    parser.add_argument("--default-branch", default="main", help="the default branch's name (default: main)")
    parser.add_argument("--sha", help="the default branch's commit the board describes (default: GITHUB_SHA)")
    return parser.parse_args(argv)


def board_dry_run(args):
    """The fixture's board, as it stands: the fixture names its own default
    branch, as the template's dry run read it."""
    text = board.render(load_json(args.dry_run))
    if args.out:
        board.write_if_changed(Path(args.out), text)
    else:
        sys.stdout.write(text)
    return 0


def board_live_run(args, environ):
    token, repo = environ.get("GITHUB_TOKEN"), environ.get("GITHUB_REPOSITORY")
    if not token or not repo:
        raise CommandError("GITHUB_TOKEN and GITHUB_REPOSITORY must be set")
    github = github_adapter()
    api = github.GitHub(token, repo)
    snap = github.snapshot(api, args.sha or environ.get("GITHUB_SHA", ""), args.default_branch)
    out = Path(args.out) if args.out else board.STATE_PATH
    changed = board.write_if_changed(out, board.render(snap))
    print(f"{out.as_posix()}: {'rewritten' if changed else 'unchanged'}")
    if not args.render_only:
        due = board.comments_due(snap, environ.get("GITHUB_EVENT_NAME", ""), environ.get("GITHUB_REF_NAME", ""),
                                 datetime.now(timezone.utc))
        for line in github.post_comments(api, due):
            print(line)
    return 0


def run_board(argv, environ=None):
    args = board_arguments(argv)
    try:
        return board_dry_run(args) if args.dry_run else board_live_run(args, environment(environ))
    except (CommandError, team.ClientError) as problem:
        print(f"board: {problem}", file=sys.stderr)
        return 1


# --- check --------------------------------------------------------------------

def read_event(environ):
    """The event GitHub hands the run, and what to say when there is none
    to check: no event, or one that is not a JSON object."""
    where = environ.get("GITHUB_EVENT_PATH")
    if not where or not Path(where).is_file():
        return None, "no pull request event to check"
    try:
        event = load_json(where)
    except ValueError:
        event = None
    if not isinstance(event, dict):
        return None, "the event is not readable; nothing checked"
    return event, ""


def run_check(argv, environ=None, root="."):
    """The protocol check: a run outside a pull request event passes, so
    the check never breaks a push build; the claim file is looked for in
    the checkout, root."""
    argparse.ArgumentParser(
        prog=f"{PROG} check",
        description="Hold the claude/<issue> pull request in GITHUB_EVENT_PATH to the claim protocol.",
    ).parse_args(argv)
    event, why_not = read_event(environment(environ))
    if event is None:
        print(why_not)
        return 0
    problems = protocol_check.problems_for(
        event.get("pull_request") or {}, lambda number: (Path(root) / claims.claim_path(number)).is_file())
    for problem in problems:
        print(problem)
    if problems:
        return 1
    print("protocol ok")
    return 0


# --- team ---------------------------------------------------------------------

def team_arguments(argv, environ=None):
    environ = environment(environ)
    parser = argparse.ArgumentParser(prog=f"{PROG} team", description="Team mode's clerical steps (docs/agents/TEAM.md).")
    parser.add_argument("--repo", default=environ.get("GITHUB_REPOSITORY"))
    sub = parser.add_subparsers(dest="command", required=True)
    seeding = sub.add_parser("seed")
    seeding.add_argument("--epic", type=int, required=True)
    seeding.add_argument("--plan", required=True, help="the lead's plan as a JSON file")
    seeding.add_argument("--team", default=environ.get("CLAUDE_TEAM", ""))
    routing = sub.add_parser("route")
    routing.add_argument("--event", required=True)
    routing.add_argument("--issue", type=int, required=True)
    routing.add_argument("--label", default=None)
    nexting = sub.add_parser("next")
    nexting.add_argument("--agent", required=True)
    nexting.add_argument("--team", default=environ.get("CLAUDE_TEAM", ""))
    nexting.add_argument("--label", default=None,
                         help="the member's label, made before a take so the session can add it")
    return parser.parse_args(argv)


def team_answer(client, args):
    """The one line the workflow reads, from the clerk's rules."""
    if args.command == "seed":
        return team.seed(client, args.epic, load_json(args.plan), team.team_members(args.team))
    if args.command == "route":
        return team.route(client, args.event, args.issue, args.label)
    return team.next_move(client, args.agent, team.team_members(args.team), args.label)


def run_team(argv, environ=None):
    environ = environment(environ)
    args = team_arguments(argv, environ)
    token = environ.get("GITHUB_TOKEN") or environ.get("GH_TOKEN")
    if not token or not args.repo:
        print("GITHUB_TOKEN and --repo (or GITHUB_REPOSITORY) are needed", file=sys.stderr)
        return 1
    try:
        print(team_answer(github_adapter().GitHub(token, args.repo), args))
    except (team.ClientError, team.TeamError, ValueError, KeyError) as problem:
        print(f"team: {problem}", file=sys.stderr)
        return 1
    return 0


# --- guard --------------------------------------------------------------------

GUARD_HELP = (f"usage: {PROG} guard < payload.json\n\n"
              "The PreToolUse hook on Bash: reads the tool call's JSON on stdin and exits 2, "
              "with the reason on stderr, to refuse a push the board rules forbid.")


def run_guard(argv, stdin=None):
    """The push guard as Claude Code's hook runs it. Its arguments are not
    parsed: exit 2 stops the tool call, so an argument error would refuse
    every Bash call of the session, and the guard fails open instead, as
    it does for a payload that is not JSON. Help is the one argument it
    answers, and it lets the call pass too."""
    if "-h" in argv or "--help" in argv:
        print(GUARD_HELP)
        return guard.PASS
    try:
        payload = json.load(sys.stdin if stdin is None else stdin)
    except ValueError:
        return guard.PASS
    code, reason = guard.run(payload)
    if reason:
        print(reason, file=sys.stderr)
    return code


# --- verify-vendored ----------------------------------------------------------

def run_verify_vendored(argv):
    """vendoring.py's check under the package's own command line, for a
    consumer that would rather name the module than the file."""
    parser = argparse.ArgumentParser(
        prog=f"{PROG} verify-vendored",
        description="Hold VENDOR_DIR/sdlc_protocol to the record VENDOR_DIR/VENDORED: exit 1 on any difference.")
    parser.add_argument("vendor_dir", metavar="VENDOR_DIR", help="the directory that holds sdlc_protocol/ and VENDORED")
    return vendoring.verify(parser.parse_args(argv).vendor_dir)


COMMANDS = {"board": run_board, "check": run_check, "team": run_team, "guard": run_guard,
            "verify-vendored": run_verify_vendored}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv[:1] in (["-h"], ["--help"]):
        print(USAGE)
        return 0
    if not argv or argv[0] not in COMMANDS:
        print(USAGE, file=sys.stderr)
        return 2
    return COMMANDS[argv[0]](argv[1:])
