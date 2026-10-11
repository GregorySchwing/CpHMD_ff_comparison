"""Refuse a git push the board rules forbid, before it runs.

Claude Code runs this as a PreToolUse hook on the Bash tool, through
cli.py's "guard" subcommand: the tool call arrives as JSON on stdin, exit
code 2 stops it, and stderr goes back into the session as the reason.
The rules are the board rules of the repository that vendors this
(its .claude/rules/board.md): the default branch takes only pull
requests, the board branch takes only the board workflow's writes, and
the one forced push the protocol allows is --force-with-lease on one's
own claude/ branch after a rebase. A refusal is the rule working; the
answer to one is a pull request, a rebase or another branch, never a way
around the hook.

The hook is the first fence and branch protection on GitHub the second,
so anything this module cannot read it lets through rather than wedging
the session: an unparseable command, a broken payload or a git that
answers nothing all pass, and GitHub still refuses what matters.

shlex reads a newline as whitespace, which would glue a line onto the one
before it, so a small pass runs first: it turns each unquoted newline into
the separator the shell reads it as, ends a comment at its newline, and
removes a backslash-newline, the shell's line continuation. It tracks
single quotes, double quotes and the backslash, and nothing else of the
shell: the quoting inside a command substitution is read as the quoting
around it. In particular a here-document body is read as commands like any
other line: a body line that starts a push to main is refused, and the
person rephrases it or puts the text in a file. Cutting bodies out was
tried and dropped: telling a here-document from an arithmetic shift, a
quoted marker or an opener inside a substitution takes more of the shell
than this hook can read, each misreading cut real commands out of sight,
and the one shape it bought is rare and answered by a rephrasing.
"""

import shlex
import subprocess

# Options of git push that take their value as a separate word; the value
# must not be read as a refspec.
VALUE_OPTIONS = {"--receive-pack", "--repo", "-o", "--push-option"}

# Global options of git itself that take a separate value, between "git"
# and the verb, as in "git -C clone push".
GIT_VALUE_OPTIONS = {"-C", "-c", "--git-dir", "--work-tree", "--namespace", "--exec-path"}

# What ends a word for the shell: a # right after one of these starts a
# comment; anywhere else it is part of the word.
WORD_BREAKS = " \t\n;&|()<>"

# The hook's exit codes: 0 lets the tool call run, 2 stops it and sends
# stderr back into the session.
PASS = 0
REFUSE = 2


def plain_shell(command):
    """The command with each unquoted newline made a ; so that shlex, which
    reads a newline as whitespace, ends a command there as the shell does.
    A comment ends at its newline, which then separates; a backslash-newline
    is the line continuation the shell removes; a newline inside quotes
    stays inside its word."""
    out, quote, word_start, i = [], None, True, 0
    while i < len(command):
        ch = command[i]
        if ch == "\\" and quote != "'":
            if command[i + 1:i + 2] != "\n":
                out.append(command[i:i + 2])
            i, word_start = i + 2, False
        elif quote:
            if ch == quote:
                quote = None
            out.append(ch)
            i, word_start = i + 1, False
        elif ch == "#" and word_start:
            end = command.find("\n", i)
            i = len(command) if end == -1 else end
        else:
            if ch in "'\"":
                quote = ch
            out.append(" ; " if ch == "\n" else ch)
            i, word_start = i + 1, ch in WORD_BREAKS
    return "".join(out)


def simple_commands(command):
    """The argument lists of the simple commands inside a shell command.

    shlex with punctuation_chars keeps the shell's operators as their own
    tokens: ;, &&, | and ( ) end one command, a redirection swallows its
    target and the file descriptor digit before it. A command shlex cannot
    read (an unclosed quote) yields nothing, and the caller lets it pass.
    Newlines and comments are read by plain_shell first, so shlex's own
    comment handling, which would read a # inside a word as a comment and
    swallow the rest of the text, is off.
    """
    lexer = shlex.shlex(plain_shell(command), posix=True, punctuation_chars=True)
    lexer.whitespace_split = True
    lexer.commenters = ""
    try:
        tokens = list(lexer)
    except ValueError:
        return []
    commands, current, eat_next = [], [], False
    for token in tokens:
        operator = all(ch in "();<>|&" for ch in token)
        if operator and (">" in token or "<" in token):
            if current and current[-1].isdigit():
                current.pop()
            eat_next = True
        elif operator:
            if current:
                commands.append(current)
            current = []
            eat_next = False
        elif eat_next:
            eat_next = False
        else:
            current.append(token)
    if current:
        commands.append(current)
    return commands


def push_arguments(argv):
    """The arguments after "push" when argv is a git push, else None."""
    if not argv or argv[0].rsplit("/", 1)[-1] != "git":
        return None
    i = 1
    while i < len(argv):
        word = argv[i]
        if word in GIT_VALUE_OPTIONS:
            i += 2
        elif word.startswith("-"):
            i += 1
        elif word == "push":
            return argv[i + 1:]
        else:
            return None
    return None


def refusal(args, current_branch, default_branch):
    """Why this push is refused, or None when it is allowed."""
    force = lease = False
    wholesale = None
    positionals = []
    skip = False
    for arg in args:
        if skip:
            skip = False
        elif arg in ("-f", "--force"):
            force = True
        elif arg == "--force-with-lease" or arg.startswith("--force-with-lease="):
            lease = True
        elif arg in ("--all", "--branches", "--mirror"):
            wholesale = arg
        elif arg in VALUE_OPTIONS:
            skip = True
        elif arg.startswith("-"):
            continue
        else:
            positionals.append(arg)
    if wholesale:
        return (f"a push with {wholesale} moves branches wholesale; "
                "push your own claim branch alone (.claude/rules/board.md)")
    if force:
        return ("--force overwrites history another account may hold; the one forced push "
                "the protocol allows is --force-with-lease on your own claude/ branch "
                "after a rebase (.claude/rules/board.md)")
    destinations = []
    for spec in positionals[1:]:
        name = spec.lstrip("+")
        if ":" in name:
            name = name.split(":", 1)[1]
        name = name.removeprefix("refs/heads/")
        if name in ("", "HEAD"):
            name = current_branch
        if name:
            destinations.append(name)
    if not positionals[1:] and current_branch:
        destinations.append(current_branch)
    for name in destinations:
        if name == default_branch:
            return (f"a push to {name} is refused: the default branch takes only "
                    "pull requests (.claude/rules/board.md)")
        if name == "board":
            return ("a push to board is refused: only the board workflow writes that "
                    "branch (.claude/rules/board.md)")
    if lease:
        strays = [name for name in destinations if not name.startswith("claude/")]
        if strays:
            return (f"--force-with-lease to {', '.join(strays)} is refused: a forced push "
                    "is for your own claude/ branch only (.claude/rules/board.md)")
    return None


def git_line(args, cwd):
    """One line from git, or None when git cannot answer."""
    try:
        done = subprocess.run(["git", *args], cwd=cwd, capture_output=True, timeout=10)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if done.returncode != 0:
        return None
    return done.stdout.decode("utf-8", "replace").strip() or None


def branches(cwd):
    """The current branch (None when detached) and the default branch."""
    current = git_line(["rev-parse", "--abbrev-ref", "HEAD"], cwd)
    if current == "HEAD":
        current = None
    head = git_line(["symbolic-ref", "--quiet", "--short", "refs/remotes/origin/HEAD"], cwd)
    default = head.removeprefix("origin/") if head else "main"
    return current, default


def command_of(payload):
    """The shell command a hook payload carries, or "" when the payload is
    not the shape Claude Code sends (an object whose tool_input is an
    object with a string command): a payload the guard cannot read passes,
    since branch protection still stands behind it."""
    tool_input = payload.get("tool_input") if isinstance(payload, dict) else None
    command = tool_input.get("command") if isinstance(tool_input, dict) else None
    return command if isinstance(command, str) else ""


def where_to_ask(payload, cwd):
    """The directory git is asked in: the session's own, from the payload,
    else the caller's, else here. A payload cwd that is not text is passed
    over, since git cannot be run in it."""
    session = payload.get("cwd") if isinstance(payload, dict) else None
    return (session if isinstance(session, str) and session else None) or cwd or "."


def run(payload, cwd=None):
    """The hook's answer to one tool call: the exit code and the reason,
    which is empty when the call passes. The payload is the hook's JSON
    already decoded; cli.py reads it from stdin and lets a payload that is
    not JSON pass, as this lets one that is not an object pass. Git is
    asked in the payload's own cwd, else in cwd, else here, and only once
    the command holds a push, so a session that never pushes never waits
    on it."""
    command = command_of(payload)
    if "push" not in command:
        return PASS, ""
    pushes = [args for argv in simple_commands(command)
              if (args := push_arguments(argv)) is not None]
    if not pushes:
        return PASS, ""
    current, default = branches(where_to_ask(payload, cwd))
    for args in pushes:
        reason = refusal(args, current, default)
        if reason:
            return REFUSE, reason
    return PASS, ""
