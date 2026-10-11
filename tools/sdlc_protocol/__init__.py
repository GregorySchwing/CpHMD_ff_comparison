"""The claim protocol for pull-request-gated Claude agents, as one package.

Accounts and agents share a repository by claiming an issue as a branch,
working in the open on a draft pull request and landing only what merges
cleanly with the default branch as it is now. The rules for that live here
once, so the places that run them cannot drift apart. Four repositories
run a copy of it, which the sdlc-protocol repository's sync keeps up to
date by pull request (ecosystem.json there says who takes which modules):
the template, for its board, its protocol check, its team clerk and its
push guard; coordbench's sdlc arm, for the merge gate and the bare-git
world; the plugin marketplace, for the push guard; and the GitHub App, for
the protocol check. The rules are plain functions and dataclasses over
what a world reports, only the modules under adapters/ call GitHub or
change a git repository, and nothing needs more than Python 3.12's
standard library. The protocol in prose is spec/PROTOCOL.md in the
sdlc-protocol repository.
"""

__version__ = "0.2.0"
