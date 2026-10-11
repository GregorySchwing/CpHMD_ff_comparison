"""The record a vendored copy of this package carries, and its check.

A repository that runs the protocol from a copy in its own tree keeps the
file VENDORED beside the copy. Its first line names the release and the
commit the copy was made from, "sdlc-protocol <version> <commit>", and
each line after it the sha256 of one copied file and the file's path
inside the copy, in the layout sha256sum writes, sorted by path. The
library's scripts/vendor.py writes the record with every copy it makes,
and the consumer's own check runs this file as a script:

  python3 <vendor_dir>/sdlc_protocol/vendoring.py verify <vendor_dir>

so a hand edit of the copy, a file lost or a file added beside it fails
there, offline, on the pull request that made it. A change to the copy is
made in sdlc-protocol and reaches the consumer by the sync's pull request.

That is why this module imports nothing from the package and nothing
outside the standard library: it runs inside a copy that may hold only a
few of the package's files, with its own directory first on sys.path.
"""

import hashlib
import re
import sys
from pathlib import Path, PurePosixPath

PACKAGE = "sdlc_protocol"
RECORD = "VENDORED"
CACHE_DIRECTORY = "__pycache__"
CACHE_SUFFIX = ".pyc"

FIRST_LINE = re.compile(r"sdlc-protocol (\S+) ([0-9a-f]{40}(?:-dirty)?|unknown)")
DIGEST_LINE = re.compile(r"([0-9a-f]{64})  (\S(?:.*\S)?)")

REVENDOR = "re-vendor with python3 <sdlc-protocol>/scripts/vendor.py <this repository>"
LET_THE_SYNC_BRING_IT = "change it in sdlc-protocol and let the sync bring it"
USAGE = "usage: python3 <vendor_dir>/sdlc_protocol/vendoring.py verify <vendor_dir>"


# --- the record -------------------------------------------------------------------

def record_text(version, commit, digests):
    """The text of VENDORED for a copy of version at commit whose files
    have these digests, a map from path inside the copy to sha256 hex."""
    first = f"sdlc-protocol {version} {commit}"
    if not FIRST_LINE.fullmatch(first):
        raise ValueError(f"cannot record version {version!r} at commit {commit!r}")
    lines = [first] + [f"{digests[path]}  {path}" for path in sorted(digests)]
    return "\n".join(lines) + "\n"


def parse_record(text):
    """The version, the commit and the digests a VENDORED text records.
    A record that cannot be read raises ValueError saying what is wrong,
    and the old record of one line, which had no digests, is one of them."""
    lines = text.splitlines()
    head = FIRST_LINE.fullmatch(lines[0]) if lines else None
    if not head:
        raise ValueError(f"the first line is not 'sdlc-protocol <version> <commit>'; {REVENDOR}")
    if len(lines) == 1:
        raise ValueError(f"it is the old record of one line, with no file hashes; {REVENDOR}")
    return head.group(1), head.group(2), digests_in(lines[1:])


def digests_in(lines):
    """The digests on the record's lines after the first, which are
    numbered from 2 in what a problem says."""
    digests = {}
    for number, line in enumerate(lines, start=2):
        found = DIGEST_LINE.fullmatch(line)
        if not found or not is_inside_the_copy(found.group(2)):
            raise ValueError(f"line {number} is not '<sha256>  <path inside the copy>': {line!r}; {REVENDOR}")
        if found.group(2) in digests:
            raise ValueError(f"line {number} lists {found.group(2)} a second time; {REVENDOR}")
        digests[found.group(2)] = found.group(1)
    return digests


def is_inside_the_copy(path):
    """Whether a recorded path names a file inside the copy, as a relative
    POSIX path that never climbs out of it."""
    parts = PurePosixPath(path).parts
    return bool(parts) and not path.startswith("/") and ".." not in parts and "\\" not in path


# --- the copy ---------------------------------------------------------------------

def is_cache(relative):
    """Whether a path inside the copy is one of Python's caches, which
    running the code writes there and which are never copied or recorded."""
    return CACHE_DIRECTORY in relative.parts or relative.suffix == CACHE_SUFFIX


def copied_files(copy_dir):
    """The files of a copy as sorted relative POSIX paths, caches left out."""
    root = Path(copy_dir)
    found = []
    for path in root.rglob("*"):
        relative = path.relative_to(root)
        if path.is_file() and not is_cache(relative):
            found.append(relative.as_posix())
    return sorted(found)


def digest_of(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digests_of(copy_dir):
    """The sha256 of every file of a copy, by its path inside the copy."""
    root = Path(copy_dir)
    return {name: digest_of(root / name) for name in copied_files(root)}


# --- the check --------------------------------------------------------------------

def problems(vendor_dir):
    """What is wrong with the copy <vendor_dir>/sdlc_protocol against the
    record <vendor_dir>/VENDORED, one plain sentence each; [] when the copy
    is exactly what the record says."""
    vendor = Path(vendor_dir)
    record, copy = vendor / RECORD, vendor / PACKAGE
    if not record.is_file():
        return [f"{record} is missing; {REVENDOR}"]
    try:
        _, _, recorded = parse_record(record.read_text(encoding="ascii"))
    except ValueError as problem:
        return [f"{record} cannot be read: {problem}"]
    if not copy.is_dir():
        return [f"{copy} is missing; {REVENDOR}"]
    return file_problems(copy, recorded, digests_of(copy))


def file_problems(copy, recorded, actual):
    """Each file the record and the copy disagree on, in path order."""
    found = []
    for name in sorted(set(recorded) | set(actual)):
        where = copy / name
        if name not in actual:
            found.append(f"{where} is recorded in {RECORD} but missing; {REVENDOR}")
        elif name not in recorded:
            found.append(f"{where} is not recorded in {RECORD}: remove it, or {LET_THE_SYNC_BRING_IT}")
        elif actual[name] != recorded[name]:
            found.append(f"{where} was edited by hand: {LET_THE_SYNC_BRING_IT}")
    return found


def short_commit(commit):
    """A commit as a person reads it: seven hex digits, keeping -dirty."""
    if commit == "unknown":
        return commit
    sha, _, dirty = commit.partition("-")
    return sha[:7] + (f"-{dirty}" if dirty else "")


def verify(vendor_dir, out=None):
    """Print each problem and return 1, or print the one line that says the
    copy matches its record and return 0."""
    out = sys.stdout if out is None else out
    found = problems(vendor_dir)
    for problem in found:
        print(problem, file=out)
    if found:
        return 1
    version, commit, digests = parse_record((Path(vendor_dir) / RECORD).read_text(encoding="ascii"))
    print(f"vendored sdlc-protocol {version} at {short_commit(commit)}: {len(digests)} files match {RECORD}",
          file=out)
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv[:1] in (["-h"], ["--help"]):
        print(USAGE)
        return 0
    if len(argv) != 2 or argv[0] != "verify":
        print(USAGE, file=sys.stderr)
        return 2
    return verify(argv[1])


if __name__ == "__main__":
    sys.exit(main())
