"""python3 -m sdlc_protocol: the package's command line, in cli.py.

A vendored copy can also be run by its path, as python3
tools/sdlc_protocol board ...: Python then runs this file with the
package's own directory first on sys.path and no package around it, so
the directory above takes that place and the package is imported by name.
"""

import sys

if not __package__:
    from pathlib import Path
    sys.path[0] = str(Path(__file__).resolve().parent.parent)
    from sdlc_protocol.cli import main
else:
    from .cli import main

if __name__ == "__main__":
    sys.exit(main())
