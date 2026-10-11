"""The worlds the merge gate decides over, one module each.

sdlc_protocol.gate decides; a world is where its decisions land and where
it reads how each branch stands. bare_git is a bare repository reached by
git, as the coordbench sdlc arm's shared remote is; memory is a world held
in memory, for tests and for simulations that set each branch's standing
and check by hand; github is a repository on GitHub, through the REST API.
Nothing is imported here, so a caller loads only the world it uses.
"""
