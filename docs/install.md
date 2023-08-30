# Installation

## For Development

1. Fork the repo into your personal github account

1. At this point, you can start a Github Codespace if you like. If you do this, 
ignore the local installation instrutions below.

1. On your local computer, install [poetry](https://python-poetry.org/) and
consider setting [the configuration such that virtual environments for a given
projects are installed in that project directory](https://python-poetry.org/docs/configuration/#local-configuration).  

1. I suggest that you do this so that poetry builds the virtual environment in 
your project directory:

```python
poetry config virtualenvs.in-project true
```

1. Next, I like working on a fork rather than the actual repository of record.
I set my [git remotes](https://git-scm.com/book/en/v2/Git-Basics-Working-with-Remotes)
so that `origin` points to my fork, and `upstream` points to the 'upstream' repository.

```bash
➜  isocomp git:(develop) ✗ git remote -v 
origin  https://github.com/cmatKhan/isocomp.git (fetch)
origin  https://github.com/cmatKhan/isocomp.git (push)
upstream        https://github.com/collaborativebioinformatics/isocomp.git (fetch)
upstream        https://github.com/collaborativebioinformatics/isocomp.git (push)
```

1. On your machine, `cd` into your local repository, `git checkout` the development 
branch, and make sure it is up-to-date with the upstream (ie the original) repository. 

__NOTE__: if you branch, in general make sure you branch off the `develop` repo, not `main`!  

1. Create a virtual environment in the current directory like so:

```bash
poetry install
```

This will install the virtual environment with the dependencies (and the dependencies' dependencies) 
listed in the [pyproject.toml](./pyproject.toml).

1. If you wish to simultaneously develop both code _and_ docs, then you can 
start the documentation development server:

```bash
mkdocs serve
```