# Installation

## For Development

Clone the repo and make a virtual environment called `.venv` 
in the repo. Activate the environment and do the following:

1. Install [poetry](https://python-poetry.org/) and consider setting [the configuration 
such that virtual environments for a given projects are installed in that project 
directory](https://python-poetry.org/docs/configuration/#local-configuration).  

2. Next, I like working on a fork rather than the actual repository of record. I set my 
[git remotes](https://git-scm.com/book/en/v2/Git-Basics-Working-with-Remotes) so that `origin` points to 
my fork, and `upstream` points to the 'upstream' repository.

```bash
➜  isocomp git:(develop) ✗ git remote -v 
origin  https://github.com/cmatKhan/isocomp.git (fetch)
origin  https://github.com/cmatKhan/isocomp.git (push)
upstream        https://github.com/collaborativebioinformatics/isocomp.git (fetch)
upstream        https://github.com/collaborativebioinformatics/isocomp.git (push)
```

3. On your machine, `cd` into your local repository, `git checkout` the development 
branch, and make sure it is up-to-date with the upstream (ie the original) repository. 

__NOTE__: if you branch, in general make sure you branch off the `develop` repo, not `main`!  

4. Create a virtual environment in the current directory 
like so:

```bash
python -m venv .venv
```

5. Activate the venv

```bash
source .venv/bin/activate
```

6. Then (assuming poetry is installed already), do:

```bash
$ poetry install
```

This will install the virtual environment with the dependencies (and the dependencies' dependencies) 
listed in the [pyproject.toml](./pyproject.toml).

7. Finally, serve the docs

```bash
mkdocs serve
```