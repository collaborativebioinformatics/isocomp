# Isocomp: comparing high-quality IsoSeq3 isoforms between samples

![](images/logo.png)

## Contributors
1. Yutong Qiu (Carnegie Mellon)
2. Chia Sin	Liew (University of Nebraska-Lincoln)
3. Chase Mateusiak (Washington University)
4. Rupesh Kesharwani (Baylor College of Medicine)
5. Bida	Gu (University of Southern California)
6. Muhammad Sohail Raza (Beijing Institute of Genomics, Chinese Academy of Sciences/China National Center for Bioinformation)
6. Evan	Biederstedt (HMS)

## Introduction
NGS targeted sequencing and WES have become routine for clinical diagnosis of Mendelian disease [CITATION]. Family sequencing (or "trio sequencing") involves sequencing a patient and parents (trio) or other relatives. This improves the diagnostic potential via the interpretation of germline mutations and enables detection of de novo mutations which underlie most Mendelian disorders. 

Transcriptomic profiling has been gaining used over the past several decades. However, this endeavor has been hampered by short-read sequencing, especially for inferring alternative splicing, allelic imbalance, and isoform variation. 

The promise of long-read sequencing has been to overcome the inherent uncertainties of short-reads. 

Something something Isoseq3: https://www.pacb.com/products-and-services/applications/rna-sequencing/

Provides high-quality, polished, assembled full isoforms. With this, we will be able to identify alternatively-spliced isoforms and detect gene fusions. 

Since the advent of HiFi reads, the error rates have plummeted. 

The goal of this project will be to extend the utility of long-read RNAseq for investigating Mendelian diseases between multiple samples. 

And what about gene fusions? We detect these in the stupidest possible way with short-read sequencing, and we think they're cancer-specific. What about the germline?


## Goals

Given high-quality assembled isoforms from 2-3 samples, we want to algorithmically (definitively) characterize the "unique" (i.e. differing) isoforms between samples.


## Description




## Flowchart
![](images/workflow.png)
### To extract sets of unique isoforms
![](images/workflow_part1.png)
### To annotate the unique isoforms
![](images/workflow_part2.png)

## Quick start
### Deployment

Eventually, `pip install isocomp`.  But not yet.

## DEPENDENCIES

python >=3.8

If you're working on `ada`, you'll need to update the old, crusty version of 
python to something more modern and exciting. 

__The easy way__ (untested, but should work):

Install miniconda and create a conda env 
with python 3.9

__The manual method ([source](https://askubuntu.com/a/1424179))__ (tested, works):

```
ssh ... # your username login to ada

mkdir /home/${USER}/.local

# use your favorite text editor. no need to be vim
vim /home/${USER}/.bashrc

# add the following to the end (or where ever)
export PATH=/home/$USER/.local/bin:$PATH

# logout of the current session and log back in
exit
ssh ... (your username, etc)

# Download a more current version of python
wget https://www.python.org/ftp/python/3.9.15/Python-3.9.15.tgz

# unpack
tar xfp Python-3.9.15.tgz 
# remove the tarball
rm Python-3.9.15.tgz 

# cd into the Python package dir, configure and make
cd Python-3.9.15/

./configure --prefix=/home/${USER}/.local --exec_prefix=/home/${USER}/.local --enable-optimizations

make # this takes some time

make altinstall

# the following should point at a python in your /home/$USER/.local/bin dir
which python3.9

# optional, but convenient
ln -s /home/$USER/.local/bin/python3.9 /home/$USER/.local/bin/python

# Download the pip installer
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
# install pip
python3.9 get-pip.py

# confirm that pip is where you think it is
which pip # location should be in your .local

# at this point, you can do:
pip install poetry

# and continue with the development install below

```
### Development

Install [poetry](https://python-poetry.org/) and consider setting [the configuration 
such that virtual environments for a given projects are installed in that project 
directory](https://python-poetry.org/docs/configuration/#local-configuration).  

Next, I like working on a fork rather than the actual repository of record. I set my 
[git remotes](https://git-scm.com/book/en/v2/Git-Basics-Working-with-Remotes) so that `origin` points to 
my fork, and `upstream` points to the 'upstream' repository.

```bash
➜  isocomp git:(develop) ✗ git remote -v 
origin  https://github.com/cmatKhan/isocomp.git (fetch)
origin  https://github.com/cmatKhan/isocomp.git (push)
upstream        https://github.com/collaborativebioinformatics/isocomp.git (fetch)
upstream        https://github.com/collaborativebioinformatics/isocomp.git (push)
```

On your machine, `cd` into your local repository, `git checkout` the development 
branch, and make sure it is up-to-date with the upstream (ie the original) repository. 

__NOTE__: if you branch, in general make sure you branch off the `develop` repo, not `main`!  

Then (assuming poetry is installed already), do:

```bash
$ poetry install
```

This will install the virtual environment with the dependencies (and the dependencies' dependencies) 
listed in the [pyproject.toml](./pyproject.toml).
### <u>Adding dependencies</u>

To add a development dependency (eg, `mkdocs` is not something a user needs), 
use `poetry add -D <dependency>` this is equivalent to `pip install`ing into your 
virtual environment with the added benefit that the dependency is tracked in the 
[pyproject.toml](./pyproject.toml).  

To add a deployment dependency, just omit the `-D` flag.  

### <u>Writing code</u>

Do this first!

```bash
$ pip install -e .
```

[This is an 'editable install'](https://stackoverflow.com/questions/35064426/when-would-the-e-editable-option-be-useful-with-pip-install) 
and means that any change you make in your code is immediately available in your environment. 
__NOTE__: If you happen to see a [Logging Error](https://github.com/pypa/pip/issues/11309) when you run the `pip install -e .` 
command, you can ignore it.  

If you use vscode, [this is a useful 
plugin](https://marketplace.visualstudio.com/items?itemName=jshaptic.autodocs-vscode-support) 
which will automatically generate docstrings for you. [Default docstring format is google](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html), 
which is what the scripts we currently have use. This is an example of what a google formatted 
docstring looks like:

```python
def get_all_windows(gene_df:pd.DataFrame, bp_df:pd.DataFrame) -> pd.DataFrame:
    """From gene boundaries and 100 bp nonzero coverage windows, produce a merged window df

    Args:
        gene_df (pd.DataFrame): one window per gene, > 0.05 avg coverage
        bp_df (pd.DataFrame): one window per 100 bp, > 0.05 avg coverage

    Returns:
        pd.DataFrame: merged windows df
    """
    ...
```

In the function definition, the [type hints](https://docs.python.org/3/library/typing.html) of the arguments (eg `gene_df:pdDataFrame`)
are *not* required, but if you include them, autoDocs will automatically 
generate the data types in the docstring skeleton, also, which is nice. The `-> <datatype>` 
at the end of the function definition is the return data type.

### <u>Tests</u>

Unit tests can be written into the [src/tests](./src/tests) directory. There is 
an example in [src/tests/test_isocomp.py](./src/tests/test_isocomp.py). There are a couple other 
examples of tests -- ie for logging and error handling -- 
[here, too](https://github.com/cmatKhan/lmdemo/blob/main/tests/test_lmdemo.py).

### <u>Build</u>

It's good to intermittently build the package as you go. To do so, use `poetry build` 
which will create a `.whl` and `.tar.gz` (`dist` is already included 
in the [gitignore](./.gitignore)). You can 'distribute' these files to others -- they 
can be installed with `pip` or `conda` -- or use them to install the software outside of 
your current virtual environment.

### <u>Documentation</u>

If you would like to write documentation (ie not docstrings, but long form letters 
to your adoring users), then this can be done in markdown 
[__or jupyter notebooks__](https://pypi.org/project/mkdocs-jupyter/) (already added as a dev dependency) in the 
[docs](./docs) directory. Add the markdown/notebook document to the `nav` section in 
the [mkdocs.yml](./mkdocs.yml) and it will be added to the menu of the documentation 
site. Use `mkdocs serve` locally to see what the documentation looks like. 
`mkdocs build` will build the site in a directory called `site`, which is in the .gitignore already. 
Like `poetry build` it is a good diea to do `mkdocs build` intermittently as you write documentation. 
Eventually, we'll use `mkdocs gh-deploy` to deploy the site to github pages. Maybe if we get fancy, we'll 
set up the github actions to build the package on mac,windows and linux OSes on every push to develop, and rebuild 
the docs and push the package to pypi on every push to `main`.


## Computational Resources / Operation

## Citations
[1] https://www.pacb.com/products-and-services/applications/rna-sequencing/

