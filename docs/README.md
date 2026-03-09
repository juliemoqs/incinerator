[![pytest](https://github.com/juliemoqs/incinerator/actions/workflows/pytest.yaml/badge.svg)](https://github.com/juliemoqs/incinerator/actions/workflows/pytest.yaml)
[![mypy](https://github.com/juliemoqs/incinerator/actions/workflows/mypy.yaml/badge.svg)](https://github.com/juliemoqs/incinerator/actions/workflows/mypy.yaml/)
[![ruff](https://github.com/juliemoqs/incinerator/actions/workflows/ruff.yaml/badge.svg)](https://github.com/juliemoqs/incinerator/actions/workflows/ruff.yaml)
[![PyPI](https://img.shields.io/pypi/v/incinerator.svg)](https://pypi.python.org/pypi/incinerator)
[![Generic badge](https://img.shields.io/badge/documentation-live-blue.svg)](https://juliemoqs.github.io/incinerator/)

# INCInERATOr
**I**ndicator of **N**onaligned **C**entroids **IN** **E**xoplanet **R**eliability **A**nalysis of **T**ransit **O**bse**R**vations

## How to install `incinerator`
The easiest way to install `incinerator` and all of its dependencies is with `pip`. 

We recommend you do this installation in a new [virtual environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

```bash
# Create new environment called incinerator with Python 3.12
conda create -n incinerator python=3.12

# Activate environment
conda activate incinerator

# Install incinerator
pip install incinerator --upgrade
```

## How to use `incinerator`

`incinerator` localizes transit signals at the pixel level by mapping transit depth across a target pixel file and fitting a model to determine the signal's ccd position. 

Check out our full documentation, which includes a tutorial, [here](https://juliemoqs.github.io/incinerator/) to get started!
