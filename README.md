# Filter structural variants against a set of genomic coordinates

Author: Bernie Pope (bjpope@unimelb.edu.au)

## License

3 Clause BSD License. See LICENSE.txt in source repository.

## Installation

#### External dependencies

`sv_filter` depends on the following programs and libraries:

 * [python](https://www.python.org/download/releases/2.7.5/) (version 2.7.5)
 * [bx-python](https://pypi.python.org/pypi/bx-python) 
 * [pyVCF](https://pypi.python.org/pypi/PyVCF)
 * [pybedtools](https://pypi.python.org/pypi/pybedtools)

I recommend using a virtual environment:

```
cd /place/to/install
virtualenv sv_filter 
source sv_filter/bin/activate
pip install -U https://github.com/bjpop/sv_filter
```

If you don't want to use a virtual environment then you can just install with pip:

```
pip install -U https://github.com/bjpop/sv_filter
```

## Worked example


## Usage

You can get a summary of the command line arguments like so:

```
sv_filter -h
```
