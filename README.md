[![Tests](https://github.com/oxfordmmm/gnomon/actions/workflows/tests.yaml/badge.svg)](https://github.com/oxfordmmm/gnomon/actions/workflows/tests.yaml)

# gnomon
Python code to integrate results of tb-pipeline and provide an antibiogram, mutations and variations

Provides a library of functions for use within scripts, as well as a CLI tool for linking the functions together to produce output

## Install
Currently there may be some issues with versions of [gumpy](https://github.com/oxfordmmm/gumpy)/[piezo](https://github.com/oxfordmmm/piezo) on pypi, so these may need to be installed from git beforehand.
```
git clone git@github.com:GlobalPathogenAnalysisService/gnomon.git
cd gnomon
pip install .
```
TODO: PyPi

## Docker
A Docker image should be built on releases. To open a shell with Gnomon installed:
```
docker run -it oxfordmmm/gnomon:latest
```
## User stories

1. As a bioinformatician, I want to be able to run `gnomon` on the command line, passing it (i) a GenBank file (or pickled `gumpy.Genome` object), (ii) a resistance catalogue and (iii) a VCF file, and get back `pandas.DataFrames` of the genetic variants, mutations, effects and predictions/antibiogram. The latter is for all the drugs described in the passed resistance catalogue.

2. As a GPAS developer, I want to be able to embed `gnomon` in a Docker image/NextFlow pipeline that consumes the outputs of [tb-pipeline](https://github.com/Pathogen-Genomics-Cymru/tb-pipeline) and emits a structured, well-designed `JSON` object describing the genetic variants, mutations, effects and predictions/antibiogram.

3. In general, I would also like the option to output fixed- and variable-length FASTA files (the latter takes into account insertions and deletions described in any input VCF file).

## Unit testing

For speed, rather than use NC_000962.3 (i.e. H37Rv *M. tuberculosis*), we shall use SARS-CoV-2 and have created a fictious drug resistance catalogue, along with some `vcf` files and the expected outputs in `tests/`.

These can be run with `pytest -vv`