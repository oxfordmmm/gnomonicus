[![Tests](https://github.com/oxfordmmm/gnomonicus/actions/workflows/tests.yaml/badge.svg)](https://github.com/oxfordmmm/gnomonicus/actions/workflows/tests.yaml)

# gnomonicus
Python code to integrate results of tb-pipeline and provide an antibiogram, mutations and variations

Provides a library of functions for use within scripts, as well as a CLI tool for linking the functions together to produce output

## Usage
```
usage: gnomonicus [-h] --vcf_file VCF_FILE --genome_object GENOME_OBJECT [--catalogue_file CATALOGUE_FILE]
              [--ignore_vcf_filter] [--progress] [--output_dir OUTPUT_DIR] [--json] [--alt_json] [--fasta FASTA]

options:
  -h, --help            show this help message and exit
  --vcf_file VCF_FILE   the path to a single VCF file
  --genome_object GENOME_OBJECT
                        the path to a compressed gumpy Genome object or a genbank file
  --catalogue_file CATALOGUE_FILE
                        the path to the resistance catalogue
  --ignore_vcf_filter   whether to ignore the FILTER field in the vcf (e.g. necessary for some versions of
                        Clockwork VCFs)
  --progress            whether to show progress using tqdm
  --output_dir OUTPUT_DIR
                        Directory to save output files to. Defaults to wherever the script is run from.
  --json                Flag to create a single JSON output as well as the CSVs
  --alt_json            Whether to produce the alternate JSON format. Requires the --json flag too
  --fasta FASTA         Use to output a FASTA file of the resultant genome. Specify either 'fixed' or 'variable'
                        for fixed length and variable length FASTA respectively.
```

## Helper usage
As the main script can utilise pickled `gumpy.Genome` objects, there is a supplied helper script. This converts a Genbank file into a pickled gumpy.Genome for significant time saving.
Due to the security implications of the pickle module, **DO NOT SEND/RECEIVE PICKLES**. This script should be used on a host VM before running nextflow to avoid reinstanciation.
Supports gzip compression to reduce file size significantly (using the `--compress` flag).
```
usage: gbkToPkl FILENAME [--compress]
```

## Install
Currently there may be some issues with versions of [gumpy](https://github.com/oxfordmmm/gumpy)/[piezo](https://github.com/oxfordmmm/piezo) on pypi, so these may need to be installed from git beforehand.
```
git clone git@github.com:GlobalPathogenAnalysisService/gnomonicus.git
cd gnomonicus
pip install .
```
TODO: PyPi

## Docker
A Docker image should be built on releases. To open a shell with gnomonicus installed:
```
docker run -it oxfordmmm/gnomonicus:latest
```
## User stories

1. As a bioinformatician, I want to be able to run `gnomonicus` on the command line, passing it (i) a GenBank file (or pickled `gumpy.Genome` object), (ii) a resistance catalogue and (iii) a VCF file, and get back `pandas.DataFrames` of the genetic variants, mutations, effects and predictions/antibiogram. The latter is for all the drugs described in the passed resistance catalogue.

2. As a GPAS developer, I want to be able to embed `gnomonicus` in a Docker image/NextFlow pipeline that consumes the outputs of [tb-pipeline](https://github.com/Pathogen-Genomics-Cymru/tb-pipeline) and emits a structured, well-designed `JSON` object describing the genetic variants, mutations, effects and predictions/antibiogram.

3. In general, I would also like the option to output fixed- and variable-length FASTA files (the latter takes into account insertions and deletions described in any input VCF file).

## Unit testing

For speed, rather than use NC_000962.3 (i.e. H37Rv *M. tuberculosis*), we shall use SARS-CoV-2 and have created a fictious drug resistance catalogue, along with some `vcf` files and the expected outputs in `tests/`.

These can be run with `pytest -vv`