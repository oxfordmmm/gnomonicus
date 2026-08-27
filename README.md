[![Tests](https://github.com/oxfordmmm/gnomonicus/actions/workflows/tests.yaml/badge.svg)](https://github.com/oxfordmmm/gnomonicus/actions/workflows/tests.yaml) 
[![Build and release Docker](https://github.com/oxfordmmm/gnomonicus/actions/workflows/build.yaml/badge.svg)](https://github.com/oxfordmmm/gnomonicus/actions/workflows/build.yaml) 
[![PyPI version](https://badge.fury.io/py/gnomonicus.svg)](https://badge.fury.io/py/gnomonicus)
[![Docs](https://github.com/oxfordmmm/gnomonicus/actions/workflows/docs.yaml/badge.svg)](https://oxfordmmm.github.io/gnomonicus/)

# gnomonicus
Python code to integrate results of tb-pipeline and provide an antibiogram, mutations and variations

Provides a library of functions for use within scripts, as well as a CLI tool for linking the functions together to produce output

## Documentation
API reference for developers, and CLI instructions can be found here: https://oxfordmmm.github.io/gnomonicus/ 
## Usage
```
usage: gnomonicus [-h] [-v] --vcf_file VCF_FILE --genome_object GENOME_OBJECT [--catalogue_file CATALOGUE_FILE] [--ignore_vcf_filter] [--output_dir OUTPUT_DIR] [--json] [--csvs CSVS [CSVS ...]] [--debug]
                  [--resistance_genes] --min_dp MIN_DP

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --vcf_file VCF_FILE   the path to a single VCF file
  --genome_object GENOME_OBJECT
                        the path to a genbank file
  --catalogue_file CATALOGUE_FILE
                        the path to the resistance catalogue
  --ignore_vcf_filter   whether to ignore the FILTER field in the vcf (e.g. necessary for some versions of Clockwork VCFs)
  --output_dir OUTPUT_DIR
                        Directory to save output files to. Defaults to wherever the script is run from.
  --json                Flag to create a single JSON output as well as the CSVs
  --csvs CSVS [CSVS ...]
                        Types of CSV to produce. Accepted values are [variants, mutations, effects, predictions, all]. `all` produces all of the CSVs
  --debug               Whether to log debugging messages to the log. Defaults to False
  --resistance_genes    Flag to filter mutations and variants to only include genes present in the resistance catalogue
  --min_dp MIN_DP       Minimum depth for a variant to be considered in the VCF. Below this value, rows are interpreted as null calls.
```

## Install
Simple install using pip for the latest release
```
pip install gnomonicus
```

Install from source
```
git clone https://github.com/oxfordmmm/gnomonicus.git
cd gnomonicus
pip install -e .
```

## Docker
A Docker image should be built on releases. To open a shell with gnomonicus installed:
```
docker run -it oxfordmmm/gnomonicus:latest
```

## Notes
When generating mutations, in cases of synonymous amino acid mutation, the nucelotides changed are also included. This can lead to a mix of nucleotides and amino acids for coding genes, but these are excluded from generating effects unless specified in the catalogue. This means that the default rule of `gene@*= --> S` is still in place regardless of the introduced `gene@*?` which would otherwise take precedence. For example:
```
  'MUTATIONS': [
      {
          'MUTATION': 'F2F',
          'GENE': 'S',
          'GENE_POSITION': 2
      },
      {
          'MUTATION': 't6c',
          'GENE': 'S',
          'GENE_POSITION': 6
      },
  ],
  'EFFECTS': {
      'AAA': [
          {
              'GENE': 'S',
              'MUTATION': 'F2F',
              'PREDICTION': 'S'
          },
          {
              'PHENOTYPE': 'S'
          }
      ],
  }
```
The nucelotide variation is included in the the `MUTATIONS`, but explictly removed from the `EFFECTS` unless it is specified within the catalogue.
In order for this variation to be included, a line in the catalogue of `S@F2F&S@t6c` would have to be present.

## Unit testing

For speed, rather than use NC_000962.3 (i.e. H37Rv *M. tuberculosis*), we shall use SARS-CoV-2 and have created a fictious drug resistance catalogue, along with some `vcf` files and the expected outputs in `tests/`.

These can be run with `pytest -vv`

## Citation
If you use `gnomonicus` in your work, please cite:
```
Westhead J, Baker CS, Brouard M, Colpus M, Constantinides B, Hall A, Knaggs J, Lopes Alves M, Spies R, Thai H, Surrall S, Govender K, Peto TEA, Crook DW, Omar SV, Turner R, Fowler PW
Characterising the performance of an antibiotic resistance prediction tool, gnomonicus, using a diverse testset of 2,663 Mycobacterium tuberculosis samples
Microbial Genomics 11:001592 doi:10.1099/mgen.0.001592
```

A BibTeX citation is included [here](CITATION.bib)