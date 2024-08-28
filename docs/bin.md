# gnomonicus
`gnomonicus` is a script which links together the functions defined in gnomonicus.py as a CLI script to produce variants, mutations and an antibiogram from a minos VCF, reference genome and a resistance catalogue.
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

# merge-vcfs
Script for merging a minos VCF with a gvcf file at specific positions for including null calls at certain positions. Required as minos will not consistently report null calls at all positions, and we want to be able to highlight occasions that we have null calls on resistance determining sites.
```
usage: merge-vcfs [-h] --minos_vcf MINOS_VCF --gvcf GVCF --resistant-positions RESISTANT_POSITIONS --output OUTPUT

Merge a minos VCF with a GVCF at certain positions (driven by the catalogue).

options:
  -h, --help            show this help message and exit
  --minos_vcf MINOS_VCF
                        The minos VCF filepath
  --gvcf GVCF           The GVCF filepath
  --resistant-positions RESISTANT_POSITIONS
                        Path to list of resistant sites
  --output OUTPUT       The output VCF file path

```

