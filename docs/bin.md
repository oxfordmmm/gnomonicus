# gnomonicus
`gnomonicus` is a script which links together the functions defined in gnomonicus.py as a CLI script to produce variants, mutations and an antibiogram from a minos VCF, reference genome and a resistance catalogue.
```
usage: gnomonicus [-h] --vcf_file VCF_FILE --genome_object GENOME_OBJECT [--catalogue_file CATALOGUE_FILE] [--ignore_vcf_filter] [--progress] [--output_dir OUTPUT_DIR] [--json] [--alt_json] [--fasta FASTA]
                  [--minor_populations MINOR_POPULATIONS]

options:
  -h, --help            show this help message and exit
  --vcf_file VCF_FILE   the path to a single VCF file
  --genome_object GENOME_OBJECT
                        the path to a compressed gumpy Genome object or a genbank file
  --catalogue_file CATALOGUE_FILE
                        the path to the resistance catalogue
  --ignore_vcf_filter   whether to ignore the FILTER field in the vcf (e.g. necessary for some versions of Clockwork VCFs)
  --progress            whether to show progress using tqdm
  --output_dir OUTPUT_DIR
                        Directory to save output files to. Defaults to wherever the script is run from.
  --json                Flag to create a single JSON output as well as the CSVs
  --alt_json            Whether to produce the alternate JSON format. Requires the --json flag too
  --fasta FASTA         Use to output a FASTA file of the resultant genome. Specify either 'fixed' or 'variable' for fixed length and variable length FASTA respectively.
  --minor_populations MINOR_POPULATIONS
                        Path to a line separated file containing genome indices of minor populations.
```

# gbkToPkl
gbkToPkl is a simple script which instanciates a gumpy.Genome object from a specified
genbank file, and saves it as a pickle. Due to the security implications of 
pickle, use at your own risk!!

Designed to be run once on the host to provide significant speed up for containerised
workflows. Resultant pickles should not be sent/recieved!!
```
usage: gbkToPkl [-h] [--compress] genbank_file

positional arguments:
  genbank_file  Path to a genbank file

options:
  -h, --help    show this help message and exit
  --compress    Whether to gzip compress the output pickle.
```


