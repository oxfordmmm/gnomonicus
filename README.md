# gnomon
Python code to integrate results of tb-pipeline and provide an antibiogram, mutations etc

## User stories

1. As a bioinformatician, I want to be able to run `gnomon` on the command line, passing it (i) a GenBank file (or pickled `gumpy.Genome` object), (ii) a resistance catalogue and (iii) a VCF file, and get back `pandas.DataFrames` of the genetic variants, mutations, effects and predictions/antibiogram. The latter is for all the drugs described in the passed resistance catalogue.

2. As a GPAS developer, I want to be able to embed `gnomon` in a Docker image/NextFlow pipeline that consumes the outputs of [tb-pipeline](https://github.com/Pathogen-Genomics-Cymru/tb-pipeline) and emits a structured, well-designed `JSON` object describing the genetic variants, mutations, effects and predictions/antibiogram.

3. In general, I would also like the option to output fixed- and variable-length FASTA files (the latter takes into account insertions and deletions described in any input VCF file).

## Unit testing

For speed, rather than use NC_000962.3 (i.e. H37Rv *M. tuberculosis*), we shall use SARS-CoV-2 and create a fictious drug resistance catalogue. PWF to do and provide `vcf` files, the catalogue and the expected results.
