# Outputs
There are a total of 7 possible outputs from `gnomonicus`. Most of these are controllable via CLI flags.

## `<guid>.gnomonicus.log`
This is a log file detailing the execution process. It is always produced.

## `<guid>.variants.csv`
This CSV is produced by default (change with the `--csvs` flag). It details the genome level variants within this sample. 

This includes fields to allow joining to the `mutations` table, including gene name. As a result, it is possible that the same variant is present in >1 gene - resulting in multiple rows within this table.

Details of the fields can be found below:

| Field | Description | 
| ----- | ----------- |
| uniqueid | Sample GUID |
| variant | Genome level variant in GARC |
| gene_name | Name of the gene which this variant affects |
| gene_pos | Position within the gene this affects (nucleotide number or codon number as appropriate) |
| codon_idx | 0-based index of the nucleotide within the codon that this variant affects (as appropriate) |
| nucleotide_index | Genome level index of this variant | 
| indel_length | Length of this indel (as appropriate) |
|indel_nucleotides | Bases inserted or deleted in this variant (as appropriate) |
| vcf_evidence | A representation of the VCF row which caused this variant |

## `<guid>.mutations.csv`
This CSV is produced by default (change with the `--csvs` flag). It details gene level mutations within this sample. Positions given refer to nucleotides (promoter positions are <0), or amino acid numbers (if the gene codes protein).

Details of the fields can be found below:

| Field | Description | 
| ----- | ----------- |
| uniqueid | Sample GUID |
| gene | Gene name |
| mutation | Mutation in GARC |
| ref | Reference base(s). For nucleotide SNPs, ref is a nucleotide. For amino acid SNPs, ref is a codon. Not given for indels |
| alt | Alternate base(s). For nucleotide SNPs, alt is a nucleotide. For amino acid SNPs, alt is a codon. In minority populations, the codon given is `zzz`. Not given for indels |
| nucleotide_number | Gene index of the mutation. Not present in amino acid SNPs due to ambiguity. |
| nucleotide_index | Genome index of the mutation. Not present in amino acid SNPs due to ambiguity. |
| gene_position | Gene position of the mutation. Given as gene index for non-coding mutations, and amino acid index for coding mutations |
| codes_protein | Boolean for whether this gene codes protein |
| indel_length | Length of this indel (as appropriate) |
| indel_nucleotides | Bases inserted or deleted (as appropriate) |
| amino_acid_number | Amino acid index of this mutation (as appropriate) |
| amino_acid_sequence | Alternate amino acid (as appropriate) |
| number_nucleotide_changes | Number of nucleotides which are changed within the ref and alt. For minor populations, this is the number of minor population SNPs within this base/codon |

## `<guid>.effects.csv`
This CSV is produced by default (change with the `--csvs` flag). It details the effects of mutations according to the supplied resistance catalogue.

Details of the fields can be found below:

| Field | Description |
| ----- | ----------- |
| uniqueid | Sample GUID |
| gene | Name of the gene this mutation is in |
| mutation | Mutation in GARC |
| drug | Drug which this refers to |
| prediction | Prediction value. Usually one of `R, U, F, S` |
| catalogue_name | Name of the catalogue used |
| catalogue_version | Version of the catalogue used |
| prediction_values | Prediction values within the catalogue. Usually `RUS` or `RUFS`. May be numeric if an MIC catalogue is used |
| evidence | Evidence to support this prediction from the catalogue |

## `<guid>.gnomonicus-out.json`
This JSON is not produced by default (use the `--json` flag to enable). It combines the outputs of all of the above CSVs in a structured manner with additional metadata.

Below is the JSON structure where fields with a preceeding '?' are not always present depending on data:

```
{
    'meta': {
        'status': If this has succeeded or not (but this isn't created in cases it doesn't succeed),
        'workflow_name': 'gnomonicus',
        'workflow_task': 'resistance_prediction' or 'virulenece prediction',
        'workflow_version': gnomonicus.__version__,
        'time_taken': Time this step took,
        'UTC_timestamp': Timestamp for the end of this run,
        'catalogue_type': discrete_values or mic,
        'catalogue_name': Name of catalogue,
        'catalogue_version': Version of the catalogue,
        'reference': Name of the reference genome used,
        'catalogue_file': Path to the catalogue,
        'reference_file': Path to the reference file,
        'vcf_file': Path to the VCF file
    },
    'data': {
        'variants': [
            {
                'variant': Genome level variant in GARC,
                'nucleotide_index': Genome index of variant,
                ?'gene_name': Name of the gene this variant lies in (if applicable),
                ?'gene_position': Gene level nucleotide number or codon number (as appropriate) for where this variant affects,
                ?'codon_idx': 0-based index of the base within the codon this affects (if applicable)
                'vcf_evidence': Parsed VCF row
            }, ...
        ],
        ?'mutations': [
            {
                'mutation': Gene level mutation in GARC,
                'gene': Gene name,
                'gene_position': Position within the gene. Amino acid or nucleotide index depending on which is appropriate,
                'vcf_evidence': Parsed VCF row,
                'ref': Ref base(s),
                'alt': Alt base(s)
            }
        ],
        ?'effects': {
            Drug name: [
                {
                    'gene': Gene name of the mutation,
                    'mutation': Gene level mutation in GARC,
                    'prediction': Prediction caused by this mutation
                }, ...,
                {
                    'phenotype': Resultant prediction for this drug based on prediciton heirarchy
                }
            ], ...
        }
        ?'antibiogram': {
            <drug> : <prediction> essentially json[data][effects][<drug>][phenotype]
        }
    }
}
```


## FASTA file
FASTA files of the sample can also be produced (using the `--fasta` flag). These are not produced by default.

### `<guid>-fixed.fasta`
FASTA file with a fixed length of the size of the reference genome.

### `<guid>-variable.fasta`
FASTA file which is not always the size of the reference genome. This includes insertions and deletions within the nucleotide sequence.
