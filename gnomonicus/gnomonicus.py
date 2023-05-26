'''gnomonicus.py is a library providing functions which pull together output VCF of the Lodestone TB pipeline
    with a reference genome and a resistance catalogue, and utilise gumpy and
    piezo to produce variants, mutations and an antibiogram.

Based on sp3predict
'''
import datetime
import gzip
import json
import logging
import os
import pickle
import re
from collections import defaultdict
from collections.abc import Iterable

import gumpy
import numpy as np
import pandas as pd
import piezo
from tqdm import tqdm

class InvalidMutationException(Exception):
    '''Custom exception raised when an invalid mutation is detected
    '''
    def __init__(self, gene: str, mutation: str):
        '''Raise this exception

        Args:
            gene (str): Name of the gene
            mutation (str): The invalid mutation
        '''
        self.message = f"{gene}@{mutation} is not a valid mutation!"
        super().__init__(self.message)

def checkGzip(path: str) -> bool:
    '''Check if a given path is a gzipped file

    Args:
        path (str): Path to the file

    Returns:
        bool: True if the file is gzipped
    '''
    try:
        with gzip.open(path) as f:
            f.read()
        return True
    except:
        return False


def loadGenome(path: str, progress: bool) -> gumpy.Genome:
    '''Load a genome from a given path. Checks if path is to a pickle dump, or if a pickle dump of the path's file exists
    Instanciates a new gumpy Genome and dumps to pickle as a last resort

    Args:
        path (str): Path to the genbank file or pickle dump. If previously run, a genbank file's Genome object is pickled and dumped for speed
        progress (bool): Boolean as whether to show progress bar for gumpy

    Returns:
        gumpy.Genome: Genome object of the reference genome
    '''
    logging.debug(f"Using file {path}")
    #Remove trailing '/' if required
    if path[-1] == '/':
        path = path[:-1]

    #Check if the file is gzipped
    gzipped = checkGzip(path)

    #Try to load as a pickle
    try:
        if gzipped:
            logging.info("Path was to a gzipped file. Decompressing...")
            f = gzip.open(path, 'rb')
        else:
            logging.info("Path was not to a gzipped file. Defaulting to normal reading")
            f = open(path, 'rb')
        return pickle.load(f)
    except Exception as e:
        logging.info(f"Genome object not a pickle, checking if pickled version exists. Error: {e}")

    #Try pickled version created by this (path+'.pkl')
    #Check if this file is gzipped
    gzipped = checkGzip(path+".pkl")
    try:
        if gzipped:
            logging.info("Path was to a gzipped file. Decompressing...")
            f = gzip.open(path+".pkl", 'rb')
        else:
            logging.info("Path was not to a gzipped file. Defaulting to normal reading")
            f = open(path+".pkl", 'rb')
        return pickle.load(f)
    except Exception as e:
        logging.info(f"No pickled version of genome object, instanciating and dumping. Error: {e}")
    
    #Create new gumpy.Genome and pickle dump for speed later
    reference = gumpy.Genome(path, show_progress_bar=progress)
    pickle.dump(reference, open(path+'.pkl', 'wb'))
    return reference

def populateVariants(vcfStem: str, outputDir: str, diff: gumpy.GenomeDifference, make_csv: bool, catalogue: piezo.ResistanceCatalogue=None) -> pd.DataFrame:
    '''Populate and save the variants DataFrame as a CSV

    Args:
        vcfStem (str): The stem of the filename for the VCF file. Used as a uniqueID
        outputDir (str): Path to the desired output directory
        diff (gumpy.GenomeDifference): GenomeDifference object between reference and the sample
        make_csv (bool): Whether to write the CSV of the dataframe
        catalogue (piezo.ResistanceCatalogue, optional): Catalogue for determining FRS or COV for minority populations. If None is given, FRS is assumed. Defaults to None
    
    Returns:
        pd.DataFrame: DataFrame of the variants
    '''
    #Populate variants table directly from GenomeDifference
    vals = {
            'variant': diff.variants, 
            'nucleotide_index': diff.nucleotide_index,
            'indel_length': diff.indel_length,
            'indel_nucleotides': diff.indel_nucleotides,
            'vcf_evidence': diff.vcf_evidences
            }
    variants = pd.DataFrame(vals)

    if diff.genome1.minor_populations or diff.genome2.minor_populations:
        variants = pd.concat([variants, minority_population_variants(diff, catalogue)])

    #If there are variants, save them to a csv
    if not variants.empty:
        #Add unique ID to each record
        variants['uniqueid'] = vcfStem

        variants = variants[['uniqueid', 'variant', 'nucleotide_index', 'indel_length', 'indel_nucleotides', 'vcf_evidence']]
        if make_csv:
            #Save CSV
            variants.to_csv(os.path.join(outputDir, f'{vcfStem}.variants.csv'), header=True, index=False)
    variants.reset_index(inplace=True)
    return variants

def get_minority_population_type(catalogue: piezo.ResistanceCatalogue) -> str:
    '''Figure out if a catalogue uses FRS or COV. If neither or both, default to FRS

    Args:
        catalogue (piezo.ResistanceCatalogue): Catalogue

    Returns:
        str: Either 'percentage' or 'reads' for FRS or COV respectively
    '''
    if catalogue is None:
        #Nothing given, so default to FRS
        return 'percentage'
    frs = 0
    cov = 0
    for minor in catalogue.catalogue.rules['MINOR']:
        for m in minor.split(","):
            if m:
                m = float(m)
                assert m > 0, f"Minor populations must be positive: {m}"
                if m < 1:
                    #FRS
                    frs += 1
                else:
                    #COV
                    cov += 1
    #We have just COV
    if cov > 0 and frs == 0:
        return 'reads'
    #We have anything else
    return 'percentage'

def populateMutations(
        vcfStem: str, outputDir: str, diff: gumpy.GenomeDifference, reference: gumpy.Genome,
        sample: gumpy.Genome, resistanceCatalogue: piezo.ResistanceCatalogue, make_csv: bool) -> (pd.DataFrame, dict):
    '''Popuate and save the mutations DataFrame as a CSV, then return it for use in predictions

    Args:
        vcfStem (str): The stem of the filename of the VCF file. Used as a uniqueID
        outputDir (str): Path to the desired output directory
        diff (gumpy.GenomeDifference): GenomeDifference object between reference and this sample
        reference (gumpy.Genome): Reference genome
        sample (gumpy.Genome): Sample genome
        resistanceCatalogue (piezo.ResistanceCatalogue): Resistance catalogue (used to find which genes to check)
        make_csv (bool): Whether to write the CSV of the dataframe

    Raises:
        MissingFieldException: Raised when the mutations DataFrame does not contain the required fields

    Returns:
        pd.DataFrame: The mutations DataFrame
        dict: Dictionary mapping gene name --> reference gumpy.Gene object
    '''
    #For the sake of consistency, moving towards including all genes with mutations (not just those in the catalogue)
    if resistanceCatalogue:
        #Find the genes which have mutations regardless of being in the catalogue
        #Still cuts back time considerably, and ensures all mutations are included in outputs
        mask = np.isin(reference.stacked_nucleotide_index, diff.nucleotide_index)
        genesWithMutations = np.unique(reference.stacked_gene_name[mask]).tolist()

        #Make sure minority population mutations are also picked up
        minor_genes = set()
        for population in sample.minor_populations:
            for gene in reference.stacked_gene_name[reference.stacked_nucleotide_index == population[0]]:
                if gene:
                    minor_genes.add(gene)
        genesWithMutations += minor_genes

        deletions = []
        #Make sure large deletions are picked up too
        for name in reference.stacked_gene_name:
            deletions += np.unique(name[sample.is_deleted]).tolist()
        genesWithMutations = set(genesWithMutations + deletions)

    else:
        #No catalogue, so just stick to genes in the sample
        genesWithMutations = sample.genes

    #Iter resistance genes with variation to produce gene level mutations - concating into a single dataframe
    mutations = None
    referenceGenes = {}
    diffs = []
    #This is where the majority of the time to process is used.
    #However, I have tried a few ways to improve this involving reducing reliance on DF concat
    #This has potential for improvement (an actual TB sample can take ~2.5mins), but likely will
    #   come from optimising the underlying gumpy.GeneDifference code...
    for gene in tqdm(genesWithMutations):
        if gene:
            logging.debug(f"Found a gene with mutation: {gene}")
            #Save the reference genes for use later in effects.csv
            refGene = reference.build_gene(gene)
            referenceGenes[gene] = refGene
            #Get gene difference
            diff = refGene - sample.build_gene(gene)
            diffs.append(diff)

            #Pull the data out of the gumpy object
            vals = {
                'mutation': diff.mutations,
                'nucleotide_number': diff.nucleotide_number,
                'nucleotide_index': diff.nucleotide_index,
                'gene_position': diff.gene_position,
                'alt': diff.alt_nucleotides,
                'ref': diff.ref_nucleotides,
                'codes_protein': diff.codes_protein,
                'indel_length': diff.indel_length,
                'indel_nucleotides': diff.indel_nucleotides,
                }
            #As diff does not populate amino acid items for non-coding genes,
            #pull out the sequence or default to None
            if refGene.codes_protein:
                vals['amino_acid_number'] = diff.amino_acid_number
                aa_seq = []
                #Pull out the amino acid sequence from the alt codons
                for idx, num in enumerate(diff.amino_acid_number):
                    if num is not None:
                        aa_seq.append(refGene.codon_to_amino_acid[diff.alt_nucleotides[idx]])
                    else:
                        aa_seq.append(None)
                vals['amino_acid_sequence'] = np.array(aa_seq)
            else:
                vals['amino_acid_number'] = None
                vals['amino_acid_sequence'] = None
            
            geneMutations = pd.DataFrame(vals)
            #Add gene name
            geneMutations['gene'] = gene

            #Add this gene's mutations to the total dataframe
            if not geneMutations.empty:
                if mutations is not None:
                    mutations = pd.concat([mutations, geneMutations])
                else:
                    mutations = geneMutations
    if mutations is not None:
        #Ensure correct datatypes
        mutations = mutations.astype({'mutation': 'str',
                                    'gene': 'str',
                                    'nucleotide_number': 'float',
                                    'nucleotide_index': 'float',
                                    'gene_position': 'float',
                                    'alt': 'str',
                                    'ref': 'str',
                                    'codes_protein': 'bool',
                                    'indel_length': 'float',
                                    'indel_nucleotides': 'str',
                                    'amino_acid_number': 'float',
                                    'amino_acid_sequence': 'str',
                                })
    #Add minor mutations (these are stored separately)
    if reference.minor_populations or sample.minor_populations:
        #Only do this if they exist
        x = minority_population_mutations(diffs, resistanceCatalogue)
        mutations = pd.concat([mutations, x])
    #If there were mutations, write them to a CSV
    if mutations is not None:

        #Add the number of mutations which occured for this mutation
        mutations['number_nucleotide_changes'] = mutations.apply(countNucleotideChanges, axis=1)

        #Add VCF stem as the uniqueID
        mutations['uniqueid'] = vcfStem

        #Reorder the columns
        mutations = mutations[['uniqueid', 'gene', 'mutation', 'ref', 'alt', 'nucleotide_number', 'nucleotide_index', 'gene_position', 'codes_protein', 'indel_length', 'indel_nucleotides', 'amino_acid_number', 'amino_acid_sequence', 'number_nucleotide_changes']]

        if make_csv:
            #Save it as CSV
            mutations.to_csv(os.path.join(outputDir, f'{vcfStem}.mutations.csv'), index=False)

        #Remove index to return
        mutations.reset_index(inplace=True)
    return mutations, referenceGenes

def minority_population_variants(diff: gumpy.GenomeDifference, catalogue: piezo.ResistanceCatalogue) -> pd.DataFrame:
    '''Handle the logic for pulling out minority population calls for genome level variants

    Args:
        diff: (gumpy.GenomeDifference): GenomeDifference object for this comparison
        catalogue: (piezo.ResistanceCatalogue): Catalogue to use. Used for determining whether FRS or COV should be used

    Returns:
        pd.DataFrame: DataFrame containing the minority population data
    '''
    #Determine if FRS or COV should be used
    minor_type = get_minority_population_type(catalogue)

    #Get the variants in GARC
    variants_ = diff.minor_populations(interpretation=minor_type)
    variants = [v.split(":")[0] for v in variants_]


    #Split to be the same format
    #Not exactly efficient, but this should be so infrequent that it shouldn't be impactful
    vals = {
        'variant': variants_, 
        'nucleotide_index': [var.split(">")[0][:-1] if ">" in var else var.split("_")[0] for var in variants],
        'indel_length': [len(var.split("_")[-1]) if "_" in var else 0 for var in variants],
        'indel_nucleotides': [var.split("_")[-1] if "_" in var else None for var in variants],
        'vcf_evidence': [diff.genome2.vcf_evidence.get(
                                                        int(var.split(">")[0][:-1])
                                                    ) 
                        if ">" in var 
                        else diff.genome2.vcf_evidence.get(
                                                        int(var.split("_")[0])
                                                        )
                        for var in variants]
        }
    #Convert everything to numpy arrays
    vals = {key: np.array(vals[key]) for key in vals.keys()}
    return pd.DataFrame(vals).astype({
                                    'vcf_evidence': 'object'
                                })


def minority_population_mutations(diffs: [gumpy.GeneDifference], catalogue: piezo.ResistanceCatalogue) -> pd.DataFrame:
    '''Handle the logic for pulling out minority population calls for gene level variants

    Args:
        diffs ([gumpy.GeneDifference]): List of GeneDifference objects for these comparisons
        catalogue: (piezo.ResistanceCatalogue): Catalogue to use. Used for determining whether FRS or COV should be used

    Returns:
        pd.DataFrame: DataFrame containing the minority population data
    '''
    #Get the mutations
    mutations_ = []
    genes = []
    gene_pos = []
    nucleotide_number = []
    nucleotide_index = []
    alt = []
    ref = []
    codes_protein = []
    indel_length = []
    indel_nucleotides = []
    is_cds = []
    is_het = []
    is_null = []
    is_promoter = []
    is_snp = []
    aa_num = []
    aa_seq = []

    #Determine if FRS or COV should be used
    minor_type = get_minority_population_type(catalogue)

    for diff in diffs:
        #As mutations returns in GARC, split into constituents for symmetry with others
        mutations = diff.minor_populations(interpretation=minor_type)
        
        #Without gene names/evidence
        muts = [mut.split(":")[0] for mut in mutations]
        #Gene numbers
        numbers = [
            int(mut.split("_")[0]) if "_" in mut #Indel index: <idx>_<type>_<bases>
            else 
                int(mut[:-1]) if "=" in mut #Synon SNP: <idx>=
                else int(mut[1:][:-1]) #SNP: <ref><idx><alt>
            for mut in muts
            ]

        #Iter these to pull out all other details from the GeneDifference objects
        for mut, num, full_mut in zip(muts, numbers, mutations):
            mutations_.append(full_mut) #Keep evidence in these
            genes.append(diff.gene1.name)
            gene_pos.append(num)
            codes_protein.append(diff.gene1.codes_protein)
            is_cds.append(num > 0 and diff.gene1.codes_protein)
            is_het.append("Z" in mut.upper())
            is_null.append("X" in mut.upper())
            is_promoter.append(num < 0)

            if "_" in mut:
                #Indel
                _, t, bases = mut.split("_")
                ref.append(None)
                alt.append(None)
                if t == "del":
                    indel_length.append(-1 * len(bases))
                else:
                    indel_length.append(len(bases))
                indel_nucleotides.append(bases)
                is_snp.append(False)
                nucleotide_number.append(num)
                nucleotide_index.append(diff.gene1.nucleotide_index[diff.gene1.nucleotide_number == num][0])
                aa_num.append(None)
                aa_seq.append(None)                
                continue
            else:
                indel_length.append(None)
                indel_nucleotides.append(None)

            if mut[0].isupper() or mut[0] == '!':
                #Protein coding SNP
                nucleotide_number.append(None)
                nucleotide_index.append(None)
                #Pull out codons for ref/alt
                ref.append(diff.gene1.codons[diff.gene1.amino_acid_number == num][0])
                alt.append(diff.gene2.codons[diff.gene2.amino_acid_number == num][0])
                is_snp.append(True)
                aa_num.append(num)
                aa_seq.append(mut[-1])
            else:
                #Other SNPs
                nucleotide_number.append(num)
                nucleotide_index.append(diff.gene1.nucleotide_index[diff.gene1.nucleotide_number == num][0])
                aa_num.append(None)
                aa_seq.append(None)
                ref.append(diff.gene1.nucleotide_sequence[diff.gene1.nucleotide_number == num][0])
                alt.append(diff.gene2.nucleotide_sequence[diff.gene2.nucleotide_number == num][0])
                is_snp.append(True)

    vals = {
        'mutation': mutations_,
        'gene': genes,
        'nucleotide_number': nucleotide_number,
        'nucleotide_index': nucleotide_index,
        'gene_position': gene_pos,
        'alt': alt,
        'ref': ref,
        'codes_protein': codes_protein,
        'indel_length': indel_length,
        'indel_nucleotides': indel_nucleotides,
        'amino_acid_number': aa_num,
        'amino_acid_sequence': aa_seq,
        }

    return pd.DataFrame(vals).astype({'mutation': 'str',
                                        'gene': 'str',
                                        'nucleotide_number': 'float',
                                        'nucleotide_index': 'float',
                                        'gene_position': 'float',
                                        'alt': 'str',
                                        'ref': 'str',
                                        'codes_protein': 'bool',
                                        'indel_length': 'float',
                                        'indel_nucleotides': 'str',
                                        'amino_acid_number': 'float',
                                        'amino_acid_sequence': 'str',
                                    })

def countNucleotideChanges(row: pd.Series) -> int:
    '''Calculate the number of nucleotide changes required for a given row's amino acid mutation

    Args:
        row (pd.Series): A row of the mutations dataframe

    Returns:
        int: The number of mutations which occured to cause this mutation
    '''
    if row['ref'] is not None and len(row['ref'])==3:
        #Numpy sum is considerably slower for this...
        return sum(i!=j for (i,j) in zip(row['ref'],row['alt'] ))
    return 0

def getMutations(mutations: pd.DataFrame, catalogue: piezo.catalogue, referenceGenes: dict) -> [[str, str]]:
    '''Get all of the mutations (including multi-mutations) from the mutations df
    Multi-mutations currently only exist within the converted WHO catalogue, and are a highly specific combination 
        of mutations which must all be present for a single resistance value.

    Args:
        mutations (pd.DataFrame): Mutations dataframe
        catalogue (piezo.catalogue): The resistance catalogue. Used to find which multi-mutations we care about
        referenceGenes (dict): Dictionary of geneName->gumpy.Gene

    Returns:
        [[str, str]]: List of [gene, mutation] or in the case of multi-mutations, [None, multi-mutation]
    '''
    mutations = list(zip(mutations['gene'], mutations['mutation']))
    #Grab the multi-mutations from the catalogue
    #By doing this, we can check a fixed sample space rather than every permutation of the mutations
    #This makes the problem tractable, but does not address a possible issue with multi-mutations not encapsulating full codons
    multis = catalogue.catalogue.rules[catalogue.catalogue.rules['MUTATION_TYPE']=='MULTI']['MUTATION']
    if len(multis) > 0:
        #We have a catalogue including multi rules, so check if any of these are present in the mutations
        joined = [gene+'@'+mut for (gene, mut) in mutations]
        for multi in multis:
            check = True
            for mutation in multi.split("&"):
                check = check and mutation in joined
            if check:
                #This exact multi mutation exists, so add it to the mutations list
                mutations.append((None, multi))
    
    #Check if the catalogue supports large deletions
    if "GENE" in set(catalogue.catalogue.rules['MUTATION_AFFECTS']):
        large_dels = True
    else:
        large_dels = False

    #Filtering out *just* nucelotide changes for cases of synon mutations
    #The important part of these should have already been found by multi-mutations
    fixed = []
    for gene, mutation in mutations:
        if gene is not None and referenceGenes[gene].codes_protein:
            #Codes protein so check for nucleotide changes
            nucleotide = re.compile(r"""
                                [acgtzx][0-9]+[acgtzx]
                                """, re.VERBOSE)
            if nucleotide.fullmatch(mutation):
                #Is a nucleotide (non-promoter) mutation in a coding gene
                #So skip it as it may cause prediction problems
                continue
        #Remove large dels if not supported
        if not large_dels:
            #Check if this is a large del
            large = re.compile(r"""
                                del_(1\.0)|(0\.[0-9][0-9])
                                """, re.VERBOSE)
            if large.fullmatch(mutation):
                continue
        fixed.append((gene, mutation))
    return fixed

def populateEffects(
        outputDir: str, resistanceCatalogue: piezo.ResistanceCatalogue,
        mutations: pd.DataFrame, referenceGenes: dict, vcfStem: str, make_csv: bool) -> (pd.DataFrame, dict):
    '''Populate and save the effects DataFrame as a CSV

    Args:
        outputDir (str): Path to the directory to save the CSV
        resistanceCatalogue (piezo.ResistanceCatalogue): Resistance catalogue for predictions
        mutations (pd.DataFrame): Mutations dataframe
        referenceGenes (dict): Dictionary mapping gene name --> reference gumpy.Gene objects
        vcfStem (str): The basename of the given VCF - used as the sample name
        make_csv (bool): Whether to write the CSV of the dataframe

    Raises:
        InvalidMutationException: Raised if an invalid mutation is detected

    Returns:
        (pd.DataFrame, dict): (DataFrame containing the effects data, A metadata dictionary mapping drugs to their predictions)
    '''
    #Assume wildtype behaviour unless otherwise specified
    phenotype = {drug: 'S' for drug in resistanceCatalogue.catalogue.drugs}

    effects = {}
    effectsCounter = 0

    #Default prediction values are RFUS but use piezo catalogue's values if existing
    values = resistanceCatalogue.catalogue.values

    for (gene, mutation) in tqdm(getMutations(mutations, resistanceCatalogue, referenceGenes)):
        #Ensure its a valid mutation
        if gene is not None and not referenceGenes[gene].valid_variant(mutation):
            logging.error(f"Not a valid mutation {gene}@{mutation}")
            raise InvalidMutationException(gene, mutation)
        
        #Get the prediction
        if gene is not None:
            prediction = resistanceCatalogue.predict(gene+'@'+mutation)
        else:
            #This is a multi-mutation so is already of required format
            prediction = resistanceCatalogue.predict(mutation)

        #If the prediction is interesting, iter through drugs to find predictions
        if prediction != 'S':
            for drug in prediction.keys():
                #Check for empty strings
                #I don't think this can be reached... Commenting out for now
                # if not drug:
                    # continue

                #Prioritise values based on order within the values list
                if values.index(prediction[drug]) < values.index(phenotype[drug]):
                    #The prediction is closer to the start of the values list, so should take priority
                    phenotype[drug] = prediction[drug]

                #Add to the dict
                effects[effectsCounter] = [vcfStem, gene, mutation, resistanceCatalogue.catalogue.name, drug, prediction[drug]]
                #Increment counter
                effectsCounter += 1
    
    #Build the DataFrame
    effects = pd.DataFrame.from_dict(effects, 
                                        orient="index", 
                                        columns=["uniqueid", "gene", "mutation", 
                                            "catalogue_name", "drug", "prediction"]
                                        )
    effects = effects[["uniqueid", "gene", "mutation", "drug", "prediction", "catalogue_name"]]
    effects['catalogue_version'] = resistanceCatalogue.catalogue.version
    effects['prediction_values'] = ''.join(resistanceCatalogue.catalogue.values)
    effects['evidence'] = "{}"
    
    #Save as CSV
    if len(effects) > 0 and make_csv:
        effects.to_csv(os.path.join(outputDir, f'{vcfStem}.effects.csv'), index=False)

    effects.reset_index(inplace=True)

    #Return  the metadata dict to log later
    return effects, {"WGS_PREDICTION_"+drug: phenotype[drug] for drug in resistanceCatalogue.catalogue.drugs}

def saveJSON(variants, mutations, effects, path: str, guid: str, catalogue: piezo.ResistanceCatalogue, gnomonicusVersion: str, time_taken: float, reference: gumpy.Genome, vcf_path: str, reference_path: str, catalogue_path: str) -> None:
    '''Create and save a single JSON output file for use within GPAS. JSON structure:
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
                    'vcf_evidence': Parsed VCF row
                }, ...
            ],
            ?'mutations': [
                {
                    'mutation': Gene level mutation in GARC,
                    'gene': Gene name,
                    'gene_position': Position within the gene. Amino acid or nucleotide index depending on which is appropriate,
                    'vcf_evidence': Parsed VCF row
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
    Where fields with a preceeding '?' are not always present depending on data

    Args:
        path (str): Path to the directory where the variant/mutation/effect CSV files are saved. Also the output dir for this.
        guid (str): Sample GUID
        catalogue (piezo.ResistanceCatalogue): Catalogue used
        gnomonicusVersion (str): Semantic versioning string for the gnomonicus module. Can be accessed by `gnomonicus.__version__`
        time_taken (float): Number of seconds taken to run this.
        reference (gumpy.Genome): Reference genome object
        vcf_path (str): Path to the VCF file used for this run
        reference_path (str): Path to the reference genome used for this run
        catalogue_path (str): Path to the catalogue used for this run
    '''
    values = catalogue.catalogue.values
    #Define some metadata for the json
    meta = {
        'status': 'success',
        'workflow_name': 'gnomonicus',
        'workflow_version': gnomonicusVersion, #gnomonicus version used
        'workflow_task': 'resistance_prediction', #TODO: Update this when we know how to detect a virulence catalogue
        'guid': guid, #Sample GUID
        'UTC-datetime-completed': datetime.datetime.utcnow().isoformat(), #ISO datetime run
        'time_taken': time_taken,
        'reference': reference.name,
        'catalogue_type': 'discrete_values', #TODO: Update this when we have MIC catalogues
        'catalogue_name': catalogue.catalogue.name,
        'catalogue_version': catalogue.catalogue.version,
        'catalogue_file': catalogue_path,
        'reference_file': reference_path,
        'vcf_file': vcf_path
        }
    data = {}
    #Variants field
    _variants = []
    for _, variant in variants.iterrows():
        row = {
            'variant': variant['variant'],
            'nucleotide_index': variant['nucleotide_index'],
            'vcf_evidence': variant['vcf_evidence']
        }
        _variants.append(row)
    data['variants'] = _variants

    #Depending on mutations/effects, populate
    _mutations = []
    if mutations is not None:
        for _, mutation in mutations.iterrows():
            row = {
                'mutation': mutation['mutation'],
                'gene': mutation['gene'],
                'gene_position': mutation['gene_position'],
            }
            _mutations.append(row)
        data['mutations'] = _mutations

    _effects = defaultdict(list)
    antibiogram = {}
    if effects is not None and len(effects) > 0:
        for _, effect in effects.iterrows():
            prediction = {
                'gene': effect['gene'],
                'mutation': effect['mutation'],
                'prediction': effect['prediction'],
                'evidence': {}
            }
            _effects[effect['drug']].append(prediction)
        
        #Get the overall predictions for each drug
        for drug, predictions in _effects.items():
            phenotype = 'S'
            for prediction in predictions:
                #Use the prediction heierarchy to use most signifiant prediction
                if values.index(prediction['prediction']) < values.index(phenotype):
                    #The prediction is closer to the start of the values list, so should take priority
                    phenotype = prediction['prediction']
            _effects[drug].append({'phenotype': phenotype})
            antibiogram[drug] = phenotype
        data['effects'] = _effects
    data['antibiogram'] = antibiogram

    #Convert fields to a list so it can be json serialised
    with open(os.path.join(path, f'{guid}.gnomonicus-out.json'), 'w') as f:
        f.write(json.dumps({'meta': meta, 'data': data}, indent=2, sort_keys=True))

