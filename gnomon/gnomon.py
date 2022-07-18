'''gnomon.py is a library providing functions which pull together output VCF of the Lodestone TB pipeline
    with a reference genome and a resistance catalogue, and utilise gumpy and
    piezo to produce variants, mutations and an antibiogram.

Based on sp3predict
'''
import logging
import os
import pickle
import datetime
import json
from collections.abc import Iterable
from collections import defaultdict

import gumpy
import numpy as np
import pandas as pd
import piezo


class MissingFieldException(Exception):
    '''Custom exception for when required fields are missing from a given table
    '''
    def __init__(self, field: str, table: str) -> None:
        '''Raise this exception

        Args:
            field (str): Field name missing
            table (str): Table which the field is missing from
        '''
        self.message = f"Field: {field} is not in the {table} table!"
        super().__init__(self.message)

class NoVariantsException(Exception):
    '''Custom exception raised when there are no variants detected
    '''
    def __init__(self):
        super().__init__("No variants were detected!")

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

    #Try to load as a pickle
    try:
        return pickle.load(open(path, "rb"))
    except Exception as e:
        logging.info(f"Genome object not a pickle, checking if pickled version exists. Error: {e}")

    #Try pickled version created by this (path+'.pkl')
    try:
        return pickle.load(open(path+'.pkl', 'rb'))
    except Exception as e:
        logging.info(f"No pickled version of genome object, instanciating and dumping. Error: {e}")
    
    #Create new gumpy.Genome and pickle dump for speed later
    reference = gumpy.Genome(path, show_progress_bar=progress)
    pickle.dump(reference, open(path+'.pkl', 'wb'))
    return reference

def populateVariants(vcfStem: str, outputDir: str, diff: gumpy.GenomeDifference) -> pd.DataFrame:
    '''Populate and save the variants DataFrame as a CSV

    Args:
        vcfStem (str): The stem of the filename for the VCF file. Used as a uniqueID
        outputDir (str): Path to the desired output directory
        diff (gumpy.GenomeDifference): GenomeDifference object between reference and the sample
    
    Returns:
        pd.DataFrame: DataFrame of the variants
    '''
    #Populate variants table directly from GenomeDifference
    vals = {
            'VARIANT': diff.variants, 
            'NUCLEOTIDE_INDEX': diff.nucleotide_index,
            'IS_INDEL': diff.is_indel,
            'IS_NULL': diff.is_null,
            'IS_HET': diff.is_het,
            'IS_SNP': diff.is_snp,
            'INDEL_LENGTH': diff.indel_length,
            'INDEL_NUCLEOTIDES': diff.indel_nucleotides
            }
    vals = handleIndels(vals)
    variants = pd.DataFrame(vals)

    #If there are variants, save them to a csv
    if not variants.empty:
        #Add unique ID to each record
        variants['UNIQUEID'] = vcfStem

        #Double check for required fields
        if 'VARIANT' not in variants.columns:
            logging.error("VARIANT not in variant table!")
            raise MissingFieldException('VARIANT', 'variant')
        if 'IS_SNP' not in variants.columns:
            logging.error("IS_SNP not in variant table!")
            raise MissingFieldException('IS_SNP', 'variant')

        #Set the index
        variants.set_index(['UNIQUEID', 'VARIANT', 'IS_SNP'], inplace=True, verify_integrity=True)
        #Save CSV
        variants.to_csv(os.path.join(outputDir, 'variants.csv'), header=True)
    variants.reset_index(inplace=True)
    return variants

def populateMutations(
        vcfStem: str, outputDir: str, diff: gumpy.GenomeDifference, reference: gumpy.Genome,
        sample: gumpy.Genome, resistanceCatalogue: piezo.ResistanceCatalogue) -> (pd.DataFrame, dict):
    '''Popuate and save the mutations DataFrame as a CSV, then return it for use in predictions

    Args:
        vcfStem (str): The stem of the filename of the VCF file. Used as a uniqueID
        outputDir (str): Path to the desired output directory
        diff (gumpy.GenomeDifference): GenomeDifference object between reference and this sample
        reference (gumpy.Genome): Reference genome
        sample (gumpy.Genome): Sample genome
        resistanceCatalogue (piezo.ResistanceCatalogue): Resistance catalogue (used to find which genes to check)

    Raises:
        MissingFieldException: Raised when the mutations DataFrame does not contain the required fields

    Returns:
        pd.DataFrame: The mutations DataFrame
        dict: Dictionary mapping gene name --> reference gumpy.Gene object
    '''
    #TODO: Determine if all mutations will be required - currently cuts back to just genes in the catalogue
    if resistanceCatalogue:
        #Find the resistance associated genes which also have mutations in this sample
        #This considerably cuts back on the number of genes which have to be explored
        mask = np.isin(reference.stacked_nucleotide_index, diff.nucleotide_index)
        genesWithMutations = np.unique(reference.stacked_gene_name[mask])
        detectedResistanceGenes = set(resistanceCatalogue.catalogue.genes).intersection(set(genesWithMutations))
    else:
        #No catalogue, so just stick to genes in the sample
        detectedResistanceGenes = sample.genes

    #Iter resistance genes with variation to produce gene level mutations - concating into a single dataframe
    mutations = None
    referenceGenes = {}
    #TODO: Might be worth optimising this to avoid concating dfs as they are expensive to work with?
    for gene in detectedResistanceGenes:
        if gene:
            logging.debug(f"Found a resistance gene with mutation: {gene}")
            #Save the reference genes for use later in effects.csv
            refGene = reference.build_gene(gene)
            referenceGenes[gene] = refGene
            #Get gene difference
            diff = refGene - sample.build_gene(gene)

            #Pull the data out of the gumpy object
            vals = {
                'MUTATION': diff.mutations,
                'NUCLEOTIDE_NUMBER': diff.nucleotide_number,
                'NUCLEOTIDE_INDEX': diff.nucleotide_index,
                'GENE_POSITION': diff.gene_position,
                'ALT': diff.alt_nucleotides,
                'REF': diff.ref_nucleotides,
                'CODES_PROTEIN': diff.codes_protein,
                'INDEL_LENGTH': diff.indel_length,
                'INDEL_NUCLEOTIDES': diff.indel_nucleotides,
                'IS_CDS': diff.is_cds,
                'IS_HET': diff.is_het,
                'IS_NULL': diff.is_null,
                'IS_PROMOTER': diff.is_promoter,
                'IS_SNP': diff.is_snp,
                'AMINO_ACID_NUMBER': diff.amino_acid_number,
                'AMINO_ACID_SEQUENCE': diff.amino_acid_sequence
                }
            vals = handleIndels(vals)
            
            geneMutations = pd.DataFrame(vals)
            #Add gene name
            geneMutations['GENE'] = gene

            #Add this gene's mutations to the total dataframe
            if not geneMutations.empty:
                if mutations is not None:
                    mutations = pd.concat([mutations, geneMutations])
                else:
                    mutations = geneMutations
    #If there were mutations, write them to a CSV
    if mutations is not None:
        #Add synonymous booleans for analysis later
        mutations[['IS_SYNONYMOUS','IS_NONSYNONYMOUS']] = mutations.apply(assignMutationBools, axis=1)

        #Add the number of mutations which occured for this mutation
        mutations['NUMBER_NUCLEOTIDE_CHANGES'] = mutations.apply(countNucleotideChanges, axis=1)

        #Add VCF stem as the uniqueID
        mutations['UNIQUEID'] = vcfStem

        #Ensure correct datatypes
        mutations = mutations.astype({'IS_SYNONYMOUS':'bool',\
                                        'IS_NONSYNONYMOUS':'bool',\
                                        'IS_HET':'bool',\
                                        'IS_NULL':'bool',\
                                        'NUMBER_NUCLEOTIDE_CHANGES':'int'})

        #Raise errors if required fields are missing
        if 'GENE' not in mutations.columns:
            logging.error("GENE is not in the mutations table")
            raise MissingFieldException('GENE', 'mutations')
        if 'MUTATION' not in mutations.columns:
            logging.error("MUTATION is not in the mutations table")
            raise MissingFieldException('MUTATION', 'mutations')

        #Set the index
        mutations.set_index(['UNIQUEID', 'GENE', 'MUTATION'], inplace=True, verify_integrity=True)

        #Save it as CSV
        mutations.to_csv(os.path.join(outputDir, 'mutations.csv'))

        #Remove index to return
        mutations.reset_index(inplace=True)
    return mutations, referenceGenes

def handleIndels(vals: dict) -> dict:
    '''Due to how piezo works, we need to add alternative representations for indels
        Gumpy returns indel mutations as <pos>_(ins|del)_<nucleotides>
        Piezo catalogues may be that specific, or may require <pos>_indel or <pos>_(ins|del)_<n bases>

    Args:
        vals (dict): Initial mutation/variant values

    Returns:
        dict: Mutation/variant values with added indels as required
    '''
    #Determine if these are values from variants or mutations (as they have different struct)
    if 'MUTATION' in vals.keys():
        access = 'MUTATION'
    elif 'VARIANT' in vals.keys():
        access = 'VARIANT'
    else:
        raise NoVariantsException()

    toAdd = {key: [] for key in vals.keys() if isinstance(vals[key], Iterable)}
    toRemove = []
    for (n, mutation) in enumerate(vals[access]):
        if 'ins' in mutation or 'del' in mutation:
            #Is an indel so prepare extras to add and remove this
            toRemove.append(n)
            pos, indel, bases = mutation.split("_")
            indels = [
                pos + "_" + indel + "_" + bases,
                pos + "_indel", 
                pos + "_" + indel + "_" + str(len(bases))
                ]
            #Add the extras to `toAdd`
            for i in indels:
                for key in toAdd.keys():
                    #Add the extra to the MUTATION/VARIANT field
                    if key == access:
                        toAdd[key].append(i)
                    #Leave the other fields unchanged
                    else:
                        toAdd[key].append(vals[key][n])
    
    #Concat the two dicts
    for key in toAdd.keys():
        vals[key] = vals[key].tolist()
        for n in toRemove:
            del vals[key][n]
        vals[key] = vals[key] + toAdd[key]
    return vals
    

def assignMutationBools(row: pd.Series) -> pd.Series:
    '''Create the appropriate values of 'IS_SYNONYMOUS' and ''IS_NONSYNONYMOUS' for 
    the given row of a mutations dataframe

    Args:
        row (pd.Series): A row of the mutations dataframe

    Returns:
        pd.Series: A pandas series corresponding to [IS_SYNONYMOUS, IS_NONSYNONYMOUS]
    '''
    #In sp3predict this adds IS_HET and IS_NULL too, but these are added from gumpy.GeneDifference

    #Synonymous mutations are within amino acids so check the mutation to see if of form X...X
    if row['MUTATION'].isupper():
        isSynon = row['MUTATION'][0] == row['MUTATION'][-1]
        isNonSynon = not isSynon
    else:
        #Not an amino acid mutation so both False
        isSynon = isNonSynon = False
    
    return pd.Series([isSynon, isNonSynon])

def countNucleotideChanges(row: pd.Series) -> int:
    '''Calculate the number of nucleotide changes required for a given row's amino acid mutation

    Args:
        row (pd.Series): A row of the mutations dataframe

    Returns:
        int: The number of mutations which occured to cause this mutation
    '''
    if row['REF'] is not None and len(row['REF'])==3:
        #Numpy sum is considerably slower for this...
        return sum(i!=j for (i,j) in zip(row['REF'],row['ALT'] ))
    return 0

def populateEffects(
        sample: gumpy.Genome, outputDir: str, resistanceCatalogue: piezo.ResistanceCatalogue,
        mutations: pd.DataFrame, referenceGenes: dict) -> (pd.DataFrame, dict):
    '''Populate and save the effects DataFrame as a CSV

    Args:
        sample (gumpy.Genome): The Genome object of the sample
        outputDir (str): Path to the directory to save the CSV
        resistanceCatalogue (piezo.ResistanceCatalogue): Resistance catalogue for predictions
        mutations (pd.DataFrame): Mutations dataframe
        referenceGenes (dict): Dictionary mapping gene name --> reference gumpy.Gene objects

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
    if not values:
        values = list("RFUS")

    for (gene, mutation) in zip(mutations['GENE'], mutations['MUTATION']):
        #Ensure its a valid mutation
        if not referenceGenes[gene].valid_variant(mutation):
            logging.error(f"Not a valid mutation {gene}@{mutation}")
            raise InvalidMutationException(gene, mutation)
        
        #Get the prediction
        prediction = resistanceCatalogue.predict(gene+'@'+mutation)

        #If the prediction is interesting, iter through drugs to find predictions
        if prediction != 'S':
            for drug in prediction.keys():
                #Check for empty strings
                if not drug:
                    continue

                #Prioritise values based on order within the values list
                if values.index(prediction[drug]) < values.index(phenotype[drug]):
                    #The prediction is closer to the start of the values list, so should take priority
                    phenotype[drug] = prediction[drug]

                #Add to the dict
                effects[effectsCounter] = [sample.name, gene, mutation, resistanceCatalogue.catalogue.name, drug, prediction[drug]]
        else:
            #Add the no effect to the dict
            effects[effectsCounter] = [sample.name, gene, mutation, resistanceCatalogue.catalogue.name, "UNK", "S"]
        #Increment counter
        effectsCounter += 1
    
    #Build the DataFrame
    effects = pd.DataFrame.from_dict(effects, 
                                        orient="index", 
                                        columns=["UNIQUEID", "GENE", "MUTATION", 
                                            "CATALOGUE_NAME", "DRUG", "PREDICTION"]
                                        )
    
    #Set the index
    effects.set_index(["UNIQUEID", "DRUG", "GENE", "MUTATION", "CATALOGUE_NAME"], inplace=True)
    
    #Save as CSV
    effects.to_csv(os.path.join(outputDir, 'effects.csv'))

    effects.reset_index(inplace=True)

    #Return  the metadata dict to log later
    return effects, {"WGS_PREDICTION_"+drug: phenotype[drug] for drug in resistanceCatalogue.catalogue.drugs}

def saveJSON(variants, mutations, effects, path: str, guid: str, values: list, gnomonVersion: str) -> None:
    '''Create and save a single JSON output file for use within GPAS. JSON structure:
    {
        'meta': {
            'version': gnomon version,
            'guid': sample GUID,
            'UTC-datetime-run': ISO formatted UTC datetime run,
            'fields': Dictionary of data fields for parsing
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': Genome level variant in GARC,
                    'NUCLEOTIDE_INDEX': Genome index of variant
                }, ...
            ],
            ?'MUTATIONS': [
                {
                    'MUTATION': Gene level mutation in GARC,
                    'GENE': Gene name,
                    'GENE_POSITION': Position within the gene. Amino acid or nucleotide index depending on which is appropriate
                }
            ],
            ?'EFFECTS': {
                Drug name: [
                    {
                        'GENE': Gene name of the mutation,
                        'MUTATION': Gene level mutation in GARC,
                        'PREDICTION': Prediction caused by this mutation
                    }, ...,
                    {
                        'PHENOTYPE': Resultant prediction for this drug based on prediciton heirarchy
                    }
                ], ...
            }
        }
    }
    Where fields with a preceeding '?' are not always present depending on data

    Args:
        path (str): Path to the directory where the variant/mutation/effect CSV files are saved. Also the output dir for this.
        guid (str): Sample GUID
        values (str): Prediction values for the resistance catalogue in priority order (values[0] is highest priority)
        gnomonVersion (str): Semantic versioning string for the gnomon module. Can be accessed by `gnomon.__version__`
    '''

    #Define some metadata for the json
    meta = {
        'version': gnomonVersion, #Gnomon version used
        'guid': guid, #Sample GUID
        'UTC-datetime-run': datetime.datetime.utcnow().isoformat(), #ISO datetime run
        'fields': dict() #Fields included. These vary according to existance of mutations/effects
        }
    data = {}
    #Variants field
    _variants = []
    meta['fields']['VARIANTS'] = ['VARIANT', 'NUCLEOTIDE_INDEX']
    for _, variant in variants.iterrows():
        row = {
            'VARIANT': variant['VARIANT'],
            'NUCLEOTIDE_INDEX': variant['NUCLEOTIDE_INDEX'],
        }
        _variants.append(row)
    data['VARIANTS'] = _variants

    #Depending on mutations/effects, populate
    _mutations = []
    if mutations is not None:
        meta['fields']['MUTATIONS'] = ['MUTATION', 'GENE', 'GENE_POSITION']
        for _, mutation in mutations.iterrows():
            row = {
                'MUTATION': mutation['MUTATION'],
                'GENE': mutation['GENE'],
                'GENE_POSITION': mutation['GENE_POSITION']
            }
            _mutations.append(row)
        data['MUTATIONS'] = _mutations

    _effects = defaultdict(list)
    if effects is not None:
        meta['fields']['EFFECTS'] = dict()
        for _, effect in effects.iterrows():
            prediction = {
                'GENE': effect['GENE'],
                'MUTATION': effect['MUTATION'],
                'PREDICTION': effect['PREDICTION']
            }
            _effects[effect['DRUG']].append(prediction)
        
        #Get the overall predictions for each drug
        for drug, predictions in _effects.items():
            phenotype = 'S'
            for prediction in predictions:
                #Use the prediction heierarchy to use most signifiant prediction
                if values.index(prediction['PREDICTION']) < values.index(phenotype):
                    #The prediction is closer to the start of the values list, so should take priority
                    phenotype = prediction['PREDICTION']
            _effects[drug].append({'PHENOTYPE': phenotype})
            meta['fields']['EFFECTS'][drug] = [['GENE', 'MUTATION', 'PREDICTION'], 'PHENOTYPE']
        data['EFFECTS'] = _effects

    #Convert fields to a list so it can be json serialised
    with open(os.path.join(path, 'gnomon-out.json'), 'w') as f:
        print(json.dumps({'meta': meta, 'data': data}, indent=2, sort_keys=True), file=f)

def toAltJSON(path: str, reference: gumpy.Genome, vcfStem: str, catalogue: str) -> None:
    '''Convert the output JSON into a similar format to the COVID workflows:
    {
        <guid>: {
            'WorkflowInformation': {
                'gnomonVersion': Version,
                'referenceIdentifier': ID of the reference,
                'sampleIdentifier': ID of the sample (VCF stem),
                'catalogueName': Name of the prediction catalogue
            },
            'Gnomon': {
                'aaDeletions': [GARC of aa deletions],
                'aaInsertions': [GARC of aa insertions],
                'aaSubsitutions': [GARC of aa SNPs],
                'deletions': [GARC of deletions], #Inclusive of AA deletions
                'insertions': [GARC of insertions], #Inclusive of AA insertions
                'substitutions': [GARC of SNPs], #Inclusive of AA SNPs
                'effects': {
                    <drug name>: [
                        {
                            'GENE': Gene name,
                            'MUTATION': Mutation in GARC,
                            'PREDICTION': Prediction caused by this mutation
                        }, ...,
                        {
                            'PHENOTYPE': Resultant prediction following catalogue heirarchy
                        }
                    ], ...
                }
            },
            'GnomonOutputJSON': Full original output JSON
        }
    }

    Args:
        path (str): Path to the directory containing the original JSON
        reference (gumpy.Genome): Reference gumpy Genome object
        vcfStem (str): Stem of the VCF file. Should be the sample GUID
        catalogue (str): Name of the catalogue
    '''
    original = json.load(open(os.path.join(path, 'gnomon-out.json'), 'r'))

    variants = [x['VARIANT'] for x in original['data']['VARIANTS']]
    #Only insertions of form <pos>_ins_<bases>
    insertions = [x for x in variants if "ins" in x and x[-1].isalpha()]
    #Only deletions of form <pos>_del_<bases>
    deletions = [x for x in variants if "del" in x and "indel" not in x and x[-1].isalpha()]
    #Should just be SNPs left
    snps = [x for x in variants if "ins" not in x and "del" not in x]

    mutations = [x['MUTATION'] for x in original['data'].get('MUTATIONS', {})]
    #Only insertions of form <gene>@<pos>_ins_<bases>
    aaInsertions = [x for x in mutations if "ins" in x and x[-1].isalpha()] 
    #Only deletions of form <gene>@<pos>_del_<bases>
    aaDeletions = [x for x in mutations if "del" in x and "indel" not in x and x[-1].isalpha()]
    #Should just be SNPs left
    aaSnps = [x for x in mutations if "ins" not in x and "del" not in x]

    effects = original['data'].get('EFFECTS', {})

    out = {
        vcfStem: {
            'WorkflowInformation': {
                'gnomonVersion': original['meta']['version'],
                'referenceIdentifier': reference.name,
                'sampleIdentifier': vcfStem,
                'catalogueName': catalogue
            },
            'Gnomon': {
                'aaDeletions': aaDeletions,
                'aaInsertions': aaInsertions,
                'aaSubsitutions': aaSnps,
                'deletions': deletions,
                'insertions': insertions,
                'substitutions': snps,
                'effects': effects
            },
            'GnomonOutputJSON': original
        }
    }

    with open(os.path.join(path, 'alt-gnomon-out.json'), 'w') as f:
        print(json.dumps(out, indent=2), file=f)
