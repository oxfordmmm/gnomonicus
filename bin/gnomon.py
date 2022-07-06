'''Gnomon is a script which pulls together output VCF of the Lodestone TB pipeline
    with a reference genome and a resistance catalogue, and utilises gumpy and 
    piezo to produce variants, mutations and an antibiogram.

Based on sp3predict
'''
import argparse
import logging
import os
import pickle

import gumpy
import pandas as pd
import numpy as np
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


def loadGenome(path: str) -> gumpy.Genome:
    '''Load a genome from a given path. Checks if path is to a pickle dump, or if a pickle dump of the path's file exists
    Instanciates a new gumpy Genome and dumps to pickle as a last resort

    Args:
        path (str): Path to the genbank file or pickle dump. If previously run, a genbank file's Genome object is pickled and dumped for speed

    Returns:
        gumpy.Genome: Genome object of the reference genome
    '''
    #Remove trailing '/' if required
    if path[-1] == '/':
        path = path[:-1]

    #Try to load as a pickle
    try:
        return pickle.load(open(path, "rb"))
    except:
        logging.info("Genome object not a pickle, checking if pickled version exists")

    #Try pickled version created by this (path+'.pkl')
    try:
        return pickle.load(open(path+'.pkl', 'rb'))
    except:
        logging.info("No pickled version of genome object, instanciating and dumping")
    
    #Create new gumpy.Genome and pickle dump for speed later
    reference = gumpy.Genome(path, show_progress_bar=True)
    pickle.dump(reference, open(path+'.pkl', 'wb'))
    return reference

def populateVariants(vcfStem: str, outputDir: str, diff: gumpy.GenomeDifference) -> None:
    '''Populate and save the variants DataFrame as a CSV

    Args:
        vcfStem (str): The stem of the filename for the VCF file. Used as a uniqueID
        outputDir (str): Path to the desired output directory
        diff (gumpy.GenomeDifference): GenomeDifference object between reference and the sample
    '''
    #Populate variants table directly from GenomeDifference
    variants = pd.DataFrame({
            'VARIANT': diff.variants, 
            'NUCLEOTIDE_INDEX': diff.nucleotide_index,
            'IS_INDEL': diff.is_indel,
            'IS_NULL': diff.is_null,
            'IS_HET': diff.is_het,
            'IS_SNP': diff.is_snp,
            'INDEL_LENGTH': diff.indel_length,
            'INDEL_NUCLEOTIDES': diff.indel_nucleotides
            })

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

def populateMutations(vcfStem: str, outputDir: str, diff: gumpy.GenomeDifference, reference: gumpy.Genome, sample: gumpy.Genome, resistanceCatalogue: piezo.ResistanceCatalogue) -> (pd.DataFrame, dict):
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
            geneMutations = pd.DataFrame({
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
                })
            #Add gene name
            geneMutations['GENE'] = gene

            #Add this gene's mutations to the total dataframe
            if not geneMutations.empty:
                if mutations:
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
        return numpy.sum([i!=j for (i,j) in zip(row['REF'],row['ALT'] ) ])
    else:
        return 0

def populateEffects(sample: gumpy.Genome, outputDir: str, resistanceCatalogue: piezo.ResistanceCatalogue, mutations: pd.DataFrame, referenceGenes: dict) -> dict:
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
        dict: A metadata dictionary mapping drugs to their predictions
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
    effects = pandas.DataFrame.from_dict(effects, 
                                        orient="index", 
                                        columns=["UNIQUEID", "GENE", "MUTATION", 
                                            "CATALOGUE_NAME", "DRUG", "PREDICTION"]
                                        )
    
    #Set the index
    effects.set_index(["UNIQUEID", "DRUG", "GENE", "MUTATION", "CATALOGUE_NAME"], inplace=True)
    
    #Save as CSV
    effects.to_csv(os.path.join(outputDir, 'effects.csv'))

    #Return  the metadata dict to log later
    return {"WGS_PREDICTION_"+drug: phenotype[drug] for drug in resistanceCatalogue.catalogue.drugs}

                

    


if __name__ == "__main__":
    #Argparser setup
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf_file",required=True,help="the path to a single VCF file")
    parser.add_argument("--genome_object",default="H37Rv_3.gbk",help="the path to a compressed gumpy Genome object or a genbank file")
    parser.add_argument("--catalogue_file",default=None,required=False,help="the path to the resistance catalogue")
    parser.add_argument("--ignore_vcf_status",action='store_true',default=False,help="whether to ignore the STATUS field in the vcf (e.g. necessary for some versions of Clockwork VCFs)")
    parser.add_argument("--ignore_vcf_filter",action='store_true',default=False,help="whether to ignore the FILTER field in the vcf (e.g. necessary for some versions of Clockwork VCFs)")
    parser.add_argument("--progress",action='store_true',default=False,help="whether to show progress using tqdm")
    parser.add_argument("--output_dir", action="store_true", required=False, default=".", help="Directory to save output files to. Defaults to wherever the script is run from.")
    options = parser.parse_args()

    #Make the output directory if it doesn't already exist
    os.makedirs(options.output_dir, exist_ok=True)

    #Logging setup
    logging.basicConfig(filename=os.path.join(options.output_dir, 'gnomon.log'), filemode='w', format='%(asctime)s -  %(levelname)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.DEBUG)
    logging.info(f"Gnomon starting with output directory {options.output_dir}")

    #Get reference genome
    reference = loadGenome(options.genome_object)
    logging.debug("Loaded reference genome")

    #Build the mutated genome using gumpy
    vcf = gumpy.VCFFile(options.vcf_file)
    sample = reference + vcf
    logging.debug("Applied the VCF to the reference")

    #Get the stem of the VCF filename for use as a unique ID
    vcfStem = os.path.split(options.vcf_file)[-1].split(".")[0]
    logging.debug("Got VCF stem")

    #Get resistance catalogue
    if options.catalogue_file:
        resistanceCatalogue = piezo.ResistanceCatalogue(options.catalogue_file, prediction_subset_only=True)
        logging.debug("Loaded resistance catalogue")
    else:
        resistanceCatalogue = None
        logging.info("No resistance catalogue provided, producing variants and mutations only")

    #Get the GenomeDifference for extracting genome level mutations
    diff = reference - sample
    logging.debug("Got the genome difference")

    #Complain if there are no variants
    if diff.variants is None:
        logging.error("No variants detected!")
        raise NoVariantsException()

    #Get the variations and mutations
    populateVariants(vcfStem, options.output_dir, diff)
    logging.debug("Populated and saved variants.csv")

    mutations, referenceGenes = populateMutations(vcfStem, options.output_dir, diff, 
                        reference, sample, resistanceCatalogue)
    if not mutations:
        logging.info("No mutations found - probably due to exclusively inter-gene variation or no variation.\n\t\t\t\t\t\t\t No effects.csv written")
    else:
        logging.debug("Populated and saved mutatons.csv")

    #Get the effects of the mutations
    if resistanceCatalogue and mutations:
        metadata = populateEffects(sample, options.output_dir, resistanceCatalogue, mutations, referenceGenes)
        logging.debug("Populated and saved effects.csv")
    else:
        metadata = {}
        logging.info("Skipped effects.csv due to lack of resistance catalogue or mutations")

    #Add data to the log
    logging.info("********** Successfully completed **********")
    
    logging.info(f"VCF file: {options.vcf_file}")
    logging.info(f"Reference genome file: {options.genome_object}")

    if resistanceCatalogue:
        logging.info(f"Catalogue reference genome: {resistanceCatalogue.catalogue.genbank_reference}")
        logging.info(f"Catalogue name: {resistanceCatalogue.catalogue.name}")
        logging.info(f"Catalogue version: {resistanceCatalogue.catalogue.version}")
        logging.info(f"Catalogue grammar: {resistanceCatalogue.catalogue.grammar}")
        logging.info(f"Catalogue values: {resistanceCatalogue.catalogue.values}")
        logging.info(f"Catalogue path: {options.catalogue_file}")
    for drug in sorted(metadata.keys()):
        logging.info(f"{drug} {metadata[drug]}")