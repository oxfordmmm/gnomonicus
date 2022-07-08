import pytest
import os
import shutil

import piezo
import gumpy
import pandas as pd
import gnomon


def setupOutput(testID: str):
    '''Ensure that the output folder exists and is empty in preparation for a test

    Args:
        testID (str): Desired folder name for this test
    '''
    path = f"tests/outputs/{testID}"
    #Make the dir if it doesn't exist
    os.makedirs(path, exist_ok=True)

    #Check for contents
    if len(os.listdir(path)) > 0:
        #Not empty, so delete and recreate
        shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)

def test_1():
    '''Input:
            NC_045512.2-S_E484K-minos.vcf
        Expect output:
            variants:    23012g>a
            mutations:   S@E484K
            predictions: {'AAA': 'R', 'BBB': 'S'}
    '''
    #Setup
    setupOutput('1')
    reference = gnomon.loadGenome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-S_E484K-minos.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-S_E484K-minos"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/1/"
    gnomon.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomon.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomon.populateEffects(sample, path, catalogue, mutations, referenceGenes)

    #Check for expected values within csvs
    variants = pd.read_csv(path+"variants.csv")
    mutations = pd.read_csv(path+"mutations.csv")
    effects = pd.read_csv(path+"effects.csv")

    assert variants['VARIANT'][0] == '23012g>a'

    assert mutations['GENE'][0] == 'S'
    assert mutations['MUTATION'][0] == 'E484K'

    assert 'AAA' in effects['DRUG'].to_list()
    # assert 'BBB' in effects['DRUG'].to_list()
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'R'
    # assert effects.loc(effects['DRUG'] == 'BBB')[0] == 'S'    

def test_2():
    '''Input:
            NC_045512.2-S_E484K-samtools.vcf
        Expect output:
            variants:    23012g>a
            mutations:   S@E484K
            predictions: {'AAA': 'R', 'BBB': 'S'}
    '''
    #Setup
    setupOutput('2')
    reference = gnomon.loadGenome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-S_E484K-samtools.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-S_E484K-samtools"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/2/"
    gnomon.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomon.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomon.populateEffects(sample, path, catalogue, mutations, referenceGenes)

    #Check for expected values within csvs
    variants = pd.read_csv(path+"variants.csv")
    mutations = pd.read_csv(path+"mutations.csv")
    effects = pd.read_csv(path+"effects.csv")

    assert variants['VARIANT'][0] == '23012g>a'

    assert mutations['GENE'][0] == 'S'
    assert mutations['MUTATION'][0] == 'E484K'

    assert 'AAA' in effects['DRUG'].to_list()
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'R'

def test_3():
    '''Input:
            NC_045512.2-S_F2F-minos.vcf
        Expect output:
            variants:    21568t>c
            mutations:   S@F2F
            predictions: {'AAA': 'S', 'BBB': 'S'}
    '''
    #Setup
    setupOutput('3')
    reference = gnomon.loadGenome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-S_F2F-minos.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-S_F2F-minos"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/3/"
    gnomon.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomon.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomon.populateEffects(sample, path, catalogue, mutations, referenceGenes)

    #Check for expected values within csvs
    variants = pd.read_csv(path+"variants.csv")
    mutations = pd.read_csv(path+"mutations.csv")
    effects = pd.read_csv(path+"effects.csv")

    assert variants['VARIANT'][0] == '21568t>c'

    assert mutations['GENE'][0] == 'S'
    assert mutations['MUTATION'][0] == 'F2F' #Known to fail with PyPi v1.0.0 gumpy due to oversight. Try with git version

    assert 'AAA' in effects['DRUG'].to_list()
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'S'


def test_4():
    '''Input:
            NC_045512.2-S_F2L-minos.vcf
        Expect output:
            variants:    21566t>c
            mutations:   S@F2L
            predictions: {'AAA': 'U', 'BBB': 'S'}
    '''
    #Setup
    setupOutput('4')
    reference = gnomon.loadGenome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-S_F2L-minos.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-S_F2L-minos"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/4/"
    gnomon.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomon.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomon.populateEffects(sample, path, catalogue, mutations, referenceGenes)

    #Check for expected values within csvs
    variants = pd.read_csv(path+"variants.csv")
    mutations = pd.read_csv(path+"mutations.csv")
    effects = pd.read_csv(path+"effects.csv")

    assert variants['VARIANT'][0] == '21566t>c'

    assert mutations['GENE'][0] == 'S'
    assert mutations['MUTATION'][0] == 'F2L'

    assert 'AAA' in effects['DRUG'].to_list()
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'U'


def test_5():
    '''Input:
            NC_045512.2-S_200_indel-minos.vcf
        Expect output:
            variants:    21762_indel, 21762_ins_1, 21762_ins_c
            mutations:   S@200_indel, S@200_ins_1, S@200_ins_c
            predictions: {'AAA': 'R', 'BBB': 'S'}
    '''
    #Setup
    setupOutput('5')
    reference = gnomon.loadGenome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-S_200_indel-minos.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-S_200_indel-minos"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/5/"
    gnomon.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomon.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomon.populateEffects(sample, path, catalogue, mutations, referenceGenes)

    #Check for expected values within csvs
    variants = pd.read_csv(path+"variants.csv")
    mutations = pd.read_csv(path+"mutations.csv")
    effects = pd.read_csv(path+"effects.csv")

    variantGARC = variants['VARIANT'].to_list()
    assert '21762_indel' in variantGARC
    assert '21762_ins_1' in variantGARC
    assert '21762_ins_c' in variantGARC

    mutationGenes = mutations['GENE'].to_list()
    for gene in mutationGenes:
        assert gene == 'S'

    mutationGARC = mutations['MUTATION'].to_list()
    assert '200_indel' in mutationGARC
    assert '200_ins_1' in mutationGARC
    assert '200_ins_c' in mutationGARC

    assert 'AAA' in effects['DRUG'].to_list()
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'R'

def test_6():
    '''Input:
            NC_045512.2-double-minos.vcf
        Expect output:
            variants:    27758g>c
            mutations:   ORF7a!122S, ORF7b@M1I
            predictions: {'AAA': 'R', 'BBB': 'R'}
    '''
    #Setup
    setupOutput('6')
    reference = gnomon.loadGenome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-double-minos.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-double-minos"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/6/"
    gnomon.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomon.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomon.populateEffects(sample, path, catalogue, mutations, referenceGenes)

    #Check for expected values within csvs
    variants = pd.read_csv(path+"variants.csv")
    mutations = pd.read_csv(path+"mutations.csv")
    effects = pd.read_csv(path+"effects.csv")

    assert variants['VARIANT'][0] == '27758g>c'


    assert 'ORF7a' in mutations['GENE'].to_list()
    assert 'ORF7b' in mutations['GENE'].to_list()

    assert mutations['MUTATION'][mutations['GENE'].to_list().index('ORF7a')] == '!122S'
    assert mutations['MUTATION'][mutations['GENE'].to_list().index('ORF7b')] == 'M1I'

    assert 'AAA' in effects['DRUG'].to_list()
    assert 'BBB' in effects['DRUG'].to_list()
    
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'R'
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('BBB')] == 'R'
