'''Suite of unit tests matching the test cases defined in tests/NC_045512.2-README.md

Run from root of git dir with `pytest -vv`
'''
import json
import os
import shutil

import gumpy
import gzip
import pandas as pd
import pickle
import piezo
import pytest
#Helpful function import for testing nested JSON equality but does not handle [{}, {}] well
from recursive_diff import recursive_eq

import gnomonicus

'''
Due to complications testing equalities of nested jsons of lists/dicts, there is a lot of 
code specificially dedicated to ensuring the lists are in the same order (they differ due to
a mixture of dictionary behaviour and different positions within files). However, we only care that
the contents of the JSON is there for these tests rather than caring about order.
'''

def setupOutput(testID: str) -> None:
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

def concatFields(d: dict) -> str:
    '''Concat the value of a dictionary in the order of the sorted keys

    Args:
        d (dict): Dictionary input

    Returns:
        str: String of values
    '''
    return ''.join([str(d[key]) for key in sorted(list(d.keys()))])

def sortValues(json: dict) -> dict:
    '''Sort the values within the VARIANTS, MUTATIONS and EFFECTS lists in a JSON.
    THis allows to test for contents equality as order is not particularly important here

    Args:
        json (dict): JSON in

    Returns:
        dict: JSON with VARIANTS, MUTATIONS and EFFECTS lists in reproducable orders for equality checks
    '''
    variants = json['data']['VARIANTS']
    mutations = json['data'].get('MUTATIONS', None)
    effects = json['data'].get('EFFECTS', None)
    
    json['data']['VARIANTS'] = sorted(variants, key=concatFields)
    if mutations is not None:
        json['data']['MUTATIONS'] = sorted(mutations, key=concatFields)
    if effects is not None:
        for drug in effects.keys():
            #Don't retain the pandas NaN values - swap to None for equality checking
            for i in range(len(effects[drug])):
                if 'GENE' in effects[drug][i].keys():
                    if pd.isnull(effects[drug][i]['GENE']):
                        effects[drug][i]['GENE'] = None
            
            #Sort and add now we have consistent values
            json['data']['EFFECTS'][drug] = sorted(effects[drug], key=concatFields)
    
    return json

def test_misc():
    '''Testing misc things which should be raised/edge case behaviour not confined neatly to a whole flow test case
    '''
    setupOutput('0')
    #Ensure that there is not an existing pickle (for clean start and testing loadGenome)
    if os.path.exists("tests/test-cases/NC_045512.2.gbk.pkl"):
        os.remove("tests/test-cases/NC_045512.2.gbk.pkl")
    #Gzipped
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk.gz", False)

    #Pickled
    if os.path.exists("tests/test-cases/NC_045512.2.gbk.pkl"):
        os.remove("tests/test-cases/NC_045512.2.gbk.pkl")
    reference_ = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)

    assert reference == reference_

    #Last edge case of loadGenome is gbk.pkl but its a gzipped file:
    r = pickle.load(open("tests/test-cases/NC_045512.2.gbk.pkl", 'rb'))
    pickle.dump(r, gzip.open("tests/test-cases/reference.gbk.pkl", "wb"))

    reference_ = gnomonicus.loadGenome("tests/test-cases/reference.gbk/", False)

    assert reference == reference_

    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-S_E484K-minos.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-S_E484K-minos"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/0/"
    gnomonicus.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    
    #Check for differences if a catalogue is not given. Should be the same mutations but different referenceGenes
    mutations_, referenceGenes_ = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, None)
    
    assert mutations.equals(mutations_)
    assert referenceGenes != referenceGenes_
    
    #Trying to raise an InvalidMutationException with malformed mutation
    should_raise_error = pd.DataFrame({
                                        'UNIQUEID': ['a'], 'GENE': ['S'], 'MUTATION': ['aa'], 'NUCLEOTIDE_NUMBER': ['a'], 
                                        'NUCLEOTIDE_INDEX': ['a'], 'GENE_POSITION': ['a'], 'ALT': ['a'], 'REF': ['a'], 
                                        'CODES_PROTEIN': ['a'], 'INDEL_LENGTH': ['a'], 'INDEL_NUCLEOTIDES': ['a'], 
                                        'IS_CDS': ['a'], 'IS_HET': ['a'], 'IS_NULL': ['a'], 'IS_PROMOTER': ['a'], 
                                        'IS_SNP': ['a'], 'AMINO_ACID_NUMBER': ['a'], 'AMINO_ACID_SEQUENCE': ['a'], 
                                        'IS_SYNONYMOUS': ['a'], 'IS_NONSYNONYMOUS': ['a'], 'NUMBER_NUCLEOTIDE_CHANGES': ['a']
        })

    with pytest.raises(gnomonicus.InvalidMutationException):
        gnomonicus.populateEffects(path, catalogue, should_raise_error, referenceGenes, vcfStem)

    
    #Quick check that invalid values in `handleIndels` raise a NoVariantsException
    with pytest.raises(gnomonicus.NoVariantsException):
        gnomonicus.handleIndels({'no': 'vars'})
    


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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-S_E484K-minos.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-S_E484K-minos"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/1/"
    gnomonicus.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")

    assert variants['VARIANT'][0] == '23012g>a'

    assert mutations['GENE'][0] == 'S'
    assert mutations['MUTATION'][0] == 'E484K'

    assert 'AAA' in effects['DRUG'].to_list()
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'R'

    gnomonicus.saveJSON(variants, mutations, effects, path, vcfStem, catalogue.catalogue.values, gnomonicus.__version__)
    gnomonicus.toAltJSON(path, reference, vcfStem, catalogue)

    expectedJSON = {
        'meta': {
            'version': gnomonicus.__version__,
            'guid': vcfStem,
            'fields': {
                "EFFECTS": {
                        "AAA": [
                        [
                            "GENE",
                            "MUTATION",
                            "PREDICTION"
                        ],
                        "PHENOTYPE"
                        ]
                    },
                "MUTATIONS": [
                    "MUTATION",
                    "GENE",
                    "GENE_POSITION"
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX"
                    ]
            }
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': '23012g>a',
                    'NUCLEOTIDE_INDEX': 23012
                }
            ],
            'MUTATIONS': [
                {
                    'MUTATION': 'E484K',
                    'GENE': 'S',
                    'GENE_POSITION':484
                }
            ],
            'EFFECTS': {
                'AAA': [
                    {
                        'GENE': 'S',
                        'MUTATION': 'E484K',
                        'PREDICTION': 'R'
                    },
                    {
                        'PHENOTYPE': 'R'
                    }
                ],
            }
        }
    }

    #Ensure the same key ordering as actual by running through json dumping and loading
    strJSON = json.dumps(expectedJSON, indent=2, sort_keys=True)
    expectedJSON_ = sortValues(json.loads(strJSON))

    actualJSON = sortValues(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON_, actualJSON)

    #Try again with alt JSON format
    expectedJSON2 = {
        vcfStem: {
            'WorkflowInformation': {
                'gnomonicusVersion': gnomonicus.__version__,
                'referenceIdentifier': reference.name,
                'sampleIdentifier': vcfStem,
                'catalogueName': 'gnomonicus_test'
            },
            'gnomonicus': {
                'aaDeletions': [],
                'aaInsertions': [],
                'aaSubsitutions': ['S@E484K'],
                'deletions': [], 
                'insertions': [], 
                'substitutions': ['23012g>a'], 
                'frameshifts': 0,
                'effects': {
                        'AAA': [
                            {
                                'GENE': 'S',
                                'MUTATION': 'E484K',
                                'PREDICTION': 'R'
                            },
                            {
                                'PHENOTYPE': 'R'
                            }
                        ],
                    }
            },
            'gnomonicusOutputJSON': expectedJSON
        }
    }
    # strJSON2 = json.dumps(expectedJSON2, indent=2, sort_keys=True)

    actualJSON2 = json.load(open(os.path.join(path, f'{vcfStem}.alt-gnomonicus-out.json'), 'r'))
    #Remove datetime as this is unreplicable
    del actualJSON2[vcfStem]['gnomonicusOutputJSON']['meta']['UTC-datetime-run']

    recursive_eq(expectedJSON2, actualJSON2)



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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-S_E484K-samtools.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-S_E484K-samtools"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/2/"
    gnomonicus.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")

    assert variants['VARIANT'][0] == '23012g>a'

    assert mutations['GENE'][0] == 'S'
    assert mutations['MUTATION'][0] == 'E484K'

    assert 'AAA' in effects['DRUG'].to_list()
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'R'

    gnomonicus.saveJSON(variants, mutations, effects, path, vcfStem, catalogue.catalogue.values, gnomonicus.__version__)

    expectedJSON = {
        'meta': {
            'version': gnomonicus.__version__,
            'guid': vcfStem,
            'fields': {
                "EFFECTS": {
                        "AAA": [
                        [
                            "GENE",
                            "MUTATION",
                            "PREDICTION"
                        ],
                        "PHENOTYPE"
                        ]
                    },
                "MUTATIONS": [
                    "MUTATION",
                    "GENE",
                    "GENE_POSITION"
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX"
                    ]
            }
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': '23012g>a',
                    'NUCLEOTIDE_INDEX': 23012
                }
            ],
            'MUTATIONS': [
                {
                    'MUTATION': 'E484K',
                    'GENE': 'S',
                    'GENE_POSITION':484
                }
            ],
            'EFFECTS': {
                'AAA': [
                    {
                        'GENE': 'S',
                        'MUTATION': 'E484K',
                        'PREDICTION': 'R'
                    },
                    {
                        'PHENOTYPE': 'R'
                    }
                ],
            }
        }
    }

    #Ensure the same key ordering as actual by running through json dumping and loading
    strJSON = json.dumps(expectedJSON, indent=2, sort_keys=True)
    expectedJSON = sortValues(json.loads(strJSON))

    actualJSON = sortValues(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON, actualJSON)

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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-S_F2F-minos.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-S_F2F-minos"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/3/"
    gnomonicus.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")

    assert variants['VARIANT'][0] == '21568t>c'

    assert mutations['GENE'][0] == 'S'
    assert mutations['MUTATION'][0] == 'F2F'
    assert mutations['GENE'][1] == 'S'
    assert mutations['MUTATION'][1] == 't6c'

    assert 'AAA' in effects['DRUG'].to_list()
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'S'

    gnomonicus.saveJSON(variants, mutations, effects, path, vcfStem, catalogue.catalogue.values, gnomonicus.__version__)

    expectedJSON = {
        'meta': {
            'version': gnomonicus.__version__,
            'guid': vcfStem,
            'fields': {
                "EFFECTS": {
                        "AAA": [
                        [
                            "GENE",
                            "MUTATION",
                            "PREDICTION"
                        ],
                        "PHENOTYPE"
                        ]
                    },
                "MUTATIONS": [
                    "MUTATION",
                    "GENE",
                    "GENE_POSITION"
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX"
                    ]
            }
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': '21568t>c',
                    'NUCLEOTIDE_INDEX': 21568
                }
            ],
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
        }
    }

    #Ensure the same key ordering as actual by running through json dumping and loading
    strJSON = json.dumps(expectedJSON, indent=2, sort_keys=True)
    expectedJSON = sortValues(json.loads(strJSON))

    actualJSON = sortValues(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON, actualJSON)


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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-S_F2L-minos.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-S_F2L-minos"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/4/"
    gnomonicus.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")

    assert variants['VARIANT'][0] == '21566t>c'

    assert mutations['GENE'][0] == 'S'
    assert mutations['MUTATION'][0] == 'F2L'

    assert 'AAA' in effects['DRUG'].to_list()
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'U'

    gnomonicus.saveJSON(variants, mutations, effects, path, vcfStem, catalogue.catalogue.values, gnomonicus.__version__)

    expectedJSON = {
        'meta': {
            'version': gnomonicus.__version__,
            'guid': vcfStem,
            'fields': {
                "EFFECTS": {
                        "AAA": [
                        [
                            "GENE",
                            "MUTATION",
                            "PREDICTION"
                        ],
                        "PHENOTYPE"
                        ]
                    },
                "MUTATIONS": [
                    "MUTATION",
                    "GENE",
                    "GENE_POSITION"
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX"
                    ]
            }
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': '21566t>c',
                    'NUCLEOTIDE_INDEX': 21566
                }
            ],
            'MUTATIONS': [
                {
                    'MUTATION': 'F2L',
                    'GENE': 'S',
                    'GENE_POSITION':2
                }
            ],
            'EFFECTS': {
                'AAA': [
                    {
                        'GENE': 'S',
                        'MUTATION': 'F2L',
                        'PREDICTION': 'U'
                    },
                    {
                        'PHENOTYPE': 'U'
                    }
                ],
            }
        }
    }

    #Ensure the same key ordering as actual by running through json dumping and loading
    strJSON = json.dumps(expectedJSON, indent=2, sort_keys=True)
    expectedJSON = sortValues(json.loads(strJSON))

    actualJSON = sortValues(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON, actualJSON)    


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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-S_200_indel-minos.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-S_200_indel-minos"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/5/"
    gnomonicus.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")

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

    gnomonicus.saveJSON(variants, mutations, effects, path, vcfStem, catalogue.catalogue.values, gnomonicus.__version__)

    expectedJSON = {
        'meta': {
            'version': gnomonicus.__version__,
            'guid': vcfStem,
            'fields': {
                "EFFECTS": {
                        "AAA": [
                        [
                            "GENE",
                            "MUTATION",
                            "PREDICTION"
                        ],
                        "PHENOTYPE"
                        ]
                    },
                "MUTATIONS": [
                    "MUTATION",
                    "GENE",
                    "GENE_POSITION"
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX"
                    ]
            }
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': '21762_indel',
                    'NUCLEOTIDE_INDEX': 21762
                },
                {
                    'VARIANT': '21762_ins_1',
                    'NUCLEOTIDE_INDEX': 21762
                },
                {
                    'VARIANT': '21762_ins_c',
                    'NUCLEOTIDE_INDEX': 21762
                },
            ],
            'MUTATIONS': [
                {
                    'MUTATION': '200_indel',
                    'GENE': 'S',
                    'GENE_POSITION':200
                },
                {
                    'MUTATION': '200_ins_1',
                    'GENE': 'S',
                    'GENE_POSITION':200
                },
                {
                    'MUTATION': '200_ins_c',
                    'GENE': 'S',
                    'GENE_POSITION':200
                },
            ],
            'EFFECTS': {
                'AAA': [
                    {
                        'GENE': 'S',
                        'MUTATION': '200_indel',
                        'PREDICTION': 'U'
                    },
                    {
                        'GENE': 'S',
                        'MUTATION': '200_ins_1',
                        'PREDICTION': 'R'
                    },
                    {
                        'GENE': 'S',
                        'MUTATION': '200_ins_c',
                        'PREDICTION': 'R'
                    },
                    {
                        'PHENOTYPE': 'R'
                    }
                ],
            }
        }
    }

    #Ensure the same key ordering as actual by running through json dumping and loading
    strJSON = json.dumps(expectedJSON, indent=2, sort_keys=True)
    expectedJSON = sortValues(json.loads(strJSON))

    actualJSON = sortValues(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON, actualJSON)

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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-double-minos.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-double-minos"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/6/"
    gnomonicus.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")

    assert variants['VARIANT'][0] == '27758g>c'


    assert 'ORF7a' in mutations['GENE'].to_list()
    assert 'ORF7b' in mutations['GENE'].to_list()

    assert mutations['MUTATION'][mutations['GENE'].to_list().index('ORF7a')] == '!122S'
    assert mutations['MUTATION'][mutations['GENE'].to_list().index('ORF7b')] == 'M1I'

    assert 'AAA' in effects['DRUG'].to_list()
    assert 'BBB' in effects['DRUG'].to_list()
    
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'R'
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('BBB')] == 'R'

    gnomonicus.saveJSON(variants, mutations, effects, path, vcfStem, catalogue.catalogue.values, gnomonicus.__version__)

    expectedJSON = {
        'meta': {
            'version': gnomonicus.__version__,
            'guid': vcfStem,
            'fields': {
                "EFFECTS": {
                        "AAA": [
                        [
                            "GENE",
                            "MUTATION",
                            "PREDICTION"
                        ],
                        "PHENOTYPE"
                        ], 
                        "BBB": [
                        [
                            "GENE",
                            "MUTATION",
                            "PREDICTION"
                        ],
                        "PHENOTYPE"
                        ], 
                    },
                "MUTATIONS": [
                    "MUTATION",
                    "GENE",
                    "GENE_POSITION"
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX"
                    ]
            }
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': '27758g>c',
                    'NUCLEOTIDE_INDEX': 27758
                }
            ],
            'MUTATIONS': [
                {
                    'MUTATION': '!122S',
                    'GENE': 'ORF7a',
                    'GENE_POSITION': 122
                },
                {
                    'MUTATION': 'M1I',
                    'GENE': 'ORF7b',
                    'GENE_POSITION': 1
                },
            ],
            'EFFECTS': {
                'AAA': [
                    {
                        'GENE': 'ORF7a',
                        'MUTATION': '!122S',
                        'PREDICTION': 'R'
                    },
                    {
                        'PHENOTYPE': 'R'
                    }
                ],
                'BBB': [
                    {
                        'GENE': 'ORF7b',
                        'MUTATION': 'M1I',
                        'PREDICTION': 'R'
                    },
                    {
                        'PHENOTYPE': 'R'
                    }
                ],
            }
        }
    }

    #Ensure the same key ordering as actual by running through json dumping and loading
    strJSON = json.dumps(expectedJSON, indent=2, sort_keys=True)
    expectedJSON = sortValues(json.loads(strJSON))

    actualJSON = sortValues(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON, actualJSON)


def test_7():
    '''Input:
            NC_045512.2-S_E484K&1450_ins_a-minos.vcf
        Expect output:
            variants:    23012g>a, 23012_ins_a, 23012_indel, 23012_ins_1
            mutations:   S@E484K, S@1450_ins_a, S@1450_ins_a&S@E484K, S@1450_indel, S@1450_ins_1
            predictions: {'AAA': 'R', 'BBB': 'R'}
    '''
    #Setup
    setupOutput('7')
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-S_E484K&1450_ins_a-minos.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-S_E484K&1450_ins_a-minos"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/7/"
    gnomonicus.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")

    #Sort the variants for comparing
    variants_ = sorted(variants['VARIANT'])
    assert variants_[0] == '23012_indel'
    assert variants_[1] == '23012_ins_1'
    assert variants_[2] == '23012_ins_a'
    assert variants_[3] == '23012g>a'

    #Sort the mutations for comparing
    mutations_ = sorted(list(zip(mutations['GENE'], mutations['MUTATION'])), key= lambda x: x[0] + x[1] if x[0] is not None else x[1])

    assert mutations_[0][0] == 'S'
    assert mutations_[0][1] == '1450_indel'

    assert mutations_[1][0] == 'S'
    assert mutations_[1][1] == '1450_ins_1'

    assert mutations_[2][0] == 'S'
    assert mutations_[2][1] == '1450_ins_a'

    assert mutations_[3][0] == 'S'
    assert mutations_[3][1] == 'E484K'


    #Sort the effects for comparing
    effects_ = [i[1] for i in sorted([(str(e), e) for _, e in effects.iterrows()], key=lambda x: x[0])]

    #Expected effects. For each row, x[0] = DRUG, x[1] = GENE, x[2] = MUTATION, x[3] = PREDICTION
    expected = [
        ['AAA', 'S', 'E484K', 'R'],
        ['AAA', 'S', '1450_indel', 'U'],
        ['AAA', 'S', '1450_ins_1', 'R'],
        ['AAA', 'S', '1450_ins_a', 'R'],
        ['BBB', None, 'S@1450_ins_a&S@E484K', 'R'],
    ]
    assert len(expected) == len(effects_)
    #Iter expected and effects to check for equality
    for row, exp in zip(effects_, expected):
        assert row['DRUG'] == exp[0]

        #Dealing with pd.nan as equality doesn't work here...
        if pd.isnull(row['GENE']):
            assert exp[1] is None
        else:
            assert row['GENE'] == exp[1]

        assert row['MUTATION'] == exp[2]
        assert row['PREDICTION'] == exp[3]

    gnomonicus.saveJSON(variants, mutations, effects, path, vcfStem, catalogue.catalogue.values, gnomonicus.__version__)

    expectedJSON = {
        'meta': {
            'version': gnomonicus.__version__,
            'guid': vcfStem,
            'fields': {
                "EFFECTS": {
                        "AAA": [
                        [
                            "GENE",
                            "MUTATION",
                            "PREDICTION"
                        ],
                        "PHENOTYPE"
                        ],
                        "BBB": [
                        [
                            "GENE",
                            "MUTATION",
                            "PREDICTION"
                        ],
                        "PHENOTYPE"
                        ]
                    },
                "MUTATIONS": [
                    "MUTATION",
                    "GENE",
                    "GENE_POSITION"
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX"
                    ]
            }
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': '23012_indel',
                    'NUCLEOTIDE_INDEX': 23012
                },
                {
                    'VARIANT': '23012_ins_1',
                    'NUCLEOTIDE_INDEX': 23012
                },
                {
                    'VARIANT': '23012_ins_a',
                    'NUCLEOTIDE_INDEX': 23012
                },
                {
                    'VARIANT': '23012g>a',
                    'NUCLEOTIDE_INDEX': 23012
                }
            ],
            'MUTATIONS': [
                {
                    'MUTATION': '1450_indel',
                    'GENE': 'S',
                    'GENE_POSITION':1450
                },
                {
                    'MUTATION': '1450_ins_1',
                    'GENE': 'S',
                    'GENE_POSITION':1450
                },
                {
                    'MUTATION': '1450_ins_a',
                    'GENE': 'S',
                    'GENE_POSITION':1450
                },
                {
                    'MUTATION': 'E484K',
                    'GENE': 'S',
                    'GENE_POSITION':484
                },
            ],
            'EFFECTS': {
                'AAA': [
                    {
                        'GENE': 'S',
                        'MUTATION': 'E484K',
                        'PREDICTION': 'R'
                    },
                    {
                        'GENE': 'S',
                        'MUTATION': '1450_indel',
                        'PREDICTION': 'U'
                    },
                    {
                        'GENE': 'S',
                        'MUTATION': '1450_ins_1',
                        'PREDICTION': 'R'
                    },
                    {
                        'GENE': 'S',
                        'MUTATION': '1450_ins_a',
                        'PREDICTION': 'R'
                    },
                    {
                        'PHENOTYPE': 'R'
                    }
                ],
                'BBB': [
                    {
                        'GENE': None,
                        'MUTATION': 'S@1450_ins_a&S@E484K',
                        'PREDICTION': 'R'
                    },
                    {
                        'PHENOTYPE': 'R'
                    }
                ],
            }
        }
    }

    #Ensure the same key ordering as actual by running through json dumping and loading
    strJSON = json.dumps(expectedJSON, indent=2, sort_keys=True)
    expectedJSON_ = sortValues(json.loads(strJSON))

    actualJSON = sortValues(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON_, actualJSON)


def test_8():
    '''Force the "S" gene to be non-coding to test non-coding behaviours
        Input:
            NC_045512.2-S_E484K-minos.vcf
        Expect output:
            variants:    23012g>a
            mutations:   g1450a
            predictions: {'AAA': 'U'}
    '''
    #Setup
    setupOutput('8')
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    #Force non coding (as all genes included are coding)
    reference.genes['S']['codes_protein'] = False

    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile("tests/test-cases/NC_045512.2-S_E484K-minos.vcf", ignore_filter=True, bypass_reference_calls=True)
    vcfStem = "NC_045512.2-S_E484K-minos"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/8/"
    gnomonicus.populateVariants(vcfStem, path, diff)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")

    assert variants['VARIANT'][0] == '23012g>a'

    assert mutations['GENE'][0] == 'S'
    assert mutations['MUTATION'][0] == 'g1450a'

    assert 'AAA' in effects['DRUG'].to_list()
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'U'

    gnomonicus.saveJSON(variants, mutations, effects, path, vcfStem, catalogue.catalogue.values, gnomonicus.__version__)
    gnomonicus.toAltJSON(path, reference, vcfStem, catalogue)

    expectedJSON = {
        'meta': {
            'version': gnomonicus.__version__,
            'guid': vcfStem,
            'fields': {
                "EFFECTS": {
                        "AAA": [
                        [
                            "GENE",
                            "MUTATION",
                            "PREDICTION"
                        ],
                        "PHENOTYPE"
                        ]
                    },
                "MUTATIONS": [
                    "MUTATION",
                    "GENE",
                    "GENE_POSITION"
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX"
                    ]
            }
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': '23012g>a',
                    'NUCLEOTIDE_INDEX': 23012
                }
            ],
            'MUTATIONS': [
                {
                    'MUTATION': 'g1450a',
                    'GENE': 'S',
                    'GENE_POSITION':1450
                }
            ],
            'EFFECTS': {
                'AAA': [
                    {
                        'GENE': 'S',
                        'MUTATION': 'g1450a',
                        'PREDICTION': 'U'
                    },
                    {
                        'PHENOTYPE': 'U'
                    }
                ],
            }
        }
    }

    #Ensure the same key ordering as actual by running through json dumping and loading
    strJSON = json.dumps(expectedJSON, indent=2, sort_keys=True)
    expectedJSON_ = sortValues(json.loads(strJSON))

    actualJSON = sortValues(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON_, actualJSON)