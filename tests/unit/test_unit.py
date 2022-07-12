'''Suite of unit tests matching the test cases defined in tests/NC_045512.2-README.md

Run from root of git dir with `pytest -vv`
'''
import os
import shutil
import json

import pandas as pd
#Helpful function import for testing nested JSON equality but does not handle [{}, {}] well
from recursive_diff import recursive_eq

import gumpy
import piezo
import pytest

import gnomon

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
            json['data']['EFFECTS'][drug] = sorted(effects[drug], key=concatFields)
    
    return json

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
    reference = gnomon.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
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
    assert effects['PREDICTION'][effects['DRUG'].to_list().index('AAA')] == 'R'

    gnomon.saveJSON(path, vcfStem, catalogue.catalogue.values, gnomon.__version__)

    expectedJSON = {
        'meta': {
            'version': gnomon.__version__,
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

    actualJSON = sortValues(json.load(open(os.path.join(path, 'gnomon-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON, actualJSON)

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
    reference = gnomon.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
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

    gnomon.saveJSON(path, vcfStem, catalogue.catalogue.values, gnomon.__version__)

    expectedJSON = {
        'meta': {
            'version': gnomon.__version__,
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

    actualJSON = sortValues(json.load(open(os.path.join(path, 'gnomon-out.json'), 'r')))
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
    reference = gnomon.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
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

    gnomon.saveJSON(path, vcfStem, catalogue.catalogue.values, gnomon.__version__)

    expectedJSON = {
        'meta': {
            'version': gnomon.__version__,
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
                    'GENE_POSITION':2
                }
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

    actualJSON = sortValues(json.load(open(os.path.join(path, 'gnomon-out.json'), 'r')))
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
    reference = gnomon.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
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

    gnomon.saveJSON(path, vcfStem, catalogue.catalogue.values, gnomon.__version__)

    expectedJSON = {
        'meta': {
            'version': gnomon.__version__,
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

    actualJSON = sortValues(json.load(open(os.path.join(path, 'gnomon-out.json'), 'r')))
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
    reference = gnomon.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
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

    gnomon.saveJSON(path, vcfStem, catalogue.catalogue.values, gnomon.__version__)

    expectedJSON = {
        'meta': {
            'version': gnomon.__version__,
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

    actualJSON = sortValues(json.load(open(os.path.join(path, 'gnomon-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    # for diff in recursive_diff(expectedJSON, actualJSON):
    #     print(diff)

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
    reference = gnomon.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
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

    gnomon.saveJSON(path, vcfStem, catalogue.catalogue.values, gnomon.__version__)

    expectedJSON = {
        'meta': {
            'version': gnomon.__version__,
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

    actualJSON = sortValues(json.load(open(os.path.join(path, 'gnomon-out.json'), 'r')))
    #Remove datetime as this is unreplicable
    del actualJSON['meta']['UTC-datetime-run']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON, actualJSON)
