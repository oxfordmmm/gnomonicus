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
#Helpful function import for testing nested JSON equality as it gives exact differences
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

def prep_json(j: dict) -> dict:
    '''Prepare a JSON for comparison by removing fields which cannot be reproduced

    Args:
        j (dict): Initial JSON

    Returns:
        dict: JSON without fields such as time and file paths
    '''
    del j['meta']['time_taken_s']
    del j['meta']['UTC-datetime-completed']
    del j['meta']['catalogue_file']
    del j['meta']['reference_file']
    del j['meta']['vcf_file']
    return j

def ordered(obj):
    '''Recursively sort a JSON for equality checking. Based on https://stackoverflow.com/questions/25851183/how-to-compare-two-json-objects-with-the-same-elements-in-a-different-order-equa

    Args:
        obj (object): Any JSON element. Probably one of dict, list, tuple, str, int, None

    Returns:
        object: Sorted JSON
    '''
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list) or isinstance(obj, tuple):
        return sorted(ordered(x) for x in obj)
    if type(obj) == type(None):
        return str(obj)
    else:
        return obj




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
    gnomonicus.populateVariants(vcfStem, path, diff, False)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue, False)
    
    #Check for differences if a catalogue is not given. Should be the same mutations but different referenceGenes
    mutations_, referenceGenes_ = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, None, False)
    
    assert mutations.equals(mutations_)
    assert referenceGenes != referenceGenes_
    
    #Trying to raise an InvalidMutationException with malformed mutation
    should_raise_error = pd.DataFrame({
                                        'uniqueid': ['a'], 'gene': ['S'], 'mutation': ['aa'], 'nucleotide_number': ['a'], 
                                        'nucleotide_index': ['a'], 'gene_position': ['a'], 'alt': ['a'], 'ref': ['a'], 
                                        'codes_protein': ['a'], 'indel_length': ['a'], 'indel_nucleotides': ['a'], 
                                        'amino_acid_number': ['a'], 'amino_acid_sequence': ['a'], 
                                        'number_nucleotide_changes': ['a']
        })

    with pytest.raises(gnomonicus.InvalidMutationException):
        gnomonicus.populateEffects(path, catalogue, should_raise_error, referenceGenes, vcfStem, False, False)

    

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
    gnomonicus.populateVariants(vcfStem, path, diff, True)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue, True)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem, True, True)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert variants['variant'][0] == '23012g>a'

    assert mutations['gene'][0] == 'S'
    assert mutations['mutation'][0] == 'E484K'

    assert 'AAA' in effects['drug'].to_list()
    assert effects['prediction'][effects['drug'].to_list().index('AAA')] == 'R'
    
    hits = []
    for _, row in predictions.iterrows():
        assert row['catalogue_name'] == 'gnomonicus_test'
        assert row['catalogue_version'] == 'v1.0'
        assert row['catalogue_values'] == 'RFUS'
        assert row['evidence'] == '{}'
        if row['drug'] == 'AAA':
            hits.append('AAA')
            assert row['prediction'] == "R"
        elif row['drug'] == 'BBB':
            hits.append('BBB')
            assert row['prediction'] == 'S'
        else:
            hits.append(None)
    assert sorted(hits) == ['AAA', 'BBB']

    gnomonicus.saveJSON(variants, mutations, effects, path, vcfStem, catalogue, gnomonicus.__version__, -1, reference, '', '', '')

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '23012g>a',
                    'nucleotide_index': 23012,
                    'gene_name': 'S',
                    'gene_position': 484,
                    'codon_idx': 0,
                    'vcf_evidence': {
                        'GT': [1, 1], 'DP': 44, 'DPF': 0.991, 'COV': [0, 44], 
                        'FRS': 1.0, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 
                        'REF': 'g', 'ALTS': ['a'], 'POS': 23012
                    },
                    'vcf_idx': 1
                }
            ],
            'mutations': [
                {
                    'mutation': 'E484K',
                    'gene': 'S',
                    'gene_position':484,
                    "ref": "gaa",
                    "alt": "aaa"
                }
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'S',
                        'mutation': 'E484K',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'R',
                'BBB': 'S'
            }
        }
    }
    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))



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
    gnomonicus.populateVariants(vcfStem, path, diff, True)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue, True)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem, True, True)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert variants['variant'][0] == '23012g>a'

    assert mutations['gene'][0] == 'S'
    assert mutations['mutation'][0] == 'E484K'

    assert 'AAA' in effects['drug'].to_list()
    assert effects['prediction'][effects['drug'].to_list().index('AAA')] == 'R'

    hits = []
    for _, row in predictions.iterrows():
        assert row['catalogue_name'] == 'gnomonicus_test'
        assert row['catalogue_version'] == 'v1.0'
        assert row['catalogue_values'] == 'RFUS'
        assert row['evidence'] == '{}'
        if row['drug'] == 'AAA':
            hits.append('AAA')
            assert row['prediction'] == "R"
        elif row['drug'] == 'BBB':
            hits.append('BBB')
            assert row['prediction'] == 'S'
        else:
            hits.append(None)
    assert sorted(hits) == ['AAA', 'BBB']

    gnomonicus.saveJSON(variants, mutations, effects, path, vcfStem, catalogue, gnomonicus.__version__, -1, reference, '', '', '')

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '23012g>a',
                    'nucleotide_index': 23012,
                    'gene_name': 'S',
                    'gene_position': 484,
                    'codon_idx': 0,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "PL": [
                            255,
                            33,
                            0
                        ],
                        "POS": 23012,
                        "REF": "g",
                        "ALTS": [
                            "a"
                        ]
                    },
                    'vcf_idx': 1
                }
            ],
            'mutations': [
                {
                    'mutation': 'E484K',
                    'gene': 'S',
                    'gene_position':484,
                    "ref": "gaa",
                    "alt": "aaa"
                }
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'S',
                        'mutation': 'E484K',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'R',
                'BBB': 'S'
            }
        }
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))

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
    gnomonicus.populateVariants(vcfStem, path, diff, True)
    mutations_, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue, True)
    gnomonicus.populateEffects(path, catalogue, mutations_, referenceGenes, vcfStem, True, True)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert variants['variant'][0] == '21568t>c'

    assert mutations['gene'][0] == 'S'
    assert mutations['mutation'][0] == 'F2F'

    assert 'AAA' in effects['drug'].to_list()
    assert effects['prediction'][effects['drug'].to_list().index('AAA')] == 'S'

    hits = []
    for _, row in predictions.iterrows():
        assert row['catalogue_name'] == 'gnomonicus_test'
        assert row['catalogue_version'] == 'v1.0'
        assert row['catalogue_values'] == 'RFUS'
        assert row['evidence'] == '{}'
        if row['drug'] == 'AAA':
            hits.append('AAA')
            assert row['prediction'] == "S"
        elif row['drug'] == 'BBB':
            hits.append('BBB')
            assert row['prediction'] == 'S'
        else:
            hits.append(None)
    assert sorted(hits) == ['AAA', 'BBB']

    gnomonicus.saveJSON(variants, mutations_, effects, path, vcfStem, catalogue, gnomonicus.__version__, -1, reference, '', '', '')

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '21568t>c',
                    'nucleotide_index': 21568,
                    'gene_name': 'S',
                    'gene_position': 2,
                    'codon_idx': 2,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            0,
                            44
                        ],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21568,
                        "REF": "t",
                        "ALTS": [
                            "c"
                        ]
                    },
                    'vcf_idx': 1
                }
            ],
            'mutations': [
                {
                    'mutation': 'F2F',
                    'gene': 'S',
                    'gene_position': 2,
                    'ref': 'ttt',
                    'alt': 'ttc'
                },
                {
                    'mutation': 't6c',
                    'gene': 'S',
                    'gene_position': 6
                },
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'S',
                        'mutation': 'F2F',
                        'prediction': 'S',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'S'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'S',
                'BBB': 'S'
            }
        }
    }


    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


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
    gnomonicus.populateVariants(vcfStem, path, diff, True)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue, True)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem, True, True)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert variants['variant'][0] == '21566t>c'

    assert mutations['gene'][0] == 'S'
    assert mutations['mutation'][0] == 'F2L'

    assert 'AAA' in effects['drug'].to_list()
    assert effects['prediction'][effects['drug'].to_list().index('AAA')] == 'U'

    hits = []
    for _, row in predictions.iterrows():
        assert row['catalogue_name'] == 'gnomonicus_test'
        assert row['catalogue_version'] == 'v1.0'
        assert row['catalogue_values'] == 'RFUS'
        assert row['evidence'] == '{}'
        if row['drug'] == 'AAA':
            hits.append('AAA')
            assert row['prediction'] == "U"
        elif row['drug'] == 'BBB':
            hits.append('BBB')
            assert row['prediction'] == 'S'
        else:
            hits.append(None)
    assert sorted(hits) == ['AAA', 'BBB']

    gnomonicus.saveJSON(variants, mutations, effects, path, vcfStem, catalogue, gnomonicus.__version__, -1, reference, '', '', '')

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '21566t>c',
                    'nucleotide_index': 21566,
                    'gene_name': 'S',
                    'gene_position': 2,
                    'codon_idx': 0,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            0,
                            44
                        ],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21566,
                        "REF": "t",
                        "ALTS": [
                            "c"
                        ]
                    },
                    'vcf_idx': 1
                }
            ],
            'mutations': [
                {
                    'mutation': 'F2L',
                    'gene': 'S',
                    'gene_position': 2,
                    'ref': 'ttt',
                    'alt': 'ctt'
                }
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'S',
                        'mutation': 'F2L',
                        'prediction': 'U',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'U'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'U',
                'BBB': 'S'
            }
        }
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))
   


def test_5():
    '''Input:
            NC_045512.2-S_200_indel-minos.vcf
        Expect output:
            variants:    21762_ins_c
            mutations:   S@200_ins_c
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
    gnomonicus.populateVariants(vcfStem, path, diff, True)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue, True)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem, True, True)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    variantGARC = variants['variant'].to_list()
    assert '21762_ins_c' in variantGARC

    mutationGenes = mutations['gene'].to_list()
    for gene in mutationGenes:
        assert gene == 'S'

    mutationGARC = mutations['mutation'].to_list()
    assert '200_ins_c' in mutationGARC

    assert 'AAA' in effects['drug'].to_list()
    assert effects['prediction'][effects['drug'].to_list().index('AAA')] == 'R'

    hits = []
    for _, row in predictions.iterrows():
        assert row['catalogue_name'] == 'gnomonicus_test'
        assert row['catalogue_version'] == 'v1.0'
        assert row['catalogue_values'] == 'RFUS'
        assert row['evidence'] == '{}'
        if row['drug'] == 'AAA':
            hits.append('AAA')
            assert row['prediction'] == "R"
        elif row['drug'] == 'BBB':
            hits.append('BBB')
            assert row['prediction'] == 'S'
        else:
            hits.append(None)
    assert sorted(hits) == ['AAA', 'BBB']

    gnomonicus.saveJSON(variants, mutations, effects, path, vcfStem, catalogue, gnomonicus.__version__, -1, reference, '', '', '')

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '21762_ins_c',
                    'nucleotide_index': 21762,
                    'gene_name': 'S',
                    'gene_position': 200,
                    'codon_idx': 1,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            0,
                            44
                        ],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21762,
                        "REF": "c",
                        "ALTS": [
                            "cc"
                        ]
                    },
                    'vcf_idx': 1
                },
            ],
            'mutations': [
                {
                    'mutation': '200_ins_c',
                    'gene': 'S',
                    'gene_position':200
                },
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'S',
                        'mutation': '200_ins_c',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'R',
                'BBB': 'S'
            }
        }
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


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
    gnomonicus.populateVariants(vcfStem, path, diff, True)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue, True)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem, True, True)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert variants['variant'][0] == '27758g>c'


    assert 'ORF7a' in mutations['gene'].to_list()
    assert 'ORF7b' in mutations['gene'].to_list()

    assert mutations['mutation'][mutations['gene'].to_list().index('ORF7a')] == '!122S'
    assert mutations['mutation'][mutations['gene'].to_list().index('ORF7b')] == 'M1I'

    assert 'AAA' in effects['drug'].to_list()
    assert 'BBB' in effects['drug'].to_list()
    
    assert effects['prediction'][effects['drug'].to_list().index('AAA')] == 'R'
    assert effects['prediction'][effects['drug'].to_list().index('BBB')] == 'R'

    hits = []
    for _, row in predictions.iterrows():
        assert row['catalogue_name'] == 'gnomonicus_test'
        assert row['catalogue_version'] == 'v1.0'
        assert row['catalogue_values'] == 'RFUS'
        assert row['evidence'] == '{}'
        if row['drug'] == 'AAA':
            hits.append('AAA')
            assert row['prediction'] == "R"
        elif row['drug'] == 'BBB':
            hits.append('BBB')
            assert row['prediction'] == 'R'
        else:
            hits.append(None)
    assert sorted(hits) == ['AAA', 'BBB']

    gnomonicus.saveJSON(variants, mutations, effects, path, vcfStem, catalogue, gnomonicus.__version__, -1, reference, '', '', '')

    expectedJSON = {
        'meta': {
            'workflow_version': gnomonicus.__version__,
            'guid': vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0"
        },
        'data': {
            'variants': [
                {
                    'variant': '27758g>c',
                    'nucleotide_index': 27758,
                    'gene_name': 'ORF7a',
                    'gene_position': 122,
                    'codon_idx': 1,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            0,
                            44
                        ],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 27758,
                        "REF": "g",
                        "ALTS": [
                            "c"
                        ]
                    },
                    'vcf_idx': 1
                },
                {
                    'variant': '27758g>c',
                    'nucleotide_index': 27758,
                    'gene_name': 'ORF7b',
                    'gene_position': 1,
                    'codon_idx': 2,
                    'vcf_evidence': {
                        "GT": [
                            1,
                            1
                        ],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [
                            0,
                            44
                        ],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 27758,
                        "REF": "g",
                        "ALTS": [
                            "c"
                        ]
                    },
                    'vcf_idx': 1
                }
            ],
            'mutations': [
                {
                    'mutation': '!122S',
                    'gene': 'ORF7a',
                    'gene_position': 122,
                    'ref': 'tga',
                    'alt': 'tca'
                },
                {
                    'mutation': 'M1I',
                    'gene': 'ORF7b',
                    'gene_position': 1,
                    'ref': 'atg',
                    'alt': 'atc'
                }
            ],
            'effects': {
                'AAA': [
                    {
                        'gene': 'ORF7a',
                        'mutation': '!122S',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
                'BBB': [
                    {
                        'gene': 'ORF7b',
                        'mutation': 'M1I',
                        'prediction': 'R',
                        'evidence': {}
                    },
                    {
                        'phenotype': 'R'
                    }
                ],
            },
            'antibiogram': {
                'AAA': 'R',
                'BBB': 'R'
            }
        }
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(json.load(open(os.path.join(path, f'{vcfStem}.gnomonicus-out.json'), 'r')))

    #assert == does work here, but gives ugly errors if mismatch
    #Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_7():
    '''Input:
            NC_045512.2-S_E484K&1450_ins_a-minos.vcf
        Expect output:
            variants:    23012g>a, 23012_ins_a
            mutations:   S@E484K, S@1450_ins_a, S@1450_ins_a&S@E484K
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
    assert variants_[0] == '23012_ins_a'
    assert variants_[1] == '23012g>a'

    #Sort the mutations for comparing
    mutations_ = sorted(list(zip(mutations['GENE'], mutations['MUTATION'])), key= lambda x: x[0] + x[1] if x[0] is not None else x[1])

    assert mutations_[0][0] == 'S'
    assert mutations_[0][1] == '1450_ins_a'

    assert mutations_[1][0] == 'S'
    assert mutations_[1][1] == 'E484K'


    #Expected effects. For each row, x[0] = DRUG, x[1] = GENE, x[2] = MUTATION, x[3] = PREDICTION
    expected = [
        ['AAA', 'S', 'E484K', 'R'],
        ['AAA', 'S', '1450_ins_a', 'R'],
        ['BBB', None, 'S@1450_ins_a&S@E484K', 'R'],
    ]

    compare_effects(effects, expected)


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
                    "GENE_POSITION",
                    'VCF_EVIDENCE'
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX",
                    'VCF_EVIDENCE'
                    ]
            }
        },
        'data': {
            'VARIANTS': [
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

    #For whatever reason, recursive_diff thinks the VCF evidence fields are different types
    #So compare separately...
    expected_vcf = str({'GT': (1, 1), 'DP': 44, 'DPF': 0.991, 'COV': (0, 44), 'FRS': 1.0, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'g', 'ALTS': ('aa',)})
    for i in range(len(actualJSON['data']['VARIANTS'])):
        assert sorted(actualJSON['data']['VARIANTS'][i]['VCF_EVIDENCE']) == sorted(expected_vcf) 

    #Clean up the VCF Evidences
    del actualJSON['data']['VARIANTS'][0]['VCF_EVIDENCE']
    del actualJSON['data']['VARIANTS'][1]['VCF_EVIDENCE']
    del actualJSON['data']['MUTATIONS'][0]['VCF_EVIDENCE']
    del actualJSON['data']['MUTATIONS'][1]['VCF_EVIDENCE']

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
                    "GENE_POSITION",
                    'VCF_EVIDENCE'
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX",
                    'VCF_EVIDENCE'
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

    #For whatever reason, recursive_diff thinks the VCF evidence fields are different types
    #So compare separately...
    expected_vcf = str({'GT': (1, 1), 'DP': 44, 'DPF': 0.991, 'COV': (0, 44), 'FRS': 1.0, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'g', 'ALTS': ('a',)})
    assert sorted(actualJSON['data']['VARIANTS'][0]['VCF_EVIDENCE']) == sorted(expected_vcf) 
    assert sorted(actualJSON['data']['MUTATIONS'][0]['VCF_EVIDENCE']) == sorted(expected_vcf) 

    #Clean up the VCF Evidences
    del actualJSON['data']['VARIANTS'][0]['VCF_EVIDENCE']
    del actualJSON['data']['MUTATIONS'][0]['VCF_EVIDENCE']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON_, actualJSON)

def test_9():
    '''Test minority populations
    Input:
        NC_045512.2-minors.vcf
    Expect output:
        variants:    25382t>c:0.045, 25283_del_g:0.045, 25252_ins_cc:0.045
        mutations:   !1274Q:0.045, 3721_del_g:0.045, 3690_ins_cc:0.045
        predictions: {'AAA': 'R'}
    '''
    #Setup
    setupOutput('9')
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)

    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-minors.vcf", 
        ignore_filter=True, 
        minor_population_indices={25382, 25283, 25252, 21558}
    )
    vcfStem = "NC_045512.2-minors"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/9/"
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
    assert variants_ == sorted(['25382t>c:0.045', '25283_del_t:0.045','25252_ins_cc:0.045', '21558g>a:0.045'])

    #Sort the mutations for comparing
    mutations_ = sorted(list(zip(mutations['GENE'], mutations['MUTATION'])), key= lambda x: x[0] + x[1] if x[0] is not None else x[1])
    assert mutations_ == sorted([('S', '!1274Q:0.045'), ('S', '3721_del_t:0.045'), ('S', '3690_ins_cc:0.045'), ('S', 'g-5a:0.045')])


    #Expected effects. For each row, x[0] = DRUG, x[1] = GENE, x[2] = MUTATION, x[3] = PREDICTION
    expected = [
        ['AAA', 'S', 'g-5a:0.045', 'U'],
        ['AAA', 'S', '!1274Q:0.045', 'R'],
        ['AAA', 'S', '3721_del_t:0.045', 'R'],
        ['AAA', 'S', '3690_ins_cc:0.045', 'R'],
    ]
    compare_effects(effects, expected)

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
                    "GENE_POSITION",
                    'VCF_EVIDENCE'
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX",
                    'VCF_EVIDENCE'
                    ]
            }
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': '25382t>c:0.045',
                    'NUCLEOTIDE_INDEX': 25382
                },
                {
                    'VARIANT': '21558g>a:0.045',
                    'NUCLEOTIDE_INDEX': 21558
                },
                {
                    'VARIANT': '25283_del_t:0.045',
                    'NUCLEOTIDE_INDEX': 25283
                },
                {
                    'VARIANT': '25252_ins_cc:0.045',
                    'NUCLEOTIDE_INDEX': 25252
                },
            ],
            'MUTATIONS': [
                {
                    'MUTATION': '!1274Q:0.045',
                    'GENE': 'S',
                    'GENE_POSITION':1274
                },
                {
                    'MUTATION': 'g-5a:0.045',
                    'GENE': 'S',
                    'GENE_POSITION':-5
                },
                {
                    'MUTATION': '3721_del_t:0.045',
                    'GENE': 'S',
                    'GENE_POSITION':3721
                },
                {
                    'MUTATION': '3690_ins_cc:0.045',
                    'GENE': 'S',
                    'GENE_POSITION':3690
                },
            ],
            'EFFECTS': {
                'AAA': [
                    {
                        'GENE': 'S',
                        'MUTATION': '!1274Q:0.045',
                        'PREDICTION': 'R'
                    },
                    {
                        'GENE': 'S',
                        'MUTATION': 'g-5a:0.045',
                        'PREDICTION': 'U'
                    },
                    {
                        'GENE': 'S',
                        'MUTATION': '3721_del_t:0.045',
                        'PREDICTION': 'R'
                    },
                    {
                        'GENE': 'S',
                        'MUTATION': '3690_ins_cc:0.045',
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

    #For whatever reason, recursive_diff thinks the VCF evidence fields are different types
    #So compare separately...
    expected_vcf = [
        "{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'g', 'ALTS': ('a',)}",
        "{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'g', 'ALTS': ('gcc',)}",
        "{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'tt', 'ALTS': ('t',)}",
        "{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 't', 'ALTS': ('c',)}",
    ]
    for i in range(len(actualJSON['data']['VARIANTS'])):
        assert sorted(actualJSON['data']['VARIANTS'][i]['VCF_EVIDENCE']) == sorted(expected_vcf[i])
    expected_vcf = [
        str([{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'g', 'ALTS': ('a',)}]),
        'nan',
        str([{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'g', 'ALTS': ('gcc',)}]),
        str([{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'tt', 'ALTS': ('t',)}]),
    ]
    for i in range(len(actualJSON['data']['MUTATIONS'])):
        assert sorted(str(actualJSON['data']['MUTATIONS'][i]['VCF_EVIDENCE'])) == sorted(expected_vcf[i])

    #Clean up the VCF Evidences
    for i in range(len(actualJSON['data']['VARIANTS'])):
        del actualJSON['data']['VARIANTS'][i]['VCF_EVIDENCE']
    for i in range(len(actualJSON['data']['MUTATIONS'])):
        del actualJSON['data']['MUTATIONS'][i]['VCF_EVIDENCE']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON_, actualJSON)


def test_10():
    '''Test minority populations
    Input:
        NC_045512.2-minors.vcf
    Expect output:
        variants:    25382t>c:2, 25283_del_g:2, 25252_ins_cc:2
        mutations:   !1274Q:2, 3721_del_g:2, 3690_ins_cc:2
        predictions: {'AAA': 'R'}
    '''
    #Setup
    setupOutput('10')
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)

    catalogue = piezo.ResistanceCatalogue("tests/test-cases/NC_045512.2-test-catalogue-COV.csv", prediction_subset_only=True)
    
    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-minors.vcf", 
        ignore_filter=True, 
        minor_population_indices={25382, 25283, 25252, 21558}
    )
    vcfStem = "NC_045512.2-minors"

    sample = reference + vcf

    diff = reference - sample

    #Populate the tables
    path = "tests/outputs/10/"
    gnomonicus.populateVariants(vcfStem, path, diff, catalogue=catalogue)
    mutations, referenceGenes = gnomonicus.populateMutations(vcfStem, path, diff, 
                                    reference, sample, catalogue)
    gnomonicus.populateEffects(path, catalogue, mutations, referenceGenes, vcfStem)

    #Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")

    #Sort the variants for comparing
    variants_ = sorted(variants['VARIANT'])
    assert variants_ == sorted(['25382t>c:2', '25283_del_t:2', '25252_ins_cc:2', '21558g>a:2'])

    #Sort the mutations for comparing
    mutations_ = sorted(list(zip(mutations['GENE'], mutations['MUTATION'])), key= lambda x: x[0] + x[1] if x[0] is not None else x[1])
    assert mutations_ == sorted([('S', '!1274Q:2'), ('S', '3721_del_t:2'), ('S', '3690_ins_cc:2'), ('S', 'g-5a:2')])


    #Expected effects. For each row, x[0] = DRUG, x[1] = GENE, x[2] = MUTATION, x[3] = PREDICTION
    expected = [
        ['AAA', 'S', 'g-5a:2', 'U'],
        ['AAA', 'S', '!1274Q:2', 'R'],
        ['AAA', 'S', '3721_del_t:2', 'R'],
        ['AAA', 'S', '3690_ins_cc:2', 'R'],
    ]
    compare_effects(effects, expected)

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
                    "GENE_POSITION",
                    'VCF_EVIDENCE'
                    ],
                "VARIANTS": [
                    "VARIANT",
                    "NUCLEOTIDE_INDEX",
                    'VCF_EVIDENCE'
                    ]
            }
        },
        'data': {
            'VARIANTS': [
                {
                    'VARIANT': '25382t>c:2',
                    'NUCLEOTIDE_INDEX': 25382
                },
                {
                    'VARIANT': '21558g>a:2',
                    'NUCLEOTIDE_INDEX': 21558
                },
                {
                    'VARIANT': '25283_del_t:2',
                    'NUCLEOTIDE_INDEX': 25283
                },
                {
                    'VARIANT': '25252_ins_cc:2',
                    'NUCLEOTIDE_INDEX': 25252
                },
            ],
            'MUTATIONS': [
                {
                    'MUTATION': '!1274Q:2',
                    'GENE': 'S',
                    'GENE_POSITION':1274
                },
                {
                    'MUTATION': 'g-5a:2',
                    'GENE': 'S',
                    'GENE_POSITION':-5
                },
                {
                    'MUTATION': '3721_del_t:2',
                    'GENE': 'S',
                    'GENE_POSITION':3721
                },
                {
                    'MUTATION': '3690_ins_cc:2',
                    'GENE': 'S',
                    'GENE_POSITION':3690
                },
            ],
            'EFFECTS': {
                'AAA': [
                    {
                        'GENE': 'S',
                        'MUTATION': '!1274Q:2',
                        'PREDICTION': 'R'
                    },
                    {
                        'GENE': 'S',
                        'MUTATION': 'g-5a:2',
                        'PREDICTION': 'U'
                    },
                    {
                        'GENE': 'S',
                        'MUTATION': '3721_del_t:2',
                        'PREDICTION': 'R'
                    },
                    {
                        'GENE': 'S',
                        'MUTATION': '3690_ins_cc:2',
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

    for i in range(len(actualJSON['data']['VARIANTS'])):
        print(actualJSON['data']['VARIANTS'][i]['VCF_EVIDENCE'])
    print()
    for i in range(len(actualJSON['data']['MUTATIONS'])):
        print(actualJSON['data']['MUTATIONS'][i]['VCF_EVIDENCE'])
    #For whatever reason, recursive_diff thinks the VCF evidence fields are different types
    #So compare separately...
    expected_vcf = [
        "{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'g', 'ALTS': ('a',)}",
        "{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'g', 'ALTS': ('gcc',)}",
        "{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'tt', 'ALTS': ('t',)}",
        "{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 't', 'ALTS': ('c',)}",
    ]
    for i in range(len(actualJSON['data']['VARIANTS'])):
        assert sorted(actualJSON['data']['VARIANTS'][i]['VCF_EVIDENCE']) == sorted(expected_vcf[i]) 
    expected_vcf = [
        "[{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'g', 'ALTS': ('a',)}]",
        "nan",
        "[{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'g', 'ALTS': ('gcc',)}]",
        "[{'GT': (0, 0), 'DP': 44, 'DPF': 0.991, 'COV': (42, 2), 'FRS': 0.045, 'GT_CONF': 300.34, 'GT_CONF_PERCENTILE': 54.73, 'REF': 'tt', 'ALTS': ('t',)}]",
    ]
    for i in range(len(actualJSON['data']['MUTATIONS'])):
        assert sorted(str(actualJSON['data']['MUTATIONS'][i]['VCF_EVIDENCE'])) == sorted(expected_vcf[i]) 

    #Clean up the VCF Evidences
    for i in range(len(actualJSON['data']['VARIANTS'])):
        del actualJSON['data']['VARIANTS'][i]['VCF_EVIDENCE']
    for i in range(len(actualJSON['data']['MUTATIONS'])):
        del actualJSON['data']['MUTATIONS'][i]['VCF_EVIDENCE']

    #This already asserts that the inputs are equal so no need for assert
    recursive_eq(expectedJSON_, actualJSON)



def compare_effects(effects: pd.DataFrame, expected: [str]) -> None:
    '''Compare an effects DataFrame with the expected values

    Args:
        effects (pd.DataFrame): Effects DataFrame
        expected ([str]): List of expected values (in order)
    '''
    #Sort the effects for comparing
    effects_ = [i[1] for i in sorted([(str(e), e) for _, e in effects.iterrows()], key=lambda x: x[0])]
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

