"""Suite of unit tests matching the test cases defined in tests/NC_045512.2-README.md

Run from root of git dir with `pytest -vv`
"""

import json
import os
import shutil

import grumpy
import pandas as pd
import piezo
import pytest

# Helpful function import for testing nested JSON equality as it gives exact differences
from recursive_diff import recursive_eq

import gnomonicus

"""
Due to complications testing equalities of nested jsons of lists/dicts, there is a lot of 
code specificially dedicated to ensuring the lists are in the same order (they differ due to
a mixture of dictionary behaviour and different positions within files). However, we only care that
the contents of the JSON is there for these tests rather than caring about order.
"""


def setupOutput(testID: str) -> None:
    """Ensure that the output folder exists and is empty in preparation for a test

    Args:
        testID (str): Desired folder name for this test
    """
    path = f"tests/outputs/{testID}"
    # Make the dir if it doesn't exist
    os.makedirs(path, exist_ok=True)

    # Check for contents
    if len(os.listdir(path)) > 0:
        # Not empty, so delete and recreate
        shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)


def prep_json(j: dict) -> dict:
    """Prepare a JSON for comparison by removing fields which cannot be reproduced

    Args:
        j (dict): Initial JSON

    Returns:
        dict: JSON without fields such as time and file paths
    """
    del j["meta"]["time_taken_s"]
    del j["meta"]["UTC-datetime-completed"]
    del j["meta"]["catalogue_file"]
    del j["meta"]["reference_file"]
    del j["meta"]["vcf_file"]
    return j


def variants_key(x):
    """Used as the sorted(key=) function for reliably sorting the variants list

    Args:
        x (list): List of the ordered values as key, value pairs

    Returns:
        str: String of the `variant+gene_name`
    """
    variant = ""
    gene = ""
    for i in x:
        if i[0] == "variant":
            variant = i[1]
        elif i[0] == "gene_name":
            gene = i[1]
    return variant + gene


def ordered(obj):
    """Recursively sort a JSON for equality checking. Based on https://stackoverflow.com/questions/25851183/how-to-compare-two-json-objects-with-the-same-elements-in-a-different-order-equa

    Args:
        obj (object): Any JSON element. Probably one of dict, list, tuple, str, int, None

    Returns:
        object: Sorted JSON
    """
    if isinstance(obj, dict):
        if "variants" in obj.keys():
            # We have the 'data' field which needs a little extra nudge
            if "effects" in obj.keys():
                # Case when we have effects populated
                return [
                    ("antibiogram", ordered(obj["antibiogram"])),
                    ("effects", ordered(obj["effects"])),
                    ("mutations", ordered(obj["mutations"])),
                    (
                        "variants",
                        sorted([ordered(x) for x in obj["variants"]], key=variants_key),
                    ),
                ]
            elif "mutations" in obj.keys():
                # Case for if we have mutations and variants but no effects
                return [
                    ("mutations", ordered(obj["mutations"])),
                    (
                        "variants",
                        sorted([ordered(x) for x in obj["variants"]], key=variants_key),
                    ),
                ]
            else:
                # Case where no effects/mutations have been populated
                return [
                    (
                        "variants",
                        sorted([ordered(x) for x in obj["variants"]], key=variants_key),
                    )
                ]
        else:
            return sorted((k, ordered(obj[k])) for k in sorted(list(obj.keys())))

    if isinstance(obj, list) or isinstance(obj, tuple):
        return sorted(ordered(x) for x in obj)

    # Nones cause issues with ordering as there is no < operator.
    # Convert to string to avoid this
    if isinstance(obj, type(None)):
        return str(obj)

    # Because nan types are helpful, `float('nan') == float('nan') -> False`
    # So check based on str value rather than direct equality
    if str(obj) == "nan":
        # Conversion to None for reproducability
        return str(None)

    if isinstance(obj, int):
        # Ints are still ordered (just not numerically) if sorted by str value, so convert to str
        # to allow sorting lists of None/int
        return str(obj)
    if isinstance(obj, float):
        # Similarly convert float, but check if they are x.0 to avoid str comparison issues
        if int(obj) == obj:
            return str(int(obj))
        else:
            return str(obj)
    else:
        return obj


def test_misc():
    """Testing misc things which should be raised/edge case behaviour not confined neatly to a whole flow test case"""
    setupOutput("0")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")

    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = grumpy.VCFFile("tests/test-cases/NC_045512.2-S_E484K-minos.vcf", True, 3)
    vcfStem = "NC_045512.2-S_E484K-minos"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/0/"
    gnomonicus.populateVariants(vcfStem, path, diff, False, False, sample)
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, False, False
    )

    # Check for differences if a catalogue is not given. Should be the same mutations
    mutations_ = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, None, False, False
    )

    assert mutations.equals(mutations_)


def test_1():
    """Input:
        NC_045512.2-S_E484K-minos.vcf
    Expect output:
        variants:    23012g>a
        mutations:   S@E484K
        predictions: {'AAA': 'R', 'BBB': 'S'}
    """
    # Setup
    setupOutput("1")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = grumpy.VCFFile("tests/test-cases/NC_045512.2-S_E484K-minos.vcf", True, 3)
    vcfStem = "NC_045512.2-S_E484K-minos"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/1/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, False, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert variants["variant"][0] == "23012g>a"

    assert mutations["gene"][0] == "S"
    assert mutations["mutation"][0] == "E484K"

    assert "AAA" in effects["drug"].to_list()
    assert effects["prediction"][effects["drug"].to_list().index("AAA")] == "R"

    hits = []
    for _, row in predictions.iterrows():
        assert row["catalogue_name"] == "gnomonicus_test"
        assert row["catalogue_version"] == "v1.0"
        assert row["catalogue_values"] == "RFUS"
        if row["drug"] == "AAA":
            hits.append("AAA")
            assert row["prediction"] == "R"
        elif row["drug"] == "BBB":
            hits.append("BBB")
            assert row["prediction"] == "S"
        else:
            hits.append(None)
    assert sorted(hits) == ["AAA", "BBB"]

    gnomonicus.saveJSON(
        variants,
        mutations,
        e,
        phenotypes,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "23012g>a",
                    "nucleotide_index": 23012,
                    "gene_name": "S",
                    "gene_position": 484,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "REF": "g",
                        "ALTS": ["a"],
                        "POS": 23012,
                    },
                    "vcf_idx": 1,
                }
            ],
            "mutations": [
                {
                    "mutation": "E484K",
                    "gene": "S",
                    "gene_position": 484,
                    "ref": "gaa",
                    "alt": "aaa",
                }
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "E484K",
                        "prediction": "R",
                        "evidence": {"row": 0},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"AAA": "R", "BBB": "S"},
        },
    }
    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))

    # Running the same test again, but with no catalogue
    # Should just produce variants and mutations sections of the JSON with empty sections for effects+antibiogram
    setupOutput("1")
    gnomonicus.populateVariants(vcfStem, path, diff, True, False, sample)
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, None, True, False
    )
    gnomonicus.populateEffects(path, None, mutations, vcfStem, True, True, reference)

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    # Neither of these should be populated
    with pytest.raises(Exception):
        effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    with pytest.raises(Exception):
        predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    gnomonicus.saveJSON(
        variants,
        mutations,
        None,
        {},
        path,
        vcfStem,
        None,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": None,
            "catalogue_name": None,
            "catalogue_version": None,
        },
        "data": {
            "variants": [
                {
                    "variant": "23012g>a",
                    "nucleotide_index": 23012,
                    "gene_name": "S",
                    "gene_position": 484,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "REF": "g",
                        "ALTS": ["a"],
                        "POS": 23012,
                    },
                    "vcf_idx": 1,
                }
            ],
            "mutations": [
                {
                    "mutation": "E484K",
                    "gene": "S",
                    "gene_position": 484,
                    "ref": "gaa",
                    "alt": "aaa",
                }
            ],
            "effects": {},
            "antibiogram": {},
        },
    }
    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_3():
    """Input:
        NC_045512.2-S_F2F-minos.vcf
    Expect output:
        variants:    21568t>c, 21763_del_t
        mutations:   S@F2F, S@201_del_t
        predictions: {'AAA': 'U', 'BBB': 'S'}
    """
    # Setup
    setupOutput("3")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = grumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_F2F-minos.vcf",
        True,
        3,
    )
    vcfStem = "NC_045512.2-S_F2F-minos"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/3/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, False, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )
    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations_csv = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert variants["variant"][0] == "21568t>c"
    assert variants["variant"][1] == "21763_del_t"

    # Sort the mutations for comparing
    mutations_ = sorted(
        list(zip(mutations_csv["gene"], mutations_csv["mutation"])),
        key=lambda x: x[0] + x[1] if x[0] is not None else x[1],
    )
    assert mutations_ == sorted(
        [
            ("S", "201_del_t"),
            ("S", "F2F"),
        ]
    )

    assert "AAA" in effects["drug"].to_list()

    hits = []
    for _, row in predictions.iterrows():
        assert row["catalogue_name"] == "gnomonicus_test"
        assert row["catalogue_version"] == "v1.0"
        assert row["catalogue_values"] == "RFUS"
        if row["drug"] == "AAA":
            hits.append("AAA")
            assert row["prediction"] == "U"
        elif row["drug"] == "BBB":
            hits.append("BBB")
            assert row["prediction"] == "S"
        else:
            hits.append(None)
    assert sorted(hits) == ["AAA", "BBB"]

    gnomonicus.saveJSON(
        variants,
        mutations_csv,
        e,
        phenotypes,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "21568t>c",
                    "nucleotide_index": 21568,
                    "gene_name": "S",
                    "gene_position": 2,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21568,
                        "REF": "t",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "21763_del_t",
                    "nucleotide_index": 21763,
                    "gene_name": "S",
                    "gene_position": 201,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21762,
                        "REF": "ct",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {
                    "mutation": "F2F",
                    "gene": "S",
                    "gene_position": 2,
                    "ref": "ttt",
                    "alt": "ttc",
                },
                {"mutation": "201_del_t", "gene": "S", "gene_position": 201},
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "F2F",
                        "prediction": "S",
                        "evidence": {"row": 6},
                    },
                    {
                        "gene": "S",
                        "mutation": "201_del_t",
                        "prediction": "U",
                        "evidence": {"row": 23},
                    },
                    {"phenotype": "U"},
                ],
            },
            "antibiogram": {"AAA": "U", "BBB": "S"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_4():
    """Input:
        NC_045512.2-S_F2L-minos.vcf
    Expect output:
        variants:    21566t>c
        mutations:   S@F2L
        predictions: {'AAA': 'U', 'BBB': 'S'}
    """
    # Setup
    setupOutput("4")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = grumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_F2L-minos.vcf",
        True,
        3,
    )
    vcfStem = "NC_045512.2-S_F2L-minos"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/4/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, False, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert variants["variant"][0] == "21566t>c"

    assert mutations["gene"][0] == "S"
    assert mutations["mutation"][0] == "F2L"

    assert "AAA" in effects["drug"].to_list()
    assert effects["prediction"][effects["drug"].to_list().index("AAA")] == "U"

    hits = []
    for _, row in predictions.iterrows():
        assert row["catalogue_name"] == "gnomonicus_test"
        assert row["catalogue_version"] == "v1.0"
        assert row["catalogue_values"] == "RFUS"
        if row["drug"] == "AAA":
            hits.append("AAA")
            assert row["prediction"] == "U"
        elif row["drug"] == "BBB":
            hits.append("BBB")
            assert row["prediction"] == "S"
        else:
            hits.append(None)
    assert sorted(hits) == ["AAA", "BBB"]

    gnomonicus.saveJSON(
        variants,
        mutations,
        e,
        phenotypes,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "21566t>c",
                    "nucleotide_index": 21566,
                    "gene_name": "S",
                    "gene_position": 2,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21566,
                        "REF": "t",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                }
            ],
            "mutations": [
                {
                    "mutation": "F2L",
                    "gene": "S",
                    "gene_position": 2,
                    "ref": "ttt",
                    "alt": "ctt",
                }
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "F2L",
                        "prediction": "U",
                        "evidence": {"row": 3},
                    },
                    {"phenotype": "U"},
                ],
            },
            "antibiogram": {"AAA": "U", "BBB": "S"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_5():
    """Input:
        NC_045512.2-S_200_indel-minos.vcf
    Expect output:
        variants:    21762_ins_c
        mutations:   S@200_ins_c
        predictions: {'AAA': 'R', 'BBB': 'S'}
    """
    # Setup
    setupOutput("5")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = grumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_200_indel-minos.vcf",
        True,
        3,
    )
    vcfStem = "NC_045512.2-S_200_indel-minos"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/5/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, False, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    variantGARC = variants["variant"].to_list()
    assert "21762_ins_c" in variantGARC

    mutationGenes = mutations["gene"].to_list()
    for gene in mutationGenes:
        assert gene == "S"

    mutationGARC = mutations["mutation"].to_list()
    assert "200_ins_c" in mutationGARC

    assert "AAA" in effects["drug"].to_list()
    assert effects["prediction"][effects["drug"].to_list().index("AAA")] == "R"

    hits = []
    for _, row in predictions.iterrows():
        assert row["catalogue_name"] == "gnomonicus_test"
        assert row["catalogue_version"] == "v1.0"
        assert row["catalogue_values"] == "RFUS"
        if row["drug"] == "AAA":
            hits.append("AAA")
            assert row["prediction"] == "R"
        elif row["drug"] == "BBB":
            hits.append("BBB")
            assert row["prediction"] == "S"
        else:
            hits.append(None)
    assert sorted(hits) == ["AAA", "BBB"]

    gnomonicus.saveJSON(
        variants,
        mutations,
        e,
        phenotypes,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "21762_ins_c",
                    "nucleotide_index": 21762,
                    "gene_name": "S",
                    "gene_position": 200,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21762,
                        "REF": "c",
                        "ALTS": ["cc"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {"mutation": "200_ins_c", "gene": "S", "gene_position": 200},
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "200_ins_c",
                        "prediction": "R",
                        "evidence": {"row": 11},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"AAA": "R", "BBB": "S"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_6():
    """Input:
        NC_045512.2-double-minos.vcf
    Expect output:
        variants:    27758g>c
        mutations:   ORF7a!122S, ORF7b@M1I
        predictions: {'AAA': 'R', 'BBB': 'R'}
    """
    # Setup
    setupOutput("6")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = grumpy.VCFFile(
        "tests/test-cases/NC_045512.2-double-minos.vcf",
        True,
        3,
    )
    vcfStem = "NC_045512.2-double-minos"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/6/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, False, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert variants["variant"][0] == "27758g>c"

    assert "ORF7a" in mutations["gene"].to_list()
    assert "ORF7b" in mutations["gene"].to_list()

    assert mutations["mutation"][mutations["gene"].to_list().index("ORF7a")] == "!122S"
    assert mutations["mutation"][mutations["gene"].to_list().index("ORF7b")] == "M1I"

    assert "AAA" in effects["drug"].to_list()
    assert "BBB" in effects["drug"].to_list()

    assert effects["prediction"][effects["drug"].to_list().index("AAA")] == "R"
    assert effects["prediction"][effects["drug"].to_list().index("BBB")] == "R"

    hits = []
    for _, row in predictions.iterrows():
        assert row["catalogue_name"] == "gnomonicus_test"
        assert row["catalogue_version"] == "v1.0"
        assert row["catalogue_values"] == "RFUS"
        if row["drug"] == "AAA":
            hits.append("AAA")
            assert row["prediction"] == "R"
        elif row["drug"] == "BBB":
            hits.append("BBB")
            assert row["prediction"] == "R"
        else:
            hits.append(None)
    assert sorted(hits) == ["AAA", "BBB"]

    gnomonicus.saveJSON(
        variants,
        mutations,
        e,
        phenotypes,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "27758g>c",
                    "nucleotide_index": 27758,
                    "gene_name": "ORF7a",
                    "gene_position": 122,
                    "codon_idx": 1,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 27758,
                        "REF": "g",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "27758g>c",
                    "nucleotide_index": 27758,
                    "gene_name": "ORF7b",
                    "gene_position": 1,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 27758,
                        "REF": "g",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {
                    "mutation": "!122S",
                    "gene": "ORF7a",
                    "gene_position": 122,
                    "ref": "tga",
                    "alt": "tca",
                },
                {
                    "mutation": "M1I",
                    "gene": "ORF7b",
                    "gene_position": 1,
                    "ref": "atg",
                    "alt": "atc",
                },
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "ORF7a",
                        "mutation": "!122S",
                        "prediction": "R",
                        "evidence": {"row": 12},
                    },
                    {"phenotype": "R"},
                ],
                "BBB": [
                    {
                        "gene": "ORF7b",
                        "mutation": "M1I",
                        "prediction": "R",
                        "evidence": {"row": 13},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"AAA": "R", "BBB": "R"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_7():
    """Input:
        NC_045512.2-S_E484K&1450_ins_a-minos.vcf
    Expect output:
        variants:    23012g>a, 23012_ins_a
        mutations:   S@E484K, S@1450_ins_a, S@1450_ins_a&S@E484K
        predictions: {'AAA': 'R', 'BBB': 'R'}
    """
    # Setup
    setupOutput("7")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = grumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_E484K&1450_ins_a-minos.vcf",
        True,
        3,
    )
    vcfStem = "NC_045512.2-S_E484K&1450_ins_a-minos"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/7/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, False, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    # Sort the variants for comparing
    variants_ = sorted(variants["variant"])
    assert variants_[0] == "23012_ins_a"
    assert variants_[1] == "23012g>a"

    # Sort the mutations for comparing
    mutations_ = sorted(
        list(zip(mutations["gene"], mutations["mutation"])),
        key=lambda x: x[0] + x[1] if x[0] is not None else x[1],
    )

    assert mutations_[0][0] == "S"
    assert mutations_[0][1] == "1450_ins_a"

    assert mutations_[1][0] == "S"
    assert mutations_[1][1] == "E484K"

    # Expected effects. For each row, x[0] = DRUG, x[1] = GENE, x[2] = MUTATION, x[3] = PREDICTION
    expected = [
        ["AAA", "S", "E484K", "R"],
        ["AAA", "S", "1450_ins_a", "R"],
        ["BBB", None, "S@1450_ins_a&S@E484K", "R"],
    ]

    compare_effects(effects, expected)

    hits = []
    for _, row in predictions.iterrows():
        assert row["catalogue_name"] == "gnomonicus_test"
        assert row["catalogue_version"] == "v1.0"
        assert row["catalogue_values"] == "RFUS"
        if row["drug"] == "AAA":
            hits.append("AAA")
            assert row["prediction"] == "R"
        elif row["drug"] == "BBB":
            hits.append("BBB")
            assert row["prediction"] == "R"
        else:
            hits.append(None)
    assert sorted(hits) == ["AAA", "BBB"]

    gnomonicus.saveJSON(
        variants,
        mutations,
        e,
        phenotypes,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "23012_ins_a",
                    "nucleotide_index": 23012,
                    "gene_name": "S",
                    "gene_position": 1450,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 23012,
                        "REF": "g",
                        "ALTS": ["aa"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "23012g>a",
                    "nucleotide_index": 23012,
                    "gene_name": "S",
                    "gene_position": 484,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 23012,
                        "REF": "g",
                        "ALTS": ["aa"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {"mutation": "1450_ins_a", "gene": "S", "gene_position": 1450},
                {
                    "mutation": "E484K",
                    "gene": "S",
                    "gene_position": 484,
                    "ref": "gaa",
                    "alt": "aaa",
                },
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "E484K",
                        "prediction": "R",
                        "evidence": {"row": 0},
                    },
                    {
                        "gene": "S",
                        "mutation": "1450_ins_a",
                        "prediction": "R",
                        "evidence": {"row": 11},
                    },
                    {"phenotype": "R"},
                ],
                "BBB": [
                    {
                        "gene": None,
                        "mutation": "S@1450_ins_a&S@E484K",
                        "prediction": "R",
                        "evidence": {"row": 1},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"AAA": "R", "BBB": "R"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_8():
    """Force the "S" gene to be non-coding to test non-coding behaviours
    Input:
        NC_045512.2-S_E484K-minos.vcf
    Expect output:
        variants:    23012g>a
        mutations:   g1450a
        predictions: {'AAA': 'U'}
    """
    # Setup
    setupOutput("8")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")
    # Force non coding (as all genes included are coding)
    genes = [gene for gene in reference.gene_definitions if gene.name != "S"]
    s = None
    for gene in reference.gene_definitions:
        if gene.name == "S":
            s = gene
            s.coding = False
            break
    assert s is not None
    genes.append(s)
    reference.gene_definitions = genes

    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = grumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_E484K-minos.vcf",
        True,
        3,
    )
    vcfStem = "NC_045512.2-S_E484K-minos"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/8/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, False, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert variants["variant"][0] == "23012g>a"

    assert mutations["gene"][0] == "S"
    assert mutations["mutation"][0] == "g1450a"

    assert "AAA" in effects["drug"].to_list()
    assert effects["prediction"][effects["drug"].to_list().index("AAA")] == "U"

    hits = []
    for _, row in predictions.iterrows():
        assert row["catalogue_name"] == "gnomonicus_test"
        assert row["catalogue_version"] == "v1.0"
        assert row["catalogue_values"] == "RFUS"
        if row["drug"] == "AAA":
            hits.append("AAA")
            assert row["prediction"] == "U"
        elif row["drug"] == "BBB":
            hits.append("BBB")
            assert row["prediction"] == "S"
        else:
            hits.append(None)
    assert sorted(hits) == ["AAA", "BBB"]

    gnomonicus.saveJSON(
        variants,
        mutations,
        e,
        phenotypes,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "23012g>a",
                    "nucleotide_index": 23012,
                    "gene_name": "S",
                    "gene_position": 1450,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 23012,
                        "REF": "g",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                }
            ],
            "mutations": [{"mutation": "g1450a", "gene": "S", "gene_position": 1450}],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "g1450a",
                        "prediction": "U",
                        "evidence": {"row": 3},
                    },
                    {"phenotype": "U"},
                ],
            },
            "antibiogram": {"AAA": "U", "BBB": "S"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_9():
    """Test minority populations
    Input:
        NC_045512.2-minors.vcf
    Expect output:
        variants:    25382t>c:0.045, 25283_del_g:0.045, 25252_ins_cc:0.045
        mutations:   !1274Q:0.045, 3721_del_g:0.045, 3690_ins_cc:0.045
        predictions: {'AAA': 'R'}
    """
    # Setup
    setupOutput("9")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")

    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = grumpy.VCFFile(
        "tests/test-cases/NC_045512.2-minors.vcf",
        True,
        2,
    )
    vcfStem = "NC_045512.2-minors"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.FRS)

    # Populate the tables
    path = "tests/outputs/9/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, True, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, True
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    # Sort the variants for comparing
    variants_ = sorted(variants["variant"])
    assert variants_ == sorted(
        ["25382t>c:0.045", "25283_del_t:0.045", "25252_ins_cc:0.045", "21558g>a:0.045"]
    )

    # Sort the mutations for comparing
    mutations_ = sorted(
        list(zip(mutations["gene"], mutations["mutation"])),
        key=lambda x: x[0] + x[1] if x[0] is not None else x[1],
    )
    assert mutations_ == sorted(
        [
            ("S", "!1274Q:0.045"),
            ("S", "3721_del_t:0.045"),
            ("S", "3690_ins_cc:0.045"),
            ("S", "g-5a:0.045"),
        ]
    )

    # Expected effects. For each row, x[0] = DRUG, x[1] = GENE, x[2] = MUTATION, x[3] = PREDICTION
    expected = [
        ["AAA", "S", "g-5a:0.045", "S"],
        ["AAA", "S", "!1274Q:0.045", "R"],
        ["AAA", "S", "3721_del_t:0.045", "R"],
        ["AAA", "S", "3690_ins_cc:0.045", "S"],
    ]
    compare_effects(effects, expected)

    hits = []
    for _, row in predictions.iterrows():
        assert row["catalogue_name"] == "gnomonicus_test"
        assert row["catalogue_version"] == "v1.0"
        assert row["catalogue_values"] == "RFUS"
        if row["drug"] == "AAA":
            hits.append("AAA")
            assert row["prediction"] == "R"
        elif row["drug"] == "BBB":
            hits.append("BBB")
            assert row["prediction"] == "S"
        else:
            hits.append(None)
    assert sorted(hits) == ["AAA", "BBB"]

    gnomonicus.saveJSON(
        variants,
        mutations,
        e,
        phenotypes,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "25382t>c:0.045",
                    "nucleotide_index": 25382,
                    "gene_name": "S",
                    "gene_position": 1274,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 2],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25382,
                        "REF": "t",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "21558g>a:0.045",
                    "nucleotide_index": 21558,
                    "gene_name": "S",
                    "gene_position": -5,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 2],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21558,
                        "REF": "g",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "25252_ins_cc:0.045",
                    "nucleotide_index": 25252,
                    "gene_name": "S",
                    "gene_position": 3690,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 2],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25252,
                        "REF": "g",
                        "ALTS": ["gcc"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "25283_del_t:0.045",
                    "nucleotide_index": 25283,
                    "gene_name": "S",
                    "gene_position": 3721,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 2],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25282,
                        "REF": "tt",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {
                    "mutation": "!1274Q:0.045",
                    "gene": "S",
                    "gene_position": 1274,
                    "ref": "taa",
                    "alt": "caa",
                },
                {"mutation": "g-5a:0.045", "gene": "S", "gene_position": -5},
                {"mutation": "3721_del_t:0.045", "gene": "S", "gene_position": 3721},
                {"mutation": "3690_ins_cc:0.045", "gene": "S", "gene_position": 3690},
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "!1274Q:0.045",
                        "prediction": "R",
                        "evidence": {"row": 20},
                    },
                    {
                        "gene": "S",
                        "mutation": "g-5a:0.045",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "S",
                        "mutation": "3721_del_t:0.045",
                        "prediction": "R",
                        "evidence": {"row": 19},
                    },
                    {
                        "gene": "S",
                        "mutation": "3690_ins_cc:0.045",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"AAA": "R", "BBB": "S"},
        },
    }
    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_10():
    """Test minority populations
    Input:
        NC_045512.2-minors.vcf
    Expect output:
        variants:    25382t>c:2, 25283_del_g:2, 25252_ins_cc:2
        mutations:   !1274Q:2, 3721_del_g:2, 3690_ins_cc:2
        predictions: {'AAA': 'R'}
    """
    # Setup
    setupOutput("10")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")

    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue-COV.csv",
        prediction_subset_only=True,
    )

    vcf = grumpy.VCFFile(
        "tests/test-cases/NC_045512.2-minors.vcf",
        True,
        2,
    )
    vcfStem = "NC_045512.2-minors"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/10/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, True, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, True
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    # Sort the variants for comparing
    variants_ = sorted(variants["variant"])
    assert variants_ == sorted(
        ["25382t>c:2", "25283_del_t:2", "25252_ins_cc:2", "21558g>a:2"]
    )

    # Sort the mutations for comparing
    mutations_ = sorted(
        list(zip(mutations["gene"], mutations["mutation"])),
        key=lambda x: x[0] + x[1] if x[0] is not None else x[1],
    )
    assert mutations_ == sorted(
        [
            ("S", "!1274Q:2"),
            ("S", "3721_del_t:2"),
            ("S", "3690_ins_cc:2"),
            ("S", "g-5a:2"),
        ]
    )

    # Expected effects. For each row, x[0] = DRUG, x[1] = GENE, x[2] = MUTATION, x[3] = PREDICTION
    expected = [
        ["AAA", "S", "g-5a:2", "S"],
        ["AAA", "S", "!1274Q:2", "R"],
        ["AAA", "S", "3721_del_t:2", "R"],
        ["AAA", "S", "3690_ins_cc:2", "S"],
    ]
    compare_effects(effects, expected)

    hits = []
    for _, row in predictions.iterrows():
        assert row["catalogue_name"] == "gnomonicus_test"
        assert row["catalogue_version"] == "v1.0"
        assert row["catalogue_values"] == "RFUS"
        if row["drug"] == "AAA":
            hits.append("AAA")
            assert row["prediction"] == "R"
        elif row["drug"] == "BBB":
            hits.append("BBB")
            assert row["prediction"] == "S"
        else:
            hits.append(None)
    assert sorted(hits) == ["AAA", "BBB"]

    gnomonicus.saveJSON(
        variants,
        mutations,
        e,
        phenotypes,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "25382t>c:2",
                    "nucleotide_index": 25382,
                    "gene_name": "S",
                    "gene_position": 1274,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 2],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25382,
                        "REF": "t",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "21558g>a:2",
                    "nucleotide_index": 21558,
                    "gene_name": "S",
                    "gene_position": -5,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 2],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21558,
                        "REF": "g",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "25252_ins_cc:2",
                    "nucleotide_index": 25252,
                    "gene_name": "S",
                    "gene_position": 3690,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 2],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25252,
                        "REF": "g",
                        "ALTS": ["gcc"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "25283_del_t:2",
                    "nucleotide_index": 25283,
                    "gene_name": "S",
                    "gene_position": 3721,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [0, 0],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [42, 2],
                        "FRS": 0.045,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 25282,
                        "REF": "tt",
                        "ALTS": ["t"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {
                    "mutation": "!1274Q:2",
                    "gene": "S",
                    "gene_position": 1274,
                    "ref": "taa",
                    "alt": "caa",
                },
                {"mutation": "g-5a:2", "gene": "S", "gene_position": -5},
                {"mutation": "3721_del_t:2", "gene": "S", "gene_position": 3721},
                {"mutation": "3690_ins_cc:2", "gene": "S", "gene_position": 3690},
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "!1274Q:2",
                        "prediction": "R",
                        "evidence": {"row": 20},
                    },
                    {
                        "gene": "S",
                        "mutation": "g-5a:2",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "S",
                        "mutation": "3721_del_t:2",
                        "prediction": "R",
                        "evidence": {"row": 19},
                    },
                    {
                        "gene": "S",
                        "mutation": "3690_ins_cc:2",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"AAA": "R", "BBB": "S"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_11():
    """Testing a catalogue and sample which have large deletions
    Input:
        TEST-DNA-large-del.vcf
    Expect output:
        variants:    3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc
        mutations:   A@-1_del_aaaaaaaaccccccccccgggggggggg, A@del_0.93, B@del_1.0, C@4_del_ggg
        predictions: {'AAA': 'R'}
    """
    # Setup
    setupOutput("11")
    reference = grumpy.Genome("tests/test-cases/TEST-DNA.gbk")

    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/TEST-DNA-catalogue.csv", prediction_subset_only=True
    )

    vcf = grumpy.VCFFile("tests/test-cases/TEST-DNA-large-del.vcf", True, 1)
    vcfStem = "TEST-DNA-large-del"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/11/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, True, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, True
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert len(variants) == 3
    assert variants["variant"].tolist() == [
        "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
        "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
        "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
    ]

    assert len(mutations) == 5
    assert sorted(mutations["mutation"].tolist()) == sorted(
        [
            "-1_del_aaaaaaaaccccccccccgggggggggg",
            "del_0.93",
            "del_1.0",
            "1_del_gggttttttttttaaaaaaaaaacccccccccc",
            "4_del_ggg",
        ]
    )

    assert len(effects) == 5
    # Expected effects. For each row, x[0] = DRUG, x[1] = GENE, x[2] = MUTATION, x[3] = PREDICTION
    expected = [
        # Odd ordering, but it works
        ["AAA", "B", "1_del_gggttttttttttaaaaaaaaaacccccccccc", "U"],
        ["AAA", "A", "-1_del_aaaaaaaaccccccccccgggggggggg", "U"],
        ["AAA", "A", "del_0.93", "R"],
        ["AAA", "B", "del_1.0", "U"],
        ["AAA", "C", "4_del_ggg", "U"],
    ]
    compare_effects(effects, expected)

    hits = []
    for _, row in predictions.iterrows():
        assert row["catalogue_name"] == "gnomonicus_test_dna"
        assert row["catalogue_version"] == "v1.0"
        assert row["catalogue_values"] == "RFUS"
        if row["drug"] == "AAA":
            hits.append("AAA")
            assert row["prediction"] == "R"
        else:
            hits.append(None)
    assert sorted(hits) == ["AAA"]

    gnomonicus.saveJSON(
        variants,
        mutations,
        e,
        phenotypes,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "TEST_DNA",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test_dna",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                    "nucleotide_index": 3,
                    "gene_name": "A",
                    "gene_position": -1,
                    "codon_idx": None,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 2,
                        "COV": [1, 1],
                        "GT_CONF": 2.05,
                        "POS": 2,
                        "REF": "aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                    "nucleotide_index": 3,
                    "gene_name": "B",
                    "gene_position": 1,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 2,
                        "COV": [1, 1],
                        "GT_CONF": 2.05,
                        "POS": 2,
                        "REF": "aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "3_del_aaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                    "nucleotide_index": 3,
                    "gene_name": "C",
                    "gene_position": 4,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 2,
                        "COV": [1, 1],
                        "GT_CONF": 2.05,
                        "POS": 2,
                        "REF": "aaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccccccccccggggggggggttttttttttaaaaaaaaaaccc",
                        "ALTS": ["a"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {
                    "mutation": "-1_del_aaaaaaaaccccccccccgggggggggg",
                    "gene": "A",
                    "gene_position": -1,
                },
                {"mutation": "del_0.93", "gene": "A", "gene_position": None},
                {"mutation": "del_1.0", "gene": "B", "gene_position": None},
                {
                    "mutation": "1_del_gggttttttttttaaaaaaaaaacccccccccc",
                    "gene": "B",
                    "gene_position": 1,
                },
                {"mutation": "4_del_ggg", "gene": "C", "gene_position": 4},
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "A",
                        "mutation": "-1_del_aaaaaaaaccccccccccgggggggggg",
                        "prediction": "U",
                        "evidence": {"row": 4},
                    },
                    {
                        "gene": "A",
                        "mutation": "del_0.93",
                        "prediction": "R",
                        "evidence": {"row": 16},
                    },
                    {
                        "gene": "B",
                        "mutation": "1_del_gggttttttttttaaaaaaaaaacccccccccc",
                        "prediction": "U",
                        "evidence": {"row": 9},
                    },
                    {
                        "gene": "B",
                        "mutation": "del_1.0",
                        "prediction": "U",
                        "evidence": {"row": 17},
                    },
                    {
                        "gene": "C",
                        "mutation": "4_del_ggg",
                        "prediction": "U",
                        "evidence": {"row": 14},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {
                "AAA": "R",
            },
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_12():
    """Sample with a few mutations + a catalogue supporting both general mutlits, and an epistasis rule
    Input:
        NC_045512.2-epistasis.vcf
    Expect output:
        variants:    23012g>a, 23012_ins_a, 21763_del_t
        mutations:   S@E484K, S@1450_ins_a, S@201_del_t, S@1450_ins_a&S@E484K, S@201_del_t&S@E484K
        predictions: {'AAA': 'S', 'BBB': 'R'}
    """
    # Setup
    setupOutput("12")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue-epistasis.csv",
        prediction_subset_only=True,
    )

    vcf = grumpy.VCFFile("tests/test-cases/NC_045512.2-epistasis.vcf", True, 2)
    vcfStem = "NC_045512.2-epistasis"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/12/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, True, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, True
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    # Sort the variants for comparing
    variants_ = sorted(variants["variant"])
    assert variants_[0] == "21763_del_t"
    assert variants_[1] == "23012_ins_a"
    assert variants_[2] == "23012g>a"

    # Sort the mutations for comparing
    mutations_ = sorted(
        list(zip(mutations["gene"], mutations["mutation"])),
        key=lambda x: x[0] + x[1] if x[0] is not None else x[1],
    )

    assert mutations_[0][0] == "S"
    assert mutations_[0][1] == "1450_ins_a"

    assert mutations_[1][0] == "S"
    assert mutations_[1][1] == "201_del_t"

    assert mutations_[2][0] == "S"
    assert mutations_[2][1] == "E484K"

    # Expected effects. For each row, x[0] = DRUG, x[1] = GENE, x[2] = MUTATION, x[3] = PREDICTION
    expected = [
        ["AAA", "S", "E484K", "R"],
        ["AAA", "S", "201_del_t", "R"],
        ["AAA", "S", "1450_ins_a", "R"],
        ["AAA", None, "S@201_del_t&S@E484K", "S"],
        ["AAA", None, "S@1450_ins_a&S@E484K", "R"],
        ["BBB", None, "S@1450_ins_a&S@E484K", "R"],
    ]

    compare_effects(effects, expected)

    hits = []
    for _, row in predictions.iterrows():
        assert row["catalogue_name"] == "gnomonicus_test"
        assert row["catalogue_version"] == "v1.0"
        assert row["catalogue_values"] == "RFUS"
        if row["drug"] == "AAA":
            hits.append("AAA")
            assert row["prediction"] == "S"
        elif row["drug"] == "BBB":
            hits.append("BBB")
            assert row["prediction"] == "R"
        else:
            hits.append(None)
    assert sorted(hits) == ["AAA", "BBB"]

    gnomonicus.saveJSON(
        variants,
        mutations,
        e,
        phenotypes,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "21763_del_t",
                    "nucleotide_index": 21763,
                    "gene_name": "S",
                    "gene_position": 201,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 21762,
                        "REF": "ct",
                        "ALTS": ["c"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "23012_ins_a",
                    "nucleotide_index": 23012,
                    "gene_name": "S",
                    "gene_position": 1450,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 23012,
                        "REF": "g",
                        "ALTS": ["aa"],
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "23012g>a",
                    "nucleotide_index": 23012,
                    "gene_name": "S",
                    "gene_position": 484,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "POS": 23012,
                        "REF": "g",
                        "ALTS": ["aa"],
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {"mutation": "201_del_t", "gene": "S", "gene_position": 201},
                {"mutation": "1450_ins_a", "gene": "S", "gene_position": 1450},
                {
                    "mutation": "E484K",
                    "gene": "S",
                    "gene_position": 484,
                    "ref": "gaa",
                    "alt": "aaa",
                },
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "E484K",
                        "prediction": "R",
                        "evidence": {"row": 0},
                    },
                    {
                        "gene": "S",
                        "mutation": "1450_ins_a",
                        "prediction": "R",
                        "evidence": {"row": 11},
                    },
                    {
                        "gene": "S",
                        "mutation": "201_del_t",
                        "prediction": "R",
                        "evidence": {"row": 11},
                    },
                    {
                        "gene": None,
                        "mutation": "S@1450_ins_a&S@E484K",
                        "prediction": "R",
                        "evidence": {"row": 22},
                    },
                    {
                        "gene": None,
                        "mutation": "S@201_del_t&S@E484K",
                        "prediction": "S",
                        "evidence": {"row": 23},
                    },
                    {"phenotype": "S"},
                ],
                "BBB": [
                    {
                        "gene": None,
                        "mutation": "S@1450_ins_a&S@E484K",
                        "prediction": "R",
                        "evidence": {"row": 1},
                    },
                    {"phenotype": "R"},
                ],
            },
            "antibiogram": {"AAA": "S", "BBB": "R"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_13():
    """Force null calls by passing a min_dp
    Input:
        NC_045512.2-S_E484K-minos.vcf
    Expect output:
        variants:    23012g>x
        mutations:   E484X
        predictions: {'AAA': 'F'}
    """
    # Setup
    setupOutput("13")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = grumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_E484K-minos.vcf",
        True,
        45,
    )
    vcfStem = "NC_045512.2-S_E484K-minos"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/13/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, True, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, True
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert variants["variant"][0] == "23012g>x"

    assert mutations["gene"][0] == "S"
    assert mutations["mutation"][0] == "E484X"

    assert "AAA" in effects["drug"].to_list()
    assert effects["prediction"][effects["drug"].to_list().index("AAA")] == "F"

    hits = []
    for _, row in predictions.iterrows():
        assert row["catalogue_name"] == "gnomonicus_test"
        assert row["catalogue_version"] == "v1.0"
        assert row["catalogue_values"] == "RFUS"
        if row["drug"] == "AAA":
            hits.append("AAA")
            assert row["prediction"] == "F"
        elif row["drug"] == "BBB":
            hits.append("BBB")
            assert row["prediction"] == "S"
        else:
            hits.append(None)
    assert sorted(hits) == ["AAA", "BBB"]

    gnomonicus.saveJSON(
        variants,
        mutations,
        e,
        phenotypes,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "variants": [
                {
                    "variant": "23012g>x",
                    "nucleotide_index": 23012,
                    "gene_name": "S",
                    "gene_position": 484,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "GT": [1, 1],
                        "DP": 44,
                        "DPF": 0.991,
                        "COV": [0, 44],
                        "FRS": 1.0,
                        "GT_CONF": 300.34,
                        "GT_CONF_PERCENTILE": 54.73,
                        "REF": "g",
                        "ALTS": ["a"],
                        "POS": 23012,
                    },
                    "vcf_idx": 1,
                }
            ],
            "mutations": [
                {
                    "mutation": "E484X",
                    "gene": "S",
                    "gene_position": 484,
                    "ref": "gaa",
                    "alt": "xaa",
                }
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "E484X",
                        "prediction": "F",
                        "evidence": {"row": 2},
                    },
                    {"phenotype": "F"},
                ],
            },
            "antibiogram": {"AAA": "F", "BBB": "S"},
        },
    }
    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_14():
    """No variants/mutations/effects/predictions
    Input:
        NC_045512.2-no-variants.vcf
    Expect output:
        variants:
        mutations:
        predictions:
    """
    # Setup
    setupOutput("14")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = grumpy.VCFFile(
        "tests/test-cases/NC_045512.2-no-variants.vcf",
        True,
        2,
    )
    vcfStem = "NC_045512.2-no-variants"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/14/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, True, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, True
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")

    assert len(variants) == 0
    assert len(mutations) == 0
    assert len(effects) == 0


def test_15():
    """VCF with a complex row which shouldn't get any evidence for the variant"""
    # Setup
    setupOutput("15")
    reference = grumpy.Genome("tests/test-cases/NC_045512.2.gbk")
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue-COV.csv",
        prediction_subset_only=True,
    )

    vcf = grumpy.VCFFile(
        "tests/test-cases/NC_045512.2-complex.vcf",
        True,
        3,
    )
    vcfStem = "NC_045512.2-complex"

    sample = grumpy.mutate(reference, vcf)

    diff = grumpy.GenomeDifference(reference, sample, grumpy.MinorType.COV)

    # Populate the tables
    path = "tests/outputs/15/"
    gnomonicus.populateVariants(
        vcfStem, path, diff, True, True, sample, catalogue=catalogue
    )
    mutations = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, True
    )
    e, phenotypes, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, vcfStem, True, True, reference
    )

    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")

    gnomonicus.saveJSON(
        variants,
        mutations,
        e,
        phenotypes,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")

    assert len(variants) == 4
    assert variants["variant"][0] == "25252g>g"
    assert (
        variants["variant"][1]
        == "25253_del_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    )
    assert variants["variant"][2] == "25252g>t:99"
    assert (
        variants["variant"][3]
        == "25253_del_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa:99"
    )

    assert len(mutations) == 3
    assert mutations["gene"][0] == "S"
    assert (
        mutations["mutation"][0]
        == "3691_del_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    )
    assert mutations["gene"][1] == "S"
    assert mutations["mutation"][1] == "V1230V:99"
    assert mutations["gene"][2] == "S"
    assert (
        mutations["mutation"][2]
        == "3691_del_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa:99"
    )

    assert len(effects) == 3

    expectedJSON = {
        "meta": {
            "workflow_version": gnomonicus.__version__,
            "guid": vcfStem,
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_task": "resistance_prediction",
            "reference": "NC_045512",
            "catalogue_type": "RFUS",
            "catalogue_name": "gnomonicus_test",
            "catalogue_version": "v1.0",
        },
        "data": {
            "antibiogram": {"AAA": "U", "BBB": "S"},
            "variants": [
                {
                    "variant": "25252g>g",
                    "nucleotide_index": 25252,
                    "gene_name": "S",
                    "gene_position": 1230,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "VCF row is complex": "Please refer to the VCF file for evidence."
                    },
                    "vcf_idx": 2,
                },
                {
                    "variant": "25252g>t:99",
                    "nucleotide_index": 25252,
                    "gene_name": "S",
                    "gene_position": 1230,
                    "codon_idx": 2,
                    "vcf_evidence": {
                        "VCF row is complex": "Please refer to the VCF file for evidence."
                    },
                    "vcf_idx": 1,
                },
                {
                    "variant": "25253_del_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                    "nucleotide_index": 25253,
                    "gene_name": "S",
                    "gene_position": 3691,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "VCF row is complex": "Please refer to the VCF file for evidence."
                    },
                    "vcf_idx": 2,
                },
                {
                    "variant": "25253_del_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa:99",
                    "nucleotide_index": 25253,
                    "gene_name": "S",
                    "gene_position": 3691,
                    "codon_idx": 0,
                    "vcf_evidence": {
                        "VCF row is complex": "Please refer to the VCF file for evidence."
                    },
                    "vcf_idx": 1,
                },
            ],
            "mutations": [
                {
                    "gene": "S",
                    "gene_position": 1230,
                    "mutation": "V1230V:99",
                    "ref": "gtg",
                    "alt": "gtt",
                },
                {
                    "gene": "S",
                    "gene_position": 3691,
                    "mutation": "3691_del_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                },
                {
                    "gene": "S",
                    "gene_position": 3691,
                    "mutation": "3691_del_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa:99",
                },
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "3691_del_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                        "prediction": "U",
                        "evidence": {"row": 9},
                    },
                    {
                        "gene": "S",
                        "mutation": "3691_del_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa:99",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {
                        "gene": "S",
                        "mutation": "V1230V:99",
                        "prediction": "S",
                        "evidence": {},
                    },
                    {"phenotype": "U"},
                ]
            },
        },
    }
    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def compare_effects(effects: pd.DataFrame, expected: [str]) -> None:
    """Compare an effects DataFrame with the expected values

    Args:
        effects (pd.DataFrame): Effects DataFrame
        expected ([str]): List of expected values (in order)
    """
    # Sort the effects for comparing
    effects_ = [
        i[1]
        for i in sorted(
            [(str(e), e) for _, e in effects.iterrows()], key=lambda x: x[0]
        )
    ]
    assert len(expected) == len(effects_)
    # Iter expected and effects to check for equality
    for row, exp in zip(effects_, expected):
        assert row["drug"] == exp[0]

        # Dealing with pd.nan as equality doesn't work here...
        if pd.isnull(row["gene"]):
            assert exp[1] is None
        else:
            assert row["gene"] == exp[1]

        assert row["mutation"] == exp[2]
        assert row["prediction"] == exp[3]
