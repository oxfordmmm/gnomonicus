"""Suite of unit tests matching the test cases defined in tests/NC_045512.2-README.md

Run from root of git dir with `pytest -vv`
"""
import gzip
import json
import os
import pickle
import shutil
import sys

import gumpy
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
    if type(obj) == type(None):
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
    # Ensure that there is not an existing pickle (for clean start and testing loadGenome)
    if os.path.exists("tests/test-cases/NC_045512.2.gbk.pkl"):
        os.remove("tests/test-cases/NC_045512.2.gbk.pkl")
    # Gzipped
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk.gz", False)

    # Pickled
    if os.path.exists("tests/test-cases/NC_045512.2.gbk.pkl"):
        os.remove("tests/test-cases/NC_045512.2.gbk.pkl")
    reference_ = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)

    assert reference == reference_

    # Last edge case of loadGenome is gbk.pkl but its a gzipped file:
    r = pickle.load(open("tests/test-cases/NC_045512.2.gbk.pkl", "rb"))
    gzip_path = "tests/test-cases/reference.gbk.pkl"
    if sys.version_info >= (3, 12):
        gzip_path += ".gz"
    pickle.dump(r, gzip.open(gzip_path, "wb"))

    reference_ = gnomonicus.loadGenome("tests/test-cases/reference.gbk", False)

    assert reference == reference_

    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_E484K-minos.vcf",
        ignore_filter=True,
        bypass_reference_calls=True,
    )
    vcfStem = "NC_045512.2-S_E484K-minos"

    sample = reference + vcf

    diff = reference - sample

    # Populate the tables
    path = "tests/outputs/0/"
    gnomonicus.populateVariants(vcfStem, path, diff, False, False)
    mutations, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, False, False
    )

    # Check for differences if a catalogue is not given. Should be the same mutations but different referenceGenes
    mutations_, referenceGenes_, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, None, False, False
    )

    assert mutations.equals(mutations_)
    assert referenceGenes != referenceGenes_

    # Trying to raise an InvalidMutationException with malformed mutation
    should_raise_error = pd.DataFrame(
        {
            "uniqueid": ["a"],
            "gene": ["S"],
            "mutation": ["aa"],
            "nucleotide_number": ["a"],
            "nucleotide_index": ["a"],
            "gene_position": ["a"],
            "alt": ["a"],
            "ref": ["a"],
            "codes_protein": ["a"],
            "indel_length": ["a"],
            "indel_nucleotides": ["a"],
            "amino_acid_number": ["a"],
            "amino_acid_sequence": ["a"],
            "number_nucleotide_changes": ["a"],
        }
    )

    with pytest.raises(gnomonicus.InvalidMutationException):
        gnomonicus.populateEffects(
            path, catalogue, should_raise_error, referenceGenes, vcfStem, False, False
        )

    # Edge case of minor variant which is in >1 gene
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk.gz", False)
    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-minors-gene-overlap.vcf",
        ignore_filter=True,
        minor_population_indices={267},
    )
    sample = reference + vcf

    diff = reference - sample

    variants = gnomonicus.minority_population_variants(diff, None, None)
    assert variants["variant"].tolist() == ["267t>c:0.045", "267t>c:0.045"]
    assert sorted(variants["gene"].tolist()) == sorted(["ORF1ab", "ORF1ab_2"])
    # These genes are weird and start in the same place (so have matching positions)
    assert variants["gene_position"].tolist() == [1, 1]
    assert variants["codon_idx"].tolist() == [2, 2]

    # Edge case of a genome which doesn't have any gene overlap
    reference = gnomonicus.loadGenome("tests/test-cases/TEST-DNA-no-overlap.gbk", False)
    vcf = gumpy.VCFFile(
        "tests/test-cases/TEST-DNA-no-overlap-snp.vcf",
        ignore_filter=True,
        minor_population_indices={78, 91, 92, 93, 94, 95},
    )
    sample = reference + vcf

    diff = reference - sample

    variants = gnomonicus.minority_population_variants(diff, None, None)
    assert variants["variant"].tolist() == ["78t>g:0.5", "93_del_cc:0.5"]
    assert variants["gene"].tolist() == [None, "C"]
    assert variants["gene_position"].tolist() == [pd.NA, 3]
    assert variants["codon_idx"].tolist() == [pd.NA, 1]


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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_E484K-minos.vcf",
        ignore_filter=True,
        bypass_reference_calls=True,
    )
    vcfStem = "NC_045512.2-S_E484K-minos"

    sample = reference + vcf

    diff = reference - sample

    # Populate the tables
    path = "tests/outputs/1/"
    gnomonicus.populateVariants(vcfStem, path, diff, True, False)
    mutations, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, _, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, referenceGenes, vcfStem, True, True
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
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
        {},
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
    gnomonicus.populateVariants(vcfStem, path, diff, True, False)
    mutations, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, None, True, False
    )
    gnomonicus.populateEffects(
        path, None, mutations, referenceGenes, vcfStem, True, True
    )

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
        path,
        vcfStem,
        None,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
        {},
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


def test_2():
    """Input:
        NC_045512.2-S_E484K-samtools.vcf
    Expect output:
        variants:    23012g>a
        mutations:   S@E484K
        predictions: {'AAA': 'R', 'BBB': 'S'}
    """
    # Setup
    setupOutput("2")
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_E484K-samtools.vcf",
        ignore_filter=True,
        bypass_reference_calls=True,
    )
    vcfStem = "NC_045512.2-S_E484K-samtools"

    sample = reference + vcf

    diff = reference - sample

    # Populate the tables
    path = "tests/outputs/2/"
    gnomonicus.populateVariants(vcfStem, path, diff, True, False)
    mutations, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, _, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, referenceGenes, vcfStem, True, True
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
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
        {},
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
                        "PL": [255, 33, 0],
                        "POS": 23012,
                        "REF": "g",
                        "ALTS": ["a"],
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


def test_3():
    """Input:
        NC_045512.2-S_F2F-minos.vcf
    Expect output:
        variants:    21568t>c
        mutations:   S@F2F
        predictions: {'AAA': 'S', 'BBB': 'S'}
    """
    # Setup
    setupOutput("3")
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_F2F-minos.vcf",
        ignore_filter=True,
        bypass_reference_calls=True,
    )
    vcfStem = "NC_045512.2-S_F2F-minos"

    sample = reference + vcf

    diff = reference - sample

    # Populate the tables
    path = "tests/outputs/3/"
    gnomonicus.populateVariants(vcfStem, path, diff, True, False)
    mutations_, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, _, _ = gnomonicus.populateEffects(
        path, catalogue, mutations_, referenceGenes, vcfStem, True, True
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert variants["variant"][0] == "21568t>c"

    assert mutations["gene"][0] == "S"
    assert mutations["mutation"][0] == "F2F"

    assert "AAA" in effects["drug"].to_list()
    assert effects["prediction"][effects["drug"].to_list().index("AAA")] == "S"

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
            assert row["prediction"] == "S"
        else:
            hits.append(None)
    assert sorted(hits) == ["AAA", "BBB"]

    gnomonicus.saveJSON(
        variants,
        mutations_,
        e,
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
        {},
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
                }
            ],
            "mutations": [
                {
                    "mutation": "F2F",
                    "gene": "S",
                    "gene_position": 2,
                    "ref": "ttt",
                    "alt": "ttc",
                },
                {"mutation": "t6c", "gene": "S", "gene_position": 6},
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "S",
                        "mutation": "F2F",
                        "prediction": "S",
                        "evidence": {"row": 6},
                    },
                    {"phenotype": "S"},
                ],
            },
            "antibiogram": {"AAA": "S", "BBB": "S"},
        },
    }

    expectedJSON = json.loads(json.dumps(expectedJSON, sort_keys=True))

    actualJSON = prep_json(
        json.load(open(os.path.join(path, f"{vcfStem}.gnomonicus-out.json"), "r"))
    )

    # assert == does work here, but gives ugly errors if mismatch
    # Recursive_eq reports neat places they differ
    recursive_eq(ordered(expectedJSON), ordered(actualJSON))


def test_3_fasta_adjudicated():
    """Input:
        NC_045512.2-S_F2F-minos.vcf
    Expect output:
        variants:    21568t>c (FASTA calling this N)
        mutations:   S@F2F
        predictions: {'AAA': 'F', 'BBB': 'S'}
    """
    # Setup
    setupOutput("3_F")
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_F2F-minos.vcf",
        ignore_filter=True,
        bypass_reference_calls=True,
    )
    vcfStem = "NC_045512.2-S_F2F-minos"

    sample = reference + vcf

    diff = reference - sample

    # Populate the tables
    path = "tests/outputs/3_F/"
    gnomonicus.populateVariants(vcfStem, path, diff, True, False)
    mutations_, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, _, mutations = gnomonicus.populateEffects(
        path,
        catalogue,
        mutations_,
        referenceGenes,
        vcfStem,
        True,
        True,
        fasta="tests/test-cases/NC_045512.2.all_n.fasta",
        reference=reference,
        make_mutations_csv=True,
    )

    # Check for expected values within csvs
    variants = pd.read_csv(path + f"{vcfStem}.variants.csv")
    mutations_csv = pd.read_csv(path + f"{vcfStem}.mutations.csv")
    effects = pd.read_csv(path + f"{vcfStem}.effects.csv")
    predictions = pd.read_csv(path + f"{vcfStem}.predictions.csv")

    assert variants["variant"][0] == "21568t>c"

    # Sort the mutations for comparing
    mutations_ = sorted(
        list(zip(mutations_csv["gene"], mutations_csv["mutation"])),
        key=lambda x: x[0] + x[1] if x[0] is not None else x[1],
    )
    assert mutations_ == sorted(
        [
            ("S", "F2F"),
            ("S", "F2X"),
            ("S", "E484X"),
            ("S", "t6c"),
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
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
        {},
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
                }
            ],
            "mutations": [
                {
                    "mutation": "F2F",
                    "gene": "S",
                    "gene_position": 2,
                    "ref": "ttt",
                    "alt": "ttc",
                },
                {
                    "mutation": "F2X",
                    "gene": "S",
                    "gene_position": 2,
                    "ref": "ttt",
                    "alt": "xxx",
                },
                {
                    "mutation": "E484X",
                    "gene": "S",
                    "gene_position": 484,
                    "ref": "gaa",
                    "alt": "xxx",
                },
                {"mutation": "t6c", "gene": "S", "gene_position": 6},
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
                        "mutation": "F2X",
                        "prediction": "F",
                        "evidence": {"row": 21, "FASTA called": "N"},
                    },
                    {
                        "gene": "S",
                        "mutation": "E484X",
                        "prediction": "F",
                        "evidence": {"row": 2, "FASTA called": "N"},
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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_F2L-minos.vcf",
        ignore_filter=True,
        bypass_reference_calls=True,
    )
    vcfStem = "NC_045512.2-S_F2L-minos"

    sample = reference + vcf

    diff = reference - sample

    # Populate the tables
    path = "tests/outputs/4/"
    gnomonicus.populateVariants(vcfStem, path, diff, True, False)
    mutations, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, _, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, referenceGenes, vcfStem, True, True
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
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
        {},
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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_200_indel-minos.vcf",
        ignore_filter=True,
        bypass_reference_calls=True,
    )
    vcfStem = "NC_045512.2-S_200_indel-minos"

    sample = reference + vcf

    diff = reference - sample

    # Populate the tables
    path = "tests/outputs/5/"
    gnomonicus.populateVariants(vcfStem, path, diff, True, False)
    mutations, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, _, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, referenceGenes, vcfStem, True, True
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
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
        {},
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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-double-minos.vcf",
        ignore_filter=True,
        bypass_reference_calls=True,
    )
    vcfStem = "NC_045512.2-double-minos"

    sample = reference + vcf

    diff = reference - sample

    # Populate the tables
    path = "tests/outputs/6/"
    gnomonicus.populateVariants(vcfStem, path, diff, True, False)
    mutations, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, _, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, referenceGenes, vcfStem, True, True
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
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
        {},
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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_E484K&1450_ins_a-minos.vcf",
        ignore_filter=True,
        bypass_reference_calls=True,
    )
    vcfStem = "NC_045512.2-S_E484K&1450_ins_a-minos"

    sample = reference + vcf

    diff = reference - sample

    # Populate the tables
    path = "tests/outputs/7/"
    gnomonicus.populateVariants(vcfStem, path, diff, True, False)
    mutations, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, _, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, referenceGenes, vcfStem, True, True
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
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
        {},
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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)
    # Force non coding (as all genes included are coding)
    reference.genes["S"]["codes_protein"] = False

    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-S_E484K-minos.vcf",
        ignore_filter=True,
        bypass_reference_calls=True,
    )
    vcfStem = "NC_045512.2-S_E484K-minos"

    sample = reference + vcf

    diff = reference - sample

    # Populate the tables
    path = "tests/outputs/8/"
    gnomonicus.populateVariants(vcfStem, path, diff, True, False)
    mutations, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, _, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, referenceGenes, vcfStem, True, True
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
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
        {},
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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)

    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue.csv", prediction_subset_only=True
    )

    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-minors.vcf",
        ignore_filter=True,
        minor_population_indices={25382, 25283, 25252, 21558},
    )
    vcfStem = "NC_045512.2-minors"

    sample = reference + vcf

    diff = reference - sample

    # Populate the tables
    path = "tests/outputs/9/"
    gnomonicus.populateVariants(vcfStem, path, diff, True, False)
    mutations, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, _, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, referenceGenes, vcfStem, True, True
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
        ["AAA", "S", "g-5a:0.045", "U"],
        ["AAA", "S", "!1274Q:0.045", "R"],
        ["AAA", "S", "3721_del_t:0.045", "R"],
        ["AAA", "S", "3690_ins_cc:0.045", "R"],
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
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
        {},
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
                    "codon_idx": 1,
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
                    "codon_idx": 0,
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
                    "codon_idx": 1,
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
                    "alt": "zzz",
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
                        "prediction": "U",
                        "evidence": {"row": 8},
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
    reference = gnomonicus.loadGenome("tests/test-cases/NC_045512.2.gbk", False)

    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/NC_045512.2-test-catalogue-COV.csv",
        prediction_subset_only=True,
    )

    vcf = gumpy.VCFFile(
        "tests/test-cases/NC_045512.2-minors.vcf",
        ignore_filter=True,
        minor_population_indices={25382, 25283, 25252, 21558},
    )
    vcfStem = "NC_045512.2-minors"

    sample = reference + vcf

    diff = reference - sample

    # Populate the tables
    path = "tests/outputs/10/"
    gnomonicus.populateVariants(vcfStem, path, diff, True, False, catalogue=catalogue)
    mutations, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, _, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, referenceGenes, vcfStem, True, True
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
        ["AAA", "S", "g-5a:2", "U"],
        ["AAA", "S", "!1274Q:2", "R"],
        ["AAA", "S", "3721_del_t:2", "R"],
        ["AAA", "S", "3690_ins_cc:2", "R"],
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
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
        {},
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
                    "codon_idx": 1,
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
                    "codon_idx": 0,
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
                    "codon_idx": 1,
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
                    "alt": "zzz",
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
                        "prediction": "U",
                        "evidence": {"row": 8},
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
    reference = gnomonicus.loadGenome("tests/test-cases/TEST-DNA.gbk", False)

    catalogue = piezo.ResistanceCatalogue(
        "tests/test-cases/TEST-DNA-catalogue.csv", prediction_subset_only=True
    )

    vcf = gumpy.VCFFile(
        "tests/test-cases/TEST-DNA-large-del.vcf",
    )
    vcfStem = "TEST-DNA-large-del"

    sample = reference + vcf

    diff = reference - sample

    # Populate the tables
    path = "tests/outputs/11/"
    gnomonicus.populateVariants(vcfStem, path, diff, True, False, catalogue=catalogue)
    mutations, referenceGenes, errors = gnomonicus.populateMutations(
        vcfStem, path, diff, reference, sample, catalogue, True, False
    )
    e, _, _ = gnomonicus.populateEffects(
        path, catalogue, mutations, referenceGenes, vcfStem, True, True
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

    assert len(mutations) == 4
    assert sorted(mutations["mutation"].tolist()) == sorted(
        ["-1_del_aaaaaaaaccccccccccgggggggggg", "del_0.93", "del_1.0", "4_del_ggg"]
    )

    assert len(effects) == 4
    # Expected effects. For each row, x[0] = DRUG, x[1] = GENE, x[2] = MUTATION, x[3] = PREDICTION
    expected = [
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
        path,
        vcfStem,
        catalogue,
        gnomonicus.__version__,
        -1,
        reference,
        "",
        "",
        "",
        {},
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
                    "codon_idx": 2,
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
                {"mutation": "4_del_ggg", "gene": "C", "gene_position": 4},
            ],
            "effects": {
                "AAA": [
                    {
                        "gene": "A",
                        "mutation": "-1_del_aaaaaaaaccccccccccgggggggggg",
                        "prediction": "U",
                        "evidence": {"row": 5},
                    },
                    {
                        "gene": "A",
                        "mutation": "del_0.93",
                        "prediction": "R",
                        "evidence": {"row": 16},
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
