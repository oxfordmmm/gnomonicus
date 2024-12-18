"""gnomonicus.py is a library providing functions which pull together output VCF of the Lodestone TB pipeline
    with a reference genome and a resistance catalogue, and utilise grumpy and
    piezo to produce variants, mutations and an antibiogram.

Based on sp3predict
"""

import copy
import datetime
import json
import logging
import os
import re
import warnings
from collections import defaultdict, OrderedDict
from typing import Dict, List, Tuple

import grumpy  # type: ignore
import pandas as pd
import piezo
from tqdm import tqdm

# Suppress warnings about concatenating empty DataFrames
warnings.simplefilter(action="ignore", category=FutureWarning)


class InvalidMutationException(Exception):
    """Custom exception raised when an invalid mutation is detected"""

    def __init__(self, gene: str, mutation: str):
        """Raise this exception

        Args:
            gene (str): Name of the gene
            mutation (str): The invalid mutation
        """
        self.message = f"{gene}@{mutation} is not a valid mutation!"
        super().__init__(self.message)


def parse_grumpy_evidence(evidence: grumpy.VCFRow) -> dict:
    """Parse the grumpy evidence into a dictionary for JSON.

    As rust doesn't have multi-type hashmaps, we need to convert all values from strings
    Should just be parsing ints and floats

    Args:
        evidence (grumpy.VCFRow): Evidence from the VCF

    Returns:
        dict: Parsed evidence
    """
    ev = {}
    for key, value in evidence.fields.items():
        item: list[int | float | None] | int | float | None = []
        if key == "GT":
            # Special case here as we need to split the string and parse int/Nones
            gt = value[0].split("/")
            item = [int(g) if g[0] != "." else None for g in gt]
        elif key == "DP":
            # Single value expected here too (most of the time)
            if len(value) == 1:
                # Odd edge cases around types though
                if value[0] == ".":
                    item = None
                else:
                    try:
                        item = int(value[0])
                    except ValueError:
                        item = float(value[0])
            else:
                item = [int(v) for v in value]
        else:
            if item is None or isinstance(item, int) or isinstance(item, float):
                # Should never happen but appease mypy
                continue
            for v in value:
                # Use duck typing to determine if it's a float or int
                try:
                    item.append(int(v))
                except ValueError:
                    try:
                        item.append(float(v))
                    except ValueError:
                        item.append(v)
        ev[key] = item
    for key, val in ev.items():
        if isinstance(val, list) and len(val) == 1:
            # Unpack single values as they probably shouldn't be lists
            ev[key] = val[0]
    # We also want to add back in some of the VCF items which aren't in the fields dict
    ev["POS"] = evidence.position
    ev["REF"] = evidence.reference
    ev["ALTS"] = evidence.alternative
    return ev


def populateVariants(
    vcfStem: str,
    outputDir: str,
    diff: grumpy.GenomeDifference,
    make_csv: bool,
    resistanceGenesOnly: bool,
    sample: grumpy.Genome,
    catalogue: piezo.ResistanceCatalogue | None = None,
) -> pd.DataFrame:
    """Populate and save the variants DataFrame as a CSV

    Args:
        vcfStem (str): The stem of the filename for the VCF file. Used as a uniqueID
        outputDir (str): Path to the desired output directory
        diff (grumpy.GenomeDifference): GenomeDifference object between reference and the sample
        make_csv (bool): Whether to write the CSV of the dataframe
        resistanceGenesOnly (bool): Whether to use just genes present in the resistance catalogue
        sample (grumpy.Genome): Sample genome object
        catalogue (piezo.ResistanceCatalogue | None, optional): Catalogue for determining FRS or COV for minority populations. If None is given, FRS is assumed. Defaults to None

    Returns:
        pd.DataFrame: DataFrame of the variants
    """
    # Populate variants table directly from GenomeDifference
    vals: dict[str, list] = {
        "variant": [],
        "nucleotide_index": [],
        "indel_length": [],
        "indel_nucleotides": [],
        "vcf_evidence": [],
        "vcf_idx": [],
        "gene": [],
        "gene_position": [],
        "codon_idx": [],
    }
    for variant in diff.variants:
        vals["variant"].append(variant.variant)
        vals["nucleotide_index"].append(variant.nucleotide_index)
        vals["indel_length"].append(variant.indel_length)
        vals["indel_nucleotides"].append(variant.indel_nucleotides)
        vals["vcf_evidence"].append(
            json.dumps(parse_grumpy_evidence(sample.get_vcf_row(variant.evidence)))
        )
        vals["vcf_idx"].append(variant.vcf_idx)
        vals["gene"].append(variant.gene_name)
        vals["gene_position"].append(variant.gene_position)
        vals["codon_idx"].append(variant.codon_idx)

    for variant in diff.minor_variants:
        vals["variant"].append(variant.variant)
        vals["nucleotide_index"].append(variant.nucleotide_index)
        vals["indel_length"].append(variant.indel_length)
        vals["indel_nucleotides"].append(variant.indel_nucleotides)
        vals["vcf_evidence"].append(
            json.dumps(parse_grumpy_evidence(sample.get_vcf_row(variant.evidence)))
        )
        vals["vcf_idx"].append(variant.vcf_idx)
        vals["gene"].append(variant.gene_name)
        vals["gene_position"].append(variant.gene_position)
        vals["codon_idx"].append(variant.codon_idx)

    # Use of Int64 rather than int is required here as pandas doesn't allow mixed int/None
    variants = pd.DataFrame(vals).astype(
        {
            "vcf_evidence": "object",
            "nucleotide_index": "Int64",
            "indel_length": "Int64",
            "vcf_idx": "Int64",
            "gene_position": "Int64",
            "codon_idx": "Int64",
        }
    )

    if catalogue is not None:
        # Figure out if we want to keep all of the variants
        genes = getGenes(sample, catalogue, resistanceGenesOnly)
        to_drop = []
        for idx, row in variants.iterrows():
            if row["gene"] not in genes:
                # Not a variant we're interested in, so remove
                to_drop.append(idx)

        variants.drop(index=to_drop, inplace=True)

    # Add unique ID to each record
    variants["uniqueid"] = vcfStem

    variants = variants[
        [
            "uniqueid",
            "variant",
            "gene",
            "gene_position",
            "codon_idx",
            "nucleotide_index",
            "indel_length",
            "indel_nucleotides",
            "vcf_evidence",
            "vcf_idx",
        ]
    ]
    variants = variants.drop_duplicates()
    if make_csv:
        # Save CSV
        variants.to_csv(
            os.path.join(outputDir, f"{vcfStem}.variants.csv"),
            header=True,
            index=False,
        )
    variants.reset_index(inplace=True)
    return variants


def get_minority_population_type(
    catalogue: piezo.ResistanceCatalogue | None,
) -> grumpy.MinorType:
    """Figure out if a catalogue uses FRS or COV. If neither or both, default to FRS

    Args:
        catalogue (piezo.ResistanceCatalogue | None): Catalogue

    Returns:
        grumpy.MinorType: Enum for FRS or COV respectively
    """
    if catalogue is None:
        # Nothing given, so default to FRS
        return grumpy.MinorType.FRS
    frs = 0
    cov = 0
    for minor in catalogue.catalogue.rules["MINOR"]:
        for m in minor.split(","):
            if m:
                m = float(m)
                assert m > 0, f"Minor populations must be positive: {m}"
                if m < 1:
                    # FRS
                    frs += 1
                else:
                    # COV
                    cov += 1
    # We have just COV
    if cov > 0 and frs == 0:
        return grumpy.MinorType.COV
    # We have anything else
    return grumpy.MinorType.FRS


def getGenes(
    sample: grumpy.Genome,
    resistanceCatalogue: piezo.ResistanceCatalogue,
    resistanceGenesOnly: bool,
) -> set:
    """Get the genes we're interested in.

    This is either just resistance genes which have variants, or all which have variants

    Args:
        sample (grumpy.Genome): Sample's genome object
        resistanceCatalogue (piezo.ResistanceCatalogue): Resistance catalogue
        resistanceGenesOnly (bool): Whether to just use genes within the catalogue

    Returns:
        set[str]: Set of gene names
    """
    if resistanceCatalogue:
        if resistanceGenesOnly:
            resistanceGenes = set(resistanceCatalogue.catalogue.genes)
            # Catch multi/epistasis rules which might not have specific instances
            multis = set()
            for _, rule in resistanceCatalogue.catalogue.rules.iterrows():
                if rule["MUTATION_TYPE"] in ["MULTI", "EPISTASIS"]:
                    mutations = rule["MUTATION"]
                    for mut in mutations.split("&"):
                        multis.add(mut.split("@")[0])
            resistanceGenes = resistanceGenes.union(multis)
        else:
            resistanceGenes = set(sample.gene_names)

        return sample.genes_with_mutations.intersection(resistanceGenes)

    else:
        # No catalogue, so just stick to genes in the sample
        return sample.genes


def count_nucleotide_changes(ref: str | None, alt: str | None) -> int | None:
    """Count number of changes between ref and alt

    Args:
        ref (str | None): Ref nucleotides
        alt (str | None): Alt nucleotides

    Returns:
        int | None: SNP distance if SNP, None if ref/alt were None
    """
    if ref is None or alt is None:
        return None
    return sum(1 for r, a in zip(ref, alt) if r != a)


def populateMutations(
    vcfStem: str,
    outputDir: str,
    genome_diff: grumpy.GenomeDifference,
    reference: grumpy.Genome,
    sample: grumpy.Genome,
    resistanceCatalogue: piezo.ResistanceCatalogue,
    make_csv: bool,
    resistanceGenesOnly: bool,
) -> pd.DataFrame | None:
    """Popuate and save the mutations DataFrame as a CSV, then return it for use in predictions

    Args:
        vcfStem (str): The stem of the filename of the VCF file. Used as a uniqueID
        outputDir (str): Path to the desired output directory
        genome_diff (grumpy.GenomeDifference): GenomeDifference object between reference and this sample
        reference (grumpy.Genome): Reference genome
        sample (grumpy.Genome): Sample genome
        resistanceCatalogue (piezo.ResistanceCatalogue): Resistance catalogue (used to find which genes to check)
        make_csv (bool): Whether to write the CSV of the dataframe
        resistanceGenesOnly (bool): Whether to use just genes present in the resistance catalogue

    Raises:
        MissingFieldException: Raised when the mutations DataFrame does not contain the required fields

    Returns:
        pd.DataFrame: The mutations DataFrame
    """
    genesWithMutations = getGenes(sample, resistanceCatalogue, resistanceGenesOnly)

    # Iter resistance genes with variation to produce gene level mutations - concating into a single dataframe
    mutations: dict[str, list] = {
        "gene": [],
        "mutation": [],
        "ref": [],
        "alt": [],
        "nucleotide_number": [],
        "nucleotide_index": [],
        "gene_position": [],
        "codes_protein": [],
        "indel_length": [],
        "indel_nucleotides": [],
        "amino_acid_number": [],
        "amino_acid_sequence": [],
        "number_nucleotide_changes": [],
    }
    for gene_name in sorted(genesWithMutations):
        gene_diff = grumpy.GeneDifference(
            reference.get_gene(gene_name),
            sample.get_gene(gene_name),
            get_minority_population_type(resistanceCatalogue),
        )
        for mutation in gene_diff.mutations:
            mutations["gene"].append(gene_name)
            mutations["mutation"].append(mutation.mutation)
            mutations["ref"].append(mutation.ref_nucleotides)
            mutations["alt"].append(mutation.alt_nucleotides)
            mutations["nucleotide_number"].append(mutation.nucleotide_number)
            mutations["nucleotide_index"].append(mutation.nucleotide_index)
            mutations["gene_position"].append(mutation.gene_position)
            mutations["codes_protein"].append(mutation.codes_protein)
            mutations["indel_length"].append(mutation.indel_length)
            mutations["indel_nucleotides"].append(mutation.indel_nucleotides)
            mutations["amino_acid_number"].append(mutation.amino_acid_number)
            mutations["amino_acid_sequence"].append(mutation.amino_acid_sequence)
            mutations["number_nucleotide_changes"].append(
                count_nucleotide_changes(
                    mutation.ref_nucleotides, mutation.alt_nucleotides
                )
            )
        for mutation in gene_diff.minor_mutations:
            mutations["gene"].append(gene_name)
            mutations["mutation"].append(mutation.mutation)
            mutations["ref"].append(mutation.ref_nucleotides)
            mutations["alt"].append(mutation.alt_nucleotides)
            mutations["nucleotide_number"].append(mutation.nucleotide_number)
            mutations["nucleotide_index"].append(mutation.nucleotide_index)
            mutations["gene_position"].append(mutation.gene_position)
            mutations["codes_protein"].append(mutation.codes_protein)
            mutations["indel_length"].append(mutation.indel_length)
            mutations["indel_nucleotides"].append(mutation.indel_nucleotides)
            mutations["amino_acid_number"].append(mutation.amino_acid_number)
            mutations["amino_acid_sequence"].append(mutation.amino_acid_sequence)
            mutations["number_nucleotide_changes"].append(
                count_nucleotide_changes(
                    mutation.ref_nucleotides, mutation.alt_nucleotides
                )
            )

    # Ensure correct datatypes
    mutations_df = pd.DataFrame(mutations).astype(
        {
            "mutation": "str",
            "gene": "str",
            "nucleotide_number": "Int64",
            "nucleotide_index": "Int64",
            "gene_position": "Int64",
            "alt": "str",
            "ref": "str",
            "codes_protein": "bool",
            "indel_length": "Int64",
            "indel_nucleotides": "str",
            "amino_acid_number": "Int64",
            "amino_acid_sequence": "str",
        }
    )
    # If there were mutations, write them to a CSV

    # Add VCF stem as the uniqueID
    mutations_df["uniqueid"] = vcfStem

    if make_csv:
        write_mutations_csv(
            mutations_df, os.path.join(outputDir, f"{vcfStem}.mutations.csv")
        )

    return mutations_df


def write_mutations_csv(
    mutations: pd.DataFrame, path: str, filter: bool = True
) -> None:
    """Prep and write the mutations CSV to the given filepath.

    Args:
        mutations (pd.DataFrame): Muations CSV
        path (str): Path to write to
        filter (bool, optional): Whether to filter nucleotide changes
    """
    # Reorder the columns
    mutations = mutations[
        [
            "uniqueid",
            "gene",
            "mutation",
            "ref",
            "alt",
            "nucleotide_number",
            "nucleotide_index",
            "gene_position",
            "codes_protein",
            "indel_length",
            "indel_nucleotides",
            "amino_acid_number",
            "amino_acid_sequence",
            "number_nucleotide_changes",
        ]
    ]

    # As we have concated several dataframes, the index is 0,1,2,0,1...
    # Reset it so that we can use it to delete
    mutations.reset_index(drop=True, inplace=True)
    mutations_ = copy.deepcopy(mutations)
    if filter:
        # Filter out nucleotide variants from synonymous mutations to avoid duplication of data
        to_drop = []
        for idx2, row in mutations_.iterrows():
            if (
                row["codes_protein"]
                and row["ref"] is not None
                and row["alt"] is not None
            ):
                # Protein coding so check if nucleotide within coding region
                if len(row["ref"]) == 1:
                    # Nucleotide SNP
                    to_drop.append(idx2)
        mutations_.drop(index=to_drop, inplace=True)
    # Save it as CSV
    mutations_.to_csv(path, index=False)


def subset_multis(
    multis: set[str],
    mutations_: list[tuple[str | None, str]],
    just_joined: bool = False,
) -> list[tuple[str | None, str]]:
    """Given a list of multis, combine existing mutations to match them

    Args:
        multis (set[str]): Set of multi mutations
        mutations_ (list[tuple[str | None, str]]): mutations in
        just_joined (bool, optional): Whether to return all mutations or just joined multis. Defaults to False

    Returns:
        list[tuple[str|None, str]]: Multi-mutations which should hit catalogue rules
    """
    mutations = [
        (gene, mut, mut.split(":")[-1] if ":" in mut else None)
        for (gene, mut) in mutations_
        if gene
        if not None
    ]
    joined = set([gene + "@" + mut for (gene, mut, _) in mutations])
    large_del_re = re.compile(
        r"""
        .*del_[01]\.[0-9]+.*
        """,
        re.VERBOSE,
    )
    snp_re = re.compile(
        r"""
        [acgtxzA-Z!]-?[0-9]+[acgtxzA-Z!]
    """,
        re.VERBOSE,
    )
    early_stop_re = re.compile(
        r"""
        [A-Z!]-?[0-9]+!
    """,
        re.VERBOSE,
    )
    indel_re = re.compile(
        r"""
        -?[0-9]+_(ins|del|indel|mixed|fs)(_(([0-9]+)|([a-z]+)))?
    """,
        re.VERBOSE,
    )

    # Find these once so they aren't fetched on every iteration
    existing_genes = set([gene for (gene, _, _) in mutations])
    large_dels = [
        (gene, mut, minor)
        for (gene, mut, minor) in mutations
        if large_del_re.fullmatch(mut)
    ]
    snps = [
        (gene, mut, minor) for (gene, mut, minor) in mutations if snp_re.fullmatch(mut)
    ]
    early_stop = [
        (gene, mut, minor)
        for (gene, mut, minor) in mutations
        if early_stop_re.fullmatch(mut)
    ]
    indels = [
        (gene, mut, minor)
        for (gene, mut, minor) in mutations
        if indel_re.fullmatch(mut)
    ]

    new_mutations: list[tuple[str | None, str]] = []
    for multi in multis:
        multi_match = []
        check = True
        if "*" in multi or "?" in multi or large_del_re.fullmatch(multi):
            # We need to be a bit more cleaver here to avoid combinatorics ruining things
            for mutation in multi.split("&"):
                gene, mut = mutation.split("@")
                rule_is_minor = ":" in mut
                this_match = []
                if gene not in existing_genes:
                    # Gene doesn't exist in this sample, so skip
                    check = False
                if "*" in mutation:
                    # A few cases here, wildcard or early-stop SNP, ins, del, indel, fs
                    if "?" in mutation:
                        # Wildcard-non-synon SNP so match any SNP in this gene
                        if "-" in mutation:
                            promoter = True
                        else:
                            promoter = False
                        matched = False
                        for g, m, minor in snps:
                            if g == gene and (
                                (minor is None and not rule_is_minor)
                                or (minor is not None and rule_is_minor)
                            ):
                                if promoter and "-" in m:
                                    matched = True
                                    if minor is not None:
                                        this_match.append(g + "@" + m + ":" + minor)
                                    else:
                                        this_match.append(g + "@" + m)
                                elif not promoter and "-" not in m:
                                    matched = True
                                    if minor is not None:
                                        this_match.append(g + "@" + m + ":" + minor)
                                    else:
                                        this_match.append(g + "@" + m)
                        check = check and matched
                    elif "!" in mutation:
                        matched = False
                        for g, m, minor in early_stop:
                            if g == gene and (
                                (minor is None and not rule_is_minor)
                                or (minor is not None and rule_is_minor)
                            ):
                                matched = True
                                if minor is not None:
                                    this_match.append(g + "@" + m + ":" + minor)
                                else:
                                    this_match.append(g + "@" + m)
                        check = check and matched
                    elif "=" in mutation:
                        matched = False
                        for g, m, minor in snps:
                            if (
                                g == gene
                                and m[0] == m[-1]
                                and (
                                    (minor is None and not rule_is_minor)
                                    or (minor is not None and rule_is_minor)
                                )
                            ):
                                matched = True
                                if minor is not None:
                                    this_match.append(g + "@" + m + ":" + minor)
                                else:
                                    this_match.append(g + "@" + m)
                        check = check and matched
                    else:
                        # Wildcard indels never have associated numbers or bases (as that wouldn't make sense)
                        # So check for matches
                        matched = False
                        if "-" in mutation:
                            promoter = True
                        else:
                            promoter = False

                        if "ins" in mutation or "del" in mutation:
                            if "ins" in mutation:
                                searching = "ins"
                            else:
                                searching = "del"
                            for g, m, minor in indels:
                                if g == gene and (
                                    (minor is None and not rule_is_minor)
                                    or (minor is not None and rule_is_minor)
                                ):
                                    if promoter and "-" not in m:
                                        continue
                                    elif not promoter and "-" in m:
                                        continue
                                    if searching in m:
                                        matched = True
                                        if minor is not None:
                                            this_match.append(g + "@" + m + ":" + minor)
                                        else:
                                            this_match.append(g + "@" + m)
                        elif "fs" in mutation:
                            # Bit more annoying here as we have to check if specific indels are framshifting
                            for g, m, minor in indels:
                                if g == gene and (
                                    (minor is None and not rule_is_minor)
                                    or (minor is not None and rule_is_minor)
                                ):
                                    # Don't need to check for promoters here as promoter framshift doesn't make sense
                                    bases_in = m.split("_")[-1]
                                    if bases_in.isnumeric():
                                        # Number of bases rather than actual bases
                                        bases = int(bases_in)
                                    else:
                                        bases = len(bases_in)
                                    if bases % 3 != 0:
                                        matched = True
                                        if minor is not None:
                                            this_match.append(g + "@" + m + ":" + minor)
                                        else:
                                            this_match.append(g + "@" + m)
                        else:
                            # Only mixed left
                            for g, m, minor in indels:
                                if g == gene and (
                                    (minor is None and not rule_is_minor)
                                    or (minor is not None and rule_is_minor)
                                ):
                                    if "mixed" in m:
                                        matched = True
                                        if minor is not None:
                                            this_match.append(g + "@" + m + ":" + minor)
                                        else:
                                            this_match.append(g + "@" + m)
                        check = check and matched

                elif "?" in mut:
                    # Specific wildcard SNP, so match on everything except the alt
                    mut = mut[:-1]
                    for g, m, minor in snps:
                        if (
                            g == gene
                            and m[:-1] == mut
                            and (
                                (minor is None and not rule_is_minor)
                                or (minor is not None and rule_is_minor)
                            )
                        ):
                            check = check and True
                            if minor is not None:
                                this_match.append(g + "@" + m + ":" + minor)
                            else:
                                this_match.append(g + "@" + m)

                elif large_del_re.fullmatch(mutation):
                    for g, m, minor in large_dels:
                        if g == gene and (
                            (minor is None and not rule_is_minor)
                            or (minor is not None and rule_is_minor)
                        ):
                            check = check and True
                            if minor is not None:
                                this_match.append(g + "@" + m + ":" + minor)
                            else:
                                this_match.append(g + "@" + m)
                else:
                    # Exact part, so check for it's existance
                    check = check and mutation in joined
                    if check:
                        # Exists, so add it to this multi
                        this_match.append(mutation)

                if check:
                    multi_match.append(this_match)
                else:
                    # Give up at first hurdle
                    break
            if check:
                # Add all if valid
                partials = multi_match[0]
                for item in multi_match[1:]:
                    these_partials = []
                    for mutation in item:
                        this_partial = []
                        for p in partials:
                            this_partial.append(p + "&" + mutation)
                        these_partials += this_partial
                    partials = these_partials

                for p in partials:
                    new_mutations.append((None, p))

        else:
            # Exact multi, so just check for existance
            for mutation in multi.split("&"):
                check = check and mutation in joined
            if check:
                # This exact multi mutation exists, so add it to the mutations list
                new_mutations.append((None, multi))
    if just_joined:
        return new_mutations
    return mutations_ + new_mutations


def getMutations(
    mutations_df: pd.DataFrame | None,
    catalogue: piezo.ResistanceCatalogue,
    reference: grumpy.Genome,
) -> List[Tuple[str | None, str]]:
    """Get all of the mutations (including multi-mutations) from the mutations df
    Multi-mutations currently only exist within the converted WHO catalogue, and are a highly specific combination
        of mutations which must all be present for a single resistance value.

    Args:
        mutations_df (pd.DataFrame): Mutations dataframe
        catalogue (piezo.ResistanceCatalogue): The resistance catalogue. Used to find which multi-mutations we care about
        reference (grumpy.Genome): Reference genome object. Used for checking if mutations are in coding regions

    Returns:
        List[Tuple[str | None, str]]: List of [gene, mutation] or in the case of multi-mutations, [None, multi-mutation]
    """
    if mutations_df is None:
        return []
    mutations: List[Tuple[str | None, str]] = list(
        zip(mutations_df["gene"], mutations_df["mutation"])
    )
    # Grab the multi-mutations from the catalogue
    # By doing this, we can check a fixed sample space rather than every permutation of the mutations
    # This makes the problem tractable, but does not address a possible issue with multi-mutations not encapsulating full codons
    multis = set(
        catalogue.catalogue.rules[
            catalogue.catalogue.rules["MUTATION_TYPE"] == "MULTI"
        ]["MUTATION"]
    )
    if len(multis) > 0:
        # We have a catalogue including multi rules, so check if any of these are present in the mutations
        mutations = subset_multis(multis, mutations)

    # Check if the catalogue supports large deletions
    if "GENE" in set(catalogue.catalogue.rules["MUTATION_AFFECTS"]):
        large_dels = True
    else:
        large_dels = False

    # Filtering out *just* nucleotide changes for cases of synon mutations
    # The important part of these should have already been found by multi-mutations
    fixed = []
    for gene, mutation in mutations:
        if gene is not None and reference.get_gene(gene).coding:
            # Codes protein so check for nucleotide changes
            nucleotide = re.compile(
                r"""
                [acgtzx][0-9]+[acgtzx]
                """,
                re.VERBOSE,
            )
            if nucleotide.fullmatch(mutation):
                # Is a nucleotide (non-promoter) mutation in a coding gene
                # So skip it as it may cause prediction problems
                continue
        # Remove large dels if not supported
        if not large_dels:
            # Check if this is a large del
            large = re.compile(
                r"""
                del_(1\.0)|(0\.[0-9][0-9]?[0-9]?)
                """,
                re.VERBOSE,
            )
            if large.fullmatch(mutation):
                continue
        fixed.append((gene, mutation))
    return sorted(fixed, key=lambda x: "".join([str(i) for i in x]))


def epistasis(
    mutations: list[tuple[str | None, str]],
    resistanceCatalogue: piezo.ResistanceCatalogue,
    phenotype: dict,
    effects: dict,
    effectsCounter: int,
    vcfStem: str,
) -> int:
    """Add epistasis rules, overriding existing predictions as required

    Args:
        mutations (list[tuple[str | None, str]]): Sample's mutations. Items are tuples of (gene, mutation)
        resistanceCatalogue (piezo.ResistanceCatalogue): Resistance catalogue to use
        phenotype (dict): Dictionary of phenotypes
        effects (dict): Dictionary of effects
        effectsCounter (int): Position to add to effects dict
        vcfStem (str): Stem of the VCF file. Used as the `uniqueid`

    Returns:
        int: Incremented effects counter

    Implicit returns:
        phenotypes: dict is updated in place
        effects: dict is updated in place
    """
    epi_rules = set(
        resistanceCatalogue.catalogue.rules[
            resistanceCatalogue.catalogue.rules["MUTATION_TYPE"] == "EPISTASIS"
        ]["MUTATION"]
    )
    if len(epi_rules) > 0:
        # We have some epistasis rules so deal with them
        mutations = subset_multis(epi_rules, mutations, just_joined=True)
        seen_multis = set([effects[key][2] for key in effects.keys()])
        for _, mutation in mutations:
            prediction = resistanceCatalogue.predict(mutation, show_evidence=True)
            if isinstance(prediction, str):
                # prediction == "S" but mypy doesn't like that
                # Default prediction so ignore (not that this should happen here)
                continue
            for drug in prediction.keys():
                pred = prediction[drug]
                if isinstance(pred, str):
                    # Shouldn't be hit but mypy complains
                    pred = pred
                    evidence = None
                else:
                    pred, evidence = pred
                if phenotype[drug] != "F":
                    # F is the only value which overrides epistasis rules
                    phenotype[drug] = pred
                # Add to the dict if not already seen
                # It's possible that this exact multi already hit an epistasis rule
                if mutation not in seen_multis:
                    effects[effectsCounter] = [
                        vcfStem,
                        None,
                        mutation,
                        resistanceCatalogue.catalogue.name,
                        drug,
                        pred,
                        evidence,
                    ]
                    # Increment counter
                    effectsCounter += 1
    return effectsCounter


def populateEffects(
    outputDir: str,
    resistanceCatalogue: piezo.ResistanceCatalogue,
    mutations: pd.DataFrame,
    vcfStem: str,
    make_csv: bool,
    make_prediction_csv: bool,
    reference: grumpy.Genome,
    make_mutations_csv: bool = False,
    append: bool = False,
) -> Tuple[pd.DataFrame, Dict, pd.DataFrame] | None:
    """Populate and save the effects DataFrame as a CSV

    Args:
        outputDir (str): Path to the directory to save the CSV
        resistanceCatalogue (piezo.ResistanceCatalogue): Resistance catalogue for predictions
        mutations (pd.DataFrame): Mutations dataframe
        vcfStem (str): The basename of the given VCF - used as the sample name
        make_csv (bool): Whether to write the CSV of the dataframe
        make_prediction_csv (bool): Whether to write the CSV of the antibiogram
        reference (grumpy.Genome | None, optional): Reference genome. Defaults to None.
        make_mutations_csv (bool, optional): Whether to write the mutations CSV to disk with new mutations. Defaults to False.
        append (bool, optional): Whether to append data to an existing df at the location (if existing).

    Raises:
        InvalidMutationException: Raised if an invalid mutation is detected

    Returns:
        (pd.DataFrame, dict): (
            DataFrame containing the effects data,
            A metadata dictionary mapping drugs to their predictions,
            DataFrame containing the mutations data
        )
    """
    if resistanceCatalogue is None:
        logging.debug("Catalogue was None, skipping effects and predictions generation")
        return None
    # Assume wildtype behaviour unless otherwise specified
    phenotype = {drug: "S" for drug in resistanceCatalogue.catalogue.drugs}

    effects = {}
    effectsCounter = 0

    # Default prediction values are RFUS but use piezo catalogue's values if existing
    values = resistanceCatalogue.catalogue.values

    # only try and build an effects table if there are mutations
    sample_mutations = getMutations(mutations, resistanceCatalogue, reference)
    for gene, mutation in tqdm(sample_mutations):
        # Ensure its a valid mutation
        # if gene is not None and not referenceGenes[gene].valid_variant(mutation):
        #     logging.error(f"Not a valid mutation {gene}@{mutation}")
        #     raise InvalidMutationException(gene, mutation)

        # Get the prediction
        if gene is not None:
            prediction = resistanceCatalogue.predict(
                gene + "@" + mutation, show_evidence=True
            )
        else:
            # This is a multi-mutation so is already of required format
            prediction = resistanceCatalogue.predict(mutation, show_evidence=True)

        # If the prediction is interesting, iter through drugs to find predictions
        if prediction != "S" and not isinstance(prediction, str):
            for drug in prediction.keys():
                drug_pred = prediction[drug]
                if isinstance(drug_pred, str):
                    # This shouldn't happen because we're showing evidence
                    # Adding to appease mypy...
                    pred: str = drug_pred
                    evidence: Dict = {}
                else:
                    pred, evidence = drug_pred
                # Prioritise values based on order within the values list
                if values.index(pred) < values.index(phenotype[drug]):
                    # The prediction is closer to the start of the values list, so should take priority
                    phenotype[drug] = pred

                # Add to the dict
                effects[effectsCounter] = [
                    vcfStem,
                    gene,
                    mutation,
                    resistanceCatalogue.catalogue.name,
                    drug,
                    pred,
                    evidence,
                ]
                # Increment counter
                effectsCounter += 1

    # Check for epistasis rules (which ignore prediction heirarchy)
    effectsCounter = epistasis(
        sample_mutations,
        resistanceCatalogue,
        phenotype,
        effects,
        effectsCounter,
        vcfStem,
    )

    # Build the DataFrame
    effects_df = pd.DataFrame.from_dict(
        effects,
        orient="index",
        columns=[
            "uniqueid",
            "gene",
            "mutation",
            "catalogue_name",
            "drug",
            "prediction",
            "evidence",
        ],
    )
    effects_df = effects_df[
        [
            "uniqueid",
            "gene",
            "mutation",
            "drug",
            "prediction",
            "catalogue_name",
            "evidence",
        ]
    ]
    effects_df["catalogue_version"] = resistanceCatalogue.catalogue.version
    effects_df["prediction_values"] = "".join(resistanceCatalogue.catalogue.values)

    # Save as CSV
    if make_csv:
        if append:
            # Check to see if there's anything there already
            try:
                old_effects = pd.read_csv(
                    os.path.join(outputDir, f"{vcfStem}.effects.csv")
                )
                effects_df = pd.concat([old_effects, effects_df])
            except FileNotFoundError:
                pass

        effects_df.to_csv(
            os.path.join(outputDir, f"{vcfStem}.effects.csv"), index=False
        )

    effects_df.reset_index(inplace=True)

    if make_prediction_csv:
        # We need to construct a simple table here
        predictions = [phenotype[drug] for drug in resistanceCatalogue.catalogue.drugs]
        vals = {
            "uniqueid": vcfStem,
            "drug": resistanceCatalogue.catalogue.drugs,
            "prediction": predictions,
            "catalogue_name": resistanceCatalogue.catalogue.name,
            "catalogue_version": resistanceCatalogue.catalogue.version,
            "catalogue_values": "".join(resistanceCatalogue.catalogue.values),
        }
        predictions_df = pd.DataFrame(vals)
        if append:
            # Check to see if there's anything there already
            try:
                old_predictions = pd.read_csv(
                    os.path.join(outputDir, f"{vcfStem}.predictions.csv")
                )
                predictions_df = pd.concat([old_predictions, predictions_df])
            except FileNotFoundError:
                pass
        predictions_df.to_csv(
            os.path.join(outputDir, f"{vcfStem}.predictions.csv"), index=False
        )
    if len(effects) == 0:
        # We have no effects to report so populate empty df
        effects_df = pd.DataFrame.from_dict(effects)

    # Return  the metadata dict to log later
    return (
        effects_df,
        {drug: phenotype[drug] for drug in resistanceCatalogue.catalogue.drugs},
        mutations,
    )


def saveJSON(
    variants,
    mutations,
    effects,
    phenotypes: dict[str, str],
    path: str,
    guid: str,
    catalogue: piezo.ResistanceCatalogue,
    gnomonicusVersion: str,
    time_taken: float,
    reference: grumpy.Genome,
    vcf_path: str,
    reference_path: str,
    catalogue_path: str,
) -> None:
    """Create and save a single JSON output file for use within GPAS. JSON structure:
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
        ?'errors': {
            <gene name>: <stack trace>
        }
        'data': {
            'variants': [
                {
                    'variant': Genome level variant in GARC,
                    'nucleotide_index': Genome index of variant,
                    'gene_name': Name of the gene this variant affects (if applicable),
                    'gene_position': Gene position which this variant affects. Nucleotide number if non coding, codon indx if coding (if applicable),
                    'codon_idx': Index of the base within the corresponding codon this affects (if applicable),
                    'vcf_evidence': Parsed VCF row,
                    'vcf_idx': Which part of the VCF row to look at for this call
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
                        'prediction': Prediction caused by this mutation,
                        'evidence': Evidence to support this prediction. Currently placeholder
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
        reference (grumpy.Genome): Reference genome object
        vcf_path (str): Path to the VCF file used for this run
        reference_path (str): Path to the reference genome used for this run
        catalogue_path (str): Path to the catalogue used for this run
        minor_errors (dict): Mapping of gene name --> stack trace of any errors occurring when parsing minor mutations
    """
    # Define some metadata for the json
    meta = OrderedDict(
        {
            "status": "success",
            "workflow_name": "gnomonicus",
            "workflow_version": gnomonicusVersion,  # gnomonicus version used
            "workflow_task": "resistance_prediction",  # TODO: Update this when we know how to detect a virulence catalogue
            "guid": guid,  # Sample GUID
            "UTC-datetime-completed": datetime.datetime.utcnow().isoformat(),  # ISO datetime run
            "time_taken_s": time_taken,
            "reference": reference.name,
            "catalogue_file": catalogue_path,
            "reference_file": reference_path,
            "vcf_file": vcf_path,
        }
    )
    if catalogue is not None:
        meta["catalogue_type"] = "".join(catalogue.catalogue.values)
        meta["catalogue_name"] = catalogue.catalogue.name
        meta["catalogue_version"] = catalogue.catalogue.version
    else:
        meta["catalogue_type"] = None
        meta["catalogue_name"] = None
        meta["catalogue_version"] = None

    # Main data collection
    data: Dict = OrderedDict()

    # Antibigram field
    data["antibiogram"] = OrderedDict(
        [(key, phenotypes[key]) for key in sorted(phenotypes.keys())]
    )

    # Variants field
    _variants = []
    for _, variant in variants.iterrows():
        row = OrderedDict(
            {
                "variant": (
                    variant["variant"] if pd.notnull(variant["variant"]) else None
                ),
                "nucleotide_index": (
                    variant["nucleotide_index"]
                    if pd.notnull(variant["nucleotide_index"])
                    else None
                ),
                "gene_name": variant["gene"] if pd.notnull(variant["gene"]) else None,
                "gene_position": (
                    variant["gene_position"]
                    if pd.notnull(variant["gene_position"])
                    else None
                ),
                "codon_idx": (
                    variant["codon_idx"] if pd.notnull(variant["codon_idx"]) else None
                ),
                "vcf_evidence": json.loads(variant["vcf_evidence"]),
                "vcf_idx": (
                    variant["vcf_idx"] if pd.notnull(variant["vcf_idx"]) else None
                ),
            }
        )
        _variants.append(row)
    _variants = sorted(
        _variants,
        key=lambda x: (
            x["gene_name"] or "z",
            x["gene_position"] or 0,
            x["variant"] or "z",
        ),
    )
    data["variants"] = _variants

    # Depending on mutations/effects, populate
    _mutations = []
    if mutations is not None:
        for _, mutation in mutations.iterrows():
            row = OrderedDict(
                {
                    "gene": mutation["gene"] if pd.notnull(mutation["gene"]) else None,
                    "gene_position": (
                        mutation["gene_position"]
                        if pd.notnull(mutation["gene_position"])
                        else None
                    ),
                    "mutation": (
                        mutation["mutation"]
                        if pd.notnull(mutation["mutation"])
                        else None
                    ),
                }
            )
            if mutation["mutation"][0].isupper() or mutation["mutation"][0] == "!":
                # Only add codon ref/alt for AA changes
                row["ref"] = mutation["ref"] if pd.notnull(mutation["ref"]) else None
                row["alt"] = mutation["alt"] if pd.notnull(mutation["alt"]) else None
            _mutations.append(row)

    _mutations = sorted(
        _mutations,
        key=lambda x: (x["gene"] or "z", x["gene_position"] or 0),
    )
    data["mutations"] = _mutations

    _effects = defaultdict(list)
    if effects is not None and len(effects) > 0:
        for _, effect in effects.iterrows():
            prediction = OrderedDict(
                {
                    "gene": effect["gene"] if pd.notnull(effect["gene"]) else None,
                    "mutation": (
                        effect["mutation"] if pd.notnull(effect["mutation"]) else None
                    ),
                    "prediction": (
                        effect["prediction"]
                        if pd.notnull(effect["prediction"])
                        else None
                    ),
                    "evidence": effect["evidence"],
                }
            )
            _effects[effect["drug"]].append(prediction)

    for drug in _effects.keys():
        _effects[drug] = sorted(
            _effects[drug],
            key=lambda x: (x["gene"] or "z", x["mutation"] or "z"),
        )
        _effects[drug].append(OrderedDict({"phenotype": phenotypes[drug]}))

    data["effects"] = OrderedDict(
        [(key, _effects[key]) for key in sorted(_effects.keys())]
    )

    # Convert fields to a list so it can be json serialised
    with open(
        os.path.join(path, f"{guid}.gnomonicus-out.json"), "w", encoding="utf-8"
    ) as f:
        f.write(json.dumps({"meta": meta, "data": data}, indent=2))
