'''gnomon.py is a library providing functions which pull together output VCF of the Lodestone TB pipeline
    with a reference genome and a resistance catalogue, and utilise gumpy and
    piezo to produce variants, mutations and an antibiogram.

Provides a CLI script (bin/gnomon) which links these functions together to produce all outputs from the inputs

Classes:
    MissingFieldException
    NoVariantsException
    InvalidMutationException

Functions:
    loadGenome
    populateVariants
    populateMutations
    populateEffects
    assignMutationBools
    countNucleotideChanges
'''

__version__ = "0.1.0"

from .gnomon import (loadGenome, populateVariants, populateMutations, populateEffects, assignMutationBools,
                        countNucleotideChanges, MissingFieldException, NoVariantsException, InvalidMutationException)
