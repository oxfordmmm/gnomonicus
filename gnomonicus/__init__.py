"""gnomonicus is a library providing functions which pull together output VCF of the Lodestone TB pipeline
    with a reference genome and a resistance catalogue, and utilise gumpy and
    piezo to produce variants, mutations and an antibiogram.

Provides a CLI script (bin/gnomonicus) which links these functions together to produce all outputs from the inputs.
Makes the assumption that VCF files are named `<GUID>.vcf`

Classes:
    * MissingFieldException
    * InvalidMutationException

Functions:
    * loadGenome
    * loadGenomeAndGenes
    * populateVariants
    * populateMutations
    * populateEffects
    * assignMutationBools
    * countNucleotideChanges
    * saveJSON
"""

import importlib.metadata

__version__ = importlib.metadata.version("gnomonicus")

from .gnomonicus_lib import (
    InvalidMutationException,  # noqa: F401
    loadGenome,  # noqa: F401
    loadGenomeAndGenes,  # noqa: F401
    minority_population_variants,  # noqa: F401
    populateEffects,  # noqa: F401
    populateMutations,  # noqa: F401
    populateVariants,  # noqa: F401
    saveJSON,  # noqa: F401
)
