#!/usr/bin/env python

"""
SELFIES: a robust representation of semantically constrained graphs with an
         example application in chemistry.

SELFIES (SELF-referencIng Embedded Strings) is a general-purpose,
sequence-based, robust representation of semantically constrained graphs.
It is based on a Chomsky type-2 grammar, augmented with two self-referencing
functions. A main objective is to use SELFIES as direct input into machine
learning models, in particular in generative models, for the generation of
outputs with high validity.

The code presented here is a concrete application of SELFIES in chemistry, for
the robust representation of molecules.

    Typical usage example:
        import selfies

        benzene = "C1=CC=CC=C1"
        selfies_benzene = selfies.encoder(benzene)
        smiles_benzene = selfies.decoder(selfies_benzene)

For comments, bug reports or feature ideas, please send an email to
mario.krenn@utoronto.ca and alan@aspuru.com.
"""

__version__ = "1.0.1"

__all__ = [
    "encoder",
    "decoder",
    "get_semantic_robust_alphabet",
    "get_semantic_constraints",
    "set_semantic_constraints",
    "len_selfies",
    "split_selfies",
    "get_alphabet_from_selfies",
    "selfies_to_encoding",
    "batch_selfies_to_flat_hot",
    "encoding_to_selfies",
    "batch_flat_hot_to_selfies",
]

from .decoder import decoder
from .encoder import encoder
from .grammar_rules import (
    get_semantic_robust_alphabet,
    get_semantic_constraints,
    set_semantic_constraints,
)
from .utils import (
    get_alphabet_from_selfies,
    len_selfies,
    split_selfies,
    selfies_to_encoding,
    batch_selfies_to_flat_hot,
    encoding_to_selfies,
    batch_flat_hot_to_selfies,
)
