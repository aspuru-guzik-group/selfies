#!/usr/bin/env python

"""
SELFIES: a robust representation of semantically constrained graphs with an
         example application in chemistry
         TODO: change this date below
         v1.0.0, 01. October 2019
by Mario Krenn, Florian Haese, AkshatKuman Nigam, Pascal Friederich,
Alan Aspuru-Guzik

SELFIES (SELF-referencing Embedded Strings) is a general-purpose,
sequence-based, robust representation of semantically constrained graphs.
It is based on a Chomsky type-2 grammar, augmented with two self-referencing
functions. A main objective is to use SELFIES as direct input into machine
learning models, in particular in generative models, for the generation of
outputs with high validity.

The code presented here is a concrete application of SELIFES in chemistry, for
the robust representation of molecule. We show the encoding and decoding of
three molecules from various databases, and the generation of a new, random
molecule with high semantic and syntactical validity.

For comments, bug reports or feature ideas, please send an email to
mario.krenn@utoronto.ca and alan@aspuru.com
"""

__version__ = "1.0.0"

from .decoder import decoder
from .encoder import encoder
from .utils import get_selfies_alphabet, get_atom_dict, set_selfies_alphabet
