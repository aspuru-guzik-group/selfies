API Reference
==================

.. currentmodule:: selfies

Core Functions
------------------------
.. autofunction:: encoder
.. autofunction:: decoder

Customization Functions
------------------------

The SELFIES grammar is derived dynamically from a set of semantic constraints,
which assign bonding capacities to various atoms.
By default, :mod:`selfies` operates under the following constraints:

.. table::
    :align: center

    +-----------+------------------------------+
    | Max Bonds | Atom(s)                      |
    +===========+==============================+
    | 1         | F, Cl, Br, I                 |
    +-----------+------------------------------+
    | 2         | O                            |
    +-----------+------------------------------+
    | 3         | B, N                         |
    +-----------+------------------------------+
    | 4         | C                            |
    +-----------+------------------------------+
    | 5         | P                            |
    +-----------+------------------------------+
    | 6         | S                            |
    +-----------+------------------------------+
    | 8         | All other atoms              |
    +-----------+------------------------------+

The +1 and -1 charged versions of O, N, C, S, and P are also constrained,
where a +1 increases the bonding capacity of the neutral atom by 1,
and a -1 decreases the bonding capacity of the neutral atom by 1.
For example, N+1 has a bonding capacity of :math:`3 + 1 = 4`,
and N-1 has a bonding capacity of :math:`3 - 1 = 2`. The charged versions
B+1 and B-1 are constrained to a capacity of 2 and 4 bonds, respectively.

However, the default constraints are inadequate for SMILES strings that violate them. For
example, nitrobenzene ``O=N(=O)C1=CC=CC=C1`` has a nitrogen with 6 bonds and
the chlorate anion ``O=Cl(=O)[O-]`` has a chlorine with 5 bonds - these
SMILES strings *cannot* be represented by SELFIES strings under the default constraints.
Additionally, users may want to specify their own custom constraints. Thus, we
provide the following methods for configuring the semantic constraints
of :mod:`selfies`.

.. warning::

    SELFIES strings may be translated differently under different semantic constraints.
    Therefore, if custom semantic constraints are used, it is recommended to report
    them for reproducibility reasons.

.. autofunction:: get_preset_constraints
.. autofunction:: get_semantic_constraints
.. autofunction:: set_semantic_constraints


Utility Functions
------------------------
.. autofunction:: len_selfies
.. autofunction:: split_selfies
.. autofunction:: get_alphabet_from_selfies
.. autofunction:: get_semantic_robust_alphabet
.. autofunction:: selfies_to_encoding
.. autofunction:: encoding_to_selfies
.. autofunction:: batch_selfies_to_flat_hot
.. autofunction:: batch_flat_hot_to_selfies


Exceptions
------------------------
.. autoexception:: EncoderError
.. autoexception:: DecoderError
