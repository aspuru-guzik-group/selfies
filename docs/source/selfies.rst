Code Documentation
==================

.. currentmodule:: selfies

Standard Functions
------------------------
.. autofunction:: encoder
.. autofunction:: decoder
.. autofunction:: len_selfies
.. autofunction:: split_selfies
.. autofunction:: get_alphabet_from_selfies
.. autofunction:: get_semantic_robust_alphabet

Advanced Functions
------------------------

By default, :mod:`selfies` operates under the following semantic constraints

.. table::
    :align: center

    +-----------+------------------------------+
    | Max Bonds | Atom(s)                      |
    +===========+==============================+
    | 1         | ``F``, ``Cl``, ``Br``, ``I`` |
    +-----------+------------------------------+
    | 2         | ``O``                        |
    +-----------+------------------------------+
    | 3         | ``N``                        |
    +-----------+------------------------------+
    | 4         | ``C``                        |
    +-----------+------------------------------+
    | 5         | ``P``                        |
    +-----------+------------------------------+
    | 6         | ``S``                        |
    +-----------+------------------------------+
    | 8         | All other atoms              |
    +-----------+------------------------------+

However, the default constraints are inadequate for SMILES that violate them. For
example, nitrobenzene ``O=N(=O)C1=CC=CC=C1`` has a nitrogen with 6 bonds and
the chlorate anion ``O=Cl(=O)[O-]`` has a chlorine with 5 bonds - these
SMILES *cannot* be represented as SELFIES under the default constraints.
Additionally, users may want to specify their own custom constraints. Thus, we
provide the following methods for configuring the semantic constraints
of :mod:`selfies`.

.. warning::

    SELFIES may be translated differently under different semantic constraints.
    Therefore, if custom semantic constraints are used, it is recommended to report
    them for reproducibility reasons. 

.. autofunction:: get_semantic_constraints
.. autofunction:: set_semantic_constraints

