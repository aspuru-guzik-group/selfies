Interpreting SELFIES
=====================

This section is an informal tutorial on interpreting SELFIES; we aim to
provide an intuition on translating SELFIES into molecules. We use the verb
*interpret* to describe the action of interpreting a SELFIES
as a molecule. SELFIES are interpreted symbol-by-symbol (from left to right),
so atoms and bonds of the molecule are derived sequentially, one-by-one.

There are three main types of SELFIES symbols: atomic symbols, branch symbols,
and ring symbols. We describe each type of symbol in each section.


Atomic Symbols
##############

Atomic symbols are of the general form ``[<B><A>]``, where
``<B> in {'', '/', '\\', '=', '#'}`` is a bond symbol and ``<A>`` is
a SMILES symbol representing an atom or ion. If the SMILES symbol is
enclosed by square brackets (e.g. ``[13C]``), then the square brackets are
dropped and ``expl`` (for "explicit brackets") is appended to obtain ``<A>``.
For example:

.. table::
    :align: center

    +---------+---------------+--------------+----------------+
    | ``<B>`` | SMILES symbol | ``<A>``      | SELFIES symbol |
    +=========+===============+==============+================+
    | ``'='`` | ``N``         | ``N``        | ``[=N]``       |
    +---------+---------------+--------------+----------------+
    | ``''``  | ``[C@@H]``    | ``C@@Hexpl`` | ``[C@@Hexpl]`` |
    +---------+---------------+--------------+----------------+
    | ``'/'`` | ``[O+]``      | ``O+expl``   | ``[/O+expl]``  |
    +---------+---------------+--------------+----------------+

An atomic symbol ``[<B><A>]`` connects atom ``<A>`` to the previously
derived atom through bond type ``<B>``. If creating this bond would violate the
bond constraints of the previous or current atom, the bond multiplicity is
reduced (minimally) such that no constraints are violated. If the previous
atom cannot make at least one bond, then SELFIES interpretation terminates, and
the current and all subsequent symbols are ignored.

In the following examples, we restrict ``C`` to 4 bonds, ``S`` to 2 bonds,
and ``I`` to 1 bond:

.. table::
    :align: center

    +---------+-----------------------------+-----------------+
    | Example | SELFIES                     | SMILES          |
    +=========+=============================+=================+
    | 1       | ``[C][=C][C][#C][13Cexpl]`` | ``C=CC#C[13C]`` |
    +---------+-----------------------------+-----------------+
    | 2       | ``[C][I][C][C][C][C]``      | ``CI``          |
    +---------+-----------------------------+-----------------+
    | 3       | ``[C][S][=C][#S][C][I]``    | ``CSC=S``       |
    +---------+-----------------------------+-----------------+

**Discussion:** In example 2, since ``I`` was constrained to 1 bond,
the subsequent ``[C][C][C][C]`` was ignored. In example 3, since ``S``
was constrained to 2 bonds, the triple bond of ``[#S]`` is first reduced
to a double bond. Then, the subsequent ``[C][I]`` is ignored because ``=S``
cannot make any more bonds.


Branch Symbols
##############

Branch symbols are of the general form ``[Branch<L>_<X>]``, where
``<L>, <X> in {1, 2, 3}``. A branch symbol specifies a branch from the
main chain, analogous to the open and closed curved brackets in SMILES.
After a branch symbol, the next ``<L>`` symbols are first read
and converted into indices according to this assignment:

.. table::
    :align: center

    +-------+-----------------+-------+-----------------+
    | Index | Symbol          | Index | Symbol          |
    +=======+=================+=======+=================+
    | 0     | ``[C]``         | 8     | ``[Branch2_3]`` |
    +-------+-----------------+-------+-----------------+
    | 1     | ``[Ring1]``     | 9     | ``[O]``         |
    +-------+-----------------+-------+-----------------+
    | 2     | ``[Ring2]``     | 10    | ``[N]``         |
    +-------+-----------------+-------+-----------------+
    | 3     | ``[Branch1_1]`` | 11    | ``[=N]``        |
    +-------+-----------------+-------+-----------------+
    | 4     | ``[Branch1_2]`` | 12    | ``[=C]``        |
    +-------+-----------------+-------+-----------------+
    | 5     | ``[Branch1_3]`` | 13    | ``[#C]``        |
    +-------+-----------------+-------+-----------------+
    | 6     | ``[Branch2_1]`` | 14    | ``[S]``         |
    +-------+-----------------+-------+-----------------+
    | 7     | ``[Branch2_2]`` | 15    | ``[P]``         |
    +-------+-----------------+-------+-----------------+
    | All other symbols assigned index 0.               |
    +-------+-----------------+-------+-----------------+

Then, the indices are read as a hexadecimal (base 16) integer ``N``.

.. table::
    :align: center

    +----------------------+---------------------------------------------+
    | SELFIES              | ``...[Branch3_1][Ring1][C][O][C][C]...``    |
    +----------------------+---------------------------------------------+
    | Next ``<L>`` symbols | ``[Ring1][C][O]``                           |
    +----------------------+---------------------------------------------+
    | Indices              | 1, 0, 9                                     |
    +----------------------+---------------------------------------------+
    | ``N``                | ``1 * (16 ** 2) + 0 * (16 ** 1) + 9 = 265`` |
    +----------------------+---------------------------------------------+

The next ``N + 1`` symbols (after the ``<L>`` symbols used to compute ``N``)
are treated as a separate SELFIES and recursively interpreted. Finally, the
branch fragment is connected to the previously derived atom. SELFIES
interpretation proceeds with the next symbol (after the ``<L> + N + 1``
symbols used to create the branch). Although we have recursively derived
new atoms that are in the branch, in practice, we still define the
"previously derived atom" as the previously derived atom
before the branch interpretation.

Work in progress... : )

Ring Symbols
############

Ring symbols are of the general form ``[Ring<L>]``, where ``<L> in {1, 2, 3}``.

Work in progress... : )

Special Symbols
###############

This subsection describes various special symbols that :mod:`selfies`
ascribes special meaning to.

+---------------+----------------------------------------------------------------------------------------------+
| Character     | Description                                                                                  |
+===============+==============================================================================================+
| ``[nop]``     | The nop (no operation) symbol is always ignored and skipped over by :func:`selfies.decoder`. |
|               |                                                                                              |
|               | Thus, it can be used as a padding symbol for SELFIES.                                        |
+---------------+----------------------------------------------------------------------------------------------+
| ``.``         | The dot symbol is used to indicate disconnected or ionic compounds, similar to how it is     |
|               |                                                                                              |
|               | used in SMILES.                                                                              |
+---------------+----------------------------------------------------------------------------------------------+
