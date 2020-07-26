Interpreting SELFIES
=====================

This section is an informal tutorial on interpreting SELFIES; we aim to
provide an intuition on translating SELFIES into molecules.

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

+---------+---------------+--------------+----------------+
| ``<B>`` | SMILES symbol | ``<A>``      | SELFIES symbol |
+=========+===============+==============+================+
| ``'='`` | ``N``         | ``N``        | ``[=N]``       |
+---------+---------------+--------------+----------------+
| ``''``  | ``[C@@H]``    | ``C@@Hexpl`` | ``[C@@Hexpl]`` |
+---------+---------------+--------------+----------------+
| ``'/'`` | ``[O+]``      | ``O+expl``   | ``[/O+expl]``  |
+---------+---------------+--------------+----------------+

An atomic symbol ``[<B><A>]`` bonds atom ``<A>`` to the previous atom
through bond type ``<B>``. If creating this bond would violate the
bond constraints of the previous or current atom, the bond multiplicity is
reduced (minimally) such that no constraints are violated. If the previous
atom cannot make any more bonds, then SELFIES derivation terminates, and
the current and all subsequent symbols are ignored.

In the following examples, we restrict ``C`` to 4 bonds, ``S`` to 2 bonds,
and ``I`` to 1 bond:

+---------+-----------------------------+-----------------+
| Example | SELFIES                     | SMILES          |
+=========+=============================+=================+
| 1       | ``[C][=C][C][#C][13Cexpl]`` | ``C=CC#C[13C]`` |
+---------+-----------------------------+-----------------+
| 2       | ``[C][I][C][C][C][C]``      | ``CI``          |
+---------+-----------------------------+-----------------+
| 3       | ``[C][S][=C][#S][C][I]``    | ``CSC=S``       |
+---------+-----------------------------+-----------------+

**Notes:** In example 2, note that since we restricted ``I`` to 1 bond,
the subsequent ``[C][C][C][C]`` is ignored. In example 3, since we
restricted ``S`` to 2 bonds, the bond multiplicity of ``[#S]`` is reduced
to a double bond. Then, all subsequent symbols are ignored because ``S``
cannot make any more bonds.

Branch Symbols
##############

Branch symbols are of the general form ``[BranchL_X]``, where
``L, X in {1, 2, 3}``.

Ring Symbols
############

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
