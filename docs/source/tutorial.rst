SELFIES Derivation
==================

This section is an informal tutorial on how molecules are derived
from a SELFIES. The SELFIES grammar has non-terminal symbols or states

.. math::

    X_0, \ldots, X_7, Q

Derivation starts with state :math:`X_0`. The SELFIES is read symbol-by-symbol,
with each symbol specifying a grammar rule. SELFIES derivation terminates
when no non-terminal symbols remain. In each subsection, we describe a type of
SELFIES symbol and the grammar rules associated with it.

Atomic Symbols
##############

Atomic symbols are of the general form ``[<B><A>]``, where
``<B> in {'', '/', '\\', '=', '#'}`` is a prefix representing a bond,
and ``<A>`` is a SMILES symbol representing an atom or ion.
If the SMILES symbol is enclosed by square brackets (e.g. ``[13C]``),
then the square brackets are dropped and ``expl`` (for "explicit brackets")
is appended to obtain ``<A>``. For example:

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

Let atomic symbol ``[<B><A>]`` be given, where ``<B>`` is a prefix
representing a bond with multiplicity :math:`\beta` and ``<A>`` is an atom
that can make :math:`\alpha` bonds maximally. The atomic symbol maps:

.. math::

    X_i \to \begin{cases}
        \texttt{<A>} X_{\alpha} & i = 0 \\
        \texttt{<B'><A>} X_{\alpha - \mu} & i \neq 0
    \end{cases}

where ``<B'>`` is a prefix representing a bond with multiplicity
:math:`\mu = \min(\beta, \alpha, i)`, and where we replace
:math:`X_{\alpha - \mu}` above with the empty string
if :math:`\alpha - \mu = 0`. Note that non-terminal states :math:`X_{1-7}`
effectively restrict the subsequent bond to a multiplicity of at most :math:`i`.
We provide an example of the derivation of the
SELFIES ``[F][=C][=C][#N]``:

.. math::

    X_0 \to \texttt{F}X_1 \to \texttt{FC}X_3 \to \texttt{FC=C}X_2 \to \texttt{FC=C=N}


**Discussion:** Intuitively, an atomic symbol ``[<B><A>]`` connects
atom ``<A>`` to the previously derived atom through bond type ``<B>``.
If creating this bond would violate the bond constraints of the previous
or current atom, the bond multiplicity is reduced (minimally) such that no
constraints are violated.

**Examples:**

.. table::
    :align: center

    +---------+-----------------------------+-----------------+
    | Example | SELFIES                     | SMILES          |
    +=========+=============================+=================+
    | 1       | ``[C][=C][C][#C][13Cexpl]`` | ``C=CC#C[13C]`` |
    +---------+-----------------------------+-----------------+
    | 2       | ``[C][F][C][C][C][C]``      | ``CF``          |
    +---------+-----------------------------+-----------------+
    | 3       | ``[C][O][=C][#O][C][F]``    | ``COC=O``       |
    +---------+-----------------------------+-----------------+

Index Symbols
#############

The state :math:`Q` is used to specify the size of branches,
and the location of ring bonds. A sequence of any SELFIES symbol(s)
:math:`s_1, \ldots, s_{\texttt{<L>}}` maps state :math:`Q` to an integer.
First, each symbol :math:`s_i` is converted to an index :math:`\text{idx}(s_i)`
according to following assignment.

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

Then :math:`Q` is mapped to the hexadecimal (base 16) integer specified
by the indices. For example, in a sequence of three SELFIES symbols,
:math:`Q` is mapped to:

.. math::

    Q \to (\text{idx}(s_1) \times 16^2) + (\text{idx}(s_2) \times 16) + \text{idx}(s_3)


Branch Symbols
##############

Branch symbols are of the general form ``[Branch<L>_<M>]``, where
``<L>, <M> in {1, 2, 3}``. A branch symbol specifies a branch from the
main chain, analogous to the open and closed curved brackets in SMILES.

A Branch symbol ``[Branch<L>_<M>]`` maps:

.. math::

    X_i \to \begin{cases}
        X_i & i \leq 1 \\
        B(Q, X_{n})X_j & i > 1
    \end{cases}

where :math:`n = \min(i - 1, \texttt{<M>})` is the initial branch
derivation state and :math:`j = i - n` is the next derivation state. In the
bottom case, the ``<L>`` symbols after the Branch symbol are read,
and used to map :math:`Q` to an index. Then
:math:`B(Q, X_{n})` takes the next :math:`Q + 1` symbols, and recursively
derives them with initial derivation state :math:`X_{n}`. The resulting
fragment is taken to be the derived branch, and derivation proceeds
with the next derivation state :math:`X_j`.

**Discussion:**  Intuitively, branch symbols are skipped for states
:math:`X_{0-1}` because the previous atom can make at most one bond
(branches require at least two bonds to be free). It is possible
that a branch is nested at the start of another branch; in SELFIES, both
branches will be connected to the same main chain atom (see Example 5).

**Examples:**

+---------+-------------------------------------------------------+---------------+-------------------------+
| Example | SELFIES                                               | :math:`Q + 1` | SMILES                  |
+=========+=======================================================+===============+=========================+
| 1       | ``[C][Branch1_1][C][F][Cl]``                          | 1             | ``C(F)Cl``              |
+---------+-------------------------------------------------------+---------------+-------------------------+
| 2       | ``[C][Branch1_2][Ring2][=C][C][C][Cl]``               | 3             | ``C(=CCC)Cl``           |
+---------+-------------------------------------------------------+---------------+-------------------------+
| 3       | ``[S][Branch1_2][C][=O][Branch1_2][C]``               | 1, 1, 1       | ``S(=O)(=O)([O-])[O-]`` |
|         |                                                       |               |                         |
|         | ``[=O][Branch1_1][C][O-expl][O-expl]``                |               |                         |
+---------+-------------------------------------------------------+---------------+-------------------------+
| 4       | ``[C][Branch2_1][Ring1][Branch1_2][C]``               | 21            | ``C(CC...CC)F``         |
|         |                                                       |               |                         |
|         | ``[C][C][C][C][C][C][C][C][C][C][C][C]``              |               |                         |
|         |                                                       |               |                         |
|         | ``[C][C][C][C][C][C][C][C][F]``                       |               |                         |
|         +-------------------------------------------------------+---------------+-------------------------+
|         | Example 4 has a single branch of 21 carbon atoms.                                               |
+---------+-------------------------------------------------------+---------------+-------------------------+
| 5       | ``[C][Branch1_2][Branch1_1][Branch1_1][C][C][Cl][F]`` | 4, 1          | ``C(C)(Cl)F``           |
+---------+-------------------------------------------------------+---------------+-------------------------+


Ring Symbols
############

Ring symbols are of the general form ``[Ring<L>]`` or ``[Expl<B>Ring<L>]``,
where ``<L> in {1, 2, 3}`` and ``<B> in {'/', '\\', '=', '#'}`` is a
prefix representing a bond. A ring symbol specifies a ring bond between two
atoms, analogous to the ring numbering digits in SMILES.

A Ring symbol ``[Ring<L>]`` maps:

.. math::

    X_i \to \begin{cases}
        X_i & i = 0 \\
        R(Q)X_i & i \neq 0
    \end{cases}

Identical to branch derivation, the ``<L>`` symbols after the Ring symbol are read,
and used to map :math:`Q` to an index. Then :math:`R(Q)` connects the current
atom to the :math:`(Q + 1)`-th preceding atom through a single bond.
More specifically, the "current" atom is the most recently derived atom,
excluding atoms derived in branches (see Example 5 below); and the ":math:`(Q + 1)`-th
preceding atom" is the atom derived :math:`Q + 1` atoms before the
current atom, counting atoms derived in branches. If the preceding atom does
not exist, then the connection is made to the 1st derived atom instead.

The Ring symbol ``[Expl<B>Ring<L>]`` has an equivalent function to
``[Ring<L>]``, except that it connects the current and :math:`(Q + 1)`-th
preceding atom through a bond of type ``<B>``.


**Discussion**: In practice, ring bonds are created during a second pass,
after all atoms and branches have been derived. The candidate ring
bonds are temporarily stored in a queue, and then made in
the order that they appear in the SELFIES. A ring bond will be made if
its connected atoms can make the ring bond without violating any
bond constraints.

It is also possible that the current atom is already bonded to the
:math:`(Q + 1)`-th preceding atom, e.g. if :math:`Q = 0`. In this case,
the multiplicity of the existing bond is increased by the minimum
of (1) the multiplicity of the ring bond candidate and (2) the number
of free bonds of both connected atoms (see Example 6). Note that
the resulting bond will be constrained to a multiplicity of at most 3.

**Examples:**

+---------+------------------------------------------------------------+---------------+------------------+
| Example | SELFIES                                                    | :math:`Q + 1` | SMILES           |
+=========+============================================================+===============+==================+
| 1       | ``[C][=C][C][=C][C][=C][Ring1][Branch1_2]``                | 5             | ``C1=CC=CC=C1``  |
+---------+------------------------------------------------------------+---------------+------------------+
| 2       | ``[C][C][=C][C][=C][C][Expl=Ring1][Branch1_2]``            | 5             | ``C=1C=CC=CC=1`` |
+---------+------------------------------------------------------------+---------------+------------------+
| 3       | ``[C][C][Expl=Ring1][C]``                                  | 1             | ``C#C``          |
+---------+------------------------------------------------------------+---------------+------------------+
| 4       | ``[C][C][C][C][C][C][C][C][C][C][C]``                      | 21            | ``C1CC...CC1``   |
|         |                                                            |               |                  |
|         | ``[C][C][C][C][C][C][C][C][C][C][C]``                      |               |                  |
|         |                                                            |               |                  |
|         | ``[Ring2][Ring1][Branch1_2]``                              |               |                  |
|         +------------------------------------------------------------+---------------+------------------+
|         | Example 4 is a single carbon ring of 22 carbon atoms.                                         |
+---------+------------------------------------------------------------+---------------+------------------+
| 5       | ``[C][C][C][C][Branch1_1][C][C][Ring1][Ring2][C][C]``      | 3             | ``C1CCC1(C)CC``  |
|         +------------------------------------------------------------+---------------+------------------+
|         | Note that the SMILES ``CC1CC(C1)CC`` is not outputted.                                        |
+---------+------------------------------------------------------------+---------------+------------------+
| 6       | ``[C][C][C][C][Expl=Ring1][Ring2][Expl#Ring1][Ring2]``     | 3, 3          | ``C#1CCC#1``     |
+---------+------------------------------------------------------------+---------------+------------------+



Special Symbols
###############

The following are symbols that have a special meaning for SELFIES:

.. _no operation: https://en.wikipedia.org/wiki/NOP_(code)

+---------------+-------------------------------------------------------------------------------------------------+
| Character     | Description                                                                                     |
+===============+=================================================================================================+
| ``[epsilon]`` | The ``[epsilon]`` symbol maps :math:`X_0 \to X_0` and :math:`X_i \to \epsilon` (the empty       |
|               | string) for all :math:`i \geq 1`.                                                               |
+---------------+-------------------------------------------------------------------------------------------------+
| ``[nop]``     | The nop (`no operation`_) symbol is always ignored and skipped over by :func:`selfies.decoder`. |
|               |                                                                                                 |
|               | Thus, it can be used as a padding symbol for SELFIES.                                           |
+---------------+-------------------------------------------------------------------------------------------------+
| ``.``         | The dot symbol is used to indicate disconnected or ionic compounds, similar to how it is        |
|               |                                                                                                 |
|               | used in SMILES.                                                                                 |
+---------------+-------------------------------------------------------------------------------------------------+
