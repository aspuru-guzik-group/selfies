Reading SELFIES
===============

This section is an informal tutorial on *reading* SELFIES. Reading a
SELFIES is the action of deriving a molecule from a SELFIES. SELFIES are
read symbol-by-symbol (from left to right), so atoms and bonds of the molecule
are derived sequentially, one-by-one. We begin by describing the various
types of SELFIES symbols, and the form they come in. Then, we delve into the
SELFIES derivation process, which involves discussion about the
SELFIES grammar rules.

----------

SELFIES Symbols
###############

There are three main types of SELFIES symbols: atomic symbols, branch symbols,
and ring symbols. Additionally, there are a few symbols that :mod:`selfies`
ascribes special meaning to.

Atomic Symbols
**************

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


Branch Symbols
**************

Branch symbols are of the general form ``[Branch<L>_<M>]``, where
``<L>, <M> in {1, 2, 3}``. A branch symbol specifies a branch from the
main chain, analogous to the open and closed curved brackets in SMILES.

Ring Symbols
************

Ring symbols are of the general form ``[Ring<L>]``, where ``<L> in {1, 2, 3}``.
A ring symbol specifies a ring bond between two atoms, analogous to the
ring numbering digits in SMILES.


Special Symbols
***************

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

----------

SELFIES Derivation
##################

The derivation of a molecule from a SELFIES occurs in two steps. First,
the main chain and projecting branches - the main scaffold - is derived.
Then, ring bonds are added between atoms of the main scaffold.

Step 1: Main Scaffold
*********************

The SELFIES grammar has non-terminal symbols or states

.. math::

    X_0, \ldots, X_7, X_{9991}, X_{9992}, X_{9993}, Q_1, Q_2, Q_3

Derivation starts with state :math:`X_0`. The SELFIES is read symbol-by-symbol,
with each symbol specifying a grammar rule. SELFIES derivation terminates
when no non-terminal symbols remain. We now describe the grammar rules
associated with each type of SELFIES symbol.

**Atomic Symbol:** Let atomic symbol ``[<B><A>]`` be given, where ``<B>`` is a prefix
representing a bond with multiplicity :math:`\beta` and ``<A>`` is an atom
that can make :math:`\alpha` bonds maximally. For :math:`i \in \{1, \ldots, 7\}`, the
atomic symbol will map:

.. math::

    \begin{align}
        X_0 &\to \texttt{<A>} X_{\alpha} \\
        X_i &\to \texttt{<B'><A><X>}
    \end{align}

where ``<B'>`` is a prefix representing a bond with multiplicity
:math:`\mu = \min(\beta, \alpha, i)`,
and ``<X>`` is the empty string if :math:`\alpha - \mu = 0` or the
non-termial symbol :math:`X_{\alpha - \mu}` otherwise. Intuitively,
non-terminal states :math:`X_i` restricts subsequent bonds to a multiplicity
of at most :math:`i`. We provide an example of the derivation of the
SELFIES ``[F][=C][=13Cexpl][#N]``:

.. math::

    X_0 \to \texttt{F}X_1 \to \texttt{FC}X_3 \to \texttt{FC=[13C]}X_2 \to \texttt{FC=C=N}


**Branch Symbol:** Let Branch symbol ``[Branch<L>_<M>]`` be given. For
:math:`i \in \{0, 1, 9991, 9992, 9993\}`, the Branch symbol maps:

.. math::

    X_i \to X_i

in other words, the Branch symbol is ignored.

in progress : )

**Ring Symbol:**

in progress : )

Step 2: Ring Formation
**********************

in progress : )
