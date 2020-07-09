Understanding SELFIES
=====================

In this section, we provide an informal tutorial on interpreting SELFIES.

Character Syntax
########################

In general, all SELFIES characters are enclosed in square brackets, with
exception of the dot-bond character. SELFIES characters share a close
connection with SMILES characters, which speeds up translation between SELFIES
and SMILES. There are four main types of SELFIES characters in the alphabet:


1.  The dot-bond symbol ``'.'`` is used to represent a SELFIES with
    disconnected parts. Similar to how it is used in SMILES, the fragments
    separated by the dot-bond symbol are treated as separate SELFIES strings.

2.  Ring symbols are of the general form ``'[Ring<N>]'`` or
    ``'[Expl<B>Ring<N>]'``, where ``N in {1, 2, 3}`` and
    ``B in {'/', '\\', '=', '#'}`` is a bond character. Ring symbols are
    used to specify ring closure bonds, analagous to the ring numbering
    in SMILES. Examples: ``'[Ring1]'`` or ``'[Expl=Ring2]'``.

3.  Branch symbols are of the general form ``'[Branch<N>_<M>]'``, where
    ``N, M in {1, 2, 3}``. Branch symbols are used to specify branches,
    analagous to the open and closed curved bracket symbols in SMILES.
    Examples: ``[Branch1_2]`` or ``[Branch2_2]``.

4.  Atom symbols are closely related to those in SMILES. In general, they
    are of the form ``'[<B><A>]'`` or ``'[<B><A>expl]'``, where
    ``B in {'/', '\\', '=', '#'}`` is a bond character and ``A`` is
