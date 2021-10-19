from selfies.utils.smiles_utils import atom_to_smiles, smiles_to_atom


def modernize_symbol(symbol):
    """Converts a SELFIES symbol from <v2 to its latest equivalent.

    :param symbol: an old SELFIES symbol.
    :return: the latest equivalent of the input symbol, or the input symbol
        itself, if no such equivalent exists.
    """

    if symbol in _SYMBOL_UPDATE_TABLE:
        return _SYMBOL_UPDATE_TABLE[symbol]

    if symbol[-5:] == "expl]":  # e.g. [XXXexpl]
        if symbol[1] in "=#/\\":
            bond_char, atom_symbol = symbol[1], symbol[2:-5]
        else:
            bond_char, atom_symbol = "", symbol[1:-5]

        atom = smiles_to_atom("[{}]".format(atom_symbol))
        if (atom is not None) and (not atom.is_aromatic):
            atom_symbol = atom_to_smiles(atom, brackets=False)  # standardize
            symbol = "[{}{}]".format(bond_char, atom_symbol)

    return symbol


def _build_update_table():
    update_table = dict()
    for L in range(1, 4):
        entries = [
            ("[Branch{}_1]", "[Branch{}]"),
            ("[Branch{}_2]", "[=Branch{}]"),
            ("[Branch{}_3]", "[#Branch{}]"),
            ("[Expl=Ring{}]", "[=Ring{}]"),
            ("[Expl#Ring{}]", "[#Ring{}]"),
            ("[Expl/Ring{}]", "[//Ring{}]"),
            ("[Expl\\Ring{}]", "[\\\\Ring{}]")
        ]

        for old, new in entries:
            update_table[old.format(L)] = new.format(L)
    return update_table


_SYMBOL_UPDATE_TABLE = _build_update_table()
