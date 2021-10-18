ELEMENTS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "Hf",
    "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Po", "At", "Rn", "Fr", "Ra", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Fl", "Lv", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
    "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
    "No", "Lr"
}

ORGANIC_SUBSET = {"B", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"}

AROMATIC_VALENCES = {
    "B": (3,), "Al": (3,),
    "C": (4,), "Si": (4,),
    "N": (3, 5), "P": (3, 5), "As": (3, 5),
    "O": (2, 4), "S": (2, 4), "Se": (2, 4), "Te": (2, 4)
}

AROMATIC_SUBSET = set(e.lower() for e in AROMATIC_VALENCES)

# =============================================================================
# SELFIES-specific constants
# =============================================================================


INDEX_ALPHABET = (
    "[C]", "[Ring1]", "[Ring2]",
    "[Branch1]", "[=Branch1]", "[#Branch1]",
    "[Branch2]", "[=Branch2]", "[#Branch2]",
    "[O]", "[N]", "[=N]", "[=C]", "[#C]", "[S]", "[P]"
)

INDEX_CODE = {c: i for i, c in enumerate(INDEX_ALPHABET)}
