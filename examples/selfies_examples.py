"""This file contains examples of how to use the ``selfies`` library.
"""

from rdkit import Chem

import selfies as sf

# 1. First we try translating SMILES --> SELFIES --> SMILES.

# Test SMILES: non-fullerene acceptor for organic solar cells.
smiles = "CN1C(=O)C2=C(c3cc4c(s3)-c3sc(-c5ncc(C#N)s5)cc3C43OCCO3)N(C)C(=O)" \
         "C2=C1c1cc2c(s1)-c1sc(-c3ncc(C#N)s3)cc1C21OCCO1"
encoded_selfies = sf.encoder(smiles)  # SMILES --> SEFLIES
decoded_smiles = sf.decoder(encoded_selfies)  # SELFIES --> SMILES

print(f"Original SMILES: {smiles}")
print(f"Translated SELFIES: {encoded_selfies}")
print(f"Translated SMILES: {decoded_smiles}")
print()

# When comparing the original and decoded SMILES, do not use == equality. Use
# RDKit to check both SMILES correspond to the same molecule.
print(f"== Equals: {smiles == decoded_smiles}")

can_smiles = Chem.CanonSmiles(smiles)
can_decoded_smiles = Chem.CanonSmiles(decoded_smiles)
print(f"RDKit Equals: {can_smiles == can_decoded_smiles}")
print()

# Advanced Usage: Customizing SELFIES

# 2. Let's view the default SELFIES semantic constraints.

default_constraints = sf.get_semantic_constraints()
print(f"Default Constraints:\n {default_constraints}")
print()

# 3. Let's modify the SELFIES constraints

# We have two compounds here, CS=CC#S and [Li]=CC in SELFIES form
c_s_compound = sf.encoder("CS=CC#S")
li_compound = sf.encoder("[Li]=CC")

# Under the default SELFIES settings, they are translated as
print("Default SELFIES:")
print(f"\t CS=CC#S --> {sf.decoder(c_s_compound)}")
print(f"\t [Li]=CC --> {sf.decoder(li_compound)}")
print()
# Since Li is not recognized by SELFIES, it is constrained to 8
# bonds by default.

# Now we add Li to the SELFIES constraints, and restrict it to 1 bond only
# Also, let's restrict S to 2 bonds (instead of its default 6).
new_constraints = default_constraints
new_constraints['Li'] = 1
new_constraints['S'] = 2

sf.set_semantic_constraints(new_constraints)  # update constraints

# Let's check on the SELFIES constraints to see if they are updated.
print(f"Updated Constraints:\n {sf.get_semantic_constraints()}")
print()

# Under our new settings, they are translated as
print("Updated SELFIES:")
print(f"\t CS=CC#S --> {sf.decoder(c_s_compound)}")
print(f"\t [Li]=CC --> {sf.decoder(li_compound)}")
# Notice that all the bond constraints are met.

# 4. Let's revert the SELFIES constraints back to the default settings
sf.set_semantic_constraints()
