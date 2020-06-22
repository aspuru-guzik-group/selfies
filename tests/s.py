from rdkit.Chem import MolFromSmiles
import selfies

print(MolFromSmiles("O=C(O)c1c(-c2cnc(-c3nc4ccccc4s3)nn2)csc1-c1cnc(-c2nc3ccccc3s2)nn1"))
print(selfies.encoder("O=C(O)c1c(-c2cnc(-c3nc4ccccc4s3)nn2)csc1-c1cnc(-c2nc3ccccc3s2)nn1"))
