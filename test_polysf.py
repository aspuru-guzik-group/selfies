import selfies as sf
from rdkit import Chem

polymer_smiles = ['*CC(*)(C)C',
                  'C1=C(SC(=C1)[*])[*]',
                  'CCCCC1=C(SC(=C1)[*])[*]',
                  'CCCCCCC1=C(SC(=C1)[*])[*]',
                  'CCCCCCCCC1=C(SC(=C1)[*])[*]',
                  'C1(=CC(=C(C=C1C=C[*])OC)[*])OCC(CC)CCCC'
                  ]

for i in polymer_smiles:
    mol = Chem.MolFromSmiles(i)
    ori_smi = Chem.MolToSmiles(mol)
    selfies = sf.encoder(ori_smi)
    de_smi = sf.decoder(selfies)
    de_smi = Chem.MolToSmiles(Chem.MolFromSmiles(de_smi))
    print('polymer smiles:', ori_smi, 'selfies:', selfies, 'decode selfies:', de_smi, 'equal?:', ori_smi == de_smi)



