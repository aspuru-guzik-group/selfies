import pytest
import selfies as sf


class EndToEndTestCase:
    def __init__(self, identifier: str, smiles: str, selfies: str = None):
        self.identifier = identifier
        self.smiles = smiles
        self.selfies = selfies


test_cases = [
    # Please use canonical SMILES to create new test cases!
    EndToEndTestCase(
        identifier='Pentylamine',
        smiles='CCCCCN',
        selfies='[C][C][C][C][C][N]'
    ),
    EndToEndTestCase(
        identifier='Phenylalanine',
        smiles='N[C@@H](Cc1ccccc1)C(=O)O',
    ),
    EndToEndTestCase(
        identifier='MDMA',
        smiles='CNC(C)Cc1ccc2c(c1)OCO2',
    ),
    EndToEndTestCase(
        identifier='DEET',
        smiles='CCN(CC)C(=O)c1cccc(C)c1',
    ),
    EndToEndTestCase(
        identifier='Paracetamol',
        smiles='CC(=O)Nc1ccc(O)cc1',
    ),
    EndToEndTestCase(
        identifier='Ibuprofen',
        smiles='CC(C)Cc1ccc(C(C)C(=O)O)cc1',
    ),
    EndToEndTestCase(
        identifier='Afatinib',
        smiles='CN(C)CC=CC(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1OC1CCOC1',
    ),
    EndToEndTestCase(
        identifier='Sorafenib',
        smiles='CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1',
    ),
    EndToEndTestCase(
        identifier='Azulene',
        smiles='c1ccc2cccc-2cc1',
    ),
    EndToEndTestCase(
        identifier='Chrysene',
        smiles='c1ccc2c(c1)ccc1c3ccccc3ccc21',
    ),
    EndToEndTestCase(
        identifier='WOJWZDJWEWKCIH-SNCNVHSMSA-N',
        smiles='Cc1c(C)c(S(=O)(=O)NC(=N)NCCC[C@H](NC(=O)[C@@H]2CCCN2C(=O)[C@H](CCC(=O)NC(c2ccccc2)(c2ccccc2)c2ccccc2)NC'
               '(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCNC(=O)OC(C)(C)C)NC(=O)[C@H](C)NC(=O)[C@@H]2CCCN2C(=O)[C@@H]2CCCN2C(=O)'
               '[C@H](CCCCNC(=O)OC(C)(C)C)NC(=O)[C@H](CCCCNC(=O)OC(C)(C)C)NC(=O)[C@H](COC(C)(C)C)NC(=O)[C@H](CCC(=O)OC('
               'C)(C)C)NC(=O)[C@H](CCCCNC(=O)OC(C)(C)C)NC(=O)[C@H](CCCNC(=N)NS(=O)(=O)c2c(C)c(C)c3c(c2C)CCC(C)(C)O3)NC('
               '=O)[C@H](CCC(=O)NC(c2ccccc2)(c2ccccc2)c2ccccc2)NC(=O)[C@H](CCC(=O)NC(c2ccccc2)(c2ccccc2)c2ccccc2)NC(=O)'
               '[C@@H](NC(=O)[C@H](CCCNC(=N)NS(=O)(=O)c2c(C)c(C)c3c(c2C)CCC(C)(C)O3)NC(=O)[C@H](CCC(=O)NC(c2ccccc2)(c2c'
               'cccc2)c2ccccc2)NC(=O)[C@H](Cc2cn(C(=O)OC(C)(C)C)cn2)NC(=O)[C@H](CCC(=O)OC(C)(C)C)NC(=O)[C@@H]2CCCN2C(=O'
               ')[C@H](COC(C)(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](COC(c2ccccc2)(c2ccccc2)c2ccccc2)'
               'NC(=O)[C@H](COC(C)(C)C)NC(=O)CNC(=O)OC(C)(C)C)C(C)C)C(=O)O)c(C)c2c1OC(C)(C)CC2',
    ),
]


@pytest.mark.parametrize("case", test_cases, ids=[tc.identifier for tc in test_cases])
def test_smiles_to_selfies(case: EndToEndTestCase):
    if not case.selfies:
        pytest.skip('Need to add an expected SELFIES result for this molecule!')
    assert sf.encoder(case.smiles) == case.selfies


@pytest.mark.parametrize("case", test_cases, ids=[tc.identifier for tc in test_cases])
def test_roundtrip_smiles_to_selfies(case: EndToEndTestCase):
    output_selfies = sf.encoder(case.smiles)
    roundtrip_smiles = sf.decoder(output_selfies)
    assert case.smiles == roundtrip_smiles
