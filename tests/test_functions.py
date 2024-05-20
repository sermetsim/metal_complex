import pytest
# Assuming functions are imported from the main script
from main_functions import *


@pytest.mark.parametrize(
    'values, expected', [
        (["C1=CC=CC=C1", "CCO", "[Cu](CC)(CC)"], ["CCCCCC", "CCO", "CuCCCC"]),
        (['c1ccccc1'], ['cccccc']),
        (['[Fe]-2'], ['Fe']),
    ]
)
def test_simplify_smiles(values, expected):
    results = simplify_smiles(values)
    assert results == expected


@pytest.mark.parametrize(
    'idx, smiles, expected', [
        ([[0,1],[0]],['CCHHCHHC','CH'],[[0,1],[0]] ),
        ([[0,4],[5]],['CHHHCHHH','CCCCHC'],[[0,1],[4]]),
    ]
)
def test_simplify_idx(idx, smiles, expected):
    results = simplify_idx(idx, smiles)
    assert results == expected


@pytest.mark.parametrize(
    'chaine, expected', [
        ('[[0,1]]',[0,1]),
        ('[[0]]',[0]),
        ('[[a]]',[]),
        ('[[a,0]]',[0]),
    ]
)
def test_extraire_nombres(chaine, expected):
    results = extraire_nombres(chaine)
    assert results == expected


@pytest.mark.parametrize(
    'ligand_list, expected', [
        (['[Cl]'],['[Cl]']),
        (['[Cl]','','','[Cl]',''],['[Cl]','[Cl]']),
        (['[H]O[H]','[H]O[H]','[Cl]',''],['[H]O[H]','[H]O[H]','[Cl]']),
    ]
)
def test_filter_my_list(ligand_list, expected):
    results = filter_my_list(ligand_list)
    assert results == expected


@pytest.mark.parametrize(
    'link_list, ligand_list, expected', [
        ([[0]],[Chem.MolFromSmiles('[Cl]')],[0]),
        ([[0],[0]],[Chem.MolFromSmiles('[Cl]'),Chem.MolFromSmiles('[Cl]')],[0,1]),
        ([[0],[0],[0]],[Chem.MolFromSmiles('[H]O[H]'), Chem.MolFromSmiles('[H]O[H]'), Chem.MolFromSmiles('[Cl]')],[0,1,2]),
    ]
)
def test_correct_link(link_list, ligand_list, expected):
    results = correct_link(link_list, ligand_list)
    assert results == expected


@pytest.mark.parametrize(
    'list_of_ligand, metal_name , expected', [
        (['[H]O[H]','[H]O[H]'], '[Fe]', 91.875),
        (['[H]O[H]','[Cl]'],'[Zr]',144.692),
        (['c1ccccc1','c1ccccc1','c1ccccc1','c1ccccc1'], '[Pt]', 507.534),
    ]
)
def test_calculate_MO(list_of_ligand, metal_name, expected):
    results = calculate_MO(list_of_ligand, metal_name)
    assert results == expected