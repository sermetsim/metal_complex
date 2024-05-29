import pytest
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src/metal_complex")))
from metal_functions import *


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


@pytest.mark.parametrize(
    "ligands,expected",
    [
        (["[Br]", "[Br]", "[Br]","[Cl]"],['ligand3-0', 'ligand3-0', 'ligand3-0', 'ligand1-0']),
        (["[C]1(C([H])(C(C([H])([H])[C](C(C([C]1[H])([H])[H])([H])[H])[H])([H])[H])[H])[H]","[Br]","S([O])([O])([O])C(F)(F)F","",], ['ligand2-0', 'ligand3-0', 'ligand8-0']),
        (["Test for error output", "[Br]","[Cl]"],"Sorry your ligand is invalid")
        
        
    ]
)
def test_smile_to_number(ligands,expected): 
    results=smile_to_number(ligands)
    assert results== expected


@pytest.mark.parametrize(
    "number_list,expected",
    [
        (['ligand1-0', 'ligand3-0', 'ligand1-0'],-3),
        (['ligand8-0', 'ligand3-0', 'ligand9-0'],-2),
        (['ligand34-0', 'ligand3-0', 'ligand1-0', 'ligand34-0'],-2),
        ("Sorry your ligand is invalid","Sorry your ligand is invalid")
        
        
    ]
)

def total_charge_of_the_ligands(number_list,expected):
    results= total_charge_of_the_ligands(number_list)
    assert results == expected





@pytest.mark.parametrize(
    "charge, total_charge_ligands, metal,expected",
    [
        (-5,-6,"Fe",1),
        (-5,7,"Pt","Impossible oxydation state of your metal. Please check your inputs"),
        (2,-2,"Zr",4),
        (4,"Sorry your ligand is invalid","Co","Sorry your ligand is invalid")
        
    ]
)

def test_metal_oxydation_state(charge, total_charge_ligands, metal, expected):
    results= metal_oxydation_state(charge, total_charge_ligands, metal)
    assert results == expected
