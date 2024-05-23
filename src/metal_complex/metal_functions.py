import copy as cp
import pandas as pd
import re
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import py3Dmol

## simplifying functions ##
def simplify_smiles(smiles):
    ''' 
    parameters: smiles - a list of SMILES strings
    returns: list of simplified SMILES strings
    usage: simplifies a list of SMILES strings by removing unnecessary characters to use these new smiles in the next function 

    Example
    ----
    >>> simplify_smiles(['[C]([H])([H])[H]'])
    ['CHHH']
    '''

    remove_chars = "[]()123456789=#-+\/:;.,!°{}"
    return [''.join([char for char in smi if char not in remove_chars]) for smi in smiles]


def simplify_idx(idx_list, smiles_list):
    ''' 
    parameters: idx_list - a list of index lists (atoms in the ligand that will link the metal)
                smiles_list - a list of SMILES string (from simplify_smiles)
    returns: list of adjusted index lists
    usage: adjusts indices in idx_list to account for removed characters (H atoms) from the corresponding SMILES strings
    
    Example
    ----
    >>> simplify_idx([[0,11],[0]],['CHHCHHHCHHHCHH','CHHH'])
    [[0,4],[0]]
    '''

    copy_list = cp.deepcopy(idx_list)
    new_list = cp.deepcopy(idx_list)
    for k in range(len(smiles_list)):
        for i in range(len(smiles_list[k])):
            if smiles_list[k][i] == 'H':
                for j in range(len(idx_list[k])):
                    if copy_list[k][j] > i:
                        new_list[k][j] -= 1
    return new_list


def extraire_nombres(chaine):
    ''' 
    parameters: chaine - a string that contain a list into a list (e.g. '[[1]]')
    returns: list of integers found in the string
    usage: extracts all list of numbers from a string and returns them as a list of integers
        
    Example
    ----
    >>> extraire_nombres([[0,1]])
    [0,1]
    '''

    return [int(nombre) for nombre in re.findall(r'\d+', chaine)]


def get_idx(ligand_list):
    ''' 
    parameters: ligand_list - a list of a ligand SMILES strings
    returns: list of index lists corresponding to each ligand
    usage:found the indices of bonding atoms for given ligand SMILES from the ligands_misc_info.csv table 
    
    Example
    ----
    >>> get_idx['[Cl]','[Cl]'])
    [[0],[0]]
    '''

    url1 = 'https://raw.githubusercontent.com/hkneiding/tmQMg-L/main/ligands_misc_info.csv'
    data_smiles = pd.read_csv(url1, sep=";")
    idx_list = []
    for j in ligand_list:
        for i in range(len(data_smiles["smiles"])):
            if j == data_smiles["smiles"][i]:
                idx = data_smiles["smiles_metal_bond_node_idx_groups"][i]
                idx_list.append(extraire_nombres(idx))
                break
    return idx_list


def filter_my_list(ligand_list):
    ''' 
    parameters: ligand_list - a list of ligands (as SMILES)
    returns: filtered list with non-empty ligands
    usage: filters out empty strings from a list 
        
    Example
    ----
    >>> filter_my_list(['CHHCHHHCHHHCHH','CHHH', '', '[Cl]'])
    ['CHHCHHHCHHHCHH','CHHH', '[Cl]']
    '''

    return [i for i in ligand_list if i != '']


def correct_link(link_list, ligand_list):
    ''' 
    parameters: link_list - a list of bond indices
                ligand_list - a list of RDKit molecule objects
    returns: 1D list of corrected indices
    usage: adjusts the atom indices to account for the atoms of the previously added ligands 
            
    Example
    ----
    >>> correct_link([[0,4],[0]], ['C([H])([H])C([H])([H])([H])C([H])([H])([H])C([H])([H])[H]','[C]([H])([H])[H]'])
    [0,4,5]
    '''

    new_link_list = cp.deepcopy(link_list)
    for i in range(len(link_list)):
        if i > 0:
            for j in range(len(link_list[i])):
                x = i - 1
                while x >= 0:
                    new_link_list[i][j] += ligand_list[x].GetNumAtoms()
                    x -= 1
    return [j for i in new_link_list for j in i]


def smiles_to_ligand(ligand_list):
    ''' 
    parameters: ligand_list - a list of SMILES strings
    returns: list of RDKit molecule objects
    usage: converts a list of SMILES strings to RDKit molecule objects 
                
    Example
    ----
    >>> correct_link(['C([H])([H])C([H])([H])([H])C([H])([H])([H])C([H])([H])[H]','[C]([H])([H])[H]'])
    [mol object of ethane, mol object of methane]
    '''

    return [Chem.MolFromSmiles(smi) for smi in ligand_list]


## complex visualisation ##

def create_molecule_in_3D(link, ligand, metal_smiles):
    ''' 
    parameters: link - a list of bond indices
                ligand - a list of RDKit molecule objects
                metal_smiles - a SMILES string of the metal
    returns: a molecule block (string) for 3D visualization
    usage: creates a 3D combined metal-ligand molecule 
    '''

    combined_molecule = Chem.MolFromSmiles(metal_smiles)
    for lig in ligand[::-1]:
        combined_molecule = Chem.CombineMols(lig, combined_molecule)

    mol = Chem.RWMol(combined_molecule)
    for atom in mol.GetAtoms():
        if f'[{atom.GetSymbol()}]' == metal_smiles:
            idx_metal = atom.GetIdx()

    for i in link:
        mol.AddBond(i, idx_metal, Chem.BondType.SINGLE)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    opt_block = Chem.MolToMolBlock(mol)
    return opt_block


def show_complex(opt_block):
    ''' 
    parameters: opt_block - a molecule block (string)
    returns: interactive 3D visualization of the molecule
    usage: visualizes a molecule in 3D using py3Dmol 
    '''

    viewer = py3Dmol.view(width=400, height=300)
    viewer.addModel(opt_block, "mol")
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    return viewer.show()


def metal_complex(list_of_ligand, metal_name):
    ''' 
    parameters: list_of_ligand - a list of ligand SMILES strings
                metal_name - a SMILES string of the metal
    returns: 3D visualization of the metal-ligand complex
    usage: handles the entire process of creating and visualizing a metal-ligand complex (using the previous functions) 
    '''

    ligand_list_not_filtered = cp.deepcopy(list_of_ligand)
    ligand_list_filtered = filter_my_list(ligand_list_not_filtered)
    not_correct_link = get_idx(ligand_list_filtered)
    easy_smiles = simplify_smiles(ligand_list_filtered)
    easy_idx = simplify_idx(not_correct_link, easy_smiles)
    ligand_mol_list = smiles_to_ligand(ligand_list_filtered)
    link_atoms = correct_link(easy_idx, ligand_mol_list)
    complex_block = create_molecule_in_3D(link_atoms, ligand_mol_list, metal_name)
    return show_complex(complex_block)


## molecular weight calculation ##

def calculate_MO(list_of_ligand, metal_name):
    ''' 
    parameters: list_of_ligand - a list of ligand SMILES strings
                metal_name - a SMILES string of the metal
    returns: molecular weight of the complex
    usage: calculates the molecular weight of the metal complex 
                    
    Example
    ----
    >>> calculate_MO(['O','O'],'[Fe]')
    91.875
    '''

    Mo = 0
    for smi in list_of_ligand:
        Mo += Descriptors.MolWt(Chem.MolFromSmiles(smi))
    Mo += Descriptors.MolWt(Chem.MolFromSmiles(metal_name))
    return Mo


## oxidation state calculation ##

def smile_to_number(ligands):
    ''' 
    parameters: ligands - a list of ligand SMILES strings
    returns: number_list - a list of strings. Each string correspond to the ligand number of each ligand
    usage: converts ligand SMILES to ligand numbers based on the ligands_misc_info.csv table  
                        
    Example
    ----
    >>> smile_to_number()
    
    '''

    url1 = 'https://raw.githubusercontent.com/hkneiding/tmQMg-L/main/ligands_misc_info.csv'
    data_smiles = pd.read_csv(url1, sep=";")
    number_list = []
    for ligand in ligands:
        if ligand == "":
            continue
        found = False
        for idx, row in data_smiles.iterrows():
            if row["smiles"] == ligand:
                number_list.append(row["name"])
                found = True
                break
        if not found:
            return "Sorry your ligand is invalid"
    return number_list


def total_charge_of_the_ligands(number_list):
    ''' 
    parameters: number_list - a list of the ligand numbers in strings
    returns: total charge of the ligands
    usage:calculates the total charge of the ligands based on their numbers determined in smile_to_number(ligands) 
                            
    Example
    ----
    >>> total_charge_of_the_ligands()
    
    '''

    for i in number_list:
        if i == "Sorry your ligand is invalid":
            return "Sorry your ligand is invalid"
    
    url2 = 'https://raw.githubusercontent.com/hkneiding/tmQMg-L/main/ligands_fingerprints.csv'
    data_number_to_charge = pd.read_csv(url2, sep=";")
    total_charge_ligands = 0
    for ligand_number in number_list:
        for index, row in data_number_to_charge.iterrows():
            if row['name'] == ligand_number:
                total_charge_ligands += row['charge']
                break
    return total_charge_ligands


def metal_oxydation_state(charge, total_charge_ligands, metal):
    ''' 
    parameters: charge - the total charge of the complex
                total_charge_ligands - the total charge of the ligands
                metal - the SMILES string of the metal
    returns: oxidation state of a metal
    usage: determines the oxidation state of the metal based on ligand charges and total complex charge 
                            
    Example
    ----
    >>> metal_oxydation_state()
    
    '''

    if total_charge_ligands == "Sorry your ligand is invalid":
        return "Sorry your ligand is invalid"
    
    url3 = 'https://raw.githubusercontent.com/sermetsim/metal_complex/main/data/oxydation%20states%20m%C3%A9taux.csv'
    data_oxydation_metal = pd.read_csv(url3, sep=";")
    oxydation_by_input = charge - total_charge_ligands
    met=[metal]
    for index, row in data_oxydation_metal.iterrows():
        if row['syllabus'] == simplify_smiles(met)[0]:
            states_str = row['oxydation states'].strip('[]')
            possibles_oxydation = list(map(int, states_str.split(',')))
            if oxydation_by_input in possibles_oxydation:
                return oxydation_by_input
            else:
                return "Impossible oxydation state of your metal. Please check your inputs"
    return "Metal not found in the data. Please check your inputs"
