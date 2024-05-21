import copy as cp
import pandas as pd
import re
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import py3Dmol

# Helper functions
def simplify_smiles(smiles):
    remove_chars = "[]()123456789=#-+\/:;.,!°{}"
    return [''.join([char for char in smi if char not in remove_chars]) for smi in smiles]


def simplify_idx(idx_list, smiles_list):
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
    return [int(nombre) for nombre in re.findall(r'\d+', chaine)]


def get_idx(ligand_list):
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
    return [i for i in ligand_list if i != '']


def correct_link(link_list, ligand_list):
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
    return [Chem.MolFromSmiles(smi) for smi in ligand_list]


def create_molecule_in_3D(link, ligand, metal_smiles):
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
    viewer = py3Dmol.view(width=400, height=300)
    viewer.addModel(opt_block, "mol")
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    return viewer.show()


def metal_complex(list_of_ligand, metal_name):
    ligand_list_not_filtered = cp.deepcopy(list_of_ligand)
    ligand_list_filtered = filter_my_list(ligand_list_not_filtered)
    not_correct_link = get_idx(ligand_list_filtered)
    easy_smiles = simplify_smiles(ligand_list_filtered)
    easy_idx = simplify_idx(not_correct_link, easy_smiles)
    ligand_mol_list = smiles_to_ligand(ligand_list_filtered)
    link_atoms = correct_link(easy_idx, ligand_mol_list)
    complex_block = create_molecule_in_3D(link_atoms, ligand_mol_list, metal_name)
    return show_complex(complex_block)


def calculate_MO(list_of_ligand, metal_name):
    Mo = 0
    for smi in list_of_ligand:
        Mo += Descriptors.MolWt(Chem.MolFromSmiles(smi))
    Mo += Descriptors.MolWt(Chem.MolFromSmiles(metal_name))
    return Mo


def smile_to_number(ligands):
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
