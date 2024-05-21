#!/usr/bin/env python
# coding: utf-8

# In[4]:


import copy as cp
import pandas as pd
import re
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import py3Dmol
import tkinter as tk
from PIL import ImageTk

#################################################
## Load data ####################################

url1 = 'https://raw.githubusercontent.com/hkneiding/tmQMg-L/main/ligands_misc_info.csv'
data_smiles = pd.read_csv(url1, sep=";")
url2 = 'https://raw.githubusercontent.com/hkneiding/tmQMg-L/main/ligands_fingerprints.csv'
data_number_to_charge = pd.read_csv(url2, sep=";")
url3 = 'https://raw.githubusercontent.com/sermetsim/metal_complex/main/data/excel%20oxydation%20states%20m%C3%A9taux.csv'
data_oxydation_metal = pd.read_csv(url3, sep=";")


#################################################
## Functions ####################################

## simplifying functions ##

def simplify_smiles(smiles):
    ''' 
    parameters: smiles - a list of SMILES strings
    returns: list of simplified SMILES strings
    usage: simplifies a list of SMILES strings by removing unnecessary characters '''

    remove_chars = "[]()123456789=#-+\/:;.,!°{}"
    return [''.join([char for char in smi if char not in remove_chars]) for smi in smiles]


def simplify_idx(idx_list, smiles_list):
    ''' 
    parameters: idx_list - a list of index lists
                smiles_list - a list of SMILES string
    returns: ist of adjusted index lists
    usage: adjusts indices in idx_list to account for removed characters from the corresponding SMILES strings '''

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
    parameters: chaine - a string
    returns: list of integers found in the string
    usage: extracts all numbers from a string and returns them as a list of integers'''

    return [int(nombre) for nombre in re.findall(r'\d+', chaine)]


## metal-ligand bond formation ##

def get_idx(ligand_list):
    ''' 
    parameters: ligand_list - a list of a ligand SMILES strings
    returns: list of index lists corresponding to each ligand
    usage: etrieves bond node indices for given ligand SMILES from the ligands_misc_info.csv table '''

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
    usage: filters out empty strings from a list '''

    return [i for i in ligand_list if i != '']


def correct_link(link_list, ligand_list):
    ''' 
    parameters: link_list - a list of bond indices
                ligand_list - a list of RDKit molecule objects
    returns: flattened list of corrected indices (only one list)
    usage: adjusts the atom indices to account for the atoms of the previously added ligands '''

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
    usage: converts a list of SMILES strings to RDKit molecule objects '''

    return [Chem.MolFromSmiles(smi) for smi in ligand_list]


## complex visualisation ##

def create_molecule_in_3D(link, ligand, metal_smiles):
    ''' 
    parameters: link - a list of bond indices
                ligand - a list of RDKit molecule objects
                metal_smiles - a SMILES string of the metal
    returns: a molecule block (string) for 3D visualization
    usage: creates a 3D combined metal-ligand molecule '''

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
    returns: 3D visualization of the molecule
    usage: visualizes a molecule in 3D using py3Dmol '''

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
    usage: handles the entire process of creating and visualizing a metal-ligand complex (using the previous functions) '''

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
    usage: calculates the molecular weight of the metal complex '''

    Mo = 0
    for smi in list_of_ligand:
        Mo += Descriptors.MolWt(Chem.MolFromSmiles(smi))
    Mo += Descriptors.MolWt(Chem.MolFromSmiles(metal_name))
    return Mo



## oxidation state calculation ##

def smile_to_number(ligands):
    ''' 
    parameters: ligands - a list of ligand SMILES strings
    returns: list of numbers corresponding to each ligand
    usage: converts ligand SMILES to numbers based on the ligands_misc_info.csv table  '''

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
    parameters: number_list - a list of ligand numbers
    returns: total charge of the ligands
    usage:calculates the total charge of the ligands based on their numbers determined in smile_to_number(ligands) '''
    
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
    usage: determines the oxidation state of the metal based on ligand charges and total complex charge '''
    
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
                return "Impossible oxidation state of your metal. Please check your inputs"
    return "Metal not found in the data. Please check your inputs"



#################################################
## Interface ####################################

## main handler ##

def handle_input():
    '''
    collects user inputs: entry_lig1, entry_lig2, ... - SMILES strings for up to six ligands 
                          entry_metal - SMILES string for the metal 
                          entry_charge - total charge of the complex converted into a float
                          
    generate and show the Metallic Complex (if prop1 is selected): calls metal_complex(all_ligands_entries, metal_entry), which generates and visualizes the metal-ligand complex in 3D
    
    calculate Molecular Mass (if prop4 is selected): calls calculate_MO(all_ligands_entries, metal_entry) to calculate the molecular mass of the complex
                                                     updates the label with the calculated molecular mass and displays it
    
    calculate Oxidation State (if prop5 is selected): converts ligand SMILES to their corresponding numbers using smile_to_number(all_ligands_entries)
                                                      calculates the total charge of the ligands using total_charge_of_the_ligands(number_list)
                                                      determines the oxidation state of the metal using metal_oxydation_state(total_charge, ligands_charge, metal_entryligand_charge,entry_charge,entry_metal)
                                                      updates the label with the calculated oxidation state and displays it
    '''

    all_ligands_entries = [entry_lig1.get(), entry_lig2.get(), entry_lig3.get(), entry_lig4.get(), entry_lig5.get(), entry_lig6.get()]
    metal_entry = entry_metal.get()
    total_charge = float(entry_charge.get())

    if prop1:
        metal_complex(all_ligands_entries, metal_entry)
    
    if prop4:
        MO_label.pack_forget()
        molar_mass = calculate_MO(all_ligands_entries, metal_entry)
        MO_label.config(text="Molecular Mass: " + str(round(molar_mass, 3)))
        MO_label.pack()
    else:
        MO_label.pack_forget()
    
    if prop5:
        ox_state_label.pack_forget()
        ox_state_label.pack()
        number_list = smile_to_number(all_ligands_entries)
        ligands_charge = total_charge_of_the_ligands(number_list)
        ox_state = metal_oxydation_state(total_charge, ligands_charge, metal_entry)
        ox_state_label.config(text="Oxidation State: " + str(ox_state))
    else:
        ox_state_label.pack_forget()

'''
The toggle functions update their corresponding global variables (prop1, prop4, prop5) wether their checkboxes are checked or not
'''
def toggle_prop1():
    global prop1
    prop1 = int(var_prop1.get())


def toggle_prop4():
    global prop4
    prop4 = int(var_prop4.get())


def toggle_prop5():
    global prop5
    prop5 = int(var_prop5.get())

'''
First all variables are considered unchecked
'''
prop1, prop4, prop5 = 0, 0, 0  


## formation of Tkinter window ##

'''
Main Window: The Tkinter window is created and titled "Metallic Complex"

Labels and Entry Fields: - a label instructs the user to insert up to six ligand SMILES strings
                         - six entry fields are provided for the ligand SMILES strings
                         - a label and entry field are provided for the metal SMILES string
                         - a label and entry field are provided for the total charge of the complex
                         
Checkboxes: - select the option to show the metallic complex
            - select the option to calculate the molecular mass
            - select the option to calculate the oxidation state
            
Submit Button: when clicked it calls handle_input() to process the inputs and display the results

Result Labels: display the calculated molecular mass and oxidation state
'''

window = tk.Tk()
window.title("Metallic Complex")

label_ligand = tk.Label(window, text="Insert all the ligands SMILES (6 max.):")
label_ligand.pack()

frame_ligands = tk.Frame(window)
frame_ligands.pack()

entry_lig1 = tk.Entry(frame_ligands, width=8)
entry_lig1.grid(row=0, column=0, padx=5, pady=5)
entry_lig2 = tk.Entry(frame_ligands, width=8)
entry_lig2.grid(row=0, column=1, padx=5, pady=5)
entry_lig3 = tk.Entry(frame_ligands, width=8)
entry_lig3.grid(row=0, column=2, padx=5, pady=5)
entry_lig4 = tk.Entry(frame_ligands, width=8)
entry_lig4.grid(row=0, column=3, padx=5, pady=5)
entry_lig5 = tk.Entry(frame_ligands, width=8)
entry_lig5.grid(row=0, column=4, padx=5, pady=5)
entry_lig6 = tk.Entry(frame_ligands, width=8)
entry_lig6.grid(row=0, column=5, padx=5, pady=5)

label_metal = tk.Label(window, text="Insert the metal SMILES:")
label_metal.pack()

entry_metal = tk.Entry(window, width=5)
entry_metal.pack()

label_charge = tk.Label(window, text="Insert the total charge of the complex:")
label_charge.pack()

entry_charge = tk.Entry(window, width=5)
entry_charge.pack()

var_prop1 = tk.IntVar()
check_prop1 = tk.Checkbutton(window, text="Show the metallic complex", variable=var_prop1, command=toggle_prop1)
check_prop1.pack()

var_prop4 = tk.IntVar()
check_prop4 = tk.Checkbutton(window, text="Molecular Mass", variable=var_prop4, command=toggle_prop4)
check_prop4.pack()

var_prop5 = tk.IntVar()
check_prop5 = tk.Checkbutton(window, text="Oxidation State", variable=var_prop5, command=toggle_prop5)
check_prop5.pack()

button = tk.Button(window, text="Submit", command=handle_input)
button.pack()

MO_label = tk.Label(window, text="")
ox_state_label = tk.Label(window, text="")

window.mainloop()


# In[ ]:




