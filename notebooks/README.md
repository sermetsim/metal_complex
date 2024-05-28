# ⚛️ Metal complex Notebook

This notebook explains the entire creation process of our program, from the objectives and problems faced to the actual code and potential improvements.

## ✅ Objectives

This project aims to develop a program that can:
- Design an interactive 3D structure of a complex from SMILES inputs.
- Determine the oxidation state of the central metal atom.
- Calculate the molecular weight of sophisticated complexes.

## 📚 Tools and Data Used

### Environment

- **Jupyter Lab**: Used for notebooks and code development.

### Packages

- **Pandas**: For data manipulation.
- **rdkit**: For molecular modeling, analysis, and design.
- **py3Dmol**: For embedding interactive 3D molecular viewers in notebooks.
- **Tkinter**: For creating and managing graphical user interfaces.
- **RegEx (re)**: For regular expression searches.

### Data

- **ligands_misc_info.csv**: Used to find the bonding atom and corresponding ligand number for each ligand.
- **ligands_fingerprints.csv**: Used to find the charge of each ligand.
- **oxidation_states_métaux.csv**: Used to verify the possible oxidation states of metals.

## 👩🏻‍💻 Creation Process

### 1. Model 3D Metal Complex

- **Initial Approach**: Attempted to use the 3Dmol Python library to draw complexes, but faced issues with complex and large structures.
- **Solution**: Developed a function to determine the coordinates of each atom in each ligand, grouping ligands and the metal in the same 2D plane before converting to 3D using RDKit.

### 2. Calculate Molecular Weight and Oxidation State

- **Molecular Weight**: Calculated using the `Chem.MolFromSmiles()` function from RDKit, summing individual SMILES molecular weights.
- **Oxidation State**: Determined using a database of ligand charges and requesting the total charge of the complex from the user.

### 3. Interface

- **Tkinter Interface**: Allows user input for ligands and metal (as SMILES), the global charge of the final complex, and choice of output (3D complex, molecular weight, oxidation state).

## 🔥 Usage

### Initialization

First, import the necessary functions:

```python
import sys
import os
notebook_path = os.getcwd()
src_path = os.path.abspath(os.path.join(notebook_path, "../src/metal_complex"))
sys.path.insert(0, src_path)
from metal_functions import *
```

### Model 3D Metal Complex

Example of creating a six-coordination number complex:

```python
list_of_ligand_SMILES = ['N([H])([H])C([H])([H])C([H])([H])N([H])[H]', 'N([H])([H])C([H])([H])C([H])([H])N([H])[H]', 'N([H])([H])C([H])([H])C([H])([H])N([H])[H]']
metal_SMILES = '[Fe]'
metal_complex(list_of_ligand_SMILES, metal_SMILES)
```

### Determine Oxidation State

Calculate oxidation state:

```python
ligands = ['N([H])([H])C([H])([H])C([H])([H])N([H])[H]', 'N([H])([H])C([H])([H])C([H])([H])N([H])[H]', 'N([H])([H])C([H])([H])C([H])([H])N([H])[H]']
metal_SMILES = '[Fe]'
charge = 3
print(metal_oxydation_state(charge, total_charge_of_the_ligands(smile_to_number(ligands)), metal_SMILES))
```

### Determine Molecular Mass

Calculate molecular mass:

```python
list_of_ligand_SMILES = ['[Cl]', '[Cl]', '[Cl]', '[Cl]', '[Cl]', '[Cl]']
metal_SMILES = '[Ni]'
print(calculate_MO(list_of_ligand_SMILES, metal_SMILES))
```

### Tkinter interface

Run the interface:

```python
%run ../src/metal_complex/interface_project.py
```


## 💪🏼 Faced Challenges

- **3D Modeling**: Initial attempts to convert SMILES directly to a 3D model faced issues. Resolved by combining ligand mol objects with metal atoms and creating 3D coordinates using RDKit.
- **Oxidation State**: Difficulty finding a suitable database, resolved by creating a custom .csv file.


## 🔮 Future Development

Future improvements could focus on:

- Determining the number of valence electrons of the central metal.
- Predicting possible reactions (e.g., ligand association/dissociation, oxidative addition).
- Enhancing model accuracy and expanding the ligand database.
- Integrating 3D visualization directly into the Tkinter interface or using a different interface for 3D visualization.

