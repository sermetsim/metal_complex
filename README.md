# ⚛️ Metal-Ligand Complex Generator

This Python script generates and visualizes metal-ligand complexes, calculates molecular properties, and determines the oxidation state of the metal in the complex. It uses the RDKit library for cheminformatics, py3Dmol for 3D visualization, and Tkinter for the graphical user interface.


## 🔗 Features

1. ***SMILES Simplification***: Simplifies SMILES strings by removing unnecessary characters.
2. ***Ligand and Index Extraction***: Extracts bond node indices from ligand SMILES.
3. ***Ligand and Link Correction***: Corrects atom indices based on the number of atoms in each ligand.
4. ***Molecule Creation***: Creates a combined metal-ligand molecule in 3D.
5. ***Visualization***: Visualizes the metal-ligand complex in 3D.
6. ***Molecular Weight Calculation***: Calculates the molecular weight of the metal complex.
7. ***Oxidation State Calculation***: Determines the oxidation state of the metal based on ligand charges and the total charge of the complex.
8. ***Graphical User Interface***: Provides a GUI for user input and displays results.

## 🛠️ Requirements

- [Python 3.x](https://www.python.org/): The programming language used to write and run the script.
- [RDKit](https://github.com/rdkit/rdkit): A collection of cheminformatics and machine learning tools used for molecular operations such as creating and manipulating molecules from SMILES strings.
- [Py3Dmol](https://github.com/3dmol/3Dmol.js): A 3D molecular visualization library used to render the metal-ligand complex.
- [Tkinter](https://wiki.python.org/moin/TkInter): Python's standard GUI package used to create the graphical user interface for input and output.
- [Pandas](https://github.com/pandas-dev/pandas): A data manipulation and analysis library used to handle and process CSV data files.



## ⚙️ Installation

1. Clone this repository:
```
git clone https://github.com/sermetsim/metal_complex
cd metal_complex
```
2. Install `RDKit`,`Py3Dmol`,`Tkinter` (usually included with Python) and `pandas`:
   ```
   pip install .
   ```


## 🔥 Usage

1. Run the script:
```
python interface_project.py
```

2. The Tkinter GUI will appear. Follow these steps:

   - Insert up to 6 ligand SMILES strings from **[ligands_misc_info.csv](https://raw.githubusercontent.com/hkneiding/tmQMg-L/main/ligands_misc_info.csv)** in the provided fields.
   - Insert the metal SMILES string.
   - Insert the total charge of the complex.
   - Select the properties you want to display:
      - Show the metallic complex.
      - Display the molecular mass.
      - Display the oxidation state.
   - Click the "Submit" button to generate the results.

Instead of run the script, you also can import the wanted functions in your own script with:  
```
from metal_functions import {the functions}
```


## 🧮 Functions

### Simplification and Extraction

- ***simplify_smiles(smiles)***: Simplifies SMILES strings by removing unnecessary characters.
- ***simplify_idx(idx_list, smiles_list)***: Simplifies indices based on the simplified SMILES.
- ***extraire_nombres(chaine)***: Extracts numbers from a string.
- ***get_idx(ligand_list)***: Gets bond node indices from ligand SMILES.
- ***filter_my_list(ligand_list)***: Removes empty entries from the ligand list.
- ***correct_link(link_list, ligand_list)***: Corrects the atom indices based on the number of atoms in each ligand.

### Molecule Creation

- ***smiles_to_ligand(ligand_list)***: Converts SMILES strings to RDKit molecule objects.
- ***create_molecule_in_3D(link, ligand, metal_smiles)***: Creates a combined metal-ligand molecule in 3D.
- ***show_complex(opt_block)***: Visualizes the molecule in 3D using Py3Dmol.

### Calculations

- ***calculate_MO(list_of_ligand, metal_name)***: Calculates the molecular weight of the metal complex.
- ***smile_to_number(ligands)***: Converts ligand SMILES to numbers.
- ***total_charge_of_the_ligands(number_list)***: Calculates the total charge of the ligands.
- ***metal_oxydation_state(charge, total_charge_ligands, metal)***: Determines the oxidation state of the metal.

### GUI Handlers

- ***handle_input()***: Handles user input from the GUI.
- ***toggle_prop1()***, ***toggle_prop4()***, ***toggle_prop5()***: Toggle functions for the GUI checkboxes.


## 📚 Data

The script uses the following CSV files for ligand and metal information:
   - [ligands_misc_info.csv](https://raw.githubusercontent.com/hkneiding/tmQMg-L/main/ligands_misc_info.csv)
   - [ligands_fingerprints.csv](https://raw.githubusercontent.com/hkneiding/tmQMg-L/main/ligands_fingerprints.csv)
   - [excel oxydation states métaux.csv](https://raw.githubusercontent.com/sermetsim/metal_complex/main/data/excel%20oxydation%20states%20m%C3%A9taux.csv)
     
These files are loaded from URLs and contain necessary data for the calculations.


## 📌 License

This project is licensed under the [MIT License](LICENSE).


## 📞 Contact

For any questions or issues, please contact [simon.sermet-magdelain@epfl.ch], [giada.foletti@epfl.ch] or [camille.graf@epfl.ch]. 


