import copy as cp
import pandas as pd
import re
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import py3Dmol
import tkinter as tk
from metal_functions import *


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
