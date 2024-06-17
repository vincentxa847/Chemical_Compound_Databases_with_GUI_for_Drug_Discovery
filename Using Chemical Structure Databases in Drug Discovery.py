from turtle import clear

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski
from rdkit.Chem import MolSurf
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Crippen
from rdkit.Chem import rdmolops
from tkinter import *
from tkinter import messagebox
from tkinter import ttk
from PIL import ImageTk,Image
import sqlite3
import re

data = Chem.SDMolSupplier('Molecules3.sdf')
# print(type(data))
mols = [x for x in data] # this mol supplies 2D coordinates

## Retrieve Name
names = []
for i in mols:
    name = i.GetProp("_Name")
    names.append(name)
print(names)

## Calculate M.W
MolecularWeights = []
for i in  mols:
    MW = Descriptors.MolWt(i)
    MolecularWeights.append(MW)
print(MolecularWeights)

## Calculate LogP
LogPs = []
for i in mols:
    LogP = Descriptors.MolLogP(i)
    LogPs.append(LogP)
print(LogPs)

## Retrieve logD directly from the mol file given
LogDs = []
with open('Molecules3.sdf', 'r') as sdf:
    contents = sdf.readlines()
for step, content in enumerate(contents):
    if content == '>  <LogD>\n':
        LogDs.append(float(contents[step+1].rstrip()))
print(LogDs)

### using rdkit.Chem.Lipinski module to fetch these properties
## Number of Hydrogen Bond Acceptors
H_Bond_Acceptors = []
for i in mols:
    H_Bond_Acceptor = Lipinski.NumHAcceptors(i)
    H_Bond_Acceptors.append(H_Bond_Acceptor)
print(H_Bond_Acceptors)

## Number of Hydrogen Bond Donors
H_Bond_Donors = []
for i in mols:
    H_Bond_Donor = Lipinski.NumHDonors(i)
    H_Bond_Donors.append(H_Bond_Donor)
print(H_Bond_Donors)

## Number of Rotatable Bonds
Rotatable_Bonds = []
for i in mols:
    Rotatable_Bond = Lipinski.NumRotatableBonds(i)
    Rotatable_Bonds.append(Rotatable_Bond)
print(Rotatable_Bonds)

## Calculate Polar Surface Area ( PSA )
PSAs = []
for i in mols:
    PSA = MolSurf.TPSA(i)
    PSAs.append(PSA)
print(PSAs)

## Number of atoms
Atoms_numbers = []
for i in mols:
    Atoms_number = Mol.GetNumAtoms(i)
    Atoms_numbers.append(Atoms_number)
print(Atoms_numbers)

## molar_refractivity
Molar_refractivitys = []
for i in mols:
    Molar_refractivity = Crippen.MolMR(i)
    Molar_refractivitys.append(Molar_refractivity)
print(Molar_refractivitys)

## formal charges
Formal_charges = []
for i in mols:
    Formal_charge = rdmolops.GetFormalCharge(i)
    Formal_charges.append(Formal_charge)
print(Formal_charges)

## Heavy_atom_count
Heavy_atoms = []
for i in mols:
    Heavy_atom = Mol.GetNumHeavyAtoms(i)
    Heavy_atoms.append(Heavy_atom)
print(Heavy_atoms)

## ------- Lipinski Rule Of five -------
Rule_of_Five_list = []
""" 
    Compounds that meet the following criteria will be marked as "Pass"
    mass <= 500 AND
    logP <= 5 AND
    donorCount <= 5 AND
    acceptorCount <= 10 
"""
for compound in range(0,len(mols)):
    if MolecularWeights[compound] <= 500 and LogPs[compound] <= 5 and H_Bond_Donors[compound] <= 5 and H_Bond_Acceptors[compound] <= 10:
        Rule_of_Five_list.append("Pass")
    else:
        Rule_of_Five_list.append("Fail")
print(Rule_of_Five_list)

## ------- Ghose filter -------
Ghose_filter_list = []
""" Molecular weight between 160 and 480
    LogP between -0.4 and +5.6
    Atom count between 20 and 70
    Molar refractivity between 40 and 130 
"""

for compound in range(0,len(mols)):
    if 160<=MolecularWeights[compound] <= 480 and -0.4<=LogPs[compound] <= 5.6 and 20<=Atoms_numbers[compound] <= 70 and 40<=Molar_refractivitys[compound] <= 130:
        Ghose_filter_list.append("Pass")
    else:
        Ghose_filter_list.append("Fail")
print(Ghose_filter_list)


## ------- REOS filter (Rapid Elimination Of Swill filter) -------
REOS_filter_list = []
""" Molecular weight between 200 and 500
    LogP between -5.0 and +5.0
    H-bond donor count between 0 and 5
    H-bond acceptor count between 0 and 10
    Formal charge between -2 and +2
    Rotatable bond count between 0 and 8
    Heavy atom count between 15 and 50 
"""

for compound in range(0,len(mols)):
    if 200 <=MolecularWeights[compound] <= 500 and -5 <=LogPs[compound] <= 5 and 0 <=H_Bond_Donors[compound] <= 5 and 0 <=H_Bond_Acceptors[compound] <= 10 and -2<= Formal_charges[compound] <=2 and 0<=Rotatable_Bonds[compound]<=8 and 15 <=Heavy_atoms[compound] <= 50:
        REOS_filter_list.append("Pass")
    else:
        REOS_filter_list.append("Fail")
print(REOS_filter_list)


## ------- Verber filter -------
Verber_filter_list = []
"""
    Rotatable bonds <= 10
    Topological polar surface area <= 140
"""

for compound in range(0,len(mols)):
    if Rotatable_Bonds[compound]<=10 and PSAs[compound]<=140:
        Verber_filter_list.append("Pass")
    else:
        Verber_filter_list.append("Fail")
print(Verber_filter_list)


## ------- Rule of three for lead-like compounds -------
Rule_of_Three_list = []
""" 
    Molecular weight <= 300
    LogP <= 3
    H-bond donor <= 3
    H-bond acceptor count <= 3
    Rotatable bond count <= 3
"""
for compound in range(0,len(mols)):
    if MolecularWeights[compound]<=300 and LogPs[compound]<=3 and H_Bond_Donors[compound]<=3 and H_Bond_Acceptors[compound] <=3 and Rotatable_Bonds[compound] <=3:
        Rule_of_Three_list.append("Pass")
    else:
        Rule_of_Three_list.append("Fail")
print(Rule_of_Three_list)

# #%%
# Create database for all compounds (100 compounds)
with open("Create.sql", 'r') as sql_file:
    sql_script = sql_file.read()
db = sqlite3.connect("CompoundDatabase.db")
try:
    cursor = db.cursor()
    cursor.executescript(sql_script)
    db.commit()
    db.close()
except sqlite3.OperationalError: # Database already exists
    pass

try:
    # Insert data into database
    database = "CompoundDatabase.db"
    conn = sqlite3.connect(database)
    cur = conn.cursor()

    compound_number = 0
    for compound in mols[:]:

        # Draw all the picture for 100 compounds
        # save the image in a folder and save the file name of the image into the database
        image = Draw.MolToImage(compound)
        image_name = "image_" + str(compound_number+1) + ".png"
        image.save("structure\\" + image_name)

        # information to store into database
        smiles = Chem.MolToSmiles(compound)
        IUPAC = compound.GetProp('_Name')
        MolecularWeight = round(Descriptors.MolWt(compound),3)
        LogP = round(Descriptors.MolLogP(compound),3)
        H_Bond_Donor = Lipinski.NumHDonors(compound)
        H_Bond_Acceptor = Lipinski.NumHAcceptors(compound)
        Rotatable_Bond = Lipinski.NumRotatableBonds(compound)
        Number_of_Atom = Mol.GetNumAtoms(compound)
        Molar_Refractivity = round(Crippen.MolMR(compound),3)
        Formal_Charge = rdmolops.GetFormalCharge(compound)
        Heavy_Atom_Count = Mol.GetNumHeavyAtoms(compound)
        PSA = round(MolSurf.TPSA(compound),3)

        Lipinski_Rule = Rule_of_Five_list[compound_number]
        Ghose_Rule = Ghose_filter_list[compound_number]
        REOS_Rule = REOS_filter_list[compound_number]
        Verber_Rule = Verber_filter_list[compound_number]


        sql = f'INSERT INTO Database VALUES ("{smiles}","{image_name}","{IUPAC}",{float(MolecularWeight)},{float(LogP)},{int(H_Bond_Donor)},{int(H_Bond_Acceptor)},{int(Rotatable_Bond)},{int(Number_of_Atom)},{float(Molar_Refractivity)},{int(Formal_Charge)},{int(Heavy_Atom_Count)},{float(PSA)},"{Lipinski_Rule}","{Ghose_Rule}","{REOS_Rule}","{Verber_Rule}")'
        # sql = "DELETE FROM Database"
        cur.execute(sql)
        compound_number += 1

    conn.commit()
    cur.close()
    conn.close()

except sqlite3.IntegrityError: # data already insert
    pass