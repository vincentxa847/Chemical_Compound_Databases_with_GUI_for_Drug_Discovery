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

#%%
## Tkinter section
root = Tk()
root.title("Chemical Database 2802815")
# root.geometry("1200x1000")


# Put the compound want to query
Smile_label = Label(root,text="Smile")
Smile_label.grid(row=1,column=0)
Smile = Entry(root,width=70)
Smile.insert(END, "CCCCC/C=C/CC(O)/C=C/C=C/C/C=C/CCCC(=O)O") # put the template  for querying
Smile.grid(row=1,column=1,padx=20)


IUPAC_label = Label(root,text="IUPAC name")
IUPAC_label.grid(row=2,column=0)
IUPAC = Entry(root,width=70)
IUPAC.insert(END, "11-dehydrothromboxane B2")  # put the template  for querying
IUPAC.grid(row=2,column=1,padx=20)

query_label = Label(root,text="") # set the default query_label to let "query_label.grid_forget()"  work


def query():
    # connect to database
    conn = sqlite3.connect("CompoundDatabase.db")
    c = conn.cursor()

    global query_label
    global image_label
    global query_result
    global query_image # image needs to be set as global variable to present normally
    global Smile
    global IUPAC

    # remove last query result
    query_label.grid_forget()

    Smile_to_be_queried = Smile.get()
    IUPAC_to_be_queried = IUPAC.get()

    # if user put a smile and IUPAC name, then a Message Boxes will pop up to avoid confusion
    if Smile.get() != "" and IUPAC.get()!= "":
        response = messagebox.showinfo("HINT", "Please input a smile OR a IUPAC name for querying")
        Label(root, text=response).pack()

    # query the database (using smile or IUPAC name to query)
    c.execute(f"SELECT * FROM Database WHERE (Smile LIKE '{Smile_to_be_queried}') OR (IUPAC_NAME LIKE '{IUPAC_to_be_queried}')")
    results = c.fetchall()


    # show the result
    query_result = ''
    for result in results: # results should be a single list
        query_result = f"Smile : {result[0]} \nIUPAC_name : {result[2]} \nMolecular Weight : {result[3]} \nLogP : {result[4]}\nH Bond Donor : {result[5]} \nH Bond Acceptor : {result[6]} \nRotatable Bonds : {result[7]} \nNumber of Atoms : {result[8]}\nMolar Refractivity : {result[9]}\nFormal Charge : {result[10]}\nLipinski Rule : {result[13]}\nGhose Rule : {result[14]}\nREOS Rule : {result[15]}\nVerber Rule : {result[16]}"
        structure = "structure/" + result[1] # file path for png file

    query_image = ImageTk.PhotoImage(Image.open(structure))

    query_label = Label(root, text=query_result)
    query_label.grid(row=5,column=1)
    image_label = Label(image=query_image)
    image_label.grid(row=5,column=0,pady=10)

    # when button press, the text in search box will be deleted
    Smile.delete(0,END)
    IUPAC.delete(0,END)

    # Set the clear button for query
    button_clear = Button(root, text="Clear", command=clear_for_query)
    button_clear.grid(row=3,column=2,pady=10)


    # Commit changes
    conn.commit()
    # close connection
    conn.close()

def Filter():
    # connect to database
    conn = sqlite3.connect("CompoundDatabase.db")
    c = conn.cursor()

    global image
    global query_label
    global filter_result

    # set as global variable to let grid.forget() in forward and back work
    global list_of_structure
    global list_of_compound_information
    global no_data_available

    global buttom_back
    global buttom_forward
    global button_save



    results = []

    # using customize filter
    # all parameter must be input to use customize filter
    if M_W_customize.get() != "" and LogP_customize.get() != "" and Hdonor_customize.get() != "" and Hacceptor_customize.get() != "" :
        M_W_customize_from, M_W_customize_to = re.split("-", M_W_customize.get())
        LogP_customize_from, LogP_customize_to = re.split("-", LogP_customize.get())
        Hdonor_customize_from, Hdonor_customize_to = re.split("-", Hdonor_customize.get())
        Hacceptor_customize_from, Hacceptor_customize_from = re.split("-", Hacceptor_customize.get())
        c.execute(f"SELECT * FROM Database WHERE Molecular_Weight BETWEEN {float(M_W_customize_from)} AND {float(M_W_customize_to)} AND LogP BETWEEN {float(LogP_customize_from)} AND {float(LogP_customize_to)} AND H_Bond_Donors BETWEEN {float(Hdonor_customize_from)} AND {float(Hdonor_customize_to)} AND H_Bond_Acceptors BETWEEN {float(Hacceptor_customize_from)} AND {float(Hdonor_customize_to)}")
        results = c.fetchall()
        print(results)


    # 1 Filter by Ghose Rule
    if Ghose_v.get() == "on" and Lipinski_Rule_v.get() == "off" and REOS_v.get() == "off" and Verber_v.get() == "off":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Pass' AND Lipinski_Rule like 'Fail' AND REOS_Rule like 'Fail' AND Verber_Rule like 'Fail'")
        results = c.fetchall()
        print(results)
    # 2 Filter by Lipinski_Rule
    if Ghose_v.get() == "off" and Lipinski_Rule_v.get() == "on" and REOS_v.get() == "off" and Verber_v.get() == "off":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Fail' AND Lipinski_Rule like 'Pass' AND REOS_Rule like 'Fail' AND Verber_Rule like 'Fail' ")
        results = c.fetchall()
        print(results)
    # 3 Filter by REOS_Rule
    if Ghose_v.get() == "off" and Lipinski_Rule_v.get() == "off" and REOS_v.get() == "on" and Verber_v.get() == "off":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Fail' AND Lipinski_Rule like 'Fail' AND REOS_Rule like 'Pass' AND Verber_Rule like 'Fail' ")
        results = c.fetchall()
        print(results)
    # 4 Filter by Verber_Rule
    if Ghose_v.get() == "off" and Lipinski_Rule_v.get() == "off" and REOS_v.get() == "off" and Verber_v.get() == "on":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Pass' AND Lipinski_Rule like 'Fail' AND REOS_Rule like 'Fail' AND Verber_Rule like 'Pass'")
        results = c.fetchall()
        print(results)
    # 5 Filter by Ghose and Lipinski_Rule
    if Ghose_v.get() == "on" and Lipinski_Rule_v.get() == "on" and REOS_v.get() == "off" and Verber_v.get() == "off":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Pass' AND Lipinski_Rule like 'Pass' AND REOS_Rule like 'Fail' AND Verber_Rule like 'Fail' ")
        results = c.fetchall()
        print(results)
    # 6 Filter by Ghose and REOS_Rule
    if Ghose_v.get() == "on" and Lipinski_Rule_v.get() == "off" and REOS_v.get() == "on" and Verber_v.get() == "off":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Pass' AND Lipinski_Rule like 'Fail' AND REOS_Rule like 'Pass' AND Verber_Rule like 'Fail'")
        results = c.fetchall()
        print(results)
    # 7 Filter by Ghose and Verber_Rule
    if Ghose_v.get() == "on" and Lipinski_Rule_v.get() == "off" and REOS_v.get() == "off" and Verber_v.get() == "on":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Pass' AND Lipinski_Rule like 'Fail' AND REOS_Rule like 'Fail' AND Verber_Rule like 'Pass'")
        results = c.fetchall()
        print(results)
    # 8 Filter by Lipinski and REOS_Rule
    if Ghose_v.get() == "off" and Lipinski_Rule_v.get() == "on" and REOS_v.get() == "on" and Verber_v.get() == "off":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Fail' AND Lipinski_Rule like 'Pass' AND REOS_Rule like 'Pass' AND Verber_Rule like 'Fail' ")
        results = c.fetchall()
        print(results)
    # 9 Filter by Lipinski and Verber_Rule
    if Ghose_v.get() == "off" and Lipinski_Rule_v.get() == "on" and REOS_v.get() == "off" and Verber_v.get() == "on":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Fail' AND Lipinski_Rule like 'Pass' AND REOS_Rule like 'Fail' AND Verber_Rule like 'Pass' ")
        results = c.fetchall()
        print(results)
    # 10 Filter by REOS and Verber_Rule
    if Ghose_v.get() == "off" and Lipinski_Rule_v.get() == "off" and REOS_v.get() == "on" and Verber_v.get() == "on":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Fail' AND Lipinski_Rule like 'Fail' AND REOS_Rule like 'Pass' AND Verber_Rule like 'Pass'")
        results = c.fetchall()
        print(results)
    # 11 Filter by Ghose, Lipinski and REOS_Rule
    if Ghose_v.get() == "on" and Lipinski_Rule_v.get() == "on" and REOS_v.get() == "on" and Verber_v.get() == "off":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Pass' AND Lipinski_Rule like 'Pass' AND REOS_Rule like 'Pass' AND Verber_Rule like 'Fail'")
        results = c.fetchall()
        print(results)
    # 12 Filter by Ghose, Lipinski and Verber_Rule
    if Ghose_v.get() == "on" and Lipinski_Rule_v.get() == "on" and REOS_v.get() == "off" and Verber_v.get() == "on":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Pass' AND Lipinski_Rule like 'Pass' AND REOS_Rule like 'Fail' AND Verber_Rule like 'Pass'")
        results = c.fetchall()
        print(results)
    # 13 Filter by Lipinski, REOS and Verber_Rule
    if Ghose_v.get() == "off" and Lipinski_Rule_v.get() == "on" and REOS_v.get() == "on" and Verber_v.get() == "on":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Fail' AND Lipinski_Rule like 'Pass' AND REOS_Rule like 'Pass' AND Verber_Rule like 'Pass' ")
        results = c.fetchall()
        print(results)
    # 14 Filter by Ghsose, REOS and Verber_Rule
    if Ghose_v.get() == "on" and Lipinski_Rule_v.get() == "off" and REOS_v.get() == "on" and Verber_v.get() == "on":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Pass' AND Lipinski_Rule like 'Fail' AND REOS_Rule like 'Pass' AND Verber_Rule like 'Pass'")
        results = c.fetchall()
        print(results)
    # 15 Filter by no_Rule
    if Ghose_v.get() == "off" and Lipinski_Rule_v.get() == "off" and REOS_v.get() == "off" and Verber_v.get() == "off":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Fail' AND Lipinski_Rule like 'Fail' AND REOS_Rule like 'Fail' AND Verber_Rule like 'Fail' ")
        results = c.fetchall()
        print(results)
    # Filter by all_Rule
    if Ghose_v.get() == "on" and Lipinski_Rule_v.get() == "on" and REOS_v.get() == "on" and Verber_v.get() == "on":
        c.execute( "SELECT * FROM Database WHERE Ghose_Rule like 'Pass' AND Lipinski_Rule like 'Pass' AND REOS_Rule like 'Pass' AND Verber_Rule like 'Pass'")
        results = c.fetchall()
        print(results)




    # put the compound information and image in a tuple and pack it in a list
    print_result = ''
    filter_result = []
    for result in results: # results should be a single list
        print_result = f"Smile : {result[0]} \nIUPAC_name : {result[2]} \nMolecular Weight : {result[3]} \nLogP : {result[4]}\nH Bond Donor : {result[5]} \nH Bond Acceptor : {result[6]} \nRotatable Bonds : {result[7]} \nNumber of Atoms : {result[8]}\nMolar Refractivity : {result[9]}\nFormal Charge : {result[10]}\nLipinski Rule : {result[13]}\nGhose Rule : {result[14]}\nREOS Rule : {result[15]}\nVerber Rule : {result[16]}"
        structure = "structure/" + result[1] # file path for png file
        image = ImageTk.PhotoImage(Image.open(structure))
        filter_result.append((print_result, image))

    ## Create table to view information of all compounds in database
    # define columns
    columns = ('IUPAC', 'Molecular Weight', 'LogP', 'H Bond Donors', 'H Bond Acceptors')

    tree = ttk.Treeview(root, columns=columns, show='headings')

    # define headings
    tree.heading('IUPAC', text='IUPAC')
    tree.heading('Molecular Weight', text='Molecular Weight')
    tree.heading('LogP', text='LogP')
    tree.heading('H Bond Donors', text='H Bond Donors')
    tree.heading('H Bond Acceptors', text='H Bond Acceptors')

    # add data to the treeview
    for row in results:
        insert = (row[2], row[3], row[4],row[5],row[6]) # IUPAC,molecular weight, LogP, Donor, Acceptor
        tree.insert('', END, values=insert)

    tree.grid(row=4, column=0, sticky='nsew', columnspan=4)

    # add a scrollbar
    scrollbar = ttk.Scrollbar(root, orient=VERTICAL, command=tree.yview)
    tree.configure(yscroll=scrollbar.set)
    scrollbar.grid(row=4, column=4, sticky='ns')



    # Commit changes
    conn.commit()
    # close connection
    conn.close()

    # set the default button status (back,forward,clear,save)
    buttom_back = Button(root, text="Back", command=back, state=DISABLED)
    buttom_forward = Button(root, text="Forward", command=lambda: forward(2))
    button_clear = Button(root, text="Clear", command=clear_for_filter)

    buttom_back.grid(row=6, column=0)
    buttom_forward.grid(row=6, column=2)
    button_clear.grid(row=3,column=2,pady=10)

    button_save = Button(root,text="Save",command=lambda : save(results))
    button_save.grid(row=3,column=3)

    if filter_result == []:
        no_data_available = Label(root, text="No compound in database match this criteria")
        no_data_available.grid(row=6,column=0,columnspan=3,sticky=W+E)
        try:
            list_of_compound_information.grid_forget()
            list_of_structure.grid_forget()
        except NameError:
            pass

        buttom_forward.grid_forget()
        buttom_back.grid_forget()
        return # exit the function without index filter_result(empty list)

    # set the first information
    list_of_structure = Label(image=filter_result[0][1])
    list_of_structure.grid(row=5,column=0)
    list_of_compound_information = Label(text=filter_result[0][0])
    list_of_compound_information.grid(row=5,column=1)



# forward and back button for viewing list of compounds when filtering
def forward(image_number):
    global list_of_structure  # renew image
    global list_of_compound_information  # Update list_of_compound_information
    global buttom_forward
    global buttom_back

    list_of_structure.grid_forget()
    list_of_compound_information.grid_forget()
    list_of_structure = Label(image=filter_result[image_number-1][1])
    list_of_compound_information = Label(text=filter_result[image_number-1][0])
    # renew button status, Here the  image_number+1 and image_number-1 is the image number for new button
    buttom_forward = Button(root,text="forward",command=lambda : forward(image_number+1))
    buttom_back = Button(root,text="back",command=lambda : back(image_number-1))

    if image_number == len(filter_result):
        buttom_forward = Button(root,text="forward",state=DISABLED)

    # put the renew button into windows
    list_of_structure.grid(row=5,column=0)
    list_of_compound_information.grid(row=5,column=1)
    buttom_back.grid(row=6,column=0)
    buttom_forward.grid(row=6,column=2)

def back(image_number):
    global list_of_structure
    global list_of_compound_information
    global buttom_forward
    global buttom_back


    list_of_structure.grid_forget()
    list_of_compound_information.grid_forget()
    list_of_structure = Label(image=filter_result[image_number-1][1])
    list_of_compound_information = Label(text=filter_result[image_number - 1][0])
    buttom_forward = Button(root,text="forward",command=lambda : forward(image_number+1))
    buttom_back = Button(root,text="back",command=lambda : back(image_number-1))

    if image_number == 1:
        buttom_back = Button(root,text="back",state=DISABLED)

    list_of_structure.grid(row=5,column=0)
    list_of_compound_information.grid(row=5,column=1)
    buttom_back.grid(row=6,column=0)
    buttom_forward.grid(row=6,column=2)

def clear_for_query():
    # clear the result of search
    query_label.grid_forget()
    image_label.grid_forget()

def clear_for_filter():
    global list_of_structure
    global list_of_compound_information
    global no_data_available

    # grid_forget would not work if put in try except block
    button_save.grid_forget()
    buttom_back.grid_forget()
    buttom_forward.grid_forget()

    list_of_structure.grid_forget()
    list_of_compound_information.grid_forget()


    try:
        # if no_data_available occur when querying, then clear button can also remove it
        no_data_available.grid_forget()
    except NameError:
        pass

    # if clean button press, then table will return to all compounds
    if __name__ == '__main__':
        create_table_for_all_compound()


# Create a function to write compounds that pass filter into a sdf file
def save(data):
    w = Chem.SDWriter("output.sdf")
    for smile in data:
        mol = Chem.MolFromSmiles(smile[0])
        AllChem.Compute2DCoords(mol)
        Chem.MolToMolBlock(mol)
        # mol = Chem.AddHs(data, addCoords=True)  # hydorgen will be added for 3D coordinates
        mol.SetProp('IUPAC_name', str(smile[2]))
        mol.SetProp('Smile', str(smile[0]))
        mol.SetProp('Molecular Weight', str(smile[3]))
        w.write(mol)
    w.close()



## Create table to view information of all compounds in database
# define columns
def create_table_for_all_compound():
    columns = ('IUPAC', 'Molecular Weight','Smile')

    tree = ttk.Treeview(root, columns=columns, show='headings')

    columns = ('IUPAC', 'Molecular Weight', 'LogP', 'H Bond Donors', 'H Bond Acceptors')

    tree = ttk.Treeview(root, columns=columns, show='headings')

    # define headings
    tree.heading('IUPAC', text='IUPAC')
    tree.heading('Molecular Weight', text='Molecular Weight')
    tree.heading('LogP', text='LogP')
    tree.heading('H Bond Donors', text='H Bond Donors')
    tree.heading('H Bond Acceptors', text='H Bond Acceptors')


    # connect to database
    conn = sqlite3.connect("CompoundDatabase.db")
    c = conn.cursor()

    # generate data
    c.execute("SELECT IUPAC_name,Molecular_Weight,LogP,H_Bond_Donors,H_Bond_Acceptors FROM Database ")
    rows = c.fetchall()


    # add data to the treeview
    for row in rows:
        tree.insert('', END, values=row)

    tree.grid(row=4, column=0, sticky='nsew',columnspan=4)

    # add a scrollbar
    scrollbar = ttk.Scrollbar(root, orient=VERTICAL, command=tree.yview)
    tree.configure(yscroll=scrollbar.set)
    scrollbar.grid(row=4, column=4, sticky='ns')

if __name__ == '__main__':
    create_table_for_all_compound()



# Create Search button
search_btn = Button(root,text="Search",command=query)
search_btn.grid(row=3,column=0,pady=10)

# Create Filter button
filter_btn = Button(root,text="Filter",command=Filter)
filter_btn.grid(row=3,column=1,pady=10)

# Create Clear button
clear_btn = Button(root,text="Clear",command=clear)
clear_btn.grid(row=3,column=2,pady=10)


# checkbox for filtering
Ghose_v = StringVar()
Ghose = Checkbutton(root, text="Ghose Filter", variable=Ghose_v, onvalue="on", offvalue="off")
Ghose.deselect()
Ghose.grid(row=1,column=2,pady=10)

Lipinski_Rule_v = StringVar()
Lipinski_Rule = Checkbutton(root, text="Lipinski Rule", variable=Lipinski_Rule_v, onvalue="on", offvalue="off")
Lipinski_Rule.deselect()
Lipinski_Rule.grid(row=1,column=3,pady=10)

REOS_v = StringVar()
REOS = Checkbutton(root, text="REOS Filter", variable=REOS_v, onvalue="on", offvalue="off")
REOS.deselect()
REOS.grid(row=2,column=2,pady=10)

Verber_v = StringVar()
Verber = Checkbutton(root, text="Verber Filter", variable=Verber_v, onvalue="on", offvalue="off")
Verber.deselect()
Verber.grid(row=2,column=3,pady=10)

# customize filter

# need a widget to hold row = 0
Empty = Label(root, text="").grid(row=0,column=0,sticky='nsew')

M_W_customize_label = Label(root, text="Molecular Weight").place(x=20 , y=10)
# M_W_customize_label.grid(row=0,column=0)
M_W_customize = Entry(root,width=10)
M_W_customize.place(x=130 , y=10) # For the entry want to get or insert,  place, grid .etc can not be done in a single line
M_W_customize.insert(0, "160-480")
# M_W_customize.grid(row=0,column=1,pady=10)
LogP_customize_label = Label(root, text="LogP").place(x=180 , y=10)
# LogP_customize_label.grid(row=0,column=2)
LogP_customize = Entry(root,width=10)
LogP_customize.place(x=220 , y=10)
# LogP_customize.grid(row=0,column=3,pady=10)
Hdonor_customize_label = Label(root, text="H Donor").place(x=270 , y=10)
# Hdonor_customize_label.grid(row=0,column=4)
Hdonor_customize = Entry(root,width=10)
Hdonor_customize.place(x=330 , y=10)
# Hdonor_customize.grid(row=0,column=5,pady=10)
Hacceptor_customize_label = Label(root, text="H Acceptor").place(x=380 , y=10)
# Hacceptor_customize_label.grid(row=0,column=6)
Hacceptor_customize = Entry(root,width=10)
Hacceptor_customize.place(x=450 , y=10)
# Hacceptor_customize.grid(row=0,column=7,pady=10)

root.mainloop()