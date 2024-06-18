# SQL Database for Chemical Compounds with Tkinter GUI for Drug Discovery
## Introduction
The development of drug discovery evolved rapidly in the past three decades due to efficient 
chemical synthesis and High-Throughput screening (HTS), which focuses on selecting small 
molecules for novel drug targets.

Although the cost of synthesis and screening per compound is low, it may become expensive when 
multiplied by millions. By focusing on “drug-like” molecule, the cost of drug discovery can 
decrease to a manageable level.

The selection criteria for drug-like compounds are diverse and encompass various aspects, including [Lipinski rule of five](https://doi.org/10.1016/s0169-409x(00)00129-0), 
[Verber](https://doi.org/10.1021/jm020017n), [Ghose](https://doi.org/10.1021/cc9800071), [Gleeson](https://doi.org/10.1021/jm701122q), 
[Rule of three for lead-like compounds](https://doi.org/10.1016/s1359-6446(03)02831-9) and 
[Physiochemical properties associated with toxicity](https://doi.org/10.1016/j.bmcl.2008.07.071).

These criteria select compounds based on properties from 1-dimensional such as logP, polar surface area (PSA), hydrogen-bonding groups (H–B
donors or acceptors) to 2-dimensional such as fingerprint similarity, up to 3-dimensional
likes the molecular electrostatic potentials. In this project, these criteria will be applied to evaluate the "drug-like" features of given compounds, demonstrating fundamental skills in the field of drug discovery. 

## Method
RDKit, a Python library, is used to parse and calculate compounds properties in the given Molecules3.sdf file. 

*Chem.Descriptors* module was utilized to calculate the LogP and Molecular weight of the compound from MOL file.
```
## Calculate M.W
MolecularWeights = []
for i in  mols:
    MW = Descriptors.MolWt(i)
    MolecularWeights.append(MW)
print(MolecularWeights)
```

*Chem.Lipinski* module was utilized to calculate the Number of hydrogen bond acceptors, Number of hydrogen bond donors and Number of
rotatable bonds of the compound.
```
H_Bond_Acceptors = []
for i in mols:
    H_Bond_Acceptor = Lipinski.NumHAcceptors(i)
    H_Bond_Acceptors.append(H_Bond_Acceptor)
print(H_Bond_Acceptors)
```

LogD of the compound was directly parsed from the MOL file given since no module available in
RDKIt that can calculate the LogD straightforward. PSA of the compound was calculated using
*Chem.MolSurf* module
```
## Calculate Polar Surface Area ( PSA )
PSAs = []
for i in mols:
    PSA = MolSurf.TPSA(i)
    PSAs.append(PSA)
print(PSAs)
```
Atom number and Heavy atom number were calculated using *Chem.rdchem*.Mol module. 
Molar refractivity was calculated using *Chem.Crippen* module. Formal charge was calculated using *Chem.Rdmolops* module. 
Then, compounds were classified by three classifiers, “Pass” will be marked on compounds that successfully
pass the criteria of classifier, otherwise it will be marked as “Fail”.

Next, the database was created to store the information of compounds (Figure 2). 
The images of compounds were produced using *Draw module*, instead of storing the images directly into database, the file
path of image was populated in the database for querying. 
Finally, a GUI interface was created using Tkinter to interact with the database.

## Result
The GUI interface can search the compound in database through IUPAC name or SMILE 
(Figure 3A). Examples of the IUPAC name and SMILE present in the search boxes to remind 
user (Figure 3A). 

This function can only take one argument at a time to avoid confusing (e.g. 
a IUPAC name or a SMILE). If the IUPAC name and SMILE are inputted at the same time, a 
message box will pop up (Figure3B). Therefore, the procedure of using this search function 
is to clear up the placeholder and enter the IUPAC name or SMILE of the compound which 
user wants to search in the database (Figure4). 

After entering the IUPAC name or SMILE, 
clicking the “Search” button to view the result of querying. “Clear” button needs to be clicked 
before further querying, the information from last querying will be cleaned and ready for 
next querying. “Filter” button can be utilized to filter the compounds in the database, there 
are four filters available for filtering, user can decide which filter to apply, the filter can be 
combined in use (e.g. user can user can apply Ghose, REOS and Verber filters or all four filters  
at the same time, Figure5A). “Forward” and “Back” button can swich the information of 
compound to another compound that passes the same set of filters (Figure5B) Also, user 
can customize molecular weight, LogP, H bond donor and H bond acceptor of the filter 
(Figure5C). But user should aware that this function require user to define all four 
parameters, and the range of each parameter need to be defined.

All compounds in the 
database are presented using the table in the middle of the panel, the information of the 
table will change when filter applied. It is worth noting that a save button will appear after 
using the filter function, this can write the coordinates and information of compounds that 
pass the filter to a sdf file named “output.sdf” in the current directory