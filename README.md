# Chemical compounds SQL database with Tkinter GUI
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
[Physiochemical properties associated with toxicity](https://doi.org/10.1016/j.bmcl.2008.07.071).\

These criteria select compounds based on properties from 1-dimensional such as logP, polar surface area (PSA), hydrogen-bonding groups (H–B
donors or acceptors) to 2-dimensional such as fingerprint similarity, up to 3-dimensional
likes the molecular electrostatic potentials. In this project, these criteria will be applied to evaluate the "drug-like" features of given compounds, demonstrating fundamental skills in the field of drug discovery. 

## Method
RDKit, a Python library, is used to parse and calculate compounds properties in the given Molecules3.sdf file. 
Compound information is stored in relational database and a GUI is configured using Tkinter for user interaction.