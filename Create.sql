CREATE TABLE Database (
	Smile VARCHAR(100) NOT NULL,
	Structure VARCHAR(100),  -- STORE THR FILE NAME OF PICTURE
	IUPAC_NAME VARCHAR(100),
	Molecular_Weight DECIMAL(3,3),
	LogP DECIMAL(2,3),
	H_Bond_Donors INT,
	H_Bond_Acceptors INT,
	Rotatable_Bonds INT,
    Number_of_Atom INT,
    Molar_Refractivity DECIMAL(3,3),
    Formal_Charge INT,
    Heavy_Atom_Count INT,
    PSA DECIMAL(3,3),
    Lipinski_Rule CHAR(4),
    Ghose_Rule CHAR(4),
    REOS_Rule CHAR(4),
    Verber_Rule CHAR(4),


	PRIMARY KEY (Smile)
	);