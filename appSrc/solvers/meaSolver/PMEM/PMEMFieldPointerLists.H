// Initilialise field pointer lists
	// Calculation Field Variables
		// Stored Gradients
		PtrList<volVectorField> potProtonGrad_PMEM(PMEMRegions.size());
		// Morphology Fields
		// Electrochemistry Fields
		// Material Property Fields
		// Nitrogen Transport Fields
		// Mole Fraction Fields
		// Fluid Property Fields
		// Binary Diffusivity Fields
		// Mixture Property Fields
		// Local Condition Fields
		// Ionomeric Material Fields
		PtrList<volScalarField> H2ODissIon_PMEM(PMEMRegions.size());
		PtrList<volScalarField> condProtonIon_PMEM(PMEMRegions.size());
		PtrList<volScalarField> condThermalIon_PMEM(PMEMRegions.size());
		// Liquid Water Fields
		// Source Term Fields
		PtrList<volScalarField> sourceEnergyOhmicHeating_PMEM(PMEMRegions.size());
	
	// Solved Field Variables	
		// Electron Transport Fields
		// Proton Transport Fields
		PtrList<volScalarField> potProton_PMEM(PMEMRegions.size());
		// Hydrogen Transport Fields
		// Water Vapour Transport Fields
		// Energy Transport Fields
		PtrList<volScalarField> T_PMEM(PMEMRegions.size());
		// Mass Transport (Pressure) Fields
	   	// Liquid Water Transport Fields

	// Coefficient Field Variables
		// Electron Transport Fields
		// Proton Transport Fields
		PtrList<volScalarField> condProtonEff_PMEM(PMEMRegions.size());
		// Hydrogen Transport Fields
		// Water Vapour Transport Fields
		// Energy Transport Fields
		PtrList<volScalarField> condThermalEff_PMEM(PMEMRegions.size());
		PtrList<volScalarField> transportCoeffEnergy_PMEM(PMEMRegions.size());
	  	// Mass Transport (Pressure) Fields
	 	// Liquid Water Transport Fields

	// Source/sink Field Variables
		// Electron Transport Fields
		// Proton Transport Fields
		PtrList<volScalarField> sourceProtonTotal_PMEM(PMEMRegions.size());
		// Hydrogen Transport Fields
		// Water Vapour Transport Fields
		// Energy Transport Fields 
		PtrList<volScalarField> sourceEnergyTotal_PMEM(PMEMRegions.size());
	  	// Mass Transport (Pressure) Fields
	 	// Liquid Water Transport Fields
