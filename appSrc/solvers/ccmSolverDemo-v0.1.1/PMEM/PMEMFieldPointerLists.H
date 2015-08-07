// Initilialise field pointer lists
//
// Calculated Field Variables
	PtrList<surfaceScalarField> protonFlux_PMEM(PMEMRegions.size());
// Solved Field Variables	
	// Electron Transport Fields
	// Proton Transport Fields
	PtrList<volScalarField> potProton_PMEM(PMEMRegions.size());

// Coefficient Field Variables
	// Electron Transport Fields
	// Proton Transport Fields
	PtrList<volScalarField> condProtonEff_PMEM(PMEMRegions.size());
// Source/sink Field Variables
	// Electron Transport Fields
	PtrList<volScalarField> sourceProtonTotal_PMEM(PMEMRegions.size());
	PtrList<volScalarField> sourceProtonSlope_PMEM(PMEMRegions.size());
