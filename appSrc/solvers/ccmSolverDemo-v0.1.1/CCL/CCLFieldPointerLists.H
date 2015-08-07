// Initilialise field pointer lists

// Calculated Field Variables
	PtrList<volScalarField> actPot_CCL(CCLRegions.size());
	PtrList<surfaceScalarField> protonFlux_CCL(CCLRegions.size());
	PtrList<surfaceScalarField> electronFlux_CCL(CCLRegions.size());

// Solved Field Variables	
	// Electron Transport Fields
	PtrList<volScalarField> potElectron_CCL(CCLRegions.size());
	// Proton Transport Fields
	PtrList<volScalarField> potProton_CCL(CCLRegions.size());
// Coefficient Field Variables
	// Electron Transport Fields
	PtrList<volScalarField> condElectronEff_CCL(CCLRegions.size());
	// Proton Transport Fields
	PtrList<volScalarField> condProtonEff_CCL(CCLRegions.size());

// Source/sink Field Variables
	// Electron Transport Fields
	PtrList<volScalarField> reactionRateBVAnodic_CCL(CCLRegions.size());
	PtrList<volScalarField> reactionRateBVCathodic_CCL(CCLRegions.size());
	PtrList<volScalarField> reactionRateBVSum_CCL(CCLRegions.size());

	PtrList<volScalarField> reactionRateBVSlopeElectron_CCL(CCLRegions.size());
	PtrList<volScalarField> reactionRateBVSlopeProton_CCL(CCLRegions.size());

	PtrList<volScalarField> sourceElectronTotal_CCL(CCLRegions.size());
	PtrList<volScalarField> sourceProtonTotal_CCL(CCLRegions.size());

	PtrList<volScalarField> sourceElectronSlope_CCL(CCLRegions.size());
	PtrList<volScalarField> sourceProtonSlope_CCL(CCLRegions.size());
