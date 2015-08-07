// Initilialise field pointer lists

// Calculated Field Variables
	PtrList<volScalarField> actPot_ACL(ACLRegions.size());
	PtrList<surfaceScalarField> protonFlux_ACL(ACLRegions.size());
	PtrList<surfaceScalarField> electronFlux_ACL(ACLRegions.size());
// Solved Field Variables	
	// Electron Transport Fields
	PtrList<volScalarField> potElectron_ACL(ACLRegions.size());
	// Proton Transport Fields
	PtrList<volScalarField> potProton_ACL(ACLRegions.size());

// Coefficient Field Variables
	// Electron Transport Fields
	PtrList<volScalarField> condElectronEff_ACL(ACLRegions.size());
	// Proton Transport Fields
	PtrList<volScalarField> condProtonEff_ACL(ACLRegions.size());

// Source/sink Field Variables
	// Electron Transport Fields
	PtrList<volScalarField> reactionRateTHVAnodic_ACL(ACLRegions.size());
	PtrList<volScalarField> reactionRateTHVCathodic_ACL(ACLRegions.size());
	PtrList<volScalarField> reactionRateTHVSum_ACL(ACLRegions.size());
	PtrList<volScalarField> reactionRateBVAnodic_ACL(ACLRegions.size());
	PtrList<volScalarField> reactionRateBVCathodic_ACL(ACLRegions.size());
	PtrList<volScalarField> reactionRateBVSum_ACL(ACLRegions.size());

	PtrList<volScalarField> reactionRateTHVSlopeElectron_ACL(ACLRegions.size());
	PtrList<volScalarField> reactionRateTHVSlopeProton_ACL(ACLRegions.size());
	PtrList<volScalarField> reactionRateBVSlopeElectron_ACL(ACLRegions.size());
	PtrList<volScalarField> reactionRateBVSlopeProton_ACL(ACLRegions.size());
	
	PtrList<volScalarField> sourceElectronTotal_ACL(ACLRegions.size());
	PtrList<volScalarField> sourceProtonTotal_ACL(ACLRegions.size());
	PtrList<volScalarField> sourceElectronSlope_ACL(ACLRegions.size());
	PtrList<volScalarField> sourceProtonSlope_ACL(ACLRegions.size());
