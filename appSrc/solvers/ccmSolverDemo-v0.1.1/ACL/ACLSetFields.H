// Create assignment to the mesh object
	fvMesh& mesh = ACLRegions[zoneID];

//	labelList NCellsI(mesh.cellZones());

	labelList NCellsI(mesh.C().size());

// Read in the transport Properties dictionary if it has been modified since last read
#include "ACLReadComponentDictionary.H"

// Assignment the value in the dictionary to the electronic conductivity
#include "ACLDimScalarAssign.H"

	// Create Assignments to Calculated Field Variables
	volScalarField& actPot = actPot_ACL[zoneID];
	surfaceScalarField& protonFlux = protonFlux_ACL[zoneID];
	surfaceScalarField& electronFlux = electronFlux_ACL[zoneID];

	// Create Assignments to Solved Field Variables
	// Electron Transport Fields
	volScalarField& potElectron = potElectron_ACL[zoneID];
	// Proton Transport Fields
	volScalarField& potProton = potProton_ACL[zoneID];

// Create Assignments to Coefficient Field Variables
	// Electron Transport Fields
	volScalarField& condElectronEff = condElectronEff_ACL[zoneID];
	// Proton Transport Fields
	volScalarField& condProtonEff = condProtonEff_ACL[zoneID];
// Create Assignments to Source/Sink Field Variables
	// Electron Transport Fields
	volScalarField& reactionRateTHVAnodic = reactionRateTHVAnodic_ACL[zoneID];
	volScalarField& reactionRateTHVCathodic = reactionRateTHVCathodic_ACL[zoneID];
	volScalarField& reactionRateTHVSum = reactionRateTHVSum_ACL[zoneID];
	volScalarField& reactionRateBVAnodic = reactionRateBVAnodic_ACL[zoneID];
	volScalarField& reactionRateBVCathodic = reactionRateBVCathodic_ACL[zoneID];
	volScalarField& reactionRateBVSum = reactionRateBVSum_ACL[zoneID];

	volScalarField& reactionRateTHVSlopeElectron = reactionRateTHVSlopeElectron_ACL[zoneID];
	volScalarField& reactionRateTHVSlopeProton = reactionRateTHVSlopeProton_ACL[zoneID];

	volScalarField& reactionRateBVSlopeElectron = reactionRateBVSlopeElectron_ACL[zoneID];
	volScalarField& reactionRateBVSlopeProton = reactionRateBVSlopeProton_ACL[zoneID];

	volScalarField& sourceElectronTotal = sourceElectronTotal_ACL[zoneID];
	volScalarField& sourceProtonTotal = sourceProtonTotal_ACL[zoneID];

	volScalarField& sourceElectronSlope = sourceElectronSlope_ACL[zoneID];
	volScalarField& sourceProtonSlope = sourceProtonSlope_ACL[zoneID];
