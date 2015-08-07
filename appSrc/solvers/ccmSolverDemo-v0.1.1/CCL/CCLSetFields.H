// Create assignment to the mesh object
	fvMesh& mesh = CCLRegions[zoneID];

//	labelList NCellsI(mesh.cellZones());

	labelList NCellsI(mesh.C().size());

// Read in the transport Properties dictionary if it has been modified since last read
#include "CCLReadComponentDictionary.H"

// Assignment the value in the dictionary to the electronic conductivity
#include "CCLDimScalarAssign.H"

	// Create Assignments to Calculated Field Variables
	volScalarField& actPot = actPot_CCL[zoneID];
	surfaceScalarField& protonFlux = protonFlux_CCL[zoneID];
	surfaceScalarField& electronFlux = electronFlux_CCL[zoneID];
	// Create Assignments to Solved Field Variables
	// Electron Transport Fields
	volScalarField& potElectron = potElectron_CCL[zoneID];
	// Proton Transport Fields
	volScalarField& potProton = potProton_CCL[zoneID];
// Create Assignments to Coefficient Field Variables
	// Electron Transport Fields
	volScalarField& condElectronEff = condElectronEff_CCL[zoneID];
	// Proton Transport Fields
	volScalarField& condProtonEff = condProtonEff_CCL[zoneID];
// Create Assignments to Source/Sink Field Variables
	// Electron Transport Fields
	volScalarField& reactionRateBVAnodic = reactionRateBVAnodic_CCL[zoneID];
	volScalarField& reactionRateBVCathodic = reactionRateBVCathodic_CCL[zoneID];
	volScalarField& reactionRateBVSum = reactionRateBVSum_CCL[zoneID];
	volScalarField& reactionRateBVSlopeElectron = reactionRateBVSlopeElectron_CCL[zoneID];
	volScalarField& reactionRateBVSlopeProton = reactionRateBVSlopeProton_CCL[zoneID];
	volScalarField& sourceElectronTotal = sourceElectronTotal_CCL[zoneID];
	volScalarField& sourceProtonTotal = sourceProtonTotal_CCL[zoneID];
	volScalarField& sourceElectronSlope = sourceElectronSlope_CCL[zoneID];
	volScalarField& sourceProtonSlope = sourceProtonSlope_CCL[zoneID];
