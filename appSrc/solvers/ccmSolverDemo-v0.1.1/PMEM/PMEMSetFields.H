// Create assignment to the mesh object
	fvMesh& mesh = PMEMRegions[zoneID];

// Read in the transport Properties dictionary if it has been modified since last read
#include "PMEMReadComponentDictionary.H"

// Assignment the value in the dictionary to the electronic conductivity
#include "PMEMDimScalarAssign.H"
	// Create Assignments to Calculated Field Variables
	surfaceScalarField& protonFlux = protonFlux_PMEM[zoneID];
	// Create Assignments to Solved Field Variables
	// Electron Transport Fields
	// Proton Transport Fields
	volScalarField& potProton = potProton_PMEM[zoneID];
// Create Assignments to Coefficient Field Variables
	// Electron Transport Fields
	// Proton Transport Fields
	volScalarField& condProtonEff = condProtonEff_PMEM[zoneID];
// Create Assignments to Source/Sink Field Variables
	// Electron Transport Fields
	volScalarField& sourceProtonTotal = sourceProtonTotal_PMEM[zoneID];
	volScalarField& sourceProtonSlope = sourceProtonSlope_PMEM[zoneID];
