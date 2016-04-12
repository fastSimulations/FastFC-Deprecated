// Create assignment to the mesh object
	fvMesh& mesh = PMEMRegions[zoneID];
	labelList NCellsI(mesh.C().size());

// Read in the transport Properties dictionary if it has been modified since last read
#include <PMEMReadComponentDictionary.H>

// Assignment the value in the dictionary to the electronic conductivity
#include <PMEMDimScalarAssign.H>

	// Creat Assignments to Calculated Field Variables
		// Stored Gradients
		volVectorField& potProtonGrad = potProtonGrad_PMEM[zoneID];
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
		volScalarField& H2ODissIon = H2ODissIon_PMEM[zoneID];
		volScalarField& condProtonIon = condProtonIon_PMEM[zoneID];
		volScalarField& condThermalIon = condThermalIon_PMEM[zoneID];
		// Liquid Water Fields
		// Source Term Fields
		volScalarField& sourceEnergyOhmicHeating = sourceEnergyOhmicHeating_PMEM[zoneID];
	
	// Solved Field Variables	
		// Electron Transport Fields
		// Proton Transport Fields
		volScalarField& potProton = potProton_PMEM[zoneID];
		// Hydrogen Transport Fields
		// Water Vapour Transport Fields
		// Energy Transport Fields
		volScalarField& T = T_PMEM[zoneID];
		// Mass Transport (Pressure) Fields
	   	// Liquid Water Transport Fields

	// Coefficient Field Variables
		// Electron Transport Fields
		// Proton Transport Fields
		volScalarField& condProtonEff = condProtonEff_PMEM[zoneID];
		// Hydrogen Transport Fields
		// Water Vapour Transport Fields
		// Energy Transport Fields
		volScalarField& condThermalEff = condThermalEff_PMEM[zoneID];
		volScalarField& transportCoeffEnergy = transportCoeffEnergy_PMEM[zoneID];
	  	// Mass Transport (Pressure) Fields
	 	// Liquid Water Transport Fields

	// Source/sink Field Variables
		// Electron Transport Fields
		// Proton Transport Fields
		volScalarField& sourceProtonTotal = sourceProtonTotal_PMEM[zoneID];
		// Hydrogen Transport Fields
		// Water Vapour Transport Fields
		// Energy Transport Fields 
		volScalarField& sourceEnergyTotal = sourceEnergyTotal_PMEM[zoneID];
	  	// Mass Transport (Pressure) Fields
	 	// Liquid Water Transport Fields
