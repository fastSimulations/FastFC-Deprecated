// Create assignment to the mesh object
	const fvMesh& mesh = CPTLRegions[zoneID];

//	Setup the variable for the cell count
	labelList NCellsI(mesh.C().size());

// Read in the transport Properties dictionary if it has been modified since last read
#include "CPTLReadComponentDictionary.H"

// Assignment the value in the dictionary to the electronic conductivity
#include "CPTLDimScalarAssign.H"

// Create Assignments to the Calculated Field Variables
	// Stored Gradients
	volVectorField& pMixGrad = pMixGrad_CPTL[zoneID];
	volVectorField& potElectronGrad = potElectronGrad_CPTL[zoneID];
	// Morphology Fields
	volScalarField& volFPore = volFPore_CPTL[zoneID];
	volScalarField& volFPoreEff = volFPoreEff_CPTL[zoneID];
	// Material Property Fields
	volScalarField& kPermDry = kPermDry_CPTL[zoneID];
	volScalarField& kPermRelLiq = kPermRelLiq_CPTL[zoneID];
	volScalarField& kPermRelGas = kPermRelGas_CPTL[zoneID];
	// Nitrogen Transport Fields
	volScalarField& wN2 = wN2_CPTL[zoneID];
	// Mole Fraction Fields
	volScalarField& xO2 = xO2_CPTL[zoneID];
	volScalarField& xN2 = xN2_CPTL[zoneID];
	volScalarField& xH2OVap = xH2OVap_CPTL[zoneID];
	// Fluid Property Fields
	volScalarField& viscO2 = viscO2_CPTL[zoneID];
	volScalarField& viscN2 = viscN2_CPTL[zoneID];
	volScalarField& viscH2OVap = viscH2OVap_CPTL[zoneID];
	// Binary Diffusivity Fields
	volScalarField& diffO2N2 = diffO2N2_CPTL[zoneID];
	volScalarField& diffO2H2O = diffO2H2O_CPTL[zoneID];
	volScalarField& diffN2H2O = diffN2H2O_CPTL[zoneID];
	// Mixture Property Fields
	volScalarField& molWeightMix = molWeightMix_CPTL[zoneID];
	volScalarField& densityMix = densityMix_CPTL[zoneID];
	volScalarField& viscMix = viscMix_CPTL[zoneID];
	volScalarField& velAvg = velAvg_CPTL[zoneID];
	// Local Condition Fields
	volScalarField& RH = RH_CPTL[zoneID];
	volScalarField& pO2 = pO2_CPTL[zoneID];
	volScalarField& pH2OVap = pH2OVap_CPTL[zoneID];
	volScalarField& pH2OSat = pH2OSat_CPTL[zoneID];
	// Liquid Water Fields
	volScalarField& jS = jS_CPTL[zoneID];
	volScalarField& pCap = pCap_CPTL[zoneID];
	volScalarField& djSDs = djSDs_CPTL[zoneID];
	volScalarField& dpCapDs = dpCapDs_CPTL[zoneID];
	volScalarField& swPhaseChangeLiqVap = swPhaseChangeLiqVap_CPTL[zoneID];
	// Source Term Fields
	volScalarField& sourceEnergyPhaseChange = sourceEnergyPhaseChange_CPTL[zoneID];
	volScalarField& sourceEnergyOhmicHeating = sourceEnergyOhmicHeating_CPTL[zoneID];
	volScalarField& sourcePhaseChangeLiqVap = sourcePhaseChangeLiqVap_CPTL[zoneID];

// Create Assignments to Solved Field Variables
	// Electron Transport Fields
	volScalarField& potElectron = potElectron_CPTL[zoneID];
	// Oxygen Transport Fields
   	volScalarField& wO2 = wO2_CPTL[zoneID];
	// Water Vapour Transport Fields
	volScalarField& wH2OVap = wH2OVap_CPTL[zoneID];
	// Energy Transport Fields
	volScalarField& T = T_CPTL[zoneID];
	// Mass Transport (Pressure) Fields
   	volScalarField& pMix = pMix_CPTL[zoneID];
   	// Liquid Water Transport Fields
	volScalarField& sLiq = sLiq_CPTL[zoneID];

// Create Assignments to Coefficient Field Variables
	// Electron Transport Fields
	volScalarField& condElectronEff = condElectronEff_CPTL[zoneID];
	// Oxygen Transport Fields
	volScalarField& diffO2MixEff = diffO2MixEff_CPTL[zoneID];
	volScalarField& transportCoeffO2 = transportCoeffO2_CPTL[zoneID];
	// Water Vapour Transport Fields
	volScalarField& diffH2OVapMixEff = diffH2OVapMixEff_CPTL[zoneID];
	volScalarField& transportCoeffH2OVap = transportCoeffH2OVap_CPTL[zoneID];
	// Energy Transport Fields
	volScalarField& condThermalEff = condThermalEff_CPTL[zoneID];
	volScalarField& transportCoeffEnergy = transportCoeffEnergy_CPTL[zoneID];
  	// Mass Transport (Pressure) Fields
	volScalarField& dynamicPerm = dynamicPerm_CPTL[zoneID];
 	// Liquid Water Transport Fields
	volScalarField& diffH2OLiq = diffH2OLiq_CPTL[zoneID];

// Create Assignments to Source/Sink Field Variables
	// Electron Transport Fields
	volScalarField& sourceElectronTotal = sourceElectronTotal_CPTL[zoneID];
	// Hydrogen Transport Fields
	volScalarField& sourceMassO2 = sourceMassO2_CPTL[zoneID];
	// Water Vapour Transport Fields
	volScalarField& sourceMassH2OVap = sourceMassH2OVap_CPTL[zoneID];
	// Energy Transport Fields 
	volScalarField& sourceEnergyTotal = sourceEnergyTotal_CPTL[zoneID];
  	// Mass Transport (Pressure) Fields
	volScalarField& sourceMassTotal = sourceMassTotal_CPTL[zoneID];
 	// Liquid Water Transport Fields
	volScalarField& sourceMassH2OLiq = sourceMassH2OLiq_CPTL[zoneID];
