// Create assignment to the mesh object
	const fvMesh& mesh = APTLRegions[zoneID];

//	Setup the variable for the cell count
	labelList NCellsI(mesh.C().size());

// Read in the transport Properties dictionary if it has been modified since last read
#include <APTLReadComponentDictionary.H>

// Assignment the value in the dictionary to the electronic conductivity
#include <APTLDimScalarAssign.H>

// Create Assignments to the Calculated Field Variables
	// Stored Gradients
	volVectorField& pMixGrad = pMixGrad_APTL[zoneID];
	volVectorField& potElectronGrad = potElectronGrad_APTL[zoneID];
	// Morphology Fields
	volScalarField& volFPore = volFPore_APTL[zoneID];
	volScalarField& volFPoreEff = volFPoreEff_APTL[zoneID];
	// Material Property Fields
	volScalarField& kPermDry = kPermDry_APTL[zoneID];
	volScalarField& kPermRelLiq = kPermRelLiq_APTL[zoneID];
	volScalarField& kPermRelGas = kPermRelGas_APTL[zoneID];
	// Nitrogen Transport Fields
	volScalarField& wN2 = wN2_APTL[zoneID];
	// Mole Fraction Fields
	volScalarField& xH2 = xH2_APTL[zoneID];
	volScalarField& xN2 = xN2_APTL[zoneID];
	volScalarField& xH2OVap = xH2OVap_APTL[zoneID];
	// Fluid Property Fields
	volScalarField& viscH2 = viscH2_APTL[zoneID];
	volScalarField& viscN2 = viscN2_APTL[zoneID];
	volScalarField& viscH2OVap = viscH2OVap_APTL[zoneID];
	// Binary Diffusivity Fields
	volScalarField& diffH2N2 = diffH2N2_APTL[zoneID];
	volScalarField& diffH2H2O = diffH2H2O_APTL[zoneID];
	volScalarField& diffN2H2O = diffN2H2O_APTL[zoneID];
	// Mixture Property Fields
	volScalarField& molWeightMix = molWeightMix_APTL[zoneID];
	volScalarField& densityMix = densityMix_APTL[zoneID];
	volScalarField& viscMix = viscMix_APTL[zoneID];
	volScalarField& velAvg = velAvg_APTL[zoneID];
	// Local Condition Fields
	volScalarField& RH = RH_APTL[zoneID];
	volScalarField& pH2 = pH2_APTL[zoneID];
	volScalarField& pH2OVap = pH2OVap_APTL[zoneID];
	volScalarField& pH2OSat = pH2OSat_APTL[zoneID];
	// Liquid Water Fields
	volScalarField& jS = jS_APTL[zoneID];
	volScalarField& pCap = pCap_APTL[zoneID];
	volScalarField& djSDs = djSDs_APTL[zoneID];
	volScalarField& dpCapDs = dpCapDs_APTL[zoneID];
	volScalarField& swPhaseChangeLiqVap = swPhaseChangeLiqVap_APTL[zoneID];
	// Source Term Fields
	volScalarField& sourceEnergyPhaseChange = sourceEnergyPhaseChange_APTL[zoneID];
	volScalarField& sourceEnergyOhmicHeating = sourceEnergyOhmicHeating_APTL[zoneID];
	volScalarField& sourcePhaseChangeLiqVap = sourcePhaseChangeLiqVap_APTL[zoneID];
	
// Create Assignments to Solved Field Variables
	// Electron Transport Fields
	volScalarField& potElectron = potElectron_APTL[zoneID];
	// Hydrogen Transport Fields
   	volScalarField& wH2 = wH2_APTL[zoneID];
	// Water Vapour Transport Fields
	volScalarField& wH2OVap = wH2OVap_APTL[zoneID];
	// Energy Transport Fields
	volScalarField& T = T_APTL[zoneID];
	// Mass Transport (Pressure) Fields
   	volScalarField& pMix = pMix_APTL[zoneID];
   	// Liquid Water Transport Fields
	volScalarField& sLiq = sLiq_APTL[zoneID];

// Create Assignments to Coefficient Field Variables
	// Electron Transport Fields
	volScalarField& condElectronEff = condElectronEff_APTL[zoneID];
	// Hydrogen Transport Fields
	volScalarField& diffH2MixEff = diffH2MixEff_APTL[zoneID];
	volScalarField& transportCoeffH2 = transportCoeffH2_APTL[zoneID];
	// Water Vapour Transport Fields
	volScalarField& diffH2OVapMixEff = diffH2OVapMixEff_APTL[zoneID];
	volScalarField& transportCoeffH2OVap = transportCoeffH2OVap_APTL[zoneID];
	// Energy Transport Fields
	volScalarField& condThermalEff = condThermalEff_APTL[zoneID];
	volScalarField& transportCoeffEnergy = transportCoeffEnergy_APTL[zoneID];
  	// Mass Transport (Pressure) Fields
	volScalarField& dynamicPerm = dynamicPerm_APTL[zoneID];
 	// Liquid Water Transport Fields
	volScalarField& diffH2OLiq = diffH2OLiq_APTL[zoneID];

// Create Assignments to Source/Sink Field Variables
	// Electron Transport Fields
	volScalarField& sourceElectronTotal = sourceElectronTotal_APTL[zoneID];
	// Hydrogen Transport Fields
	volScalarField& sourceMassH2 = sourceMassH2_APTL[zoneID];
	// Water Vapour Transport Fields
	volScalarField& sourceMassH2OVap = sourceMassH2OVap_APTL[zoneID];
	// Energy Transport Fields 
	volScalarField& sourceEnergyTotal = sourceEnergyTotal_APTL[zoneID];
  	// Mass Transport (Pressure) Fields
	volScalarField& sourceMassTotal = sourceMassTotal_APTL[zoneID];
 	// Liquid Water Transport Fields
	volScalarField& sourceMassH2OLiq = sourceMassH2OLiq_APTL[zoneID];
