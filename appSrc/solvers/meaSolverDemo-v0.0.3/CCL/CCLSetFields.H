// Create assignment to the mesh object
	fvMesh& mesh = CCLRegions[zoneID];

//	labelList NCellsI(mesh.cellZones());

	labelList NCellsI(mesh.C().size());

// Read in the transport Properties dictionary if it has been modified since last read
#include <CCLReadComponentDictionary.H>

// Assignment the value in the dictionary to the electronic conductivity
#include <CCLDimScalarAssign.H>

	// Create Assignments to Calculated Field Variables
		// Stored Gradients
		volVectorField& pMixGrad = pMixGrad_CCL[zoneID];
		volVectorField& potElectronGrad = potElectronGrad_CCL[zoneID];
		volVectorField& potProtonGrad = potProtonGrad_CCL[zoneID];
		// Morphology Fields
		volScalarField& massSupp = massSupp_CCL[zoneID];
		volScalarField& massIon = massIon_CCL[zoneID];
		volScalarField& sAC = sAC_CCL[zoneID];
		volScalarField& sVAgg = sVAgg_CCL[zoneID];
		volScalarField& ratioPtSupp = ratioPtSupp_CCL[zoneID];
		volScalarField& ratioCatIon = ratioCatIon_CCL[zoneID];
		volScalarField& volFSupp = volFSupp_CCL[zoneID];
		volScalarField& volFPt = volFPt_CCL[zoneID];
		volScalarField& volFIon = volFIon_CCL[zoneID];
		volScalarField& volFPore = volFPore_CCL[zoneID];
		volScalarField& volFPoreEff = volFPoreEff_CCL[zoneID];
		volScalarField& tIonFilm = tIonFilm_CCL[zoneID];
		volScalarField& numAgg = numAgg_CCL[zoneID];
		volScalarField& volFSuppAgg = volFSuppAgg_CCL[zoneID];
		volScalarField& volFPtAgg = volFPtAgg_CCL[zoneID];
		volScalarField& volFIonAgg = volFIonAgg_CCL[zoneID];
		// Electrochemistry Fields
		volScalarField& chTransCoeffAnodic = chTransCoeffAnodic_CCL[zoneID];
		volScalarField& chTransCoeffCathodic = chTransCoeffCathodic_CCL[zoneID];
		volScalarField& theileModulus = theileModulus_CCL[zoneID];
		volScalarField& effectiveFactor = effectiveFactor_CCL[zoneID];
		volScalarField& kCAgg = kCAgg_CCL[zoneID];
		volScalarField& exchCurrRefCorr = exchCurrRefCorr_CCL[zoneID];
		// Material Property Fields
		volScalarField& kPermDry = kPermDry_CCL[zoneID];
		volScalarField& kPermRelLiq = kPermRelLiq_CCL[zoneID];
		volScalarField& kPermRelGas = kPermRelGas_CCL[zoneID];
		volScalarField& henrySolO2 = henrySolO2_CCL[zoneID];
		volScalarField& diffO2Ion = diffO2Ion_CCL[zoneID];
		volScalarField& diffO2IonAgg = diffO2IonAgg_CCL[zoneID];
		// Nitrogen Transport Fields
		volScalarField& wN2 = wN2_CCL[zoneID];
		// Mole Fraction Fields
		volScalarField& xO2 = xO2_CCL[zoneID];
		volScalarField& xN2 = xN2_CCL[zoneID];
		volScalarField& xH2OVap = xH2OVap_CCL[zoneID];
		// Fluid Property Fields
		volScalarField& viscO2 = viscO2_CCL[zoneID];
		volScalarField& viscN2 = viscN2_CCL[zoneID];
		volScalarField& viscH2OVap = viscH2OVap_CCL[zoneID];
		// Binary Diffusivity Fields
		volScalarField& diffO2N2 = diffO2N2_CCL[zoneID];
		volScalarField& diffO2H2O = diffO2H2O_CCL[zoneID];
		volScalarField& diffN2H2O = diffN2H2O_CCL[zoneID];
		// Mixture Property Fields
		volScalarField& molWeightMix = molWeightMix_CCL[zoneID];
		volScalarField& densityMix = densityMix_CCL[zoneID];
		volScalarField& viscMix = viscMix_CCL[zoneID];
		volScalarField& velAvg = velAvg_CCL[zoneID];
		// Local Condition Fields
		volScalarField& RH = RH_CCL[zoneID];
		volScalarField& pO2 = pO2_CCL[zoneID];
		volScalarField& pH2OVap = pH2OVap_CCL[zoneID];
		volScalarField& pH2OSat = pH2OSat_CCL[zoneID];
		volScalarField& actPot = actPot_CCL[zoneID];
		// Ionomeric Material Fields
		volScalarField& H2ODissIon = H2ODissIon_CCL[zoneID];
		volScalarField& condProtonIon = condProtonIon_CCL[zoneID];
		// Liquid Water Fields
		volScalarField& jS = jS_CCL[zoneID];
		volScalarField& pCap = pCap_CCL[zoneID];
		volScalarField& djSDs = djSDs_CCL[zoneID];
		volScalarField& dpCapDs = dpCapDs_CCL[zoneID];
		volScalarField& swPhaseChangeLiqVap = swPhaseChangeLiqVap_CCL[zoneID];
		// Intermediate Source Term Fields
		volScalarField& sourceEnergyPhaseChange = sourceEnergyPhaseChange_CCL[zoneID];
		volScalarField& sourceEnergyOhmicHeating = sourceEnergyOhmicHeating_CCL[zoneID];
		volScalarField& sourceEnergyEntropicHeating = sourceEnergyEntropicHeating_CCL[zoneID];
		volScalarField& sourcePhaseChangeLiqVap = sourcePhaseChangeLiqVap_CCL[zoneID];
		volScalarField& reactionRateAnodic = reactionRateAnodic_CCL[zoneID];
		volScalarField& reactionRateCathodic = reactionRateCathodic_CCL[zoneID];
		volScalarField& sourceVolCurrent = sourceVolCurrent_CCL[zoneID];

	// Solved Field Variables	
		// Electron Transport Fields
		volScalarField& potElectron = potElectron_CCL[zoneID];
		// Proton Transport Fields
		volScalarField& potProton = potProton_CCL[zoneID];
		// Oxygen Transport Fields
	   	volScalarField& wO2 = wO2_CCL[zoneID];
		// Water Vapour Transport Fields
		volScalarField& wH2OVap = wH2OVap_CCL[zoneID];
		// Energy Transport Fields
		volScalarField& T = T_CCL[zoneID];
		// Mass Transport (Pressure) Fields
	   	volScalarField& pMix = pMix_CCL[zoneID];
	   	// Liquid Water Transport Fields
		volScalarField& sLiq = sLiq_CCL[zoneID];
	
	// Coefficient Field Variables
		// Electron Transport Fields
		volScalarField& condElectronEff = condElectronEff_CCL[zoneID];
		// Proton Transport Fields
		volScalarField& condProtonEff = condProtonEff_CCL[zoneID];
		// Oxygen Transport Fields
		volScalarField& diffO2MixEff = diffO2MixEff_CCL[zoneID];
		volScalarField& transportCoeffO2 = transportCoeffO2_CCL[zoneID];
		// Water Vapour Transport Fields
		volScalarField& diffH2OVapMixEff = diffH2OVapMixEff_CCL[zoneID];
		volScalarField& transportCoeffH2OVap = transportCoeffH2OVap_CCL[zoneID];
		// Energy Transport Fields
		volScalarField& condThermalEff = condThermalEff_CCL[zoneID];
		volScalarField& transportCoeffEnergy = transportCoeffEnergy_CCL[zoneID];
	  	// Mass Transport (Pressure) Fields
		volScalarField& dynamicPerm = dynamicPerm_CCL[zoneID];
	 	// Liquid Water Transport Fields
		volScalarField& diffH2OLiq = diffH2OLiq_CCL[zoneID];

	// Source/sink Field Variables
		// Electron Transport Fields
		volScalarField& sourceElectronTotal = sourceElectronTotal_CCL[zoneID];
		// Proton Transport Fields
		volScalarField& sourceProtonTotal = sourceProtonTotal_CCL[zoneID];
		// Oxygen Transport Fields
		volScalarField& sourceMassO2 = sourceMassO2_CCL[zoneID];
		// Water Vapour Transport Fields
		volScalarField& sourceMassH2OVap = sourceMassH2OVap_CCL[zoneID];
		// Energy Transport Fields 
		volScalarField& sourceEnergyTotal = sourceEnergyTotal_CCL[zoneID];
	  	// Mass Transport (Pressure) Fields
		volScalarField& sourceMassTotal = sourceMassTotal_CCL[zoneID];
	 	// Liquid Water Transport Fields
		volScalarField& sourceMassH2OLiq = sourceMassH2OLiq_CCL[zoneID];
