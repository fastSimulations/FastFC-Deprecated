// Create assignment to the mesh object
	fvMesh& mesh = ACLRegions[zoneID];
	labelList nCellsI(mesh.C().size());

// Read in the transport Properties dictionary if it has been modified since last read
#include <ACLReadComponentDictionary.H>

// Assignment the value in the dictionary to the electronic conductivity
#include <ACLDimScalarAssign.H>

	// Create Assignments to Calculated Field Variables
		// Stored Gradients
		volVectorField& pMixGrad = pMixGrad_ACL[zoneID];
		volVectorField& potElectronGrad = potElectronGrad_ACL[zoneID];
		volVectorField& potProtonGrad = potProtonGrad_ACL[zoneID];
		// Morphology Fields
		volScalarField& massSupp = massSupp_ACL[zoneID];
		volScalarField& massIon = massIon_ACL[zoneID];
		volScalarField& sAC = sAC_ACL[zoneID];
		volScalarField& sVAgg = sVAgg_ACL[zoneID];
		volScalarField& ratioPtSupp = ratioPtSupp_ACL[zoneID];
		volScalarField& ratioCatIon = ratioCatIon_ACL[zoneID];
		volScalarField& volFSupp = volFSupp_ACL[zoneID];
		volScalarField& volFPt = volFPt_ACL[zoneID];
		volScalarField& volFIon = volFIon_ACL[zoneID];
		volScalarField& volFPore = volFPore_ACL[zoneID];
		volScalarField& volFPoreEff = volFPoreEff_ACL[zoneID];
		volScalarField& tIonFilm = tIonFilm_ACL[zoneID];
		volScalarField& numAgg = numAgg_ACL[zoneID];
		volScalarField& volFSuppAgg = volFSuppAgg_ACL[zoneID];
		volScalarField& volFPtAgg = volFPtAgg_ACL[zoneID];
		volScalarField& volFIonAgg = volFIonAgg_ACL[zoneID];

		// Material Property Fields
		volScalarField& kPermDry = kPermDry_ACL[zoneID];
		volScalarField& kPermRelLiq = kPermRelLiq_ACL[zoneID];
		volScalarField& kPermRelGas = kPermRelGas_ACL[zoneID];
		volScalarField& henrySolH2 = henrySolH2_ACL[zoneID];
		volScalarField& diffH2Ion = diffH2Ion_ACL[zoneID];
		volScalarField& diffH2IonAgg = diffH2IonAgg_ACL[zoneID];

		// Nitrogen Transport Fields
		volScalarField& wN2 = wN2_ACL[zoneID];

		// Mole Fraction Fields
		volScalarField& xH2 = xH2_ACL[zoneID];
		volScalarField& xN2 = xN2_ACL[zoneID];
		volScalarField& xH2OVap = xH2OVap_ACL[zoneID];

		// Fluid Property Fields
		volScalarField& viscH2 = viscH2_ACL[zoneID];
		volScalarField& viscN2 = viscN2_ACL[zoneID];
		volScalarField& viscH2OVap = viscH2OVap_ACL[zoneID];

		// Binary Diffusivity Fields
		volScalarField& diffH2N2 = diffH2N2_ACL[zoneID];
		volScalarField& diffH2H2O = diffH2H2O_ACL[zoneID];
		volScalarField& diffN2H2O = diffN2H2O_ACL[zoneID];

		// Mixture Property Fields
		volScalarField& molWeightMix = molWeightMix_ACL[zoneID];
		volScalarField& densityMix = densityMix_ACL[zoneID];
		volScalarField& viscMix = viscMix_ACL[zoneID];
		volScalarField& velAvg = velAvg_ACL[zoneID];

		// Local Condition Fields
		volScalarField& RH = RH_ACL[zoneID];
		volScalarField& pH2 = pH2_ACL[zoneID];
		volScalarField& pH2OVap = pH2OVap_ACL[zoneID];
		volScalarField& pH2OSat = pH2OSat_ACL[zoneID];

		// Ionomeric Material Fields
		volScalarField& H2ODissIon = H2ODissIon_ACL[zoneID];
		volScalarField& condProtonIon = condProtonIon_ACL[zoneID];

		// Liquid Water Fields
		volScalarField& jS = jS_ACL[zoneID];
		volScalarField& pCap = pCap_ACL[zoneID];
		volScalarField& djSDs = djSDs_ACL[zoneID];
		volScalarField& dpCapDs = dpCapDs_ACL[zoneID];
		volScalarField& swPhaseChangeLiqVap = swPhaseChangeLiqVap_ACL[zoneID];

		// Electrochemistry Fields
		volScalarField& chTransCoeffAnodic = chTransCoeffAnodic_ACL[zoneID];
		volScalarField& chTransCoeffCathodic = chTransCoeffCathodic_ACL[zoneID];
		volScalarField& theileModulus = theileModulus_ACL[zoneID];
		volScalarField& effectiveFactor = effectiveFactor_ACL[zoneID];
		volScalarField& kCAgg = kCAgg_ACL[zoneID];
		volScalarField& exchCurrRefCorr = exchCurrRefCorr_ACL[zoneID];
		volScalarField& actPot = actPot_ACL[zoneID];


	// Solved Field Variables	
		// Electron Transport Fields
		volScalarField& potElectron = potElectron_ACL[zoneID];
		// Proton Transport Fields
		volScalarField& potProton = potProton_ACL[zoneID];
		// Hydrogen Transport Fields
	   	volScalarField& wH2 = wH2_ACL[zoneID];
		// Water Vapour Transport Fields
		volScalarField& wH2OVap = wH2OVap_ACL[zoneID];
		// Energy Transport Fields
		volScalarField& T = T_ACL[zoneID];
		// Mass Transport (Pressure) Fields
	   	volScalarField& pMix = pMix_ACL[zoneID];
	   	// Liquid Water Transport Fields
		volScalarField& sLiq = sLiq_ACL[zoneID];
	
	// Coefficient Field Variables
		// Electron Transport Fields
		volScalarField& condElectronEff = condElectronEff_ACL[zoneID];
		// Proton Transport Fields
		volScalarField& condProtonEff = condProtonEff_ACL[zoneID];
		// Hydrogen Transport Fields
		volScalarField& diffH2MixEff = diffH2MixEff_ACL[zoneID];
		volScalarField& transportCoeffH2 = transportCoeffH2_ACL[zoneID];
		// Water Vapour Transport Fields
		volScalarField& diffH2OVapMixEff = diffH2OVapMixEff_ACL[zoneID];
		volScalarField& transportCoeffH2OVap = transportCoeffH2OVap_ACL[zoneID];
		// Energy Transport Fields
		volScalarField& condThermalEff = condThermalEff_ACL[zoneID];
		volScalarField& transportCoeffEnergy = transportCoeffEnergy_ACL[zoneID];
	  	// Mass Transport (Pressure) Fields
		volScalarField& dynamicPerm = dynamicPerm_ACL[zoneID];
	 	// Liquid Water Transport Fields
		volScalarField& diffH2OLiq = diffH2OLiq_ACL[zoneID];

	// Source/sink Field Variables
		// Intermediate Source Term Fields
		volScalarField& sourceEnergyPhaseChange = sourceEnergyPhaseChange_ACL[zoneID];
		volScalarField& sourceEnergyEntropicHeating = sourceEnergyEntropicHeating_ACL[zoneID];
		volScalarField& sourceEnergyOhmicHeating = sourceEnergyOhmicHeating_ACL[zoneID];
		volScalarField& sourcePhaseChangeLiqVap = sourcePhaseChangeLiqVap_ACL[zoneID];
		volScalarField& reactionRateAnodic = reactionRateAnodic_ACL[zoneID];
		volScalarField& reactionRateCathodic = reactionRateCathodic_ACL[zoneID];
		volScalarField& sourceVolCurrent = sourceVolCurrent_ACL[zoneID];
		// Electron Transport Fields
		volScalarField& sourceElectronTotal = sourceElectronTotal_ACL[zoneID];
		// Proton Transport Fields
		volScalarField& sourceProtonTotal = sourceProtonTotal_ACL[zoneID];
		// Hydrogen Transport Fields
		volScalarField& sourceMassH2 = sourceMassH2_ACL[zoneID];
		// Water Vapour Transport Fields
		volScalarField& sourceMassH2OVap = sourceMassH2OVap_ACL[zoneID];
		// Energy Transport Fields 
		volScalarField& sourceEnergyTotal = sourceEnergyTotal_ACL[zoneID];
	  	// Mass Transport (Pressure) Fields
		volScalarField& sourceMassTotal = sourceMassTotal_ACL[zoneID];
	 	// Liquid Water Transport Fields
		volScalarField& sourceMassH2OLiq = sourceMassH2OLiq_ACL[zoneID];
