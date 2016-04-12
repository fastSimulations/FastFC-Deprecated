// Initilialise field pointer lists
	// Calculated Field Variables
		// Stored Gradients
		PtrList<volVectorField> pMixGrad_CCL(CCLRegions.size());
		PtrList<volVectorField> potElectronGrad_CCL(CCLRegions.size());
		PtrList<volVectorField> potProtonGrad_CCL(CCLRegions.size());
		// Morphology Fields
		PtrList<volScalarField> massSupp_CCL(CCLRegions.size());
		PtrList<volScalarField> massIon_CCL(CCLRegions.size());
		PtrList<volScalarField> sAC_CCL(CCLRegions.size());
		PtrList<volScalarField> sVAgg_CCL(CCLRegions.size());
		PtrList<volScalarField> ratioPtSupp_CCL(CCLRegions.size());
		PtrList<volScalarField> ratioCatIon_CCL(CCLRegions.size());
		PtrList<volScalarField> volFSupp_CCL(CCLRegions.size());
		PtrList<volScalarField> volFPt_CCL(CCLRegions.size());
		PtrList<volScalarField> volFIon_CCL(CCLRegions.size());
		PtrList<volScalarField> volFPore_CCL(CCLRegions.size());
		PtrList<volScalarField> volFPoreEff_CCL(CCLRegions.size());
		PtrList<volScalarField> tIonFilm_CCL(CCLRegions.size());
		PtrList<volScalarField> numAgg_CCL(CCLRegions.size());
		PtrList<volScalarField> volFSuppAgg_CCL(CCLRegions.size());
		PtrList<volScalarField> volFPtAgg_CCL(CCLRegions.size());
		PtrList<volScalarField> volFIonAgg_CCL(CCLRegions.size());
		// Electrochemistry Fields
		PtrList<volScalarField> chTransCoeffAnodic_CCL(CCLRegions.size());
		PtrList<volScalarField> chTransCoeffCathodic_CCL(CCLRegions.size());
		PtrList<volScalarField> theileModulus_CCL(CCLRegions.size());
		PtrList<volScalarField> effectiveFactor_CCL(CCLRegions.size());
		PtrList<volScalarField> kCAgg_CCL(CCLRegions.size());
		PtrList<volScalarField> exchCurrRefCorr_CCL(CCLRegions.size());
		// Material Property Fields
		PtrList<volScalarField> kPermDry_CCL(CCLRegions.size());
		PtrList<volScalarField> kPermRelLiq_CCL(CCLRegions.size());
		PtrList<volScalarField> kPermRelGas_CCL(CCLRegions.size());
		PtrList<volScalarField> henrySolO2_CCL(CCLRegions.size());
		PtrList<volScalarField> diffO2Ion_CCL(CCLRegions.size());
		PtrList<volScalarField> diffO2IonAgg_CCL(CCLRegions.size());
		// Nitrogen Transport Fields
		PtrList<volScalarField> wN2_CCL(CCLRegions.size());
		// Mole Fraction Fields
		PtrList<volScalarField> xO2_CCL(CCLRegions.size());
		PtrList<volScalarField> xN2_CCL(CCLRegions.size());
		PtrList<volScalarField> xH2OVap_CCL(CCLRegions.size());
		// Fluid Property Fields
		PtrList<volScalarField> viscO2_CCL(CCLRegions.size());
		PtrList<volScalarField> viscN2_CCL(CCLRegions.size());
		PtrList<volScalarField> viscH2OVap_CCL(CCLRegions.size());
		// Binary Diffusivity Fields
		PtrList<volScalarField> diffO2N2_CCL(CCLRegions.size());
		PtrList<volScalarField> diffO2H2O_CCL(CCLRegions.size());
		PtrList<volScalarField> diffN2H2O_CCL(CCLRegions.size());
		// Mixture Property Fields
		PtrList<volScalarField> molWeightMix_CCL(CCLRegions.size());
		PtrList<volScalarField> densityMix_CCL(CCLRegions.size());
		PtrList<volScalarField> viscMix_CCL(CCLRegions.size());
		PtrList<volScalarField> velAvg_CCL(CCLRegions.size());
		// Local Condition Fields
		PtrList<volScalarField> RH_CCL(CCLRegions.size());
		PtrList<volScalarField> pO2_CCL(CCLRegions.size());
		PtrList<volScalarField> pH2OVap_CCL(CCLRegions.size());
		PtrList<volScalarField> pH2OSat_CCL(CCLRegions.size());
		PtrList<volScalarField> actPot_CCL(CCLRegions.size());
		// Ionomeric Material Fields
		PtrList<volScalarField> H2ODissIon_CCL(CCLRegions.size());
		PtrList<volScalarField> condProtonIon_CCL(CCLRegions.size());
		// Liquid Water Fields
		PtrList<volScalarField> jS_CCL(CCLRegions.size());
		PtrList<volScalarField> pCap_CCL(CCLRegions.size());
		PtrList<volScalarField> djSDs_CCL(CCLRegions.size());
		PtrList<volScalarField> dpCapDs_CCL(CCLRegions.size());
		PtrList<volScalarField> swPhaseChangeLiqVap_CCL(CCLRegions.size());
		// Source Term Fields
		PtrList<volScalarField> sourceEnergyPhaseChange_CCL(CCLRegions.size());
		PtrList<volScalarField>	sourceEnergyEntropicHeating_CCL(CCLRegions.size());
		PtrList<volScalarField> sourceEnergyOhmicHeating_CCL(CCLRegions.size());
		PtrList<volScalarField> sourcePhaseChangeLiqVap_CCL(CCLRegions.size());
		PtrList<volScalarField> reactionRateAnodic_CCL(CCLRegions.size());
		PtrList<volScalarField> reactionRateCathodic_CCL(CCLRegions.size());
		PtrList<volScalarField> sourceVolCurrent_CCL(CCLRegions.size());

	// Solved Field Variables	
		// Electron Transport Fields
		PtrList<volScalarField> potElectron_CCL(CCLRegions.size());
		// Proton Transport Fields
		PtrList<volScalarField> potProton_CCL(CCLRegions.size());
		// Oxygen Transport Fields
	   	PtrList<volScalarField> wO2_CCL(CCLRegions.size());
		// Water Vapour Transport Fields
		PtrList<volScalarField> wH2OVap_CCL(CCLRegions.size());
		// Energy Transport Fields
		PtrList<volScalarField> T_CCL(CCLRegions.size());
		// Mass Transport (Pressure) Fields
	   	PtrList<volScalarField> pMix_CCL(CCLRegions.size());
	   	// Liquid Water Transport Fields
		PtrList<volScalarField> sLiq_CCL(CCLRegions.size());
	
	// Coefficient Field Variables
		// Electron Transport Fields
		PtrList<volScalarField> condElectronEff_CCL(CCLRegions.size());
		// Proton Transport Fields
		PtrList<volScalarField> condProtonEff_CCL(CCLRegions.size());
		// Oxygen Transport Fields
		PtrList<volScalarField> diffO2MixEff_CCL(CCLRegions.size());
		PtrList<volScalarField> transportCoeffO2_CCL(CCLRegions.size());
		// Water Vapour Transport Fields
		PtrList<volScalarField> diffH2OVapMixEff_CCL(CCLRegions.size());
		PtrList<volScalarField> transportCoeffH2OVap_CCL(CCLRegions.size());
		// Energy Transport Fields
		PtrList<volScalarField> condThermalEff_CCL(CCLRegions.size());
		PtrList<volScalarField> transportCoeffEnergy_CCL(CCLRegions.size());
	  	// Mass Transport (Pressure) Fields
		PtrList<volScalarField> dynamicPerm_CCL(CCLRegions.size());
	 	// Liquid Water Transport Fields
		PtrList<volScalarField> diffH2OLiq_CCL(CCLRegions.size());

	// Source/sink Field Variables
		// Electron Transport Fields
		PtrList<volScalarField> sourceElectronTotal_CCL(CCLRegions.size());
		// Proton Transport Fields
		PtrList<volScalarField> sourceProtonTotal_CCL(CCLRegions.size());
		// Oxygen Transport Fields
		PtrList<volScalarField> sourceMassO2_CCL(CCLRegions.size());
		// Water Vapour Transport Fields
		PtrList<volScalarField> sourceMassH2OVap_CCL(CCLRegions.size());
		// Energy Transport Fields 
		PtrList<volScalarField> sourceEnergyTotal_CCL(CCLRegions.size());
	  	// Mass Transport (Pressure) Fields
		PtrList<volScalarField> sourceMassTotal_CCL(CCLRegions.size());
	 	// Liquid Water Transport Fields
		PtrList<volScalarField> sourceMassH2OLiq_CCL(CCLRegions.size());
