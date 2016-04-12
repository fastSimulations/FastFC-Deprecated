// Initilialise field pointer lists
	// Calculated Field Variables
		// Stored Gradients
		PtrList<volVectorField> pMixGrad_ACL(ACLRegions.size());
		PtrList<volVectorField> potElectronGrad_ACL(ACLRegions.size());
		PtrList<volVectorField> potProtonGrad_ACL(ACLRegions.size());
		// Morphology Fields
		PtrList<volScalarField> massSupp_ACL(ACLRegions.size());
		PtrList<volScalarField> massIon_ACL(ACLRegions.size());
		PtrList<volScalarField> sAC_ACL(ACLRegions.size());
		PtrList<volScalarField> sVAgg_ACL(ACLRegions.size());
		PtrList<volScalarField> ratioPtSupp_ACL(ACLRegions.size());
		PtrList<volScalarField> ratioCatIon_ACL(ACLRegions.size());
		PtrList<volScalarField> volFSupp_ACL(ACLRegions.size());
		PtrList<volScalarField> volFPt_ACL(ACLRegions.size());
		PtrList<volScalarField> volFIon_ACL(ACLRegions.size());
		PtrList<volScalarField> volFPore_ACL(ACLRegions.size());
		PtrList<volScalarField> volFPoreEff_ACL(ACLRegions.size());
		PtrList<volScalarField> tIonFilm_ACL(ACLRegions.size());
		PtrList<volScalarField> numAgg_ACL(ACLRegions.size());
		PtrList<volScalarField> volFSuppAgg_ACL(ACLRegions.size());
		PtrList<volScalarField> volFPtAgg_ACL(ACLRegions.size());
		PtrList<volScalarField> volFIonAgg_ACL(ACLRegions.size());
		// Electrochemistry Fields
		PtrList<volScalarField> chTransCoeffAnodic_ACL(ACLRegions.size());
		PtrList<volScalarField> chTransCoeffCathodic_ACL(ACLRegions.size());
		PtrList<volScalarField> theileModulus_ACL(ACLRegions.size());
		PtrList<volScalarField> effectiveFactor_ACL(ACLRegions.size());
		PtrList<volScalarField> kCAgg_ACL(ACLRegions.size());
		PtrList<volScalarField> exchCurrRefCorr_ACL(ACLRegions.size());
		// Material Property Fields
		PtrList<volScalarField> kPermDry_ACL(ACLRegions.size());
		PtrList<volScalarField> kPermRelLiq_ACL(ACLRegions.size());
		PtrList<volScalarField> kPermRelGas_ACL(ACLRegions.size());
		PtrList<volScalarField> henrySolH2_ACL(ACLRegions.size());
		PtrList<volScalarField> diffH2Ion_ACL(ACLRegions.size());
		PtrList<volScalarField> diffH2IonAgg_ACL(ACLRegions.size());
		// Nitrogen Transport Fields
		PtrList<volScalarField> wN2_ACL(ACLRegions.size());
		// Mole Fraction Fields
		PtrList<volScalarField> xH2_ACL(ACLRegions.size());
		PtrList<volScalarField> xN2_ACL(ACLRegions.size());
		PtrList<volScalarField> xH2OVap_ACL(ACLRegions.size());
		// Fluid Property Fields
		PtrList<volScalarField> viscH2_ACL(ACLRegions.size());
		PtrList<volScalarField> viscN2_ACL(ACLRegions.size());
		PtrList<volScalarField> viscH2OVap_ACL(ACLRegions.size());
		// Binary Diffusivity Fields
		PtrList<volScalarField> diffH2N2_ACL(ACLRegions.size());
		PtrList<volScalarField> diffH2H2O_ACL(ACLRegions.size());
		PtrList<volScalarField> diffN2H2O_ACL(ACLRegions.size());
		// Mixture Property Fields
		PtrList<volScalarField> molWeightMix_ACL(ACLRegions.size());
		PtrList<volScalarField> densityMix_ACL(ACLRegions.size());
		PtrList<volScalarField> viscMix_ACL(ACLRegions.size());
		PtrList<volScalarField> velAvg_ACL(ACLRegions.size());
		// Local Condition Fields
		PtrList<volScalarField> RH_ACL(ACLRegions.size());
		PtrList<volScalarField> pH2_ACL(ACLRegions.size());
		PtrList<volScalarField> pH2OVap_ACL(ACLRegions.size());
		PtrList<volScalarField> pH2OSat_ACL(ACLRegions.size());
		PtrList<volScalarField> actPot_ACL(ACLRegions.size());
		// Ionomeric Material Fields
		PtrList<volScalarField> H2ODissIon_ACL(ACLRegions.size());
		PtrList<volScalarField> condProtonIon_ACL(ACLRegions.size());
		// Liquid Water Fields
		PtrList<volScalarField> jS_ACL(ACLRegions.size());
		PtrList<volScalarField> pCap_ACL(ACLRegions.size());
		PtrList<volScalarField> djSDs_ACL(ACLRegions.size());
		PtrList<volScalarField> dpCapDs_ACL(ACLRegions.size());
		PtrList<volScalarField> swPhaseChangeLiqVap_ACL(ACLRegions.size());
		// Source Term Fields
		PtrList<volScalarField> sourceEnergyPhaseChange_ACL(ACLRegions.size());
		PtrList<volScalarField>	sourceEnergyEntropicHeating_ACL(ACLRegions.size());
		PtrList<volScalarField> sourceEnergyOhmicHeating_ACL(ACLRegions.size());
		PtrList<volScalarField> sourcePhaseChangeLiqVap_ACL(ACLRegions.size());
		PtrList<volScalarField> reactionRateAnodic_ACL(ACLRegions.size());
		PtrList<volScalarField> reactionRateCathodic_ACL(ACLRegions.size());
		PtrList<volScalarField> sourceVolCurrent_ACL(ACLRegions.size());

	// Solved Field Variables	
		// Electron Transport Fields
		PtrList<volScalarField> potElectron_ACL(ACLRegions.size());
		// Proton Transport Fields
		PtrList<volScalarField> potProton_ACL(ACLRegions.size());
		// Hydrogen Transport Fields
	   	PtrList<volScalarField> wH2_ACL(ACLRegions.size());
		// Water Vapour Transport Fields
		PtrList<volScalarField> wH2OVap_ACL(ACLRegions.size());
		// Energy Transport Fields
		PtrList<volScalarField> T_ACL(ACLRegions.size());
		// Mass Transport (Pressure) Fields
	   	PtrList<volScalarField> pMix_ACL(ACLRegions.size());
	   	// Liquid Water Transport Fields
		PtrList<volScalarField> sLiq_ACL(ACLRegions.size());
	
	// Coefficient Field Variables
		// Electron Transport Fields
		PtrList<volScalarField> condElectronEff_ACL(ACLRegions.size());
		// Proton Transport Fields
		PtrList<volScalarField> condProtonEff_ACL(ACLRegions.size());
		// Hydrogen Transport Fields
		PtrList<volScalarField> diffH2MixEff_ACL(ACLRegions.size());
		PtrList<volScalarField> transportCoeffH2_ACL(ACLRegions.size());
		// Water Vapour Transport Fields
		PtrList<volScalarField> diffH2OVapMixEff_ACL(ACLRegions.size());
		PtrList<volScalarField> transportCoeffH2OVap_ACL(ACLRegions.size());
		// Energy Transport Fields
		PtrList<volScalarField> condThermalEff_ACL(ACLRegions.size());
		PtrList<volScalarField> transportCoeffEnergy_ACL(ACLRegions.size());
	  	// Mass Transport (Pressure) Fields
		PtrList<volScalarField> dynamicPerm_ACL(ACLRegions.size());
	 	// Liquid Water Transport Fields
		PtrList<volScalarField> diffH2OLiq_ACL(ACLRegions.size());

	// Source/sink Field Variables
		// Electron Transport Fields
		PtrList<volScalarField> sourceElectronTotal_ACL(ACLRegions.size());
		// Proton Transport Fields
		PtrList<volScalarField> sourceProtonTotal_ACL(ACLRegions.size());
		// Hydrogen Transport Fields
		PtrList<volScalarField> sourceMassH2_ACL(ACLRegions.size());
		// Water Vapour Transport Fields
		PtrList<volScalarField> sourceMassH2OVap_ACL(ACLRegions.size());
		// Energy Transport Fields 
		PtrList<volScalarField> sourceEnergyTotal_ACL(ACLRegions.size());
	  	// Mass Transport (Pressure) Fields
		PtrList<volScalarField> sourceMassTotal_ACL(ACLRegions.size());
	 	// Liquid Water Transport Fields
		PtrList<volScalarField> sourceMassH2OLiq_ACL(ACLRegions.size());
