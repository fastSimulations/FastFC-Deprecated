// Initilialise field pointer lists
	// Calculated Field Variables
		// Stored Gradients
		PtrList<volVectorField> pMixGrad_APTL(APTLRegions.size());
		PtrList<volVectorField> potElectronGrad_APTL(APTLRegions.size());
		// Morphology Fields
		PtrList<volScalarField> volFPore_APTL(APTLRegions.size());
		PtrList<volScalarField> volFPoreEff_APTL(APTLRegions.size());
		// Material Property Fields
		PtrList<volScalarField> kPermDry_APTL(APTLRegions.size());
		PtrList<volScalarField> kPermRelLiq_APTL(APTLRegions.size());
		PtrList<volScalarField> kPermRelGas_APTL(APTLRegions.size());
		// Nitrogen Transport Fields
		PtrList<volScalarField> wN2_APTL(APTLRegions.size());
		// Mole Fraction Fields
		PtrList<volScalarField> xH2_APTL(APTLRegions.size());
		PtrList<volScalarField> xN2_APTL(APTLRegions.size());
		PtrList<volScalarField> xH2OVap_APTL(APTLRegions.size());
		// Fluid Property Fields
		PtrList<volScalarField> viscH2_APTL(APTLRegions.size());
		PtrList<volScalarField> viscN2_APTL(APTLRegions.size());
		PtrList<volScalarField> viscH2OVap_APTL(APTLRegions.size());
		// Binary Diffusivity Fields
		PtrList<volScalarField> diffH2N2_APTL(APTLRegions.size());
		PtrList<volScalarField> diffH2H2O_APTL(APTLRegions.size());
		PtrList<volScalarField> diffN2H2O_APTL(APTLRegions.size());
		// Mixture Property Fields
		PtrList<volScalarField> molWeightMix_APTL(APTLRegions.size());
		PtrList<volScalarField> densityMix_APTL(APTLRegions.size());
		PtrList<volScalarField> viscMix_APTL(APTLRegions.size());
		PtrList<volScalarField> velAvg_APTL(APTLRegions.size());
		// Local Condition Fields
		PtrList<volScalarField> RH_APTL(APTLRegions.size());
		PtrList<volScalarField> pH2_APTL(APTLRegions.size());
		PtrList<volScalarField> pH2OVap_APTL(APTLRegions.size());
		PtrList<volScalarField> pH2OSat_APTL(APTLRegions.size());
		// Liquid Water Fields
		PtrList<volScalarField> jS_APTL(APTLRegions.size());
		PtrList<volScalarField> pCap_APTL(APTLRegions.size());
		PtrList<volScalarField> djSDs_APTL(APTLRegions.size());
		PtrList<volScalarField> dpCapDs_APTL(APTLRegions.size());
		PtrList<volScalarField> swPhaseChangeLiqVap_APTL(APTLRegions.size());
		// Source Term Fields
		PtrList<volScalarField> sourceEnergyPhaseChange_APTL(APTLRegions.size());
		PtrList<volScalarField> sourceEnergyOhmicHeating_APTL(APTLRegions.size());
		PtrList<volScalarField> sourcePhaseChangeLiqVap_APTL(APTLRegions.size());
	
	// Solved Field Variables	
		// Electron Transport Fields
		PtrList<volScalarField> potElectron_APTL(APTLRegions.size());
		// Hydrogen Transport Fields
	   	PtrList<volScalarField> wH2_APTL(APTLRegions.size());
		// Water Vapour Transport Fields
		PtrList<volScalarField> wH2OVap_APTL(APTLRegions.size());
		// Energy Transport Fields
		PtrList<volScalarField> T_APTL(APTLRegions.size());
		// Mass Transport (Pressure) Fields
	   	PtrList<volScalarField> pMix_APTL(APTLRegions.size());
	   	// Liquid Water Transport Fields
		PtrList<volScalarField> sLiq_APTL(APTLRegions.size());


	// Coefficient Field Variables
		// Electron Transport Fields
		PtrList<volScalarField> condElectronEff_APTL(APTLRegions.size());
		// Hydrogen Transport Fields
		PtrList<volScalarField> diffH2MixEff_APTL(APTLRegions.size());
		PtrList<volScalarField> transportCoeffH2_APTL(APTLRegions.size());
		// Water Vapour Transport Fields
		PtrList<volScalarField> diffH2OVapMixEff_APTL(APTLRegions.size());
		PtrList<volScalarField> transportCoeffH2OVap_APTL(APTLRegions.size());
		// Energy Transport Fields
		PtrList<volScalarField> condThermalEff_APTL(APTLRegions.size());
		PtrList<volScalarField> transportCoeffEnergy_APTL(APTLRegions.size());
	  	// Mass Transport (Pressure) Fields
		PtrList<volScalarField> dynamicPerm_APTL(APTLRegions.size());
	 	// Liquid Water Transport Fields
		PtrList<volScalarField> diffH2OLiq_APTL(APTLRegions.size());

	// Source/sink Field Variables
		// Electron Transport Fields
		PtrList<volScalarField> sourceElectronTotal_APTL(APTLRegions.size());
		// Hydrogen Transport Fields
		PtrList<volScalarField> sourceMassH2_APTL(APTLRegions.size());
		// Water Vapour Transport Fields
		PtrList<volScalarField> sourceMassH2OVap_APTL(APTLRegions.size());
		// Energy Transport Fields 
		PtrList<volScalarField> sourceEnergyTotal_APTL(APTLRegions.size());
	  	// Mass Transport (Pressure) Fields
		PtrList<volScalarField> sourceMassTotal_APTL(APTLRegions.size());
	 	// Liquid Water Transport Fields
		PtrList<volScalarField> sourceMassH2OLiq_APTL(APTLRegions.size());
