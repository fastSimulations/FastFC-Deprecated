// Initilialise field pointer lists
	// Calculated Field Variables
		// Stored Gradients
		PtrList<volVectorField> pMixGrad_CPTL(CPTLRegions.size());
		PtrList<volVectorField> potElectronGrad_CPTL(CPTLRegions.size());
		// Morphology Fields
		PtrList<volScalarField> volFPore_CPTL(CPTLRegions.size());
		PtrList<volScalarField> volFPoreEff_CPTL(CPTLRegions.size());
		// Material Property Fields
		PtrList<volScalarField> kPermDry_CPTL(CPTLRegions.size());
		PtrList<volScalarField> kPermRelLiq_CPTL(CPTLRegions.size());
		PtrList<volScalarField> kPermRelGas_CPTL(CPTLRegions.size());
		// Nitrogen Transport Fields
		PtrList<volScalarField> wN2_CPTL(CPTLRegions.size());
		// Mole Fraction Fields
		PtrList<volScalarField> xO2_CPTL(CPTLRegions.size());
		PtrList<volScalarField> xN2_CPTL(CPTLRegions.size());
		PtrList<volScalarField> xH2OVap_CPTL(CPTLRegions.size());
		// Fluid Property Fields
		PtrList<volScalarField> viscO2_CPTL(CPTLRegions.size());
		PtrList<volScalarField> viscN2_CPTL(CPTLRegions.size());
		PtrList<volScalarField> viscH2OVap_CPTL(CPTLRegions.size());
		// Binary Diffusivity Fields
		PtrList<volScalarField> diffO2N2_CPTL(CPTLRegions.size());
		PtrList<volScalarField> diffO2H2O_CPTL(CPTLRegions.size());
		PtrList<volScalarField> diffN2H2O_CPTL(CPTLRegions.size());
		// Mixture Property Fields
		PtrList<volScalarField> molWeightMix_CPTL(CPTLRegions.size());
		PtrList<volScalarField> densityMix_CPTL(CPTLRegions.size());
		PtrList<volScalarField> viscMix_CPTL(CPTLRegions.size());
		PtrList<volScalarField> velAvg_CPTL(CPTLRegions.size());
		// Local Condition Fields
		PtrList<volScalarField> RH_CPTL(CPTLRegions.size());
		PtrList<volScalarField> pO2_CPTL(CPTLRegions.size());
		PtrList<volScalarField> pH2OVap_CPTL(CPTLRegions.size());
		PtrList<volScalarField> pH2OSat_CPTL(CPTLRegions.size());
		// Liquid Water Fields
		PtrList<volScalarField> jS_CPTL(CPTLRegions.size());
		PtrList<volScalarField> pCap_CPTL(CPTLRegions.size());
		PtrList<volScalarField> djSDs_CPTL(CPTLRegions.size());
		PtrList<volScalarField> dpCapDs_CPTL(CPTLRegions.size());
		PtrList<volScalarField> swPhaseChangeLiqVap_CPTL(CPTLRegions.size());
		// Source Term Fields
		PtrList<volScalarField> sourceEnergyPhaseChange_CPTL(CPTLRegions.size());
		PtrList<volScalarField> sourceEnergyOhmicHeating_CPTL(CPTLRegions.size());
		PtrList<volScalarField> sourcePhaseChangeLiqVap_CPTL(CPTLRegions.size());
	
	// Solved Field Variables	
		// Electron Transport Fields
		PtrList<volScalarField> potElectron_CPTL(CPTLRegions.size());
		// Hydrogen Transport Fields
	   	PtrList<volScalarField> wO2_CPTL(CPTLRegions.size());
		// Water Vapour Transport Fields
		PtrList<volScalarField> wH2OVap_CPTL(CPTLRegions.size());
		// Energy Transport Fields
		PtrList<volScalarField> T_CPTL(CPTLRegions.size());
		// Mass Transport (Pressure) Fields
	   	PtrList<volScalarField> pMix_CPTL(CPTLRegions.size());
	   	// Liquid Water Transport Fields
		PtrList<volScalarField> sLiq_CPTL(CPTLRegions.size());


	// Coefficient Field Variables
		// Electron Transport Fields
		PtrList<volScalarField> condElectronEff_CPTL(CPTLRegions.size());
		// Oxygen Transport Fields
		PtrList<volScalarField> diffO2MixEff_CPTL(CPTLRegions.size());
		PtrList<volScalarField> transportCoeffO2_CPTL(CPTLRegions.size());
		// Water Vapour Transport Fields
		PtrList<volScalarField> diffH2OVapMixEff_CPTL(CPTLRegions.size());
		PtrList<volScalarField> transportCoeffH2OVap_CPTL(CPTLRegions.size());
		// Energy Transport Fields
		PtrList<volScalarField> condThermalEff_CPTL(CPTLRegions.size());
		PtrList<volScalarField> transportCoeffEnergy_CPTL(CPTLRegions.size());
	  	// Mass Transport (Pressure) Fields
		PtrList<volScalarField> dynamicPerm_CPTL(CPTLRegions.size());
	 	// Liquid Water Transport Fields
		PtrList<volScalarField> diffH2OLiq_CPTL(CPTLRegions.size());

	// Source/sink Field Variables
		// Electron Transport Fields
		PtrList<volScalarField> sourceElectronTotal_CPTL(CPTLRegions.size());
		// Hydrogen Transport Fields
		PtrList<volScalarField> sourceMassO2_CPTL(CPTLRegions.size());
		// Water Vapour Transport Fields
		PtrList<volScalarField> sourceMassH2OVap_CPTL(CPTLRegions.size());
		// Energy Transport Fields 
		PtrList<volScalarField> sourceEnergyTotal_CPTL(CPTLRegions.size());
	  	// Mass Transport (Pressure) Fields
		PtrList<volScalarField> sourceMassTotal_CPTL(CPTLRegions.size());
	 	// Liquid Water Transport Fields
		PtrList<volScalarField> sourceMassH2OLiq_CPTL(CPTLRegions.size());
