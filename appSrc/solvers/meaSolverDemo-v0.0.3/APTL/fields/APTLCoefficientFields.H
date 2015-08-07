// Start of New file
//
//
// Calculate the Coefficient Field Variables
	// Electron Transport Fields (1D)
	condElectronEff = condElectronCarbFibInput_APTL;
			// = condElectronCarbFibInput_APTL*Foam::pow((1.-volFPore), torCarbFibXInput_APTL);
	// Hydrogen Transport Fields
	diffH2MixEff = ((1.-xH2)/((xN2/diffH2N2)+(xH2OVap/diffH2H2O)))*Foam::pow(volFPoreEff,torPoreYInput_APTL);
	// Water Vapour Transport Fields
	diffH2OVapMixEff = (1.-xH2OVap)/((xN2/diffN2H2O)+(xH2/diffH2H2O))*Foam::pow(volFPoreEff,torPoreYInput_APTL);
	// Energy Transport Fields
	condThermalEff = condThermalCarbFibInput_APTL;
			// = condThermalCarbFib_APTL*Foam::pow((1-volFPore),torCarbFibYInput_APTL);
  	// Mass Transport (Pressure) Fields
	dynamicPerm = densityMix*(kPermRelGas/viscMix);
 	// Liquid Water Transport Fields
	diffH2OLiq = (Fast::Constants::densityH2OLiq*kPermRelLiq/Fast::Constants::viscH2OLiqRef);
