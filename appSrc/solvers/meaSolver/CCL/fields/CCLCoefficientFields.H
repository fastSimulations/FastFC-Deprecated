// Coefficient Field Variables
	// Electron Transport Fields
	{
		condElectronEff = condElectronSuppInput_CCL*Foam::pow(volFSupp,torSuppYInput_CCL);
	}
	// Proton Transport Fields
	{
		condProtonEff = condProtonIon*Foam::pow(volFIon,torIonYInput_CCL);
	}
	// Oxygen Transport Fields
	{
		diffO2MixEff = ((1.-xO2)/((xN2/diffO2N2)+(xH2OVap/diffO2H2O)))
						*Foam::pow(volFPoreEff,torPoreYInput_CCL);
	}
	// Water Vapour Transport Fields
	{
		diffH2OVapMixEff = (1.-xH2OVap)/((xN2/diffN2H2O)+(xO2/diffO2H2O))
							*Foam::pow(volFPoreEff,torPoreYInput_CCL);
	}
	// Energy Transport Fields
	{
		condThermalEff =  condThermalSuppInput_CCL;
	}
	
	// Volume Averaged thermal conductivity - test after implementation of main code
	//((volFPt*condThermalPtInput_CCL)+(volFSupp*condThermalSuppInput_CCL)+(volFIon*condThermalIonInput_CCL))*Foam::pow((volFPt+volFSupp+volFIon),(torSuppYInput_CCL*(condThermalSuppInput_CCL/(condThermalSuppInput_CCL+condThermalIonInput_CCL)))+(torIonYInput_CCL*(condThermalIonInput_CCL/(condThermalSuppInput_CCL+condThermalIonInput_CCL))));
  		
	// Mass Transport (Pressure) Fields
	{	
		dynamicPerm = densityMix*(kPermRelGas/viscMix);
	}
 	// Liquid Water Transport Fields
	{	
		diffH2OLiq = (Fast::Constants::densityH2OLiq*kPermRelLiq/Fast::Constants::viscH2OLiqRef);
	}
