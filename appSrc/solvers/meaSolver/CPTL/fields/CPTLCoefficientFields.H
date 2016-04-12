// Coefficient Field Variables
	// Electron Transport Fields
	{
		condElectronEff = condElectronCarbFibInput_CPTL;
				// = condElectronCarbFibInput_CPTL*Foam::pow((1.-volFPore), torCarbFibXInput_CPTL);
	}
	// Oxygen Transport Fields
	{
		diffO2MixEff = ((1.-xO2)/((xN2/diffO2N2)+(xH2OVap/diffO2H2O)))
						*Foam::pow(volFPoreEff,torPoreYInput_CPTL);
	}
	// Water Vapour Transport Fields
	{
		diffH2OVapMixEff = (1.-xH2OVap)/((xN2/diffN2H2O)+(xO2/diffO2H2O))
							*Foam::pow(volFPoreEff,torPoreYInput_CPTL);
	}
	// Energy Transport Fields
	{
		condThermalEff = condThermalCarbFibInput_CPTL;
				// = condThermalCarbFibInput_CPTL*Foam::pow((1-volFPore),torCarbFibYInput_CPTL);
	}
	
	// Volume Averaged thermal conductivity - test after implementation of main code
	//((volFPt*condThermalPtInput_CPTL)+(volFSupp*condThermalSuppInput_CPTL)+(volFIon*condThermalIonInput_CPTL))*Foam::pow((volFPt+volFSupp+volFIon),(torSuppYInput_CPTL*(condThermalSuppInput_CPTL/(condThermalSuppInput_CPTL+condThermalIonInput_CPTL)))+(torIonYInput_CPTL*(condThermalIonInput_CPTL/(condThermalSuppInput_CPTL+condThermalIonInput_CPTL))));
  		
	// Mass Transport (Pressure) Fields
	{	
		dynamicPerm = densityMix*(kPermRelGas/viscMix);
	}
 	// Liquid Water Transport Fields
	{	
		diffH2OLiq = (Fast::Constants::densityH2OLiq*kPermRelLiq/Fast::Constants::viscH2OLiqRef);
	}
