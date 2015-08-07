// Coefficient Field Variables
	// Electron Transport Fields
	{
		condElectronEff = condElectronSuppInput_ACL*Foam::pow(volFSupp,torSuppYInput_ACL);
	}
	// Proton Transport Fields
	{
		condProtonEff = condProtonIon*Foam::pow(volFIon,torIonYInput_ACL);
	}
	// Hydrogen Transport Fields
	{
		diffH2MixEff = ((1.-xH2)/((xN2/diffH2N2)+(xH2OVap/diffH2H2O)))
						*Foam::pow(volFPoreEff,torPoreYInput_ACL);
	}
	// Water Vapour Transport Fields
	{
		diffH2OVapMixEff = (1.-xH2OVap)/((xN2/diffN2H2O)+(xH2/diffH2H2O))
							*Foam::pow(volFPoreEff,torPoreYInput_ACL);
	}
	// Energy Transport Fields
	{
		condThermalEff =  condThermalSuppInput_ACL;
	}
	
	// Volume Averaged thermal conductivity - test after implementation of main code
	//((volFPt*condThermalPtInput_ACL)+(volFSupp*condThermalSuppInput_ACL)+(volFIon*condThermalIonInput_ACL))*Foam::pow((volFPt+volFSupp+volFIon),(torSuppYInput_ACL*(condThermalSuppInput_ACL/(condThermalSuppInput_ACL+condThermalIonInput_ACL)))+(torIonYInput_ACL*(condThermalIonInput_ACL/(condThermalSuppInput_ACL+condThermalIonInput_ACL))));
  		
	// Mass Transport (Pressure) Fields
	{	
		dynamicPerm = densityMix*(kPermRelGas/viscMix);
	}
 	// Liquid Water Transport Fields
	{	
		diffH2OLiq = (Fast::Constants::densityH2OLiq*kPermRelLiq/Fast::Constants::viscH2OLiqRef);
	}
