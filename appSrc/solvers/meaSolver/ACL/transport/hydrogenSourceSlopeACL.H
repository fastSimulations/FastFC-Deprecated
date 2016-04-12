// Calculate the volumetric current production
if(1==anodeKineticsTHV)
{
	sourceVolHydrogenLeft_ = -(Fast::Constants::molWeightH2/(2.*Fast::Constants::faraday))*
	 	tafelHeyrovskyVolmer
		(
			sourceVolCurrent,
			pMix,
			xH2Left_,
			T,
			sAC,
			actPot,
			exchCurrTVRefInput_ACL,
			exchCurrHVRefInput_ACL,
			concH2RefInput_ACL,
			reactionOrderH2Input_ACL,
			kappaTHVInput_ACL
		);	
	
	sourceVolHydrogenRight_ = -(Fast::Constants::molWeightH2/(2.*Fast::Constants::faraday))*
	 	tafelHeyrovskyVolmer
		(
			sourceVolCurrent,
			pMix,
			xH2Right_,
			T,
			sAC,
			actPot,
			exchCurrTVRefInput_ACL,
			exchCurrHVRefInput_ACL,
			concH2RefInput_ACL,
			reactionOrderH2Input_ACL,
			kappaTHVInput_ACL
		);
}
else
{
	sourceVolHydrogenLeft_ = -(Fast::Constants::molWeightH2/(2.*Fast::Constants::faraday))*
		butlerVolmerDiscrete
		(
 			sourceVolCurrent,
			T,
			sAC,
			actPot,
			pMix,
			xH2Left_,
			henrySolH2,
			nElectronsInput_ACL,
			symmetryFactorBVInput_ACL,
			exchCurrBVRefInput_ACL,
			actEnergyBVInput_ACL,
			tempBVRefInput_ACL,
			concH2RefInput_ACL
		);

	sourceVolHydrogenRight_ = -(Fast::Constants::molWeightH2/(2.*Fast::Constants::faraday))*
		butlerVolmerDiscrete
		(
 			sourceVolCurrent,
			T,
			sAC,
			actPot,
			pMix,
			xH2Right_,
			henrySolH2,
			nElectronsInput_ACL,
			symmetryFactorBVInput_ACL,
			exchCurrBVRefInput_ACL,
			actEnergyBVInput_ACL,
			tempBVRefInput_ACL,
			concH2RefInput_ACL
		);
}
