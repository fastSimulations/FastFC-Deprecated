// Calculate the volumetric current production
if(1==anodeKineticsTHV)
{
	sourceVolCurrentLeft_ =
	 	tafelHeyrovskyVolmer
		(
			sourceVolCurrent,
			pMix,
			xH2,
			T,
			sAC,
			actPotLeft_,
			exchCurrTVRefInput_ACL,
			exchCurrHVRefInput_ACL,
			concH2RefInput_ACL,
			reactionOrderH2Input_ACL,
			kappaTHVInput_ACL
		);	
	
	sourceVolCurrentRight_ =
	 	tafelHeyrovskyVolmer
		(
			sourceVolCurrent,
			pMix,
			xH2,
			T,
			sAC,
			actPotRight_,
			exchCurrTVRefInput_ACL,
			exchCurrHVRefInput_ACL,
			concH2RefInput_ACL,
			reactionOrderH2Input_ACL,
			kappaTHVInput_ACL
		);
}
else
{
	sourceVolCurrentLeft_ =
		butlerVolmerDiscrete
		(
 			sourceVolCurrent,
			T,
			sAC,
			actPotLeft_,
			pMix,
			xH2,
			henrySolH2,
			nElectronsInput_ACL,
			symmetryFactorBVInput_ACL,
			exchCurrBVRefInput_ACL,
			actEnergyBVInput_ACL,
			tempBVRefInput_ACL,
			concH2RefInput_ACL
		);

	sourceVolCurrentRight_ =
		butlerVolmerDiscrete
		(
 			sourceVolCurrent,
			T,
			sAC,
			actPotRight_,
			pMix,
			xH2,
			henrySolH2,
			nElectronsInput_ACL,
			symmetryFactorBVInput_ACL,
			exchCurrBVRefInput_ACL,
			actEnergyBVInput_ACL,
			tempBVRefInput_ACL,
			concH2RefInput_ACL
		);
}
