// Calculate the volumetric current production
if(1==anodeKineticsTHV)
{
	volScalarField sourceVolCurrentOld_ = sourceVolCurrent;
	sourceVolCurrent = 
 		tafelHeyrovskyVolmer
		(
			sourceVolCurrent,
			pMix,
			xH2,
			T,
			sAC,
			actPot,
			exchCurrTVRefInput_ACL,
			exchCurrHVRefInput_ACL,
			concH2RefInput_ACL,
			reactionOrderH2Input_ACL,
			kappaTHVInput_ACL
		);	
	sourceVolCurrent = sourceVolCurrentOld_ + 0.001*(sourceVolCurrent-sourceVolCurrentOld_);
}
else
{
	sourceVolCurrent = 
		butlerVolmerDiscrete
		(
 			sourceVolCurrent,
			T,
			sAC,
			actPot,
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
