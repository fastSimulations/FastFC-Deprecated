// Calculate the volumetric current production
if(1==cathodeKineticsAgglomerate)
{
	volScalarField sourceVolCurrentOld_ = sourceVolCurrent;
 
 	sourceVolCurrent = 
 		butlerVolmerAgglomerate
		(
			sourceVolCurrent,
			T,
			sAC,
			volFPore,
			actPot,
			diffO2Ion,
			pMix,
			xO2,
			henrySolO2,
			tIonFilm,
			sVAgg,
			diffO2IonAgg,
			nElectronsInput_CCL,
			symmetryFactorInput_CCL,
			exchCurrBVRefInput_CCL,
			actEnergyBVInput_CCL,
			tempBVRefInput_CCL,
			concO2BVRefInput_CCL,
			radiusAggInput_CCL 
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
			xO2,
			henrySolO2,
			nElectronsInput_CCL,
			symmetryFactorInput_CCL,
			exchCurrBVRefInput_CCL,
			actEnergyBVInput_CCL,
			tempBVRefInput_CCL,
			concO2BVRefInput_CCL
		);
}
