// Calculate the H2O Liquid Source/Sink Field Variables
{
	// Source Term Intermediate Fields
	{
		if(0==anodeSinglePhase)
		{
			volScalarField ratePhaseChangeH2OEvap_ = 
			(
			 	((condRateInput_ACL*volFPoreEff)*xH2OVap)
				/(2.*Fast::Constants::gasConstantR*T)
			)*(1+swPhaseChangeLiqVap)*Foam::mag(swPhaseChangeLiqVap);

			volScalarField ratePhaseChangeH2OCond_ = 
			(
			 	((evapRateInput_ACL*volFPore*sLiq)*Fast::Constants::densityH2OLiq)
				/(2.*Fast::Constants::molWeightH2O)
			)*(1-swPhaseChangeLiqVap)*Foam::mag(swPhaseChangeLiqVap);

			// sourceEnergyPhaseChangeLiqVap
			sourcePhaseChangeLiqVap = Fast::Constants::molWeightH2O
									*(ratePhaseChangeH2OEvap_+ratePhaseChangeH2OCond_)*(pMix*xH2OVap-pH2OSat);
		}
		else
		{
			dimensionedScalar sourceNullPhaseChange_ 
			(
			 	"sourceNullPhaseChange_", 
				dimensionSet(sourcePhaseChangeLiqVap.dimensions()), 
				0.
			);
			
			// sourceEnergyPhaseChangeLiqVap
			sourcePhaseChangeLiqVap = sourceNullPhaseChange_;
		}
	}
}
