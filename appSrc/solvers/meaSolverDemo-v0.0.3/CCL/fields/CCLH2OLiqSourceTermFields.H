// Calculate the Liquid Water Source/Sink Field Variables
{
	// Set equal to the phase change source/sink between Liquid and Vapour
	if(0==cathodeSinglePhase) // 1 == means single phase
	{
		sourceMassH2OLiq = sourcePhaseChangeLiqVap
				-((Fast::Constants::molWeightH2O/(2.*Fast::Constants::faraday))*sourceVolCurrent);

	}
	else
	{
		sourceMassH2OLiq = sourcePhaseChangeLiqVap; // note that it is set to zero
	}
}
