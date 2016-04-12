// Calculate the H2O Vapour Transport Source/Sink Field Variables
{
	// Water Vapour Transport Fields
	if(0==cathodeSinglePhase)
	{
		sourceMassH2OVap = -sourcePhaseChangeLiqVap; 	
	}
	else
	{
		sourceMassH2OVap = -((Fast::Constants::molWeightH2O/(2.*Fast::Constants::faraday))*sourceVolCurrent);
	}
}
