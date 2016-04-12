// Calculate the O2 Mass Transport Source/Sink Field Variables
{
	// Determine the total O2 Mass Source Term Fields
	sourceMassO2 = (Fast::Constants::molWeightO2/(4.*Fast::Constants::faraday))*sourceVolCurrent;
}
