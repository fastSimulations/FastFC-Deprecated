// Calculate the Energy Transport Source/Sink Field Variables
{
	// sourceEnergyOhmicHeating
	sourceEnergyOhmicHeating = Foam::pow(potProtonGrad.component(1),2.)*condProtonEff;

	// Determine the total energy source/sink term
	sourceEnergyTotal = sourceEnergyOhmicHeating;
}
