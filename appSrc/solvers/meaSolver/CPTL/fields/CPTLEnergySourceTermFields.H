// Calculate the Energy Transport Source/Sink Field Variables
{
	// sourceEnergyPhaseChange
	sourceEnergyPhaseChange = Fast::Constants::latentHeatH2OVap*sourcePhaseChangeLiqVap;
	// sourceEnergyOhmicHeating
	sourceEnergyOhmicHeating = Foam::pow(potElectronGrad.component(1),2.)*condElectronEff;
	// Determine the total energy source/sink term
	sourceEnergyTotal = sourceEnergyPhaseChange + sourceEnergyOhmicHeating;
}
