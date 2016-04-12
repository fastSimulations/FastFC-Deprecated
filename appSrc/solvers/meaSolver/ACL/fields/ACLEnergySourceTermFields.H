// Calculate the Energy Transport Source/Sink Field Variables
{
	// sourceEnergyPhaseChange
	sourceEnergyPhaseChange = Fast::Constants::latentHeatH2OVap*sourcePhaseChangeLiqVap;
	// sourceEnergyOhmicHeating
	sourceEnergyOhmicHeating = Foam::pow(potElectronGrad.component(1),2.)*condElectronEff 
								+Foam::pow(potProtonGrad.component(1),2.)*condProtonEff;
	// sourceEnergyEntropicHeating
	sourceEnergyEntropicHeating.internalField() = 0.;
	sourceEnergyEntropicHeating.boundaryField() = 0.;

	// Determine the total energy source/sink term
	sourceEnergyTotal = sourceEnergyPhaseChange + sourceEnergyOhmicHeating + sourceEnergyEntropicHeating;
}
