//	Set up necessary Gradients (These are not stored)

volVectorField potElectronGrad = fvc::grad(potElectron);
volVectorField potProtonGrad = fvc::grad(potProton);

// Set the scalar for electronic conductivity
{
	condElectronEff = Foam::pow(condElectron_CCL, torSuppY_CCL);
}
// Set the scalar for the proton conductivity
{
	condProtonEff = Foam::pow(condProton_CCL, torIonY_CCL);
}
