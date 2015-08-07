//	Set up necessary Gradients (These are not stored)

volVectorField potElectronGrad = fvc::grad(potElectron);
volVectorField potProtonGrad = fvc::grad(potProton);

// Set the scalar for electronic conductivity
{
	condElectronEff = Foam::pow(condElectron_ACL, torSuppY_ACL);
}
// Set the scalar for the proton conductivity
{
	condProtonEff = Foam::pow(condProton_ACL, torIonY_ACL); 
}
