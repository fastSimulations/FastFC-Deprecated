//	Set up necessary Gradients (These are not stored)

volVectorField potProtonGrad = fvc::grad(potProton);

// Set the proton conductivity
{
	condProtonEff = Foam::pow(condProton_PMEM,torIonY_PMEM);
}
