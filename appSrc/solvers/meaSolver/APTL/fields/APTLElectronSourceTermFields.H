// Calculate the Electrony Transport Source/Sink Field Variables
{
	// Determine the total electron source/sink term
	dimensionedScalar sourceNullElectron_
	(
	 	"sourceNullElectron_", 
		dimensionSet(sourceElectronTotal.dimensions()), 
		0.
	);

	sourceElectronTotal = sourceNullElectron_;
}
