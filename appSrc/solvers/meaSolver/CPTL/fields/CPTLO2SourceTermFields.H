// Calculate the O2 Mass Transport Source/Sink Field Variables
{
	// Determine the total O2 Mass Source Term Fields
	dimensionedScalar sourceNullMass_
	(
	 	"sourceNullMass_", 
		dimensionSet(sourceMassO2.dimensions()), 
		0.
	);

	sourceMassO2 = sourceNullMass_;
}
