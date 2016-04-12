// Calculate the H2 Mass Transport Source/Sink Field Variables
{
	// Determine the total H2 Mass Source Term Fields
	dimensionedScalar sourceNullMass_
	(
	 	"sourceNullMass_", 
		dimensionSet(sourceMassH2.dimensions()), 
		0.
	);

	sourceMassH2 = sourceNullMass_;
}
