// Initialise electron field pointer lists
#include "CCLFieldPointerLists.H"

// Populate CCL field pointer lists
forAll(CCLRegions, zoneID)
{
	Info<< nl << "    Reading transportProperties file in constants/" << CCLRegionNames[zoneID] << endl;
    // Read in all related dictionaries
    #include "CCLReadComponentDictionary.H"
	
    // Assign any Dimensioned Scalars for input values or field initialisation
    #include "CCLDimScalarAssign.H"
	
    // Create Field Objects for solved variables
    #include "CCLSolvedFieldObjects.H"

    // Create Field Objects for transport coefficient variables
    #include "CCLCoefficientFieldObjects.H"

    // Create Field Objects for source term variables
    #include "CCLSourceTermFieldObjects.H"

	// Create Field Objects for calculated variables
	#include "CCLCalculatedFieldObjects.H"
}
