// Initialise electron field pointer lists
#include <APTLFieldPointerLists.H>

// Populate APTL field pointer lists
forAll(APTLRegions, zoneID)
{
	Info<< nl << "    Reading transportProperties file in constants/" << APTLRegionNames[zoneID] << endl;
    // Read in all related dictionaries
    #include <APTLReadComponentDictionary.H>
	
    // Assign any Dimensioned Scalars for input values or field initialisation
    #include <APTLDimScalarAssign.H>
	
    // Create Field Objects for solved variables
    #include <APTLSolvedFieldObjects.H>

    // Create Field Objects for transport coefficient variables
    #include <APTLCoefficientFieldObjects.H>

    // Create Field Objects for source term variables
    #include <APTLSourceTermFieldObjects.H>

	// Create Field Objects for the calculated variables
	#include <APTLGeneralFieldObjects.H>
}
