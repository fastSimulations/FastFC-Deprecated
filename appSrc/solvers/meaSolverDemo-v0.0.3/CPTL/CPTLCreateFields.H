// Initialise electron field pointer lists
#include "CPTLFieldPointerLists.H"

// Populate CPTL field pointer lists
forAll(CPTLRegions, zoneID)
{
	Info<< nl << "    Reading transportProperties file in constants/" << CPTLRegionNames[zoneID] << endl;
    // Read in all related dictionaries
    #include "CPTLReadComponentDictionary.H"
	
    // Assign any Dimensioned Scalars for input values or field initialisation
    #include "CPTLDimScalarAssign.H"
	
    // Create Field Objects for solved variables
    #include "CPTLSolvedFieldObjects.H"

    // Create Field Objects for transport coefficient variables
    #include "CPTLCoefficientFieldObjects.H"

    // Create Field Objects for source term variables
    #include "CPTLSourceTermFieldObjects.H"

	// Create Field Objects for the calculated variables
	#include "CPTLGeneralFieldObjects.H"
}
